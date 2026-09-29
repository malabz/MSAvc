#define _POSIX_C_SOURCE 200809L
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>

#if defined(__linux__)
#  include <sys/sendfile.h>
#endif

static std::vector<std::string>        contig_lines;
static std::unordered_set<std::string> contig_seen;

struct FileInfo {
    int fd;
    off_t data_start;
    size_t file_size;
    const char* data;
    std::vector<std::string> sample_names;
    bool has_samples;
    std::string filename;
};

static inline void trim_newline(std::string& s) {
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r'))
        s.pop_back();
}

// 在已 mmap 的内存上扫头部：无固定缓冲、不会截断、O(1) 分配
static bool read_header(FileInfo& info) {
    info.data_start = 0;
    info.has_samples = false;

    const char* const data = info.data;
    const size_t size      = info.file_size;

    size_t pos = 0;
    while (pos < size) {
        const char* const line_start = data + pos;

        const char* line_end = (const char*)memchr(line_start, '\n', size - pos);
        size_t line_len;
        bool has_nl;
        if (line_end) {
            line_len = (size_t)(line_end - line_start);
            has_nl   = true;
        } else {
            line_end = data + size;
            line_len = size - pos;
            has_nl   = false;
        }

        if (line_len == 0 || line_start[0] != '#') {
            info.data_start = pos;
            break;
        }

        if (line_len >= 8 && memcmp(line_start, "##contig", 8) == 0) {
            std::string ctg(line_start, line_len + (has_nl ? 1 : 0));
            if (contig_seen.insert(ctg).second)
                contig_lines.push_back(std::move(ctg));
        }
        else if (line_len >= 6 && memcmp(line_start, "#CHROM", 6) == 0) {
            const char* p   = line_start;
            const char* end = line_start + line_len;
            int tabs = 0;
            while (p < end && tabs < 9) {
                if (*p == '\t') ++tabs;
                ++p;
            }

            while (p < end) {
                const char* q = (const char*)memchr(p, '\t', end - p);
                if (!q) q = end;

                const char* q_trim = q;
                while (q_trim > p && (q_trim[-1] == '\r' || q_trim[-1] == '\n'))
                    --q_trim;

                if (q_trim > p)
                    info.sample_names.emplace_back(p, (size_t)(q_trim - p));

                if (q == end) break;
                p = q + 1;
            }

            info.has_samples = !info.sample_names.empty();
            info.data_start = has_nl ? pos + line_len + 1 : size;
            break;
        }

        if (!has_nl) break;
        pos += line_len + 1;
    }

    return true;
}

static void output_header(const std::vector<std::string>& sample_names) {
    printf("##fileformat=VCFv4.1\n");
    for (const auto& ctg : contig_lines) printf("%s", ctg.c_str());
    printf("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate allele count, for each ALT allele, in the same order as listed\">\n");
    printf("##INFO=<ID=VT,Number=1,Type=String,Description=\"Type of small variant\">\n");
    printf("##INFO=<ID=VLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n");
    printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for (const auto& name : sample_names) printf("\t%s", name.c_str());
    printf("\n");
    fflush(stdout);
}

static void merge_fast(std::vector<FileInfo>& files,
                       const std::vector<std::string>& sample_names) {
    output_header(sample_names);
    int out_fd = fileno(stdout);

    for (auto& info : files) {
        off_t  off    = info.data_start;
        size_t remain = info.file_size - info.data_start;

#if defined(__linux__)
        // Linux：文件→文件真零拷贝
        while (remain > 0) {
            ssize_t sent = sendfile(out_fd, info.fd, &off, remain);
            if (sent < 0) { if (errno == EINTR) continue; break; }
            remain -= sent;
        }
#else
        // macOS 及其他 POSIX：mmap 页已是 page cache，write 即近似零拷贝
        const char* p = info.data + off;
        while (remain > 0) {
            ssize_t w = write(out_fd, p, remain);
            if (w < 0) { if (errno == EINTR) continue; break; }
            p      += w;
            remain -= (size_t)w;
        }
#endif
    }
}

static void merge_reorder(std::vector<FileInfo>& files) {
    std::unordered_map<std::string, int> global_sample_map;
    std::vector<std::string> global_sample_names;
    for (const auto& info : files) {
        for (const auto& sample : info.sample_names) {
            if (global_sample_map.find(sample) == global_sample_map.end()) {
                global_sample_map[sample] = global_sample_names.size();
                global_sample_names.push_back(sample);
            }
        }
    }
    std::sort(global_sample_names.begin(), global_sample_names.end());
    for (size_t i = 0; i < global_sample_names.size(); ++i)
        global_sample_map[global_sample_names[i]] = i;

    output_header(global_sample_names);

    const int global_n = global_sample_names.size();

    struct Mapping {
        std::vector<int> global_to_local;
    };
    std::vector<Mapping> mappings(files.size());
    for (size_t fi = 0; fi < files.size(); ++fi) {
        auto& info = files[fi];
        int local_n = info.sample_names.size();
        auto& map = mappings[fi];
        map.global_to_local.assign(global_n, -1);
        for (int i = 0; i < local_n; ++i) {
            int gidx = global_sample_map[info.sample_names[i]];
            map.global_to_local[gidx] = i;
        }
    }

    const size_t OUT_BUF_SIZE = 1 << 20;
    std::vector<char> out_buf(OUT_BUF_SIZE);
    size_t out_pos = 0;

    auto flush_buffer = [&]() {
        if (out_pos) {
            fwrite(out_buf.data(), 1, out_pos, stdout);
            out_pos = 0;
        }
    };

    auto append_string = [&](const char* str, size_t len) {
        if (out_pos + len + 1 > OUT_BUF_SIZE) flush_buffer();
        memcpy(out_buf.data() + out_pos, str, len);
        out_pos += len;
    };

    auto append_char = [&](char c) {
        if (out_pos + 1 > OUT_BUF_SIZE) flush_buffer();
        out_buf[out_pos++] = c;
    };

    for (size_t fi = 0; fi < files.size(); ++fi) {
        const auto& info = files[fi];
        const auto& map = mappings[fi];
        const char* data = info.data;
        size_t len = info.file_size;
        const char* end = data + len;
        const char* line_start = data + info.data_start;

        struct Field { const char* ptr; size_t len; };
        std::vector<Field> fixed_fields(9);
        std::vector<Field> gt_fields(info.sample_names.size());

        while (line_start < end) {
            const char* line_end = (const char*)memchr(line_start, '\n', end - line_start);
            if (!line_end) line_end = end;
            if (*line_start == '#') { line_start = line_end + 1; continue; }

            const char* p = line_start;
            for (int i = 0; i < 9; ++i) {
                const char* q = (const char*)memchr(p, '\t', line_end - p);
                if (!q) q = line_end;
                fixed_fields[i] = {p, (size_t)(q - p)};
                p = q + 1;
                if (p > line_end) break;
            }

            for (size_t i = 0; i < info.sample_names.size(); ++i) {
                const char* q = (const char*)memchr(p, '\t', line_end - p);
                if (!q) q = line_end;
                gt_fields[i] = {p, (size_t)(q - p)};
                p = q + 1;
                if (p > line_end) break;
            }

            for (int i = 0; i < 9; ++i) {
                append_string(fixed_fields[i].ptr, fixed_fields[i].len);
                if (i < 8) append_char('\t');
            }

            for (int g = 0; g < global_n; ++g) {
                int local = map.global_to_local[g];
                append_char('\t');
                if (local < 0) {
                    append_char('.');
                } else {
                    const auto& f = gt_fields[local];
                    append_string(f.ptr, f.len);
                }
            }
            append_char('\n');

            line_start = line_end + 1;
        }
        flush_buffer();
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "\033[31mUsage: %s <vcf1> [vcf2 ...]\033[0m\n", argv[0]);
        return 1;
    }

    std::vector<FileInfo> files;
    files.reserve(argc - 1);

    for (int i = 1; i < argc; ++i) {
        FileInfo info;
        info.filename = argv[i];
        info.fd = open(argv[i], O_RDONLY);
        if (info.fd < 0) {
            fprintf(stderr, "\033[31mopen: %s\033[0m\n", strerror(errno));
            return 1;
        }
        struct stat st;
        if (fstat(info.fd, &st) < 0) {
            fprintf(stderr, "\033[31mfstat: %s\033[0m\n", strerror(errno));
            return 1;
        }
        info.file_size = st.st_size;
        if (info.file_size == 0) {
            fprintf(stderr, "\033[33mWarning: skipping empty file %s\033[0m\n", argv[i]);
            close(info.fd);
            continue;
        }
        info.data = (const char*)mmap(NULL, info.file_size, PROT_READ, MAP_PRIVATE, info.fd, 0);
        if (info.data == MAP_FAILED) {
            fprintf(stderr, "\033[31mmmap: %s\033[0m\n", strerror(errno));
            return 1;
        }
        if (!read_header(info)) {
            fprintf(stderr, "\033[31mFailed to read header of %s\033[0m\n", argv[i]);
            return 1;
        }
        files.push_back(std::move(info));
    }

    if (files.empty()) {
        fprintf(stderr, "\033[31mNo non-empty input files.\033[0m\n");
        return 1;
    }

    bool any_samples = false;
    for (const auto& info : files) if (info.has_samples) { any_samples = true; break; }

    if (!any_samples) {
        merge_fast(files, {});
    } else {
        bool same_order = true;
        const auto& first = files[0].sample_names;
        for (size_t i = 1; i < files.size(); ++i) {
            if (files[i].sample_names != first) { same_order = false; break; }
        }
        if (same_order) {
            merge_fast(files, first);
        } else {
            merge_reorder(files);
        }
    }

    for (auto& info : files) {
        munmap((void*)info.data, info.file_size);
        close(info.fd);
    }

    return 0;
}
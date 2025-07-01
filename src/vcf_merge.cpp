#include <iostream>
#include <cstring>
#include <filesystem>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <unordered_map>
#include <sstream>
#include <vector>
#include <algorithm>


std::unordered_map<std::string, int> mp;
std::vector<int> this_order;
std::vector<std::string> contigs, var_data;
int map_new_item = 0;
char line[1024000];

void find_CHROM_line_and_detect_sample(FILE *file)
{
    while (fgets(line, sizeof(line), file))
    {
        if (strncmp(line, "##contig", 8) == 0)
        {
            contigs.emplace_back(line);
        }
        else if (strncmp(line, "#CHROM", 6) == 0)
        {
            std::stringstream ss(line);
            std::string col;

            for (int i = 0; i <= 8; i++)
            { // skip #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
                getline(ss, col, '\t');
            }

            while (getline(ss, col, '\t'))
            {
                col.erase(std::remove(col.begin(), col.end(), '\n'), col.end());
                if(mp.find(col) == mp.end())
                {
                    mp[col] = map_new_item ++;
                }
            }
            break;
        }
    }
    rewind(file);
}

void print_header()
{
    printf("##fileformat=VCFv4.1\n");
    for(auto &contig: contigs) printf("%s", contig.c_str());
    printf("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate allele count, for each ALT allele, in the same order as listed\">\n");
    printf("##INFO=<ID=VT,Number=1,Type=String,Description=\"Type of small variant\">\n");
    printf("##INFO=<ID=VLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n");
    printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    fflush(stdout);
    for(auto &item: var_data) printf("\t%s", item.c_str());
    putchar('\n');
}

void process_middle()
{
    var_data.resize(mp.size());
    for(auto &item: mp) var_data[item.second] = item.first;
    std::sort(var_data.begin(), var_data.end());
    for(int i = 0; i < var_data.size(); ++ i) mp[var_data[i]] = i;
}

void print_line_by_mp(FILE *file)
{
    this_order.clear();
    while (fgets(line, sizeof(line), file))
    {
        if (strncmp(line, "#CHROM", 6) == 0)
        {
            std::stringstream ss(line);
            std::string col;

            for (int i = 0; i <= 8; i++)
            { // skip #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
                getline(ss, col, '\t');
            }

            while (getline(ss, col, '\t')) // take order
            {
                col.erase(std::remove(col.begin(), col.end(), '\n'), col.end());
                this_order.emplace_back(mp[col]);
            }
            continue;
        }
        else if(line[0] == '#') continue;
        // split string        
        std::fill(var_data.begin(), var_data.end(), ".");
        std::stringstream ss(line);
        std::string col;
        
        // fprintf(stderr, "%s\n", line);
        for (int i = 0; i <= 8; i++)
        {
            getline(ss, col, '\t');
            if (i) putchar('\t');
            col.erase(std::remove(col.begin(), col.end(), '\n'), col.end());
            printf("%s", col.c_str());
        }
        auto now_order = this_order.begin();
        while (getline(ss, col, '\t'))
        {
            col.erase(std::remove(col.begin(), col.end(), '\n'), col.end());
            var_data[*now_order] = col;
            now_order ++;
        }
        for(auto &item: var_data) printf("\t%s", item.c_str());
        printf("\n");
    }
    fclose(file);
}

int main(int argc, char **argv)
{
    if(argc == 1) { fprintf(stderr, "Error: No vcf file found.\nWill print to stdout\nUsage: %s vcf1 vcf2 ...\n", argv[0]); exit(1); }
    std::ios::sync_with_stdio(false);
    // found CHROM line
    contigs.reserve(argc);
    FILE **files = (FILE**)malloc(sizeof(FILE*) * argc);
    if(files == nullptr) { fprintf(stderr, "Error: No memory. Program will exit.\n"); exit(1); }
    for(int argi = 1; argi < argc; ++ argi)
    {
        files[argi] = fopen(argv[argi], "r");
        if(files[argi] == nullptr) { fprintf(stderr, "Error: can not open file %s. Program will exit.\n", argv[argi]); exit(1); }
        find_CHROM_line_and_detect_sample(files[argi]);
    }
    // must seperate these two progresses
    process_middle();
    print_header();
    for(int argi = 1; argi < argc; ++ argi)
        print_line_by_mp(files[argi]);
    return 0;
}
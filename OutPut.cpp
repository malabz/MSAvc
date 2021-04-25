#include "OutPut.hpp"

#include <fstream>
#include <algorithm>
#include <iterator>
#include <cstring>
#include <sstream>

void output(const Fasta &infile, std::unordered_map<Mutation, std::vector<size_t>, hash> &mutations,
            const std::string &outfile_name, const std::string &auxfile_name)
{
    static constexpr char separator[] = ", ";

    const std::vector<std::string> &matrix = infile.sequences;
    const std::string &centre = matrix[0];
    const size_t col = centre.size();
    const size_t row = matrix.size();

    bool *mutated_sequence_mark = new bool[row]();
    for (const auto &mutation : mutations)
        for (auto sequence_index : mutation.second)
            mutated_sequence_mark[sequence_index] = true;
    unsigned mutated_sequence_number = std::count(mutated_sequence_mark, mutated_sequence_mark + row, true);

    size_t *map_to_centre_site = new size_t[col];
    for (size_t i = 0, index = 0; i != col; ++i)
        if (centre[i] != '-') map_to_centre_site[i] = index++;

    size_t *map_to_mutated_sequence = new size_t[row];
    for (size_t i = 0, index = 0; i != row; ++i)
        if (mutated_sequence_mark[i]) map_to_mutated_sequence[i] = index++;

    std::ofstream ofs(outfile_name);
    std::ofstream ofs_aux(auxfile_name);
    if (!ofs || !ofs_aux)
    {
        std::cout << "cannot write\n";
        exit(0);
    }

    std::ostringstream oss;
    oss << "#CHROM\tPOS\tpos\tID\tREF\tALT\tTYPE\tNUM\tQUAL\tFILTER\tINFO\tFORMAT\t";
    for (size_t i = 0; i != row; ++i)
        if (mutated_sequence_mark[i]) oss << infile.identifications[i] << '\t';

    std::string line = oss.str();
    if (line.back() == '\t') line.pop_back();
    ofs << line << '\n';

    ofs_aux << '\n';

    using pair_type = std::pair<Mutation, std::vector<size_t>>;
    std::vector<pair_type> mutation_vector(mutations.cbegin(), mutations.cend());
    std::sort(mutation_vector.begin(), mutation_vector.end(),
            [](const pair_type &lhs, const pair_type &rhs)
            {
                if (lhs.first.first == rhs.first.first)
                    return lhs.second.front() < rhs.second.front();
                return lhs.first.first < rhs.first.first;
            });

    bool *current_mutated_sequence_mark = new bool[mutated_sequence_number];
    for (const auto &current : mutation_vector)
    {
        const auto &mutation = current.first;
        const auto &occurrence = current.second;

        ofs << "1\t"
            << map_to_centre_site[mutation.first] << "\t"
            << mutation.first << "\t"
               ".\t"
            << mutation.l << '\t'
            << mutation.r << '\t'
            << to_string(mutation.flag) << '\t'
            << occurrence.size() << '\t'
            << ".\t"
               ".\t"
               ".\t"
               "GT\t";

        memset(current_mutated_sequence_mark, 0, sizeof(bool) * mutated_sequence_number);
        for (auto i : occurrence)
            current_mutated_sequence_mark[map_to_mutated_sequence[i]] = true;

        std::copy(current_mutated_sequence_mark,
                  current_mutated_sequence_mark + mutated_sequence_number - 1,
                  std::ostream_iterator<bool>(ofs, "\t"));
        ofs << current_mutated_sequence_mark[mutated_sequence_number - 1] << '\n';

        ofs_aux << mutation.first << "p;";
        for (size_t i = 0; i != occurrence.size(); ++i)
            ofs_aux << occurrence[i] + 1 << "p;";
        ofs_aux << '\n';
    }

    delete[] map_to_mutated_sequence;
    delete[] map_to_centre_site;
    delete[] current_mutated_sequence_mark;
    delete[] mutated_sequence_mark;
}

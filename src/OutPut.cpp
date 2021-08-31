#include "Arguments.hpp"
#include "OutPut.hpp"

#include <fstream>
#include <algorithm>
#include <iterator>
#include <cstring>
#include <sstream>

using mutation_and_occurrences_iterator = std::map<mutation::Mutation, std::vector<size_t>>::const_iterator;
using string_iterator = std::string::const_iterator;

static void output_dot_if_empty(std::ostream &os, string_iterator first, string_iterator last)
{
    if (*first == '^')
        ++first;

    if (first == last)
        os << '.';
    else
        for (; first != last; ++first)
            os << *first;
}

static void output_reference(std::ostream &os, string_iterator first, string_iterator last)
{
    if (*first == '^')
        ++first;

    if (first == last)
        os << '.';
    else
        for (; first != last; ++first)
            if (*first != '-')
                os << *first;
}

static void output_alterations(std::ostream &os, mutation_and_occurrences_iterator first, mutation_and_occurrences_iterator last)
{
    output_dot_if_empty(os, first->first.counterpart.cbegin(), first->first.counterpart.cend());

    for (++first; first != last; ++first)
    {
        os << ',';
        output_dot_if_empty(os, first->first.counterpart.cbegin(), first->first.counterpart.cend());
    }
}

static size_t compute_variation_length(const mutation::Mutation &m, const std::string &lhs)
{
    size_t ref_length;
    size_t alt_length;

    switch (m.variation_type)
    {
    case mutation::SUB:
        return m.last - m.first;

    case mutation::DEL:
        return m.last - m.first - 1;

    case mutation::INS:
        return m.counterpart.size() - 1;

    // the max length of ref and alt excluding the prefixing character
    case mutation::REP:
        ref_length = 0;
        for (size_t i = m.first + 1; i != m.last; ++i)
            if (lhs[i] != '-') ++ref_length;

        alt_length = m.counterpart.size() - 1;

        return std::max(ref_length, alt_length);
        break;

    default:
        throw std::logic_error("unexpected variation_type");
        break;
    }

    return 0;
}

static void output_info(std::ostream &os, mutation_and_occurrences_iterator first, mutation_and_occurrences_iterator last, size_t variation_length)
{
    const mutation::Mutation &current_mutation = first->first;

    size_t count = 0;
    for (auto i = first; i != last; ++i) count += i->second.size();

    os << "AC=" << count
       << ";VT=" << mutation::to_string(first->first.variation_type)
       << ";VELN=" << variation_length
          ;
}

void output(const utils::Fasta &infile, std::map<mutation::Mutation, std::vector<size_t>> &mutations)
{
    const std::string &reference = infile.sequences[arguments::reference_index];
    const std::string &reference_identifier = infile.identifications[arguments::reference_index];
    const size_t col = reference.size();
    const size_t row = infile.sequences.size();

    bool *does_the_sequence_contain_mutations = new bool[row]();
    for (const auto &mutation : mutations)
        for (auto sequence_index : mutation.second)
            does_the_sequence_contain_mutations[sequence_index] = true;
    const unsigned mutated_sequence_number = std::count(does_the_sequence_contain_mutations, does_the_sequence_contain_mutations + row, true);

    size_t *map_to_ref_site = new size_t[col];
    for (size_t i = 0, index = 0; i != col; ++i)
        if (reference[i] != '-') map_to_ref_site[i] = index++;

    size_t *map_to_mutated_sequence = new size_t[row];
    for (size_t i = 0, index = 0; i != row; ++i)
        if (does_the_sequence_contain_mutations[i]) map_to_mutated_sequence[i] = index++;

    std::ofstream ofs(arguments::outfile_name);
    if (!ofs) { std::cerr << "cannot open " << arguments::outfile_name << '\n'; exit(0); }

    static constexpr char *const header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";

    ofs << header << reference_identifier;
    for (size_t i = 0; i != row; ++i)
        if (i != arguments::reference_index) ofs << '\t' << infile.identifications[i];

    unsigned *current_mutated_sequence_mark = new unsigned[row];
    for (auto i = mutations.cbegin(), j = i; i != mutations.cend(); i = j)
    {
        for (++j; j != mutations.cend(); ++j)
            if (within_the_same_line(i->first, j->first) == false) break;
        const size_t count_of_mutations_of_this_line = std::distance(i, j);

        const auto &first_mutation = i->first;

        ofs << reference_identifier
            << '\t' << map_to_ref_site[first_mutation.first]
            << "\t."
               "\t"; output_reference(ofs, reference.cbegin() + first_mutation.first, reference.cbegin() + first_mutation.last);
        ofs << '\t'; output_alterations(ofs, i, j);
        ofs << "\t."
               "\t."
               "\t"; output_info(ofs, i, j, compute_variation_length(first_mutation, reference));
        ofs << "\tGT"
               "\t0"
            ;

        memset(current_mutated_sequence_mark, 0, sizeof(unsigned) * row);
        for (auto k = i; k != j; ++k)
        {
            const size_t current_alteration_number = std::distance(i, k) + 1;
            for (const auto occurrence : k->second)
                current_mutated_sequence_mark[occurrence] = current_alteration_number;
        }

        for (size_t k = 0; k != row; ++k)
            if (k != arguments::reference_index)
                ofs << '\t' << current_mutated_sequence_mark[k];
        ofs << '\n';

        // ofs_aux << first_mutation.first << "p;";
        // for (size_t i = 0; i != occurrence.size(); ++i)
        //     ofs_aux << occurrence[i] + 1 << "p;";
        // ofs_aux << '\n';
    }

    delete[] map_to_mutated_sequence;
    delete[] map_to_ref_site;
    delete[] current_mutated_sequence_mark;
    delete[] does_the_sequence_contain_mutations;
}

bool within_the_same_line(const mutation::Mutation &lhs, const mutation::Mutation &rhs) noexcept
{
    return lhs.variation_type == mutation::SUB
            && rhs.variation_type == mutation::SUB
            && lhs.first == rhs.first
            && lhs.last == rhs.last
            ;
}

#include "Mutation.hpp"
#include "Fasta.hpp"

#include <algorithm>
#include <iterator>

using matrix_type = std::vector<std::string>;

mutation::Mutation extract_mutation(const std::string &lhs, const std::string &rhs, size_t &position);

std::unordered_map<mutation::Mutation, std::vector<size_t>, mutation::hash> search_in(const utils::Fasta &infile, size_t reference_index)
{
    const matrix_type &matrix = infile.sequences;
    const std::string &reference = matrix[reference_index];
    const size_t col = reference.size();
    const size_t row = matrix.size();

    std::unordered_map<mutation::Mutation, std::vector<size_t>, mutation::hash> mutations;

    for (size_t i = 0; i != row; ++i)
    {
        if (i == reference_index) continue;

        const std::string &sequence = matrix[i];

        for (size_t j = 0; j != col; ++j)
            if (reference[j] != sequence[j])
                mutations[extract_mutation(reference, sequence, j)].push_back(i);
    }

    return mutations;
}

// assume lhs != rhs
unsigned get_flag(char lhs, char rhs)
{
    if (lhs == '-') return mutation::INS;

    else if (rhs == '-') return mutation::DEL;

    return mutation::SUB;
}

// assume *lhs_first != *rhs_first
template<typename lhs_iter, typename rhs_iter>
unsigned get_flag(lhs_iter lhs_first, lhs_iter lhs_last,
                  rhs_iter rhs_first)
{
    unsigned first_flag = get_flag(*lhs_first, *rhs_first);

    for (; lhs_first != lhs_last; ++lhs_first, ++rhs_first)
        if (*lhs_first != *rhs_first && get_flag(*lhs_first, *rhs_first) != first_flag)
            return MIX;

    return first_flag;
}

mutation::Mutation extract_mutation(const std::string &lhs, const std::string &rhs, size_t &position)
{
    mutation::Mutation mutation;
    mutation.first = position;

    size_t last  = position + 1;
    while (last != lhs.size() &&
           (lhs[last] != rhs[last] || lhs[last] == '-' && rhs[last] == '-'))
        ++last;

    while (lhs[last] == '-' && rhs[last] == '-') --last;

    mutation.flag = get_flag(lhs.cbegin() + mutation.first, lhs.cbegin() + last,
                             rhs.cbegin() + mutation.first);

    while (lhs[mutation.first] == '-' || rhs[mutation.first] == '-')
        --mutation.first;

    for (size_t i = mutation.first; i != last; ++i)
        if (lhs[i] != '-' || rhs[i] != '-')
        {
            mutation.l.push_back(lhs[i]);
            mutation.r.push_back(rhs[i]);
        }

    mutation.l = mutation::remove_gap(mutation.l.cbegin(), mutation.l.cend());
    mutation.r = mutation::remove_gap(mutation.r.cbegin(), mutation.r.cend());

    position = last - 1;
    return mutation;
}

bool contain_del(const matrix_type &sequences, size_t column)
{
    for (size_t i = 1; i != sequences.size(); ++i)
        if (sequences[i][column] == '-') return true;

    return false;
}

bool contain_snp(const matrix_type &sequences, size_t column)
{
    const char c = sequences[0][column];
    for (size_t i = 1; i != sequences.size(); ++i)
        if (sequences[i][column] != c) return true;

    return false;
}

size_t mutation::hash::operator()(const Mutation &mutation) const noexcept
{
    size_t hash_value = 0;
    hash_value ^= mutation.first;
    hash_value ^= mutation.flag;

    std::hash<std::string> hash;
    hash_value ^= hash(mutation.l);
    hash_value ^= hash(mutation.r);

    return hash_value;
}

std::string to_string(unsigned flag)
{
    static const std::string mix = "MIX";
    static const std::string ins = "INS";
    static const std::string del = "DEL";
    static const std::string snp = "SUB";

    if (flag & mutation::MIX) return mix;
    if (flag & mutation::INS) return ins;
    if (flag & mutation::DEL) return del;
    return snp;
}

bool mutation::Mutation::operator==(const Mutation &rhs) const noexcept
{
    return first == rhs.first &&
           l     == rhs.l     &&
           r     == rhs.r;
}

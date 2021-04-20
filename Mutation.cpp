#include "Mutation.hpp"
#include "Fasta.hpp"

#include <algorithm>
#include <iterator>

void to_lower(Fasta &file);

using matrix_type = std::vector<std::string>;

unsigned is_mutation(const std::string &lhs, const std::string &rhs, size_t &position);
Mutation extract_mutation(const std::string &lhs, const std::string &rhs, size_t &position);

std::unordered_map<Mutation, std::vector<size_t>, hash> search_in(const Fasta &infile)
{
    const matrix_type &matrix = infile.sequences;
    const std::string &centre = matrix[0];
    const size_t col = centre.size();
    const size_t row = matrix.size();

    std::unordered_map<Mutation, std::vector<size_t>, hash> mutations;

    for (size_t i = 1; i != row; ++i)
    {
        const std::string &curr_sequence = matrix[i];

        for (size_t j = 0; j != col; ++j)
            if (centre[j] != curr_sequence[j])
                mutations[extract_mutation(centre, curr_sequence, j)].push_back(i);
    }

    return mutations;
}

unsigned is_mutation(const std::string &lhs, const std::string &rhs, size_t &position)
{
    unsigned flag = 0;

    if (position != lhs.size() - 1 && lhs[position + 1] == '-')
    {
        for (size_t i = position + 1; i != lhs.size() && lhs[i] == '-'; ++i)
            if (rhs[i] != '-')
            {
                flag |= INS;
                break;
            }
    }

    if (lhs[position] != rhs[position])
    {
        if (rhs[position] == '-') flag |= DEL;
        else flag |= SNP;
    }

    return flag;
}

Mutation extract_mutation(const std::string &lhs, const std::string &rhs, size_t &position)
{
    Mutation mutation;
    mutation.first = position;
    size_t last  = position + 1;

    while (lhs[last] != rhs[last] || lhs[last] == '-' && rhs[last] == '-') ++last;

    if (last - mutation.first != 1)
        mutation.flag = MIX;
    else if (lhs[position] == '-')
        mutation.flag = INS;
    else if (rhs[position] == '-')
        mutation.flag = DEL;
    else
        mutation.flag = SNP;

    while (lhs[mutation.first] == '-' || rhs[mutation.first] == '-')
        --mutation.first;

    for (size_t i = mutation.first; i != last; ++i)
        if (lhs[i] != '-' || rhs[i] != '-')
        {
            mutation.l.push_back(lhs[i]);
            mutation.r.push_back(rhs[i]);
        }

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

size_t hash::operator()(const Mutation &mutation) const noexcept
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
    static const std::string snp = "SNP";

    if (flag & MIX) return mix;

    if (flag & INS)
        return ins;

    if (flag & DEL)
        return del;

    return snp;
}

bool Mutation::operator== (const Mutation &rhs) const noexcept
{
    return first == rhs.first &&
           flag  == rhs.flag &&
           l     == rhs.l &&
           r     == rhs.r;
}

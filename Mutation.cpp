#include "Mutation.hpp"
#include "Fasta.hpp"

#include <algorithm>
#include <iterator>

void to_lower(Fasta &file);

using matrix_type = std::vector<std::string>;

unsigned is_mutation(const std::string &lhs, const std::string &rhs, size_t &position);
Mutation extract_mutation(const std::string &lhs, const std::string &rhs, size_t &position, unsigned flag);

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
            if (unsigned flag = is_mutation(centre, curr_sequence, j))
                mutations[extract_mutation(centre, curr_sequence, j, flag)].push_back(i);
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
        if (rhs[position] == '-') flag |= DEL;
        else flag |= SNP;

    return flag;
}

Mutation extract_mutation(const std::string &lhs, const std::string &rhs, size_t &position, unsigned flag)
{
    Mutation mutation;
    mutation.first = position;
    mutation.last  = position + 1;
    mutation.flag  = flag;

    if (flag & INS)
    {
        while (lhs[mutation.last] == '-') ++mutation.last;
        mutation.snp = rhs[position];
        mutation.str = remove_gap(rhs.substr(position + 1, mutation.last - mutation.first - 1));
    }
    else if (flag & DEL)
    {
        while (mutation.last != lhs.size() && rhs[mutation.last] == '-' && lhs[mutation.last] != '-')
            ++mutation.last;
    }
    else if (flag & SNP)
    {
        while (mutation.last != lhs.size() && lhs[mutation.last] != rhs[mutation.last] && lhs[mutation.last] != '-')
            ++mutation.last;
        mutation.str = remove_gap(rhs.substr(position, mutation.last - mutation.first));
    }

    position = mutation.last - 1;
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

std::string remove_gap(const std::string &str)
{
    std::string gap_removed;
    gap_removed.reserve(str.size());

    for (size_t i = 0; i != str.size(); ++i)
        if (str[i] != '-') gap_removed.push_back(str[i]);

    return gap_removed;
}

size_t hash::operator()(const Mutation &mutation) const noexcept
{
    size_t hash_value = 0;
    hash_value ^= mutation.first;
    hash_value ^= mutation.last;
    hash_value ^= mutation.flag;

    if (mutation.flag & INS)
    {
        hash_value ^= std::hash<std::string>()(mutation.str);
        if (mutation.flag & SNP)
            hash_value ^= mutation.snp;
    }
    else if (mutation.flag & SNP)
    {
        hash_value ^= std::hash<std::string>()(mutation.str);
    }

    return hash_value;
}

bool Mutation::operator==(const Mutation &rhs) const noexcept
{
    if (first != rhs.first || last != rhs.last || flag != rhs.flag)
        return false;

    if (flag & INS)
    {
        if (str != rhs.str) return false;
        if (flag & SNP && snp != rhs.snp) return false;
        return true;
    }

    if (flag & SNP)
        return str == rhs.str;

    return true;
}

std::string Mutation::to_string() const
{
    if (flag & INS)
        return flag & DEL ? str : snp + str;

    if (flag & DEL)
        return "*";

    return str;
}

std::string to_string(unsigned flag)
{
    if (flag & INS)
        return flag & ~INS ? "MIX" : "INS";

    if (flag & DEL)
        return "DEL";

    return "SNP";
}

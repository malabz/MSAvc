#include "Mutation.hpp"
#include "Fasta.hpp"
#include "Arguments.hpp"

using matrix_type = std::vector<std::string>;

static unsigned get_flag(char lhs, char rhs)
{
    // assert lhs != rhs

    if (lhs == '-')
        return mutation::INS;
    else if (rhs == '-')
        return mutation::DEL;
    else
        return mutation::SUB;
}

template<typename lhs_iter, typename rhs_iter>
static unsigned get_flag(lhs_iter lhs_first, lhs_iter lhs_last, rhs_iter rhs_first)
{
    // assert *lhs_first != *rhs_first

    const unsigned flag = get_flag(*lhs_first++, *rhs_first++);
    for (; lhs_first != lhs_last; ++lhs_first, ++rhs_first)
        if (*lhs_first != *rhs_first && flag != get_flag(*lhs_first, *rhs_first))
            return mutation::REP;

    return flag;
}

mutation::Mutation extract_mutation(const std::string &lhs, const std::string &rhs, size_t &position)
{
    // assert lhs.size() == rhs.size()
    // assert postition < lhs.size()
    // assert lhs[position] != rhs[position]

    mutation::Mutation mutation;
    mutation.first = position;

    for (++position; position != lhs.size() && (lhs[position] != rhs[position] || lhs[position] == '-'); ++position) ;

    // position will never reach 0 in this loop
    while (lhs[position - 1] == '-' && rhs[position - 1] == '-')
        --position;

    mutation.variation_type = get_flag(lhs.cbegin() + mutation.first, lhs.cbegin() + position, rhs.cbegin() + mutation.first);

    if (mutation.variation_type != mutation::SUB)
        for (--mutation.first; lhs[mutation.first] == '-'; --mutation.first) ;

    for (size_t i = mutation.first; i != position; ++i)
        if (rhs[i] != '-') mutation.counterpart.push_back(rhs[i]);

    mutation.last = position;
    for (; lhs[mutation.last - 1] == '-'; --mutation.last) ;

    return mutation;
}

std::map<mutation::Mutation, std::vector<size_t>> mutation::search_in(const utils::Fasta &infile)
{
    const std::string &reference = infile.sequences[arguments::reference_index];
    const size_t col = reference.size();
    const size_t row = infile.sequences.size();

    std::map<mutation::Mutation, std::vector<size_t>> mutations;
    for (size_t i = 0; i != row; ++i)
        if (infile.sequences[i].size() == col)
            for (size_t j = 0; j != col; )
                if (reference[j] != infile.sequences[i][j])
                    mutations[extract_mutation(reference, infile.sequences[i], j)].push_back(i);
                else
                    ++j;

    return mutations;
}

size_t mutation::hash::operator()(const Mutation &mutation) const noexcept
{
    static constexpr std::hash<size_t> int_hash;
    static constexpr std::hash<std::string> string_hash;

    size_t hash_value = 0;

    // ...

    return hash_value;
}

bool mutation::Mutation::operator<(const Mutation &rhs) const noexcept
{
    if (first != rhs.first) return first < rhs.first;

    if (variation_type != rhs.variation_type) return variation_type < rhs.variation_type;

    if (last != rhs.last) return last < rhs.last;

    return counterpart < rhs.counterpart;
}

std::string mutation::to_string(unsigned flag)
{
    static const std::string rep = "REP";
    static const std::string ins = "INS";
    static const std::string del = "DEL";
    static const std::string sub = "SUB";

    switch (flag)
    {
    case mutation::SUB: return sub;
    case mutation::DEL: return del;
    case mutation::INS: return ins;
    case mutation::REP: return rep;

    default:
        throw std::logic_error("unexpected variation type");
    }
}

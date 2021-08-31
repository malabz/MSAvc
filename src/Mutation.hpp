#pragma once

#include "Fasta.hpp"

#include <string>
#include <vector>
#include <map>

namespace mutation
{

    enum flag_enum : unsigned
    {
        SUB = 1u << 0,
        DEL = 1u << 1,
        INS = 1u << 2,
        REP = 1u << 3
    };

    struct Mutation
    {
        unsigned variation_type;
        size_t first; // indexes the first character of the ref field
        size_t last;  // indexes the character past the last character of the ref field
        std::string counterpart;

        bool operator<(const Mutation &rhs) const noexcept;
    };

    // std::ostream &operator <<(std::ostream &os, const Mutation &m);

    class hash
    {
    public:
        size_t operator()(const Mutation &mutation) const noexcept;
    };

    std::map<Mutation, std::vector<size_t>> search_in(const utils::Fasta &infile);

    std::string to_string(unsigned flag);

}

#pragma once

#include "Fasta.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>

namespace mutation
{

    enum flag_enum : unsigned
    {
        SUB = 1u << 0,
        DEL = 1u << 1,
        INS = 1u << 2,
        MIX = 1u << 3
    };

    struct Mutation
    {
        unsigned flag;
        size_t first;
        std::string l;    // insertion or consecutive snp
        std::string r;

        bool operator== (const Mutation &rhs) const noexcept;
    };

    // std::ostream &operator <<(std::ostream &os, const Mutation &m);

    class hash
    {
    public:
        size_t operator()(const Mutation &mutation) const noexcept;
    };

    std::unordered_map<Mutation, std::vector<size_t>, hash> search_in(const utils::Fasta &infile, size_t reference_index);

    std::string to_string(unsigned flag);

    template<typename InputIterator>
    std::string remove_gap(InputIterator first, InputIterator last)
    {
        std::string str;
        for (; first != last; ++first)
            if (*first != '-')
                str.push_back(*first);

        return str;
    }
}

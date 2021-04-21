#pragma once

#include "Fasta.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>

constexpr unsigned SNP = 0x1;
constexpr unsigned DEL = 0x2;
constexpr unsigned INS = 0x4;
constexpr unsigned MIX = 0x8;

struct Mutation
{
    unsigned flag;
    size_t first;
    std::string l;    // insertion or consecutive snp
    std::string r;

    bool operator== (const Mutation &rhs) const noexcept;
    std::string to_string() const;
};

std::ostream &operator <<(std::ostream &os, const Mutation &m);

class hash
{
public:
    size_t operator()(const Mutation &mutation) const noexcept;
};

std::unordered_map<Mutation, std::vector<size_t>, hash> search_in(const Fasta &infile);

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

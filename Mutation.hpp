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
    size_t first, last;
    unsigned flag;
    char snp;           // valid if (flag & (INS | SNP))
    std::string str;    // insertion or consecutive snp

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

std::string remove_gap(const std::string &str);

std::string to_string(unsigned flag);

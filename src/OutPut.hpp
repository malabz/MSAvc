#pragma once

#include <vector>

#include "Mutation.hpp"
#include "Fasta.hpp"

void output(const utils::Fasta &infile, std::map<mutation::Mutation, std::vector<size_t>> &mutations);

bool within_the_same_line(const mutation::Mutation &lhs, const mutation::Mutation &rhs) noexcept;

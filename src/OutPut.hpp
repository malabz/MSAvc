#pragma once

#include <cstdint>
#include <vector>

#include "Mutation.hpp"
#include "Fasta.hpp"
#include "MultipleAlignmentFormat.hpp"

void output(utils::MultipleAlignmentFormat const &maf, mut::MutationContainer const &mutations);

// 0-based [begin, end)
void output_sub_block(utils::MultipleAlignmentFormat const &infile,std::uint64_t begin, std::uint64_t end);
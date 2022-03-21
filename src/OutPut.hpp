#pragma once

#include <vector>

#include "Mutation.hpp"
#include "Fasta.hpp"
#include "MultipleAlignmentFormat.hpp"

void output(utils::MultipleAlignmentFormat const &maf, mut::MutationContainer const &mutations);

// 0-based [begin, end)
void output_sub_block(utils::MultipleAlignmentFormat const &infile, unsigned begin, unsigned end);

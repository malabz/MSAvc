#pragma once

#include "Mutation.hpp"
#include "Fasta.hpp"

#include <string>
#include <vector>

void output(const utils::Fasta &infile, std::unordered_map<mutation::Mutation, std::vector<size_t>, mutation::hash> &mutations,
            const std::string &outfile_name, const std::string &autfile_name);

#pragma once

#include "Mutation.hpp"
#include "Fasta.hpp"

#include <string>
#include <vector>

void output(const Fasta &infile, std::unordered_map<Mutation, std::vector<size_t>, hash> &mutations,
            const std::string &outfile_name, const std::string &autfile_name);

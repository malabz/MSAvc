// snp_indel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "OutPut.hpp"

#include <fstream>
#include <algorithm>

void to_lower(Fasta &file) noexcept;

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::cout << argv[0] << " infile outfile\n";
        exit(0);
    }

    Fasta infile(argv[1]);
    to_lower(infile);

    auto results = search_in(infile);
    output(infile, results, argv[2]);
}

void to_lower(Fasta &file) noexcept
{
    for (auto &sequence : file.sequences)
        transform(sequence.cbegin(), sequence.cend(), sequence.begin(), tolower);
}

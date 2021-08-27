// snp_indel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "OutPut.hpp"

#include <fstream>

int main(int argc, const char **argv)
{
    if (argc != 4)
    {
        std::cout << argv[0] << " infile outfile auxfile\n";
        exit(0);
    }

    std::ifstream ifs(argv[1]);
    if (!ifs)
    {
        std::cout << "cannot access file " << argv[1] << '\n';
        exit(0);
    }

    utils::Fasta infile(ifs);
    infile.transform_tolower();

    auto results = search_in(infile, 0);
    output(infile, results, argv[2], argv[3]);
}

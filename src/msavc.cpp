// snp_indel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "OutPut.hpp"
#include "Arguments.hpp"

#include <fstream>

int main(int argc, const char **argv)
{
    arguments::parse_arguments(argc, argv);

    std::ifstream ifs(arguments::infile_name);
    if (!ifs)
    {
        std::cout << "cannot access file " << argv[1] << '\n';
        exit(0);
    }

    utils::Fasta infile(ifs);
    infile.transform_tolower();
    infile.prefix_caret_to_sequences();

    auto results = mutation::search_in(infile);
    output(infile, results);
}

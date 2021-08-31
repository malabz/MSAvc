#include "Arguments.hpp"

#include <iostream>

std::string arguments::infile_name;
std::string arguments::outfile_name;
std::string arguments::auxfile_name;

size_t arguments::reference_index;

bool arguments::parse_arguments(unsigned argc, const char *const *const argv)
{
    if (argc < 3) argument_error();
    printf("   1st argument: %s\n", argv[0]);
    printf("fasta file name: %s\n", argv[1]);
    printf("  out file name: %s\n", argv[2]);

    infile_name = argv[1];
    outfile_name = argv[2];
    auxfile_name = "tmp";

    reference_index = 0;

    return true;
}

void arguments::argument_error()
{
    std::cerr << "argument error\n";
    print_usage();
    exit(0);
}

void arguments::print_usage()
{
}

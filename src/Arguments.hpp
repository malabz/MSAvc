#pragma once

#include <string>

namespace arguments
{

    extern std::string infile_name;
    extern std::string outfile_name;
    extern std::string auxfile_name;

    extern size_t reference_index;


    bool parse_arguments(unsigned argc, const char *const *const argv);

    void argument_error();

    void print_usage();

}

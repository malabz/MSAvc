// This file contains the 'main' function. Program execution begins and ends there.
//

#include <fstream>
#include <algorithm>

#include "OutPut.hpp"
#include "Arguments.hpp"

static void standardize(std::string &str)
{
    str.resize(str.size() + 1);
    for (unsigned i = str.size() - 1; i != 0; --i)
        str[i] = std::tolower(str[i - 1]);
    str.front() = '^';
}

static void standardize(std::vector<std::string> &vs)
{
    for (std::string &str : vs)
        standardize(str);
}

int main(int argc, char **argv)
{
    std::ios::sync_with_stdio(false);

    arguments::parse_arguments(argc, argv);

    std::ifstream ifs(arguments::infile_path);
    if (!ifs) { std::cerr << "cannot access file " << arguments::infile_path << '\n'; exit(1); }

    utils::MultipleAlignmentFormat infile;
 
    if (arguments::infile_in_fasta)
    {
        utils::Fasta fasta(ifs);
        arguments::check_arguments(fasta);
        infile.read(std::move(fasta));
    }
    else
    {
        infile.read(ifs);
        arguments::check_arguments(infile);
        infile.reverse_record_if_necessary(arguments::reference_index);
    }

    for (auto &record : infile.records)
    {
        standardize(record.sequences);
        record.build_map_if_necessary(arguments::reference_index);
    }

    arguments::print_arguments();

    auto const mutations_and_corresponding_occurrences = mut::search_in(infile, arguments::reference_index);
    output(infile, mutations_and_corresponding_occurrences);

    if (arguments::sub_block)
        output_sub_block(infile, arguments::lpos, arguments::rpos);

    std::cerr << "success!\n";

    return 0;
}

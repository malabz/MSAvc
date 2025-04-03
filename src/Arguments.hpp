#pragma once

#include <string>

#include "Fasta.hpp"
#include "MultipleAlignmentFormat.hpp"

namespace arguments
{

    extern std::string infile_path;
    extern std::string outfile_path;
    extern bool infile_in_fasta; // false if infile in maf

    extern std::string reference_name;
    extern std::string reference_genome_prefix;
    extern unsigned reference_index;

    extern bool genotype_matrix;
    extern bool combine_like_substitutions;

    // [lpos : rpos)
    extern unsigned lpos; // open
    extern unsigned rpos; // close

    extern unsigned minimum_alternative_allele_count_acceptable;
    extern unsigned maximum_alternative_allele_count_acceptable;
    extern unsigned minimum_variation_length_acceptable;
    extern unsigned maximum_variation_length_acceptable;
    extern bool variation_type_acceptable[];

    extern bool force;

    extern bool sub_block;
    extern std::string sub_block_outfile_path;

    extern bool check_duplicate;
    extern bool compress_bgz;

    extern unsigned buffer_size;

    void parse_arguments(unsigned argc, const char *const *argv);

    void check_arguments(utils::Fasta const &infile);
    void check_arguments(utils::MultipleAlignmentFormat const &infile);
    void check_arguments(); // some common operations

    void infile_format_unexpected();
    void argument_error(std::string const &message);

    void produce_version_message();
    void produce_help_message();

    void specify_infile_format();
    void deduce_subblock_file_path();

    void print_arguments();
    std::ostream &print_variation_type(std::ostream &os);

}

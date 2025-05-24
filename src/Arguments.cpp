#include <iostream>
#include <cstring>
#include <filesystem>

#include <boost/program_options.hpp>

#include "Arguments.hpp"
#include "Mutation.hpp"

static char constexpr version[]                                     = "v0.1.20250523.beta";
static char constexpr help_description[]                            = "";
static char constexpr version_description[]                         = "";
static char constexpr infile_description[]                          = "";
static char constexpr outfile_description[]                         = "";
static char constexpr reference_description[]                       = "";
static char constexpr genotype_matrix_description[]                 = "";
static char constexpr no_merge_sub_description[]                    = "";
static char constexpr position_filter_begin_description[]           = "";
static char constexpr position_filter_end_description[]             = "";
static char constexpr alternative_allele_count_filter_description[] = "";
static char constexpr variation_type_filter_description[]           = "";
static char constexpr variation_length_filter_description[]         = "";
static char constexpr force_descrition[]                            = "";
static char constexpr sub_block_description[]                       = "";
static char constexpr no_check_duplicate_description[]              = "";
static char constexpr compress_bgz_description[]                    = "";
static char constexpr buffer_size_description[]                     = "";

std::string arguments::infile_path;
std::string arguments::outfile_path;
bool arguments::infile_in_fasta;

std::string arguments::reference_name;
std::string arguments::reference_genome_prefix;
unsigned arguments::reference_index;

bool arguments::genotype_matrix;
bool arguments::combine_like_substitutions;

unsigned arguments::lpos;
unsigned arguments::rpos;

unsigned arguments::minimum_alternative_allele_count_acceptable;
unsigned arguments::maximum_alternative_allele_count_acceptable;
unsigned arguments::minimum_variation_length_acceptable;
unsigned arguments::maximum_variation_length_acceptable;
bool arguments::variation_type_acceptable[4];

bool arguments::force;

bool arguments::sub_block;
std::string arguments::sub_block_outfile_path;
bool arguments::check_duplicate;
bool arguments::compress_bgz;
unsigned arguments::buffer_size;

void arguments::specify_infile_format()
{
    std::filesystem::path const path(infile_path);
    if (path.has_extension() == false) infile_format_unexpected();

    std::string const extension = path.extension().string();
    if (extension.size() == 1) infile_format_unexpected(); // only a '.' in the end

    char const flag = std::tolower(extension[1]);
    if (flag != 'f' && flag != 'm') infile_format_unexpected();

    infile_in_fasta = flag == 'f';
}

void arguments::deduce_subblock_file_path()
{
    std::filesystem::path const path(outfile_path);
    if (path.has_stem() == false) {
        std::cerr << "output file path illegal: " << outfile_path << '\n';
        exit(1);
    }

    std::filesystem::path rhs = path.parent_path();
    rhs.append(path.stem().string() + "_subblock*.fasta");
    sub_block_outfile_path = rhs.string();
}

void arguments::parse_arguments(unsigned argc, const char *const *argv)
{
    // default to true
    memset(variation_type_acceptable, true, sizeof(variation_type_acceptable));

    auto const infile_option = boost::program_options::value<std::string>(&infile_path);
    infile_option->required();

    auto const outfile_option = boost::program_options::value<std::string>(&outfile_path);
    outfile_option->required();

    auto const reference_option = boost::program_options::value<std::string>(&reference_name);

    auto const position_filter_begin_option = boost::program_options::value<unsigned>(&lpos);
    auto const position_filter_end_option = boost::program_options::value<unsigned>(&rpos);
    position_filter_begin_option->default_value(1); // 1-based
    position_filter_end_option->default_value(std::numeric_limits<unsigned>::max() - 1); // -1 in case of overflow

    auto const minimum_alternative_allele_count_filter_option = boost::program_options::value<unsigned>(&minimum_alternative_allele_count_acceptable);
    minimum_alternative_allele_count_filter_option->default_value(1);

    auto const maximum_alternative_allele_count_filter_option = boost::program_options::value<unsigned>(&maximum_alternative_allele_count_acceptable);
    maximum_alternative_allele_count_filter_option->default_value(std::numeric_limits<unsigned>::max() - 1);

    auto const variation_type_filter_option = boost::program_options::value<std::vector<std::string>>();

    auto const minimum_variation_length_filter_option = boost::program_options::value<unsigned>(&minimum_variation_length_acceptable);
    minimum_variation_length_filter_option->default_value(1);

    auto const maximum_variation_length_filter_option = boost::program_options::value<unsigned>(&maximum_variation_length_acceptable);
    maximum_variation_length_filter_option->default_value(std::numeric_limits<unsigned>::max() - 1);

    auto const buffer_size_option = boost::program_options::value<unsigned>(&buffer_size);
    buffer_size_option->default_value(1 << 20);

    boost::program_options::options_description desc("allowed options");
    desc.add_options()
        ("help,h",                                                              help_description)
        ("version,v",                                                           version_description)
        ("in,i",                infile_option,                                  infile_description)
        ("out,o",               outfile_option,                                 outfile_description)
        ("reference,r",         reference_option,                               reference_description)
        ("genotype-matrix,g",                                                   genotype_matrix_description)
        ("nomerge-sub,n",                                                       no_merge_sub_description)
        ("filter-begin,b",      position_filter_begin_option,                   position_filter_begin_description)
        ("filter-end,e",        position_filter_end_option,                     position_filter_end_description)
        ("ac-greater,c",        minimum_alternative_allele_count_filter_option, alternative_allele_count_filter_description)
        ("ac-less,d",           maximum_alternative_allele_count_filter_option, alternative_allele_count_filter_description)
        ("filter-vt,t",         variation_type_filter_option,                   variation_length_filter_description)
        ("vl-greater,l",        minimum_variation_length_filter_option,         variation_length_filter_description)
        ("vl-less,m",           maximum_variation_length_filter_option,         variation_length_filter_description)
        ("force,f",                                                             force_descrition)
        ("sub-block,s",                                                         sub_block_description)
        ("no-duplicate-name,N",                                                 no_check_duplicate_description)
        ("compress-bgz,C",                                                      compress_bgz_description)
        ("buffer-size,B",       buffer_size_option,                             buffer_size_description)
    ;

    try
    {
        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);

        bool const help_flag = vm.find("help") != vm.cend();
        bool const version_flag = vm.find("version") != vm.cend();
        if (help_flag) produce_help_message();
        if (version_flag) produce_version_message();
        if (help_flag || version_flag) exit(0);

        boost::program_options::notify(vm);

        genotype_matrix = vm.find("genotype-matrix") != vm.cend();
        combine_like_substitutions = vm.find("nomerge-sub") != vm.cend();
        force = vm.find("force") != vm.cend();
        specify_infile_format();
        sub_block = vm.find("sub-block") != vm.cend();
        check_duplicate = vm.find("no-duplicate-name") == vm.cend();
        compress_bgz = vm.find("compress-bgz") != vm.cend();
        if (sub_block)
            deduce_subblock_file_path();

        auto const found = vm.find("filter-vt");
        if (found != vm.cend())
        {
            memset(variation_type_acceptable, false, sizeof(variation_type_acceptable));

            auto const &variation_type_filter_arguments = found->second.as<std::vector<std::string>>();
            if (variation_type_filter_arguments.empty())
                argument_error("no enough valid variation types provided\n");

            for (auto const &variation_type : variation_type_filter_arguments)
            {
                // auto const variation_type_specified = tolower(variation_type.front());
                if (false) ;
                else if (variation_type == "sub") variation_type_acceptable[mut::SUB] = true;
                else if (variation_type == "del") variation_type_acceptable[mut::DEL] = true;
                else if (variation_type == "ins") variation_type_acceptable[mut::INS] = true;
                else if (variation_type == "rep") variation_type_acceptable[mut::REP] = true;
                else argument_error("unexpected variation type: " + variation_type);
            }
        }
    }
    catch(std::exception const &e)
    {
        std::cerr << e.what() << ", see " << argv[0] << " --help for more information\n";
        exit(1);
    }

    // closed to open
    ++rpos; // TODO: determine this value
}

void arguments::check_arguments(utils::Fasta const &infile)
{
    unsigned const row = infile.sequences.size();

    if (row == 0) {
        std::cerr << "fasta format error\n";
        exit(1);
    }

    if (row == 1) {
        std::cerr << "only one sequence found\n";
        exit(1);
    }

    unsigned const col = infile.sequences[0].size();
    for (unsigned i = 1; i != row; ++i)
        if (infile.sequences[i].size() != col) {
            std::cerr << "sequences not with the same length\n";
            exit(1);
        }

    if (reference_name.size())
    {
        bool reference_found = false;
        for (unsigned i = 0; i != infile.sequences.size(); ++i)
            if (infile.names[i] == reference_name)
            {
                reference_index = i;
                reference_found = true;
                break;
            }

        if (reference_found == false)
            argument_error("reference name not found");
    }
    else
    {
        reference_index = 0;
        reference_name = infile.names[0];
    }

    if (rpos > col)
        rpos = col;

    if (lpos >= col)
        argument_error("position index out of bounds");

    check_arguments();
}

void arguments::check_arguments(utils::MultipleAlignmentFormat const &infile)
{
    if (infile.records.size() == 0) {
        std::cerr << "maf format error\n";
        exit(1);
    }

    for (auto const &record : infile.records)
    {
        unsigned const row = record.sequences.size();
        if (row < 2) continue;

        unsigned const col = record.sequences[0].size();
        for (unsigned i = 1; i != row; ++i)
            if (record.sequences[i].size() != col) {
                std::cerr << "sequences not with the same length\n";
                exit(1);
            }
    }

    if (reference_name.size())
    {
        auto const found = infile.name_to_index.find(reference_name);
        if (found == infile.name_to_index.cend())
            argument_error("reference name not found");

        reference_index = found->second;
    }
    else if(reference_genome_prefix.size())
    {
        reference_genome_prefix += '.';
        reference_index = std::numeric_limits<unsigned>::max() - 1;
    }
    else
    {
        reference_index = 0;
        reference_name = infile.names[0];
    }

    check_arguments();
}

void arguments::check_arguments()
{
    if (lpos >= rpos)
        argument_error("position index out of bounds");

    if (force == false)
    {
        if(! compress_bgz)
        {
            if (std::filesystem::exists(outfile_path))
                argument_error("file " + outfile_path + " exists; use --force or -f to demand an overwrite");

            if (sub_block && std::filesystem::exists(sub_block_outfile_path))
                argument_error("file " + sub_block_outfile_path + " exists; use --force or -f to demand an overwrite");
        }
        else
        {
            if (std::filesystem::exists(outfile_path + ".gz"))
                argument_error("file " + outfile_path + ".gz exists; use --force or -f to demand an overwrite");

            if (sub_block && std::filesystem::exists(sub_block_outfile_path + ".gz"))
                argument_error("file " + sub_block_outfile_path + ".gz exists; use --force or -f to demand an overwrite");
        }
    }
}

void arguments::infile_format_unexpected()
{
    argument_error("cannot specify the format of " + infile_path);
}

void arguments::argument_error(std::string const &message)
{
    std::cerr << "argument error: " << message << '\n';
    exit(1);
}

void arguments::produce_version_message()
{
    std::cerr << version << '\n';
}

void arguments::produce_help_message()
{
    std::cerr << "msavc " << version <<
        "\n   MSAvc: variation calling for genome-scale multiple sequence alignments"
        "\n   See https://github.com/malabz/msavc for the most up-to-date documentation."
        "\n"
        "\nUsage:"
        "\n   msavc -i <inputfile> -o <outputfile> -r <seqname> [options]"
        "\n"
        "\nOptions:"
        "\n   -i, --in <inputfile>              Specify the multi-FASTA/MAF input file"
        "\n   -o, --out <outputfile>            Specify the output VCF file name"
        "\n"
        "\n   -r, --reference <seqname>         Specify the reference sequence name during extracting variations"
        "\n"
        "\n   -R, --reference-genome <genome>   Specify the reference genome prefix during extracting variations"
        "\n                                     (deprecated, please use msavc_genome)"
        "\n"
        "\n   -g, --genotype-matrix             Output genotype matrix (default=off)"
        "\n"
        "\n   -n, --nomerge-sub                 Do not merge the SUB variations with the same \"POS\" and \"REF\" "
        "\n                                     into one row (default=off)"
        "\n"
        "\n   -b, --filter-begin <integer>      Filtration via the POS column by specifying an integer in terms "
        "\n                                     of the reference genome: \"-b 24\" keep the variations with POS>=24 "
        "\n                                     (default=1, 1-based index)"
        "\n"
        "\n   -e, --filter-end <integer>        Filtration via the POS column by specifying an integer: \"-e 1000\"" 
        "\n                                     keep the variations with POS<=1000 (default=the last base index, "
        "\n                                     1-based index)"
        "\n"
        "\n   -c, --ac-greater <integer>        Filtration via the AC tag in the INFO column by specifying an "
        "\n                                     integer: \"-c 10\" output the variations with AC>=10 (default=0)"
        "\n"
        "\n   -d, --ac-less <integer>           Filtration via the AC tag: \"-d 100\" output the variations with "
        "\n                                     AC<=100 (default=the total number of sequences)"
        "\n"
        "\n   -t, --filter-vt <variationtype>   Filtration via the VT tag in the INFO column by specifying one of "
        "\n                                     sub/ins/del/rep (lowercase) flags: \"-t sub\" only output the "
        "\n                                     substitution variations (default=off)"
        "\n"
        "\n   -l, --vl-greater <integer>        Filtration via the VLEN tag in the INFO column by specifying an "
        "\n                                     integer: \"-l 5\" output the variations with VLEN>=5bp (default=0)"
        "\n"
        "\n   -m, --vl-less <integer>           Filtration via the VLEN tag: \"-m 10\" output the variations with "
        "\n                                     VLEN<=10bp (default=length of reference genome)"
        "\n"
        "\n   -s, --sub-block                   Output MSA sub-block into FASTA file, \"-s\" option works only when "
        "\n                                     \"-b\" and \"-e\" are both specified: \"-b 24 -e 1000 -s\" produce a "
        "\n                                     sub MSA block; the slice interval is 24=<POS<=1000 in terms of "
        "\n                                     the reference genome (default=off)"
        "\n"
        "\n   -N, --no-duplicate-name           Do not check duplicate name in file (default=off)"
        "\n   -C, --compress-bgz                Compress the VCF output file. As the number of sequences and"
        "\n                                     variations increases, the VCF file with \"-g\" becomes super large."
        "\n"
        "\n   -B, --buffer-size                 String buffer size: For large files, set this size smaller in order to"
        "\n                                     lower the memory usage, only used in FASTA file (default=1048576)"
        "\n   -f, --force-overwrite             Overwrite existing file (default=off)"
        "\n   -h, --help                        Help message"
        "\n   -v, --version                     Version"
        "\n"
    ;
}

static std::string const &yes_or_no(bool which) noexcept
{
    static std::string const ans[2] = { "no", "yes" };
    return ans[which];
}

std::ostream &arguments::print_variation_type(std::ostream &os)
{
    for (unsigned i = 0, cnt = 0; i != 4; ++i)
        if (variation_type_acceptable[i])
        {
            if (cnt++) std::cerr << "; ";
            std::cerr << mut::mutation_types[i];
        }
    return os;
}

void arguments::print_arguments()
{
    std::cerr
        << "\ninput file path: "
        << "\n\t" << infile_path
        << "\noutput file path: "
        << "\n\t" << outfile_path << (compress_bgz ? ".gz" : "")
        << "\ninput file format: "
        << "\n\t" << (infile_in_fasta ? "fasta" : "maf")
        << "\ncheck duplicate: "
        << "\n\t" << yes_or_no(check_duplicate)
        << "\noutput compress gz:"
        << "\n\t" << yes_or_no(compress_bgz)
        << "\n"
        << "reference name: \n\t" << reference_name
        << "\ngenotype matrix output: "
        << "\n\t" << yes_or_no(genotype_matrix)
        << "\ncombine like substitutions: "
        << "\n\t" << yes_or_no(combine_like_substitutions)
        << "\nleft limit: "
        << "\n\t" << lpos
        << "\nright limit: "
        << "\n\t" << rpos - 1 // closed to open
        << "\nminimum alternative allele count acceptable: "
        << "\n\t" << minimum_alternative_allele_count_acceptable
        << "\nmaximum alternative allele count acceptable: "
        << "\n\t" << maximum_alternative_allele_count_acceptable
        << "\nminimum variation length acceptable: "
        << "\n\t" << minimum_variation_length_acceptable
        << "\nmaximum variation length acceptable: "
        << "\n\t" << maximum_variation_length_acceptable
        << "\nvariation type acceptable: "
        << "\n\t"; print_variation_type(std::cerr)
        << "\noverwrite the file if existing file path provided: "
        << "\n\t" << yes_or_no(force)
        << "\noutput sub-block: "
        << "\n\t" << (sub_block ? sub_block_outfile_path + (compress_bgz ? ".gz" : "") : "no")
        << '\n'
    ;
}

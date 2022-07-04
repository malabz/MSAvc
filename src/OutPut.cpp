#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cstring>
#include <sstream>

#include "Arguments.hpp"
#include "OutPut.hpp"

static void output_dot_if_empty(std::ostream &os, char const *str)
{
    // assert *str != '\0'

    if (*str == '^')
        ++str;

    // if (*str == '\0')
    //     os << '.';
    // else
        os << str;
}

using mutation_and_occurrences_iterator = mut::MutationContainer::const_iterator;

static void output_alterations(std::ostream &os, mutation_and_occurrences_iterator first, mutation_and_occurrences_iterator last)
{
    for (bool head = true; first != last; ++first)
    {
        if (first->second.size() < arguments::minimum_alternative_allele_count_acceptable)
            continue;

        if (head == false) os << ',';

        output_dot_if_empty(os, first->first.counterpart_segment.data());
        head = false;
    }
}

static void output_info(std::ostream &os, mutation_and_occurrences_iterator first, mutation_and_occurrences_iterator last, unsigned variation_length)
{
    os << "\tAC=";

    bool head = true;
    for (auto i = first; i != last; ++i)
    {
        if (i->second.size() < arguments::minimum_alternative_allele_count_acceptable)
            continue;

        if (head == false) os << ',';

        os << i->second.size();
        head = false;
    }

    os << ";VT=" << mut::abbreviated_mutation_types[first->first.variation_type] << ";VLEN=" << variation_length;
}

static unsigned calculate_variation_length(mut::Mutation const &mutation)
{
    // unsigned ref_length, alt_length;

    switch (mutation.variation_type)
    {
    case mut::SUB:
        return mutation.last - mutation.first;

    case mut::DEL:
        if (mutation.front_anchored)
            return mutation.reference_segment.size() - 1;
        else
            return mutation.reference_segment.size() - 2;

    case mut::INS:
        if (mutation.front_anchored)
            return mutation.counterpart_segment.size() - 1;
        else
            return mutation.counterpart_segment.size() - 2;

    // the max length of ref and alt excluding the prefixing character
    case mut::REP:
        if (mutation.front_anchored)
            return std::max(mutation.reference_segment.size() - 1, mutation.counterpart_segment.size() - 1);
        else
            return std::max(mutation.reference_segment.size() - 2, mutation.counterpart_segment.size() - 2);

    default:
        throw std::logic_error("unexpected variation_type");
        break;
    }
}

static unsigned get_pos(mut::Mutation const &mutation) noexcept
{
    unsigned const pos = mutation.first;
    if (mutation.reference_segment.front() == '^')
		return pos + 1;
    return pos;
}

static bool within_the_same_line(mut::Mutation const &lhs, mut::Mutation const &rhs) noexcept
{
    return arguments::combine_like_substitutions == false
            && lhs.variation_type == mut::SUB
            && rhs.variation_type == mut::SUB
            && lhs.first == rhs.first
            && lhs.last == rhs.last
            ;
}

void output_file_head(std::ofstream &ofs, utils::MultipleAlignmentFormat const &infile)
{
    ofs <<   "##fileformat=VCFv4.1"
        << "\n##contig=<ID=" << arguments::reference_name << ",length=" << infile.lengths[arguments::reference_index]
        << "\n##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate allele count, for each ALT allele, in the same order as listed\">"
           "\n##INFO=<ID=VT,Number=1,Type=String,Description=\"Type of small variant\">"
           "\n##INFO=<ID=VLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
           "\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    ;
}

void output(utils::MultipleAlignmentFormat const &maf, mut::MutationContainer const &mutations)
{
    std::ofstream ofs(arguments::outfile_path);
    if (!ofs) { std::cerr << "cannot open " << arguments::outfile_path << '\n'; exit(0); }

    auto const &names = maf.names;
    unsigned const n = names.size();

    output_file_head(ofs, maf);
    ofs << "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    if (arguments::genotype_matrix)
    {
        ofs << '\t' << arguments::reference_name;
        for (unsigned i = 0; i != arguments::reference_index; ++i) ofs << '\t' << names[i];
        for (unsigned i = arguments::reference_index + 1; i != n; ++i) ofs << '\t' << names[i];
    }
    ofs << '\n';

    // used only for arguments::matrix
    std::vector<unsigned> alteration_numbers;
    static constexpr unsigned not_visited = std::numeric_limits<unsigned>::max();

    static constexpr auto acceptable = [](mut::MutationContainer::value_type const &i) -> bool
            { return i.second.size() >= arguments::minimum_alternative_allele_count_acceptable; };

    for (auto i = mutations.cbegin(), j = i; i != mutations.cend(); i = j)
    {
        for (++j; j != mutations.cend(); ++j)
            if (within_the_same_line(i->first, j->first) == false)
                break;

        auto const &mutation_delegate = i->first;

        if (arguments::variation_type_acceptable[mutation_delegate.variation_type] == false)
            continue;

        if (mutation_delegate.first < arguments::lpos && !(mutation_delegate.first == 0 && arguments::lpos == 1))
            continue;
        if (mutation_delegate.first >= arguments::rpos)
            continue;

        unsigned const variation_length = calculate_variation_length(mutation_delegate);
        if (variation_length < arguments::minimum_variation_length_acceptable)
            continue;

        if (std::count_if(i, j, acceptable) == 0)
            continue;

        ofs << arguments::reference_name;
        ofs << '\t' << get_pos(mutation_delegate);
        ofs << "\t."
               "\t"; output_dot_if_empty(ofs, mutation_delegate.reference_segment.data());
        ofs << '\t'; output_alterations(ofs, i, j);
        ofs << "\t."
               "\t.";
        output_info(ofs, i, j, variation_length);
        ofs << "\tGT";

        if (arguments::genotype_matrix)
        {
            alteration_numbers.clear();
            alteration_numbers.resize(n, not_visited);

            // assumes that where.record is the same in the mutation
            unsigned const which_record = i->second.front().record;
            utils::Record const &record = maf.records[which_record];

            unsigned cnt = 0;
            for (auto k = i; k != j; ++k)
            {
                if (k->second.size() < arguments::minimum_alternative_allele_count_acceptable)
                    continue;

                ++cnt;
                for (auto const where : k->second)
                    alteration_numbers[record.belongs[where.sequence]] = cnt;
            }

            unsigned const reference_index_in_this_record = record.where_is(arguments::reference_index);
            std::string const &reference = record.sequences[reference_index_in_this_record];
            unsigned const reference_offset = record.begins[reference_index_in_this_record];

            unsigned l, r;
            if (mutation_delegate.first == reference_offset)
            {
                l = record.map_from_source_site[1];
                // for (auto const [_, which_sequence] : i->second)
                // {
                //     unsigned candidate = 1;
                //     for (; reference[candidate] == record.sequences[which_sequence][candidate]; ++candidate) ;
                //     if (l > candidate) l = candidate;
                // }

                if (mutation_delegate.last == record.map_to_source_site.back() + reference_offset)
                    r = record.map_from_source_site[mutation_delegate.last - reference_offset];
                    // equals to:
                    // r = reference.size();
                else
                    r = record.map_from_source_site[mutation_delegate.last - reference_offset + 1];
            }
            else
            {
                l = record.map_from_source_site[mutation_delegate.first - reference_offset];
                r = record.map_from_source_site[mutation_delegate.last - reference_offset - 1] + 1;

                // for (auto const [_, which_sequence] : i->second)
                // {
                //     unsigned candidate = record.map_from_source_site[mutation_delegate.last - reference_offset];
                //     for (; reference[candidate - 1] == record.sequences[which_sequence][candidate - 1]; --candidate) ;
                //     if (r > candidate) r = candidate;
                // }
            }

            for (unsigned index_based_on_record = 0; index_based_on_record != record.sequences.size(); ++index_based_on_record)
            {
                unsigned const index_of_current_sequence = record.belongs[index_based_on_record];
                if (alteration_numbers[index_of_current_sequence] != not_visited)
                    continue;

                alteration_numbers[index_of_current_sequence] = 0;
                if (memcmp(record.sequences[index_based_on_record].data() + l, reference.data() + l, r - l))
                    alteration_numbers[index_of_current_sequence] = not_visited;
            }

            ofs << "\t0";
            for (unsigned k = 0; k != n; ++k)
                if (k != arguments::reference_index)
                {
                    ofs << '\t';
                    if (alteration_numbers[k] == not_visited)
                        ofs << '.';
                    else
                        ofs << alteration_numbers[k];
                }
        }

        ofs << '\n';
    }
}

// 0-based [begin, end)
void output_sub_block(utils::MultipleAlignmentFormat const &infile, unsigned begin, unsigned end)
{
    // assert begin < end
    --begin;
    --end;

    for (auto const &record : infile.records)
    {
        unsigned const reference_index_of_record = record.where_is(arguments::reference_index);
        if (reference_index_of_record == record.sequences.size())
            continue;

        std::string const &reference = record.sequences[reference_index_of_record];
        unsigned const reference_offset = record.begins[reference_index_of_record];
        if (!(begin >= reference_offset && end <= reference_offset + record.map_to_source_site.back()))
            continue;

        unsigned const l = record.map_from_source_site[begin - reference_offset] + 1;
        unsigned const r = record.map_from_source_site[end - reference_offset] + 1;

        utils::Fasta fasta;
        fasta.names = infile.names;
        fasta.sequences.reserve(record.sequences.size());
        for (auto const &sequence : record.sequences)
            fasta.sequences.push_back(sequence.substr(l, r - l));

        std::ofstream ofs(arguments::sub_block_outfile_path);
        if (!ofs) { std::cerr << "cannot open " << arguments::outfile_path << '\n'; return; }
        fasta.write_to(ofs);

        return;
    }

    std::cerr << "error found when output subblock\n";
}

/**
 * deprecated
 * 
 * static void output_reference(std::ostream &os, std::string::const_iterator first, std::string::const_iterator last)
 * {
 *     if (*first == '^')
 *         ++first;
 * 
 *     if (first == last)
 *         os << '.';
 *     else
 *         for (; first != last; ++first)
 *             if (*first != '-')
 *                 os << *first;
 * }
 * 
 * static unsigned calculate_alternative_allele_count(mutation_and_occurrences_iterator first, mutation_and_occurrences_iterator last) noexcept
 * {
 *     unsigned cnt = 0;
 *     for (auto i = first; i != last; ++i)
 *         cnt += i->second.size();
 *     return cnt;
 * }
 * 
 */

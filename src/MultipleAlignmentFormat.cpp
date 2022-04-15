#include <algorithm>
#include <unordered_set>
#include <numeric>
#include <array>

#include "MultipleAlignmentFormat.hpp"

static std::array<std::vector<unsigned>, 2> split(std::string const &str)
{
    std::array<std::vector<unsigned>, 2> avi;
    unsigned i = std::find_if(str.cbegin(), str.cend(), 
            [](unsigned char c) -> bool { return isspace(c) == 0; }) - str.cbegin();
    avi[0].push_back(i);

    unsigned const n = str.size();
    for (; i != n; )
    {
        while (i != n && isspace(str[i]) == 0) ++i;
        avi[1].push_back(i);
        while (i != n && isspace(str[i]) != 0) ++i;
        avi[0].push_back(i);
    }

    avi[0].pop_back();
    return avi;
}

static unsigned string_to_unsigned(std::string const &num) noexcept
{
    unsigned u;
    try
    {
        u = std::stoul(num);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        exit(0);
    }

    return u;
}

static char first_not_space(std::string const &str)
{
    auto const first = std::find_if(str.cbegin(), str.cend(),
            [](unsigned char c) -> bool { return isspace(c) == 0; });

    if (first == str.cend())
        return '\0';
    else
        return *first;
}

void utils::MultipleAlignmentFormat::read(std::istream &is)
{
    std::string line;
    for (; std::getline(is, line); )
    {
        if (first_not_space(line) != 'a')
            continue;

        Record &record = records.emplace_back();
        std::unordered_set<std::string> names_of_current_record;

        for (; std::getline(is, line); )
        {
            if (line.empty() || first_not_space(line) != 's')
                break; // to do: this line should not be dropped because it might be a leading line of a record

            auto const [begins, ends] = split(line);
            unsigned const n = begins.size();
            if (n < 7) format_error();

            unsigned token = n - 1;
            record.sequences.emplace_back(line.cbegin() + begins.back(), line.cbegin() + ends.back());

            // length of the sequence the segment belongs to
            token = n - 2;
            unsigned const length_of_parent_sequence = string_to_unsigned(line.substr(begins[token], ends[token] - begins[token]));

            // '-' / '+'
            token = n - 3;
            char const forward_flag = line[begins[token]];
            if (ends[token] != begins[token] + 1 || forward_flag != '+' && forward_flag != '-')
                format_error();
            record.forward.push_back(forward_flag == '+');

            // count of residues in this segment
            token = n - 4;
            unsigned const number_provided = string_to_unsigned(line.substr(begins[token], ends[token] - begins[token]));
            unsigned const count_of_residues = record.sequences.back().size() -
                    std::count(record.sequences.back().cbegin(), record.sequences.back().cend(), '-');
            if (number_provided != count_of_residues)
                format_error();

            // a 0-based index, which indexes the segment position in its parent sequence
            // if forward_flag == '-', the number indexes the segment position in the reversed sequence of its parent sequence
            token = n - 5;
            unsigned const position = string_to_unsigned(line.substr(begins[token], ends[token] - begins[token]));
            record.begins.push_back(position);

            // name of the sequence, which might be seperated by spaces
            token = n - 6;
            std::string name = line.substr(begins[1], ends[token] - begins[1]);

            if (names_of_current_record.contains(name))
                Fasta::duplicate_name();
            names_of_current_record.insert(name);

            // leading flag
            token = 0;
            if (begins[token] + 1 != ends[token])
                format_error();

            auto found = name_to_index.find(name);
            if (found == name_to_index.cend())
            {
                record.belongs.push_back(names.size());

                name_to_index.emplace(name, names.size());
                names.push_back(std::move(name));
                lengths.push_back(length_of_parent_sequence);
            }
            else
            {
                record.belongs.push_back(found->second);
                if (lengths[found->second] != length_of_parent_sequence)
                    format_error();
            }
        }
    }
}

void utils::MultipleAlignmentFormat::read(utils::Fasta &&fasta)
{
    unsigned const n = fasta.sequences.size();

    names = std::move(fasta.names);
    for (unsigned i = 0; i != n; ++i)
    {
        name_to_index[names[i]] = i;
        lengths.push_back(fasta.sequences[i].size());
    }

    Record &record = records.emplace_back();
    record.begins.resize(n, 0);
    record.belongs.resize(n);
    std::iota(record.belongs.begin(), record.belongs.end(), 0);
    record.sequences = std::move(fasta.sequences);
    record.forward.resize(n, true);
}

void utils::MultipleAlignmentFormat::format_error()
{
    std::cerr << "maf format error\n";
    exit(0);
}

unsigned utils::Record::where_is(unsigned index) const
{
    return std::find(belongs.cbegin(), belongs.cend(), index) - belongs.cbegin();
}

void utils::Record::build_map_if_necessary(unsigned reference_index)
{
    unsigned const reference_index_in_this_record = where_is(reference_index);
    if (reference_index_in_this_record == sequences.size()) // one record might not contain the reference sequence
        return;

    // map.empty() == true <=> the record does not contain reference

    std::string const &reference = sequences[reference_index_in_this_record];
    unsigned const col = reference.size();

    auto &map_to = map_to_source_site; map_to.reserve(col + 1);
    auto &map_from = map_from_source_site; map_from.reserve(col + 1);

    for (unsigned i = 0; i != col; ++i)
        if (reference[i] != '-') {
            map_to.push_back(map_from.size());
            map_from.push_back(i);
        } else {
            map_to.push_back(map_from.size());
        }

    map_to.push_back(map_from.size());
    map_from.push_back(col);
}

void utils::Record::reverse(std::vector<unsigned> const &lengths_of_parent_sequences, std::array<char, 128> const &map_to_complemented) noexcept
{
    for (auto &sequence : sequences)
    {
        std::reverse(sequence.begin(), sequence.end());
        for (char &site : sequence) site = map_to_complemented[site];
    }

    unsigned const n = sequences.size();
    for (unsigned i = 0; i != n; ++i)
        begins[i] = lengths_of_parent_sequences[belongs[i]] - begins[i]
                - sequences[i].size() + std::count(sequences[i].cbegin(), sequences[i].cend(), '-');

    forward.flip();
}

void utils::MultipleAlignmentFormat::reverse_record_if_necessary(unsigned reference_index) noexcept
{
    for (auto &record : records)
    {
        unsigned const reference_index_in_this_record = record.where_is(reference_index);
        if (reference_index_in_this_record == record.sequences.size()) continue;
        if (record.forward[reference_index_in_this_record]) continue;

        record.reverse(lengths, map_to_complemented);
    }
}

std::array<char, 128> utils::MultipleAlignmentFormat::map_to_complemented = utils::MultipleAlignmentFormat::construct_map_to_complemented();

std::array<char, 128> utils::MultipleAlignmentFormat::construct_map_to_complemented()
{
    std::array<char, 128> map;
    std::iota(map.begin(), map.end(), 0);

    map['a'] = 't';
    map['A'] = 'T';

    map['t'] = 'a';
    map['T'] = 'A';

    // rna sequences could only be in fasta format, which does not require any reversing or complementing
    // map['u'] = 'a';
    // map['U'] = 'A';

    map['g'] = 'c';
    map['G'] = 'C';

    map['c'] = 'g';
    map['C'] = 'G';

    return map;
}

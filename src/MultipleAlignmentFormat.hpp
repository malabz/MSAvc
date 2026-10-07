#pragma once

#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <array>
#include <string>
#include <vector>

#include "Fasta.hpp"

namespace utils
{

    struct Record
    {
        std::vector<std::string> sequences;
        std::vector<bool> forward;
        std::vector<std::uint64_t> begins;          

        std::vector<unsigned> belongs;              
        unsigned where_is(unsigned index) const;

        std::vector<std::uint64_t> map_to_source_site;    
        std::vector<std::uint64_t> map_from_source_site;  
        void build_map_if_necessary(unsigned reference_index);

        void reverse(std::vector<std::uint64_t> const &lengths_of_parent_sequences,
                     std::array<char, 128> const &map_to_complemented) noexcept;
    };

    struct MultipleAlignmentFormat
    {
        MultipleAlignmentFormat() = default;

        void read(std::istream &is, std::string &ref);
        void read(utils::Fasta &&fasta);

        static void format_error();

        std::vector<std::uint64_t> lengths;         
        std::vector<std::string> names;
        std::vector<bool> is_prefix;
        std::unordered_map<std::string, unsigned> name_to_index;

        std::vector<Record> records;

        static std::array<char, 128> map_to_complemented;

        void reverse_record_if_necessary(unsigned reference_index) noexcept;

        static std::array<char, 128> construct_map_to_complemented();
    };

}
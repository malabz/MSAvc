#pragma once

#include <map>

#include "Fasta.hpp"
#include "MultipleAlignmentFormat.hpp"

namespace mut
{

    enum flag_enum : unsigned
    {
        SUB = 0,
        DEL = 1,
        INS = 2,
        REP = 3
    };

    struct Mutation
    {
        // ******* important members *******

        unsigned variation_type;

        // zero-based indexes of the source of the reference sequence
        unsigned first, last;
        unsigned seq_id;

        std::string counterpart_segment;


        // ******* unimportant members *******

        // determined by first and last, logged only for convenience
        std::string reference_segment;

        // makes sense if variation_type is not SUB
        // true if there is an anchor before the mutation
        bool front_anchored;

        bool operator<(Mutation const &rhs) const noexcept;
    };

    struct WhereAbout
    {
        unsigned record;
        unsigned sequence;

        WhereAbout() = default;
        WhereAbout(unsigned which_record, unsigned which_sequence);
    };

    using MutationContainer = std::map<Mutation, std::vector<WhereAbout>>;

    // std::ostream &operator<<(std::ostream &os, Mutation const &m);

    void extract_mutation(MutationContainer &mutations, std::vector<std::string> const &sequences, unsigned reference_index, unsigned counterpart_index, unsigned &position);
    MutationContainer search_in(utils::MultipleAlignmentFormat const &infile, unsigned reference_index);
    MutationContainer search_in(std::vector<std::string> const &sequences, unsigned reference_index);
    MutationContainer search_in(utils::MultipleAlignmentFormat &infile, const std::string &reference_prefix);
    MutationContainer search_in(std::vector<std::string> const &sequences, std::vector<std::string> const &names, const std::string &reference_prefix);

    template<typename lhs_iter, typename rhs_iter>
    unsigned deduce_variation_type(lhs_iter lhs_first, lhs_iter lhs_last, rhs_iter rhs_first) noexcept;
    unsigned deduce_variation_type(char lhs, char rhs) noexcept;

    char constexpr abbreviated_mutation_types[4][4] = { "SUB", "DEL", "INS", "REP" };
    std::string const mutation_types[4] = { "substitution", "deletion", "insertion", "replace" };

}

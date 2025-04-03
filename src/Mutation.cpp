#include <iostream>

#include "Mutation.hpp"
#include "Arguments.hpp"

mut::MutationContainer mut::search_in(std::vector<std::string> const &sequences, unsigned reference_index)
{
    MutationContainer mutations;

    std::string const &reference = sequences[reference_index];
    unsigned const col = reference.size();
    unsigned const row = sequences.size();

    for (unsigned i = 0; i != row; ++i)
    {
        if (i == reference_index)
            continue;

        for (unsigned j = 0; j != col; )
            if (reference[j] != sequences[i][j])
                extract_mutation(mutations, sequences, reference_index, i, j);
            else
                ++j;
    }

    return mutations;
}

mut::MutationContainer mut::search_in(utils::MultipleAlignmentFormat const &infile, unsigned reference_index)
{
    MutationContainer mutations;

    for (unsigned i = 0; i != infile.records.size(); ++i)
    {
        auto const &record = infile.records[i];
        unsigned const row = record.sequences.size();

        unsigned const reference_index_in_this_record = record.where_is(reference_index);
        if (reference_index_in_this_record == row) // one record might not contain the reference sequence
            continue;

        // check this block is in [lpos, rpos)
        unsigned const offset = record.begins[reference_index_in_this_record];
        // assert the length of sequences is same
        /*
        if(offset >= arguments::rpos || offset + record.sequences[0].size() < arguments::lpos && 
           !(offset == 0 && arguments::lpos == 1))
            continue;
        */

        std::string const &reference = record.sequences[reference_index_in_this_record];
        unsigned const col = reference.size();

        auto const &map_to_source_site = record.map_to_source_site;

        MutationContainer const raw_mutations = search_in(record.sequences, reference_index_in_this_record);

        for (auto const &raw_mutation : raw_mutations)
        {
            auto mutation = raw_mutation.first;
            /*
            // check this mutations is in [lpos, rpos)
            if(map_to_source_site[mutation.first] + offset < arguments::lpos ||
               map_to_source_site[mutation.first] + offset >= arguments::rpos &&
               !(map_to_source_site[mutation.first] + offset == 0 && arguments::lpos == 1))
                continue;
            */
            // assert reference[mutation.first] != '-'
            mutation.first = map_to_source_site[mutation.first] + offset;

            for (; mutation.last != col && reference[mutation.last] == '-'; ++mutation.last)
                ;
            mutation.last = map_to_source_site[mutation.last] + offset;

            auto &positions = mutations[std::move(mutation)];
            for (auto const raw_where : raw_mutation.second)
                positions.emplace_back(i, raw_where.sequence);
        }
    }

    return mutations;
}

mut::MutationContainer mut::search_in(std::vector<std::string> const &sequences, std::vector<std::string> const &names, const std::string &reference_prefix)
{
    MutationContainer mutations;
    // may contain multiple references
    unsigned const row = sequences.size();
    std::vector<bool> is_ref(row, false);
    for (unsigned ref_id = 0; ref_id != row; ++ ref_id)
        if(names[ref_id].starts_with(reference_prefix))
            is_ref[ref_id] = true;

    for (unsigned ref_id = 0; ref_id != row; ++ ref_id)
    {
        if(is_ref[ref_id])
        {
            std::string const &reference = sequences[ref_id];
            unsigned const col = reference.size();
            for (unsigned i = 0; i != row; ++i)
            {
                if (is_ref[i]) continue;

                for (unsigned j = 0; j != col; )
                    if (reference[j] != sequences[i][j])
                        extract_mutation(mutations, sequences, ref_id, i, j);
                    else
                        ++j;
            }

        }
    }
    return mutations;
}

mut::MutationContainer mut::search_in(utils::MultipleAlignmentFormat &infile, const std::string &reference_prefix)
{
    MutationContainer mutations;

    for (unsigned i = 0; i != infile.records.size(); ++i)
    {
        auto &record = infile.records[i];
        unsigned const row = record.sequences.size();

        unsigned reference_index_in_this_record = row;
        bool once = false;
        for (unsigned j = 0; j != row; ++ j)
        {
            if (infile.is_prefix[record.belongs[j]])
            {
                reference_index_in_this_record = j;
                // assume we found once
                if (once) { std::cerr << "Warning: found two fitted prefix. Please check the block " << i << "in maf. Program will ignore the reference\n"; }
                once = true;
            }
        }
        if (reference_index_in_this_record == row) // one record might not contain the reference sequence
            continue;
        // must need re-bulit this map
        {
            std::string const &reference = record.sequences[reference_index_in_this_record];
            unsigned const col = reference.size();

            auto &map_to = record.map_to_source_site; map_to.reserve(col + 1);
            auto &map_from = record.map_from_source_site; map_from.reserve(col + 1);

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

        // check this block is in [lpos, rpos)
        unsigned const offset = record.begins[reference_index_in_this_record];
        // assert the length of sequences is same
        /*
        if(offset >= arguments::rpos || offset + record.sequences[0].size() < arguments::lpos && 
           !(offset == 0 && arguments::lpos == 1))
            continue;
        */

        std::string const &reference = record.sequences[reference_index_in_this_record];
        unsigned const col = reference.size();

        auto const &map_to_source_site = record.map_to_source_site;

        MutationContainer const raw_mutations = search_in(record.sequences, reference_index_in_this_record);

        for (auto const &raw_mutation : raw_mutations)
        {
            auto mutation = raw_mutation.first;
            /*
            // check this mutations is in [lpos, rpos)
            if(map_to_source_site[mutation.first] + offset < arguments::lpos ||
               map_to_source_site[mutation.first] + offset >= arguments::rpos &&
               !(map_to_source_site[mutation.first] + offset == 0 && arguments::lpos == 1))
                continue;
            */
            // assert reference[mutation.first] != '-'
            mutation.first = map_to_source_site[mutation.first] + offset;

            for (; mutation.last != col && reference[mutation.last] == '-'; ++mutation.last)
                ;
            mutation.last = map_to_source_site[mutation.last] + offset;

            auto &positions = mutations[std::move(mutation)];
            for (auto const raw_where : raw_mutation.second)
                positions.emplace_back(i, raw_where.sequence);
        }
    }

    return mutations;
}

unsigned mut::deduce_variation_type(char lhs, char rhs) noexcept
{
    // assert lhs != rhs

    if (lhs == '-')
        return INS;
    else if (rhs == '-')
        return DEL;
    else
        return SUB;
}

template<typename lhs_iter, typename rhs_iter>
unsigned mut::deduce_variation_type(lhs_iter lhs_first, lhs_iter lhs_last, rhs_iter rhs_first) noexcept
{
    // assert *lhs_first != *rhs_first
    unsigned const flag = deduce_variation_type(*lhs_first++, *rhs_first++);

    for (; lhs_first != lhs_last; ++lhs_first, ++rhs_first)
        // if *lhs_first == *rhs_first
        //     assert *lhs_first == '-'
        if (*lhs_first != *rhs_first && flag != deduce_variation_type(*lhs_first, *rhs_first))
            return REP;

    return flag;
}

void mut::extract_mutation(mut::MutationContainer &mutations, std::vector<std::string> const &sequences, unsigned reference_index, unsigned which_sequence, unsigned &position)
{
    // assert postition < until
    // assert reference_index < sequences.size()
    // assert which_sequence < sequences.size()

    std::string const &lhs = sequences[reference_index];
    std::string const &rhs = sequences[which_sequence];
    // assert lhs.size() == rhs.size()
    // assert lhs[position] != rhs[position]

    Mutation mutation;
    mutation.first = position;

    // str[str.size()] returns charT() since c++11
    // position will end no after than lhs.size()
    for (++position; lhs[position] != rhs[position] || lhs[position] == '-'; ++position)
        ;

    mutation.last = position;

    // position is anchored here, i.e., the next site pending check
    if (position != lhs.size())
        ++position;

    // step back while the end is double '-'
    // mutation.last will never reach 0 in this loop because of the caret prefix
    while (lhs[mutation.last - 1] == '-' && rhs[mutation.last - 1] == '-')
        --mutation.last;

    mutation.variation_type = deduce_variation_type(lhs.cbegin() + mutation.first, lhs.cbegin() + mutation.last, rhs.cbegin() + mutation.first);

    // step back first if the variation type is not substitution
    // while making sure that reference[mutation.first] is not '-'
    if (mutation.variation_type != SUB)
        // assert lhs[mutation.first - 1] == rhs[mutation.first - 1]
        for (--mutation.first; lhs[mutation.first] == '-'; --mutation.first)
            ;
    // assert lhs[mutation.first] == rhs[mutation.first]

    for (unsigned i = mutation.first; i != mutation.last; ++i)
    {
        if (lhs[i] != '-') mutation.reference_segment.push_back(lhs[i]);
        if (rhs[i] != '-') mutation.counterpart_segment.push_back(rhs[i]);
    }

    // if (mutation.reference_segment.size() == 1 && mutation.reference_segment[0] == '^'
    //     || mutation.counterpart_segment.size() == 1 && mutation.counterpart_segment[0] == '^')
    if (mutation.reference_segment.front() == '^')
    {
        unsigned i = mutation.last;
        while (lhs[i] == '-') ++i;
        if (lhs[i] == '\0')
        {
            std::cerr << "cannot cope with the condition where no aligned nucleotides are the same. Skip this mutation\n";
            // exit(0);
        }
        else
        {
            mutation.reference_segment.push_back(lhs[i]);
            mutation.counterpart_segment.push_back(lhs[i]);
            mutation.front_anchored = false;
        }
    }
    else
    {
        mutation.front_anchored = true;
    }
    mutation.seq_id = reference_index;


    // .record will be valued in the calling function
    mutations[std::move(mutation)].emplace_back().sequence = which_sequence;
}

bool mut::Mutation::operator<(const Mutation &rhs) const noexcept
{
    if (first != rhs.first)
        return first < rhs.first;

    if (variation_type != rhs.variation_type)
        return variation_type < rhs.variation_type;

    if (last != rhs.last)
        return last < rhs.last;

    return counterpart_segment < rhs.counterpart_segment;
}

static void helper_intersection_of(unsigned ll, unsigned lr, unsigned rl, unsigned rr, unsigned &l, unsigned &r) noexcept
{
    // assert ll <= rl
    if (rl < lr) { l = rl; r = std::min(lr, rr); }
    else { l = r = 0; }
}

static void intersection_of(unsigned ll, unsigned lr, unsigned rl, unsigned rr, unsigned &l, unsigned &r) noexcept
{
    if (ll > rl) helper_intersection_of(rl, rr, ll, lr, l, r);
    else helper_intersection_of(ll, lr, rl, rr, l, r);
}

mut::WhereAbout::WhereAbout(unsigned which_record, unsigned which_sequence)
    : record(which_record)
    , sequence(which_sequence)
{
}

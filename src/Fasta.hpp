#pragma once

#include <string>
#include <vector>
#include <iostream>

namespace utils
{

    class Fasta
    {
    private:
        void _read(std::istream &is);

    public:
        static constexpr unsigned max_line_length = 80;

        std::vector<std::string> sequences;
        std::vector<std::string> names;

        explicit Fasta(std::istream &is);
        Fasta() = default;

        void write_to(std::ostream &os, bool output_name = true) const;

        static void write_with_wrapping(std::ostream &os, const std::string &sequence);

        template<typename InputIterator>
        static void write_to(std::ostream &os, InputIterator sequence_first, InputIterator sequence_last)
        {
            if (sequence_first == sequence_last) return;

            using difference_type = decltype(std::distance(sequence_first, sequence_last));
            const difference_type len = std::distance(sequence_first, sequence_last);

            for (difference_type i = 0; i != len; ++sequence_first, ++i)
            {
                os << *sequence_first;
                if (i != len - 1) os << '\n';
            }
        }

        template<typename InputIterator1, typename InputIterator2>
        static void write_to(std::ostream &os, InputIterator1 sequence_first, InputIterator1 sequence_last,
                InputIterator2 name_first)
        {
            if (sequence_first == sequence_last) return;

            using difference_type = decltype(std::distance(sequence_first, sequence_last));
            const difference_type len = std::distance(sequence_first, sequence_last);

            for (difference_type i = 0; i != len; ++sequence_first, ++name_first, ++i)
            {
                os << '>' << *name_first << '\n';
                write_with_wrapping(os, *sequence_first);
                if (i != len - 1) os << '\n';
            }
        }

        static void duplicate_name();

    };

}

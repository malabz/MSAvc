#include <cstring>
#include <unordered_set>

#include "Fasta.hpp"

utils::Fasta::Fasta(std::istream &is)
{
    _read(is);
}

void utils::Fasta::write_to(std::ostream &os, bool output_name) const
{
    if (output_name)
        write_to(os, sequences.cbegin(), sequences.cend(), names.cbegin());
    else
        write_to(os, sequences.cbegin(), sequences.cend());
}

void utils::Fasta::_read(std::istream &is)
{
    std::string each_line;
    std::string each_sequence;
    std::unordered_set<std::string> names_of_current_file;

    for (bool flag = false; std::getline(is, each_line); )
    {
        if (each_line.size() == 0)
            continue;

        if (each_line[0] == '>')
        {
            std::string name = each_line.substr(1);
            if (names_of_current_file.contains(name))
                duplicate_name();
            names_of_current_file.insert(name);

            names.push_back(std::move(name));
            if (flag)
                sequences.push_back(std::move(each_sequence));
            flag = true;
        }
        else if (flag)
        {
            each_sequence += each_line;
        }
    }
    sequences.push_back(each_sequence);
}

void utils::Fasta::write_with_wrapping(std::ostream &os, const std::string &sequence)
{
    const unsigned sequence_length = sequence.size();

    char *cut_sequence = new char[sequence_length + sequence_length / max_line_length + 1];
    unsigned des_index = 0;
    for (unsigned src_index = 0; src_index < sequence_length; src_index += max_line_length)
    {
        if (src_index) cut_sequence[des_index++] = '\n';

        unsigned write_length = sequence_length - src_index;
        if (write_length > max_line_length) write_length = max_line_length;

        memcpy(cut_sequence + des_index, sequence.data() + src_index, write_length);
        des_index += write_length;
    }
    cut_sequence[des_index] = 0;

    os << cut_sequence;
    delete[] cut_sequence;
}

void utils::Fasta::duplicate_name()
{
    std::cerr << "duplicate sequence name found\n";
    exit(0);
}

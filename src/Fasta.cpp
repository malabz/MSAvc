#include <cstring>
#include <unordered_set>

#include "Fasta.hpp"
#include "Arguments.hpp"

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

void utils::Fasta::write_to(std::stringstream &ss, bool output_name) const
{
    if (output_name)
        write_to(ss, sequences.cbegin(), sequences.cend(), names.cbegin());
    else
        write_to(ss, sequences.cbegin(), sequences.cend());
}

void utils::Fasta::_read(std::istream &is)
{
    std::unordered_set<std::string> names_of_current_file;
    std::string tmp;

    char * buffer = new char [arguments::buffer_size + 1];
    int blocks = 0;
    bool is_name = false, flag = false;
    // read data as so many blocks
    while(true)
    {
        is.read(buffer, arguments::buffer_size);
        for(int bi = 0; bi < is.gcount(); ++ bi)
        {
            if(buffer[bi] == '\r') continue;
            if(buffer[bi] == '\n')
            {
                if(! is_name) continue; // this \n is in the sequence
                else
                {
                    is_name = false;
                    if (arguments::check_duplicate)
                    {
                        // std::cerr << arguments::check_duplicate << std::endl;
                        if (names_of_current_file.contains(tmp))
                            duplicate_name();
                        names_of_current_file.insert(tmp);
                    }
                    names.emplace_back(std::move(tmp));
                }
            }
            else if(buffer[bi] == '>')
            {
                // here, this string will become sequence name
                is_name = true;
                if(flag) sequences.emplace_back(std::move(tmp));
                flag = true;
            }
            else tmp += buffer[bi];
        }
        if (! is) break;
    }
    sequences.emplace_back(std::move(tmp));
    delete[] buffer;
}



void utils::Fasta::write_with_wrapping(std::stringstream &ss, const std::string &sequence)
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

    ss << cut_sequence;
    delete[] cut_sequence;
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
    exit(1);
}

void utils::Fasta::pre_length_error()
{
    std::cerr << "sequence length is different\n";
    exit(1);
}

#include "Fasta.hpp"

#include <fstream>

Fasta::Fasta(const std::string& infile) :
    size(_read(infile))
{}

void Fasta::write_to(std::ostream& os, bool with_identification) const
{
    if (with_identification)
        write_to(os, sequences.cbegin(), sequences.cend(), identifications.cbegin());
    else
        write_to(os, sequences.cbegin(), sequences.cend());
}

size_t Fasta::_read(const std::string& infile)
{
    std::ifstream ifs(infile);

    if (!ifs) exit(0);
    std::cout << "reading " << infile << "..." << std::flush;

    std::string line, current = "^";

    while (std::getline(ifs, line))
        if (line.size() && line[0] == '>')
        {
            identifications.push_back(line.substr(1));
            break;
        }

    while (std::getline(ifs, line))
    {
        if (line.size() == 0) continue;

        if (line[0] == '>')
        {
            identifications.push_back(line.substr(1));
            sequences.push_back(std::move(current));
            current = "^";
        }
        else
        {
            current += line;
        }
    }
    sequences.push_back(current);

    std::cout << "finished, " << sequences.size() << " sequences found" << std::endl;
    return sequences.size();
}

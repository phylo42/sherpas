#include <iostream>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <io/fasta.h>

using namespace rappas::io;
using std::vector;
using std::string, std::string_view;
using std::cout, std::endl;
using std::move;


fasta::fasta(string&& header, string&& sequence) noexcept
    : _header(move(header)), _sequence(move(sequence))
{}

string_view fasta::header() const noexcept
{
    return _header;
}

string_view fasta::sequence() const noexcept
{
    return _sequence;
}

//------------------------------------------------------------------------------------
namespace bio = boost::iostreams;

vector<fasta> rappas::io::read_fasta(const string& file_name)
{
    cout << "Loading fasta: " + file_name << endl;
    bio::mapped_file_source mmap(file_name);
    bio::stream<bio::mapped_file_source> is(mmap, std::ios::in);

    vector<fasta> fasta_records;
    string line, header, sequence;
    while (std::getline(is, line))
    {
        if (boost::starts_with(line, ">"))
        {
            if (sequence.size())
            {
                fasta_records.emplace_back(move(header), move(sequence));
                sequence = "";
                sequence.reserve(1024);
            }

            header = { line.c_str() + 1, line.size() - 1 };
        }
        else if (!line.empty())
        {
            sequence.append(line);
        }
    }
    fasta_records.emplace_back(move(header), move(sequence));
    cout << "Loaded " << fasta_records.size() << " sequences.\n\n" << std::flush;
    return fasta_records;
}

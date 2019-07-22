#ifndef RAPPAS_CORE_FASTA_H
#define RAPPAS_CORE_FASTA_H

#include <string>
#include <vector>

namespace rappas
{
    namespace io
    {
        class fasta
        {
        public:
            fasta() = delete;
            fasta(std::string&& header, std::string&& sequence) noexcept;
            fasta(fasta&& other) noexcept = default;
            fasta(const fasta&) = delete;
            ~fasta() = default;

            std::string_view header() const noexcept;
            std::string_view sequence() const noexcept;
        private:
            std::string _header;
            std::string _sequence;
        };

        std::vector<fasta> read_fasta(const std::string& filename);
    }
}

#endif //RAPPAS_CORE_FASTA_H

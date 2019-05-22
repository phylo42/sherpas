#ifndef RAPPAS_CORE_KMER_ITERATOR_H
#define RAPPAS_CORE_KMER_ITERATOR_H

#include <string_view>
#include "phylo_kmer.h"

namespace core
{
    /// \brief An iterator class for k-mers of an input sequence.
    /// \details Iterates over all the k-mers of a given sequence and encodes them
    /// in a rolling fashion. Returns a lightweight std::string_view to a current k-mer
    /// of a given sequence and its code.
    class kmer_iterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;

        /// Value of a kmer_iterator, which is a pair of a string view to a k-mer and its code
        using value_type = std::pair<std::string_view, phylo_kmer::key_type>;

        kmer_iterator(std::string_view sequence_view, size_t kmer_size) noexcept;
        kmer_iterator(const kmer_iterator&) = delete;
        kmer_iterator(kmer_iterator&&) = delete;
        kmer_iterator& operator=(const kmer_iterator&) = delete;
        kmer_iterator& operator=(kmer_iterator&&) = delete;
        ~kmer_iterator() noexcept = default;

        kmer_iterator& operator++();

        bool operator==(const kmer_iterator& rhs) const noexcept;
        bool operator!=(const kmer_iterator& rhs) const noexcept;

        value_type operator*() const noexcept;
    private:
        std::pair<std::pair<phylo_kmer::key_type, size_t>, bool> encode_from(size_t position) const;
        std::pair<std::pair<phylo_kmer::key_type, size_t>, bool> next_kmer() const;

        /// String view of the original sequence
        std::string_view _sequence_view;
        size_t _kmer_position;
        size_t _kmer_size;
        phylo_kmer::key_type _current_key;
    };

    /// \brief A little proxy class that just creates a k-mer iterator over a sequence.
    /// Usage: for (const auto& [kmer, code] : to_kmers(sequence_view, kmer_size)) { ... }
    class to_kmers
    {
    public:
        using const_iterator = kmer_iterator;

        to_kmers(std::string_view sequence_view, size_t kmer_size) noexcept;

        const_iterator begin() const;
        const_iterator end() const;
    private:
        std::string_view _sequence_view;
        size_t _kmer_size;
    };
}

#endif //RAPPAS_CORE_KMER_ITERATOR_H

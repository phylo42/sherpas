#include <tuple>
#include "core/kmer_iterator.h"
#include "core/phylo_kmer.h"

using namespace core;

kmer_iterator::kmer_iterator(std::string_view sequence_view, size_t kmer_size) noexcept
    : _sequence_view{ sequence_view }
    , _kmer_position{ 0 }
    , _kmer_size{ kmer_size }
{
    /// _kmer_size == 0 means npos
    if (_kmer_size == 0)
    {
        _kmer_position = std::string_view::npos;
    }
    else
    {
        /// calculate the code of the first valid k-mer
        bool stop = false;
        while (!stop)
        {
            const auto [key_and_pos, valid] = encode_from(_kmer_position);
            std::tie(_current_key, _kmer_position) = key_and_pos;
            stop = valid || (_kmer_position == std::string_view::npos);
        }
    }
}

kmer_iterator& kmer_iterator::operator++()
{
    bool last_kmer_is_valid = true;
    bool stop = false;
    while (!stop)
    {
        /// if the last k-mer is valid, we can reuse its code in a rolling fashion
        if (last_kmer_is_valid)
        {
            const auto [key_and_pos, valid] = next_kmer();
            std::tie(_current_key, _kmer_position) = key_and_pos;
            stop = valid || (_kmer_position == std::string_view::npos);
            last_kmer_is_valid = valid;
        }
        /// otherwise we need to calculate the code from scratch
        else
        {
            const auto [key_and_pos, valid] = encode_from(_kmer_position);
            std::tie(_current_key, _kmer_position) = key_and_pos;
            stop = valid || (_kmer_position == std::string_view::npos);
        }
    }
    return *this;
}

bool kmer_iterator::operator==(const kmer_iterator& rhs) const noexcept
{
    return (_kmer_position == rhs._kmer_position) && (_sequence_view == rhs._sequence_view);
}

bool kmer_iterator::operator!=(const kmer_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

kmer_iterator::value_type kmer_iterator::operator*() const noexcept
{
    return { _sequence_view.substr(_kmer_position, _kmer_size), _current_key };
}

std::pair<std::pair<phylo_kmer::key_type, size_t>, bool> kmer_iterator::next_kmer() const
{
    auto key = _current_key;
    auto position = _kmer_position;
    bool success = false;

    if (position <= _sequence_view.size() - _kmer_size)
    {
        /// get the last base of the new k-mer
        ++position;
        const auto new_base = _sequence_view[position + _kmer_size - 1];

        /// check for non-valid symbols (gaps, '*', '.' etc.)
        if (auto new_code = core::encode(new_base); new_code)
        {
            /// shift bits and remove the first base of the current k-mer
            key <<= core::bit_length<core::seq_type>();
            key &= ((1u << (core::bit_length<core::seq_type>() * _kmer_size)) - 1);

            /// add the new base
            key |= *new_code;
            success = true;
        }
        /// if we found a gap or ambiguity, skip it and go to the next possible k-mer
        else
        {
            position = position + _kmer_size;
        }
    }
    /// if we reached the end of sequence
    else
    {
        position = std::string_view::npos;
    }
    return { { key, position }, success };
}

std::pair<std::pair<phylo_kmer::key_type, size_t>, bool> kmer_iterator::encode_from(size_t position) const
{
    phylo_kmer::key_type key = phylo_kmer::nan_key;
    bool success = false;

    /// string is too small to start from _kmer_position
    if (position > _sequence_view.size() - _kmer_size)
    {
        position = std::string_view::npos;
    }
    else
    {
        const auto kmer_view = _sequence_view.substr(_kmer_position, _kmer_size);

        /// correct k-mer with gaps
        if (auto code = core::encode_kmer(kmer_view); code)
        {
            key = *code;
            success = true;
        }
        /// skip a k-mer if we found a gap or ambiguity
        else
        {
            /// WARNING:
            /// this is not effective, because we can skip all the bases up to the one
            /// that is not valid, not just only one. But this is not clear how to return
            /// this position from core::encode_kmer, still having the API of the core clear and simple.
            /// TODO: find a solution
            position += 1;
        }
    }
    return { { key, position }, success };
}

to_kmers::to_kmers(std::string_view sequence_view, size_t kmer_size) noexcept
    : _sequence_view{ sequence_view }, _kmer_size{ kmer_size }
{
    /// If the input string is too small for a given k, we force begin() == end()
    /// to represent an empty set of k-mers.
    if (_sequence_view.size() < _kmer_size)
    {
        _kmer_size = 0;
    }
}

to_kmers::const_iterator to_kmers::begin() const
{
    return { _sequence_view, _kmer_size };
}

to_kmers::const_iterator to_kmers::end() const
{
    return { _sequence_view, 0 };
}
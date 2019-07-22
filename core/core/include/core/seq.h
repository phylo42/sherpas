#ifndef RAPPAS_CORE_SEQ_H
#define RAPPAS_CORE_SEQ_H

#include <cstdint>
#include <cstddef>
#include <optional>

namespace core
{
    /// \brief Sequence traits type.
    /// \details Describes all the properties based on the size of an alphabet, the alphabet itself etc.
    ///
    /// Must declare:
    ///     - size_t alphabet_size
    ///     - size_t max_kmer_length
    ///     - char_type decode(size_t)
    ///     - size_t encode(char_type)
    template<typename SeqType>
    struct seq_traits_impl;

#ifdef SEQ_TYPE_DNA

    /// \brief Auxiliary structure to define the seq_type.
    /// \sa seq_type
    struct dna
    {};

    /// \brief A compile-time constant for a current sequence type.
    /// \details Sequence type is determined at compile time for efficiency reasons.
    using seq_type = dna;

    template<>
    struct seq_traits_impl<dna>
    {
        /// \brief The type used to store one base of a sequence.
        using char_type = char;

        /// \brief The type used to store a k-mer value.
        /// \details K-mers are not stored as strings or char* values, but as values of this type instead.
        /// For example for DNA and k==3 the k-mer "AAA" == 000ul if kmer_t is unsigned long.
        using key_type = uint32_t;

        static constexpr std::optional<uint8_t> key_to_code(char_type base)
        {
            switch (base)
            {
                case 'A':
                    [[fallthrough]];
                case 'a':
                    return { 0 };
                case 'C':
                    [[fallthrough]];
                case 'c':
                    return { 1 };
                case 'G':
                    [[fallthrough]];
                case 'g':
                    return { 2 };
                case 'T':
                    [[fallthrough]];
                case 't':
                    [[fallthrough]];
                case 'U':
                    [[fallthrough]];
                case 'u':
                    return { 3 };
                case 'N':
                    [[fallthrough]];
                case 'n':
                    [[fallthrough]];
                case '-':
                    [[fallthrough]];
                case '.':
                    [[fallthrough]];
                default:
                    return std::nullopt;
            }
        };

        static constexpr char_type code_to_key[] = {'A', 'C', 'G', 'T'};
        static constexpr char_type ambiguous_chars[] =  {'N', '.', '-', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'};

        /// \brief Alphabet size
        /// \details The number of different *codes* of the alphabet. For DNA, 'T' and 'U' have the
        /// same code, which counts only once.
        static constexpr size_t alphabet_size = 4;
        static constexpr size_t max_kmer_length = 16;
    };

#elif SEQ_TYPE_AA
    static_assert(false, """SEQ_TYPE_AA is not supported yet. Supported types:\n"""
                         """SEQ_TYPE_DNA""");

    struct aa
    {};

    using seq_type = aa;

    template<>
    struct seq_traits_impl<aa>
    {
        using char_type = uint8_t;
        using key_type = uint64_t;

        static constexpr char_type char_set[] = { 'A' };

        static constexpr char_type decode(size_t /*code*/)
        {
            return 0;
        }

        static constexpr size_t encode(char_type /*base*/)
        {
            return 0;
        }

        static constexpr char_type ambiguous_chars[] =  { 'A' };

        static constexpr size_t alphabet_size = sizeof(char_set);
        static constexpr size_t max_kmer_length = 0;
    };
#else

    static_assert(false, """Please define a sequence type to compile core. Supported types:\n"""
                         """SEQ_TYPE_DNA""");
#endif

    using seq_traits = seq_traits_impl<seq_type>;

    /// \brief Returns amount of bits used to store one base of a sequence of given type.
    template<typename SeqType>
    constexpr seq_traits::key_type bit_length();

    /// \brief Returns a value of kmer_t type that can mask the rightmost base of a kmer_t.
    /// \details E.g. for DNA 0b1111...1100
    template<typename SeqType>
    constexpr seq_traits::key_type rightest_symbol_mask();


#ifdef SEQ_TYPE_DNA
    template<>
    constexpr seq_traits::key_type bit_length<dna>()
    {
        return seq_traits::key_type{ 2u };
    }

    template<>
    constexpr seq_traits::key_type rightest_symbol_mask<dna>()
    {
        return seq_traits::key_type{ ~0b11u };
    }

#elif SEQ_TYPE_AA
    /// ...
#else
    /// ...
#endif

    template<typename SeqTraits>
    constexpr typename SeqTraits::char_type decode_impl(uint8_t code)
    {
        return SeqTraits::code_to_key[code];
    }

    template<typename SeqTraits>
    constexpr std::optional<uint8_t> encode_impl(typename SeqTraits::char_type key)
    {
        return SeqTraits::key_to_code(key);
    }

    const auto decode = decode_impl<seq_traits>;
    const auto encode = encode_impl<seq_traits>;
}

#endif
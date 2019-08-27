#include "core/phylo_kmer.h"
#include "core/seq.h"
#include <cmath>
#include <vector>

using namespace core;

/// Compares two floats for almost equality.
/// From: https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type almost_equal(T x, T y, int ulp = 1) noexcept
{
    // The machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x - y) <= std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
           // unless the result is subnormal
           || std::abs(x - y) < std::numeric_limits<T>::min();
}

bool phylo_kmer::is_nan() const
{
    return (key == nan_key) && (score == nan_score);
}

bool core::operator==(const phylo_kmer& lhs, const phylo_kmer& rhs) noexcept
{
    if (lhs.is_nan() || rhs.is_nan())
    {
        return false;
    }
    else
    {
        return (lhs.key == rhs.key) && (almost_equal<phylo_kmer::score_type>(lhs.score, rhs.score));
    }
}

phylo_kmer::score_type core::score_threshold(phylo_kmer::score_type omega, size_t kmer_size)
{
    return std::log10(powf(omega / seq_traits::alphabet_size, phylo_kmer::score_type(kmer_size)));
}

std::optional<phylo_kmer::key_type> core::encode_kmer(std::string_view kmer)
{
    phylo_kmer::key_type key = 0;
    for (const auto base : kmer)
    {
        if (const auto& base_code = core::encode(base); base_code)
        {
            key <<= core::bit_length<core::seq_type>();
            key |= *base_code;
        }
        else
        {
            return std::nullopt;
        }
    }
    return key;
}

std::optional<phylo_kmer::key_type> core::encode_kmer(const std::string& kmer)
{
    return encode_kmer(std::string_view{ kmer });
}

std::string core::decode_kmer(phylo_kmer::key_type key, size_t kmer_size)
{
    std::vector<uint8_t> result;
    result.reserve(kmer_size);

    while (key > 0)
    {
        result.push_back(core::decode(key & ~core::rightest_symbol_mask<seq_type>()));
        key >>= bit_length<seq_type>();
    }

    for (size_t i = result.size(); i < kmer_size; ++i)
    {
        result.push_back(core::decode(0));
    }

    return { result.rbegin(), result.rend() };
}


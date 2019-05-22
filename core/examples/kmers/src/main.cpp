#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>
#include <core/seq.h>
#include <core/phylo_kmer.h>
#include <core/kmer_iterator.h>

/// Iterate over all the combinations with repetition
/// http://shoaib-ahmed.com/2018/for-each-combination-with-repetetion-c++/
template<typename V, typename Callable>
void for_each_combination(V &v, size_t gp_sz, Callable f) {
    V gp(gp_sz);
    auto total_n = std::pow(v.size(), gp.size());
    for (auto i = 0; i < total_n; ++i) {
        auto n = i;
        for (auto j = 0ul; j < gp.size(); ++j) {
            gp[gp.size() - j - 1] = v[n % v.size()];
            n /= v.size();
        }
        f(gp);
    }
}

void encode_string()
{
    auto alphabet = std::vector<char>{ 'A', 'C', 'G', 'T', '-' };
    const size_t kmer_size = 3;

    for_each_combination(alphabet, kmer_size,
                         [&](std::vector<char>& bases) {
                             const auto kmer = std::string{ bases.begin(), bases.end() };
                             if (const auto key = core::encode_kmer(kmer); key)
                             {
                                 std::cout << kmer << ": " << *key << std::endl;
                                 assert(kmer == core::decode_kmer(*key, kmer.size()));
                             }
                             else
                             {
                                 std::cout << kmer << " skipped" << std::endl;
                             }
                         });
}

void encode_string_views()
{
    std::cout << "\nIteration: " << std::endl;

    const auto long_read = std::string{ "--TTTAT-AAATGNNNN-CAAAN.NNTTTT---" };
    const size_t kmer_size = 4;

    /// Iterate over k-mers of a string_view (implicitly created from long_read).
    /// This is more effective than calling core::encode_kmer for every k-mer of a string,
    /// because core::to_kmers calculates codes in a rolling fashion, if possible.
    for (const auto& [kmer, code] : core::to_kmers(long_read, kmer_size))
    {
        std::cout << kmer << ": " << code << " " << std::endl;
        assert(core::encode_kmer(kmer) == code);
    }
}

int main()
{
    /// An example of core::encode_kmer and core::decode_kmer for std::string as input
    encode_string();

    /// An example of core::encode_kmer for std::string_view as input, and a iteration
    /// over an input seqeuence
    encode_string_views();
}
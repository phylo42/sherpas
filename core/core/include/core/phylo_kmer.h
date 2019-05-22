#ifndef RAPPAS_CORE_PHYLO_KMER_H
#define RAPPAS_CORE_PHYLO_KMER_H

#include "seq.h"
#include <limits>
#include <string>
#include <optional>

namespace core
{
    /// \brief A phylo k-mer structure.
    /// \details A key-value pair for a phylo-kmer, where key is a key_type value of a k-mer, and value is
    /// posterior probability score of this k-mer. Branch node id, position etc. omitted here, because these values
    /// are shared among multiple phylo k-mers and can be stored more effectively.
    struct phylo_kmer
    {
        /// The same type used to store a k-mer value
        /// (essentially we do not distinguish between "k-mer value" and "phylo-kmer value")
        using key_type = seq_traits::key_type;

        /// The type of a "posterior probability" score of a phylokmer
        using score_type = float;

        /// The type of a branch node id, that a phylokmer is mapped to.
        using branch_type = uint16_t;

        /// The type of a phylokmer's position in the alignment
        using pos_type = uint32_t;

        static constexpr key_type nan_key = std::numeric_limits<phylo_kmer::key_type>::max();
        static constexpr score_type nan_score = std::numeric_limits<phylo_kmer::score_type>::quiet_NaN();
        static constexpr branch_type nan_branch = std::numeric_limits<phylo_kmer::branch_type>::max();

        bool is_nan() const;

        key_type key;
        score_type score;
    };

    bool operator==(const phylo_kmer& lhs, const phylo_kmer& rhs) noexcept;

    /// Returns a phylo_kmer with special values, considered as NotAPhyloKmer. This phylokmer
    /// can not be equeal to any other phylo kmer (including itself)
    phylo_kmer make_napk();

    /// Returns a minumum score
    phylo_kmer::score_type score_threshold(size_t kmer_size);

    /// \brief Returns a code of input k-mer.
    /// \details Assumes that the size of input sequence equals k
    std::optional<phylo_kmer::key_type> encode_kmer(const std::string& kmer);
    std::optional<phylo_kmer::key_type> encode_kmer(std::string_view kmer);

    /// Creates a string of size kmer_size by given key
    std::string decode_kmer(phylo_kmer::key_type key, size_t kmer_size);
}

#endif
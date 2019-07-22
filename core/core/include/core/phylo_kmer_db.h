#ifndef RAPPAS_CORE_PHYLO_KMER_DB_H
#define RAPPAS_CORE_PHYLO_KMER_DB_H


#ifdef USE_SKA_FLAT_HASH_MAP
#include <flat_hash_map/flat_hash_map.hpp>
#elif USE_SKA_BYTELL_HASH_MAP
#include <flat_hash_map/bytell_hash_map.hpp>
#elif USE_ABSL_FLAT_HASH_MAP
#include <absl/container/flat_hash_map.h>
#elif USE_FOLLY_F14_FAST_MAP
#include <folly/container/F14Map.h>
#include <folly/hash/Hash.h>
#elif USE_PHMAP_FLAT_HASH_MAP
#include <parallel_hashmap/phmap.h>
#elif USE_TSL_ROBIN_MAP
#include <tsl/robin_map.h>
#elif USE_TSL_HOPSCOTCH_MAP
#include <tsl/hopscotch_map.h>
#endif

#include "phylo_kmer.h"

namespace core
{
    class phylo_kmer_db;
}

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void load(Archive& ar, ::core::phylo_kmer_db& db, const unsigned int /* file_version */);
    }
}

namespace core
{
    template<typename... Args>
#ifdef USE_SKA_FLAT_HASH_MAP
    using hash_map = ska::flat_hash_map<Args...>;
#elif USE_SKA_BYTELL_HASH_MAP
    using hash_map = ska::bytell_hash_map<Args...>;
#elif USE_ABSL_FLAT_HASH_MAP
    using hash_map = absl::flat_hash_map<Args...>;
#elif USE_FOLLY_F14_FAST_MAP
    using hash_map = folly::F14FastMap<Args...>;
#elif USE_PHMAP_FLAT_HASH_MAP
    using hash_map = phmap::flat_hash_map<Args...>;
#elif USE_TSL_ROBIN_MAP
    using hash_map = tsl::robin_map<Args...>;
#elif USE_TSL_HOPSCOTCH_MAP
    using hash_map = tsl::hopscotch_map<Args...>;
#endif
    /// requires manual intervention (does not support noexcept destructor and can not be defined this way)
    ///using hash_map = robin_hood::unordered_flat_map<Args...>;
    ///using hash_map = robin_hood::unordered_map<Args...>;

    namespace impl {
        class search_result;
    }

    /// \brief A value of phylo-kmer stored in phylo_kmer_db.
    /// \details One k-mer can correspond to multiple values, so the key of a k-mer is not stored here.
    struct pkdb_value
    {
        phylo_kmer::branch_type branch;
        phylo_kmer::score_type score;
        ///phylo_kmer::pos_type position;
    };

    /// \brief A Phylo-kmer database, that stores all the phylo-kmers.
    class phylo_kmer_db
    {
        /// We can get rid of this friend declaration if we provide a public set_kmer_size method.
        /// I thinks it is better to inject an invasive dependency here than to provide a public access to this
        /// variable.
        template<class Archive>friend void boost::serialization::load(Archive& ar,
            ::core::phylo_kmer_db& db, const unsigned int /* file_version */);
    public:

        /// Member types
        using key_type = phylo_kmer::key_type;
        using value_type = std::vector<pkdb_value>;

        /// \brief A storage of a phylo-kmer information.
        /// \details Note that phylo-kmers are not stored as objects of phylo_kmer,
        /// which is just a temporary storage for a phylo-kmer information.
        using storage = hash_map<key_type, value_type>;
        using const_iterator = storage::const_iterator;


        /// Ctors, dtor and operator=
        phylo_kmer_db(size_t kmer_size, const std::string& tree);
        phylo_kmer_db(const phylo_kmer_db&) noexcept = delete;
        phylo_kmer_db(phylo_kmer_db&&) = default;
        phylo_kmer_db& operator=(const phylo_kmer_db&) = delete;
        phylo_kmer_db& operator=(phylo_kmer_db&&) = default;
        ~phylo_kmer_db() noexcept = default;


        /// Access
        /// \brief Searches for a key against the database.
        /// \details WARNING: This method does not know how the key was calculated. It is required
        /// to provide keys of substrings of size _kmer_size to get correct results.
        /// \sa _kmer_size
        std::optional<impl::search_result> search(key_type key) const noexcept;


        /// Iterators
        /// \brief Returns an iterator to the beginning
        const_iterator begin() const noexcept;
        /// \brief Returns an iterator to the end
        const_iterator end() const noexcept;


        /// Capacity
        /// \brief Returns the number of keys
        size_t size() const noexcept;
        /// \brief Returns the k-mer size.
        size_t kmer_size() const noexcept;
        /// \brief Returns omega (score threshold parameter)
        core::phylo_kmer::score_type omega() const noexcept;
        /// \brief Returns a view to the newick formatted phylogenetic tree
        std::string_view tree() const noexcept;


        /// Modifiers
        /// \brief Puts a phylo-kmer information in the database.
        /// \details This method is unsafe, which means it does not control if the value has
        /// a correct branch id, the score is maximal etc. All of this checks must be done before calling
        /// this method. It just puts the value in a hash map.
        /// WARNING: This method does not know how the key was calculated. Here we assume it represents
        /// a string of size _kmer_size.
        /// \sa _kmer_size
        void insert(key_type key, const pkdb_value& value);

    private:
        storage _map;

        /// \brief K-mer size.
        /// \details This number is given by user to the constructor. We can not guarantee
        /// that the keys stored in hash tables actually correspond to substrings of the size _kmer_size.
        /// Example: DNA ('A', 'C', 'G', 'T")
        ///     key('AAA') == key('AA') == 0
        /// e.g. putting 0 in hashtable, we assume it corresponds to 'AAA' having _kmer_size == 3,
        /// but we can not guarantee that it was not calculated for another k-mer size by mistake.
        size_t _kmer_size;

        /// \brief Newick formatted phylogenetic tree
        std::string _tree;
    };

    namespace impl
    {
        /// \brief A search result wrapper around a collection of pairs [branch, score]
        /// to iterate over.
        class search_result
        {
        public:
            using const_iterator = phylo_kmer_db::value_type::const_iterator;

            search_result() noexcept;
            search_result(const_iterator begin, const_iterator end) noexcept;
            search_result(const search_result&) noexcept = default;
            ~search_result() noexcept = default;

            const_iterator begin() const noexcept;
            const_iterator end() const noexcept;

        private:
            const_iterator _begin;
            const_iterator _end;
        };
    }
}


#endif

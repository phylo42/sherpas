#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()

#include <boost/filesystem.hpp>
#include <catch2/catch.hpp>
#include <core/phylo_kmer_db.h>
#include <core/serialization.h>
#include <core/kmer_iterator.h>

namespace fs = boost::filesystem;

auto create_test_map()
{
    std::unordered_map<core::phylo_kmer::key_type,
        std::unordered_map<core::phylo_kmer::branch_type, core::phylo_kmer::score_type>> values =
        {
            {
                0, { { 0, 0.00f } }
            },
            {
                1, { { 0, 0.10f }, { 1, 0.11f } }
            },
            {
                2, { { 0, 0.20f }, { 1, 0.21f }, { 2, 0.22f } }
            },
            {
                3, { { 1, 0.31f }, { 2, 0.32f } }
            },
            {
                4, { { 2, 0.42f } }
            }
        };
    return values;
}

template<typename MapType>
core::phylo_kmer_db create_db_from_map(const MapType& values, size_t kmer_size)
{
    core::phylo_kmer_db db { kmer_size };
    for (const auto& [key, entries] : values)
    {
        for (const auto& [branch, score] : entries)
        {
            db.put(key, branch, score);
        }
    }
    return db;
}

TEST_CASE("Database size", "[database]")
{
    {
        const auto values = create_test_map();
        const auto db = create_db_from_map(values, 3);
        REQUIRE(db.size() == values.size());
    }

    {
        const core::phylo_kmer_db db { 3 };
        REQUIRE(db.size() == 0);
    }
}

TEST_CASE("K-mer size", "[database]")
{
    const size_t kmer_size = 5;
    const core::phylo_kmer_db db { kmer_size };

    REQUIRE(db.kmer_size() == kmer_size);
}

TEST_CASE("Duplicate key-value pair", "[database]")
{
    core::phylo_kmer_db db { 3 };

    db.put(0, 0, 1.0);
    REQUIRE(db.size() == 1);
    db.put(0, 0, 0.99);
    REQUIRE(db.size() == 1);
    db.put(0, 0, 1.01);
    REQUIRE(db.size() == 1);

    db.put(0, 1, 1.0);
    REQUIRE(db.size() == 1);
    db.put(0, 1, 1.001);
    REQUIRE(db.size() == 1);
    db.put(0, 1, 1.0);
    REQUIRE(db.size() == 1);
    db.put(0, 1, 0.99);
    REQUIRE(db.size() == 1);

    db.put(1, 0, 1.0);
    REQUIRE(db.size() == 2);
    db.put(1, 0, 0.99);
    REQUIRE(db.size() == 2);
    db.put(1, 0, 1.01);
    REQUIRE(db.size() == 2);
}

template<typename MapType>
void compare_db(const MapType& values, const core::phylo_kmer_db& db)
{
    for (const auto& [key, entries] : values)
    {
        auto db_entries = db.search(key);
        REQUIRE((bool)db_entries);
        for (const auto&[branch, score] : *db_entries)
        {
            REQUIRE(entries.find(branch) != entries.end());
            REQUIRE(entries.find(branch)->second == Approx(score));
        }
    }
}

TEST_CASE("Database search", "[database]")
{
    const auto filename = fs::unique_path().string();
    const auto values = create_test_map();
    const auto db = create_db_from_map(values, 3);
    compare_db(values, db);
}


TEST_CASE("(De-)serialization", "[database]")
{
    const auto filename = fs::unique_path().string();
    const auto values = create_test_map();
    const size_t kmer_size = 3;

    {
        const auto db = create_db_from_map(values, kmer_size);
        core::save(db, filename);
    }

    {
        const auto db = core::load(filename);
        REQUIRE(db.size() == values.size());
        REQUIRE(db.kmer_size() == kmer_size);
    }
}

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

TEST_CASE("Encoding and decoding k-mers", "[kmers]")
{
    auto alphabet = std::vector<char>{ 'A', 'C', 'G', 'T', '-', 'N' };
    const size_t kmer_size = 3;
    size_t count = 0;
    for_each_combination(alphabet, kmer_size,
                         [&count](std::vector<char>& bases) {
                             const auto kmer = std::string{ bases.begin(), bases.end() };
                             if (const auto key = core::encode_kmer(kmer); key)
                             {
                                 REQUIRE(kmer == core::decode_kmer(*key, kmer.size()));
                                 REQUIRE(*key == count);
                                 ++count;
                             }
                         });
}

TEST_CASE("core::to_kmers iteration", "[kmers]")
{
    /// Simple iteration
    const auto long_read = std::string{ "--TTTAT-AAATGNNNN-CAAAN.NNTTTT---" };
    const size_t kmer_size = 4;

    size_t count = 0;
    for (const auto& [kmer, code] : core::to_kmers(long_read, kmer_size))
    {
        REQUIRE(core::encode_kmer(kmer) == code);
        ++count;
    }
    REQUIRE(count == 6);
}

TEST_CASE("core::to_kmers empty set", "[kmers]")
{
    const size_t kmer_size = 3;

    for (const auto& [kmer, code] : core::to_kmers("", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("-", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("---", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("----", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("NNN", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("AAN", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("ANA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("ANA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("NAA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("-AA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("A-A", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : core::to_kmers("AA-", kmer_size))
    {
        REQUIRE_FALSE(true);
    }
}
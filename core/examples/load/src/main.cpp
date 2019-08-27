#include <iostream>
#include <boost/filesystem.hpp>
#include <core/phylo_kmer_db.h>
#include <core/serialization.h>
#include <iomanip>

core::phylo_kmer_db create_db()
{
    const size_t kmer_size = 3;
    const core::phylo_kmer::score_type omega = 1.0;
    const std::string tree;

    core::phylo_kmer_db db { kmer_size, omega, tree };

    /// branch 0
    db.insert(0, { 0, 0.00f });
    db.insert(1, { 0, 0.10f });
    db.insert(2, { 0, 0.20f });

    /// branch 1
    db.insert(1, { 1, 0.11f });
    db.insert(2, { 1, 0.21f });
    db.insert(3, { 1, 0.31f });

    /// branch 2
    db.insert(2, { 2, 0.22f });
    db.insert(3, { 2, 0.32f });
    db.insert(4, { 2, 0.42f });

    return db;
}

std::ostream& operator<<(std::ostream& out, const core::phylo_kmer_db& db)
{
    for (const auto& [key, entries] : db)
    {
        out << key << ":\n";
        for (const auto& [branch, score] : entries)
        {
            out << '\t' << branch << ": " << score << '\n';
        }
    }
    return out;
}

int main()
{
    const auto filename = boost::filesystem::unique_path().string();
    core::save(create_db(), filename);

    const auto db = core::load(filename);
    std::cout << "K-mer size: " << db.kmer_size() << std::endl;
    std::cout << "Omega: " << std::setprecision(2) << std::fixed << db.omega() << std::endl;
    std::cout << db;
}
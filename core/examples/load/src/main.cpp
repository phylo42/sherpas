#include <iostream>
#include <boost/filesystem.hpp>
#include <core/phylo_kmer_db.h>
#include <core/serialization.h>

core::phylo_kmer_db create_db()
{
    core::phylo_kmer_db db { 3 };

    /// branch 0
    db.put(0, 0, 0.00f);
    db.put(1, 0, 0.10f);
    db.put(2, 0, 0.20f);

    /// branch 1
    db.put(1, 1, 0.11f);
    db.put(2, 1, 0.21f);
    db.put(3, 1, 0.31f);

    /// branch 2
    db.put(2, 2, 0.22f);
    db.put(3, 2, 0.32f);
    db.put(4, 2, 0.42f);

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
    std::cout << db;
}
#include <iostream>
#include <core/phylo_kmer_db.h>

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

void search(const core::phylo_kmer_db& db, core::phylo_kmer_db::key_type key)
{
    if (auto entries = db.search(key); entries)
    {
        std::cout << "Found " << key << ":\n";
        for (const auto& [branch, score] : *entries)
        {
            std::cout << "\tbranch " << branch << ": " << score << '\n';
        }
    }
    else
    {
        std::cout << "Key " << key << " not found.\n";
    }
}

int main()
{
    const auto db = create_db();
    search(db, 0);
    search(db, 2);
    search(db, 42);
}
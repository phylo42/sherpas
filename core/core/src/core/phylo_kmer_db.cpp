#include "core/phylo_kmer_db.h"

using namespace core;

phylo_kmer_db::phylo_kmer_db(size_t kmer_size, const std::string& tree)
    : _kmer_size{ kmer_size }, _tree { tree }
{}

void phylo_kmer_db::insert(key_type key, const pkdb_value& value)
{
    _map[key].push_back(value);
}

std::optional<impl::search_result> phylo_kmer_db::search(key_type key) const noexcept
{
    if (auto it = _map.find(key); it != _map.end())
    {
        return impl::search_result{ it->second.begin(), it->second.end() };
    }
    else
    {
        return std::nullopt;
    }
}

phylo_kmer_db::const_iterator phylo_kmer_db::begin() const noexcept
{
    return std::begin(_map);
}

phylo_kmer_db::const_iterator phylo_kmer_db::end() const noexcept
{
    return std::end(_map);
}

size_t phylo_kmer_db::size() const noexcept
{
    return _map.size();
}

size_t phylo_kmer_db::kmer_size() const noexcept
{
    return _kmer_size;
}

core::phylo_kmer::score_type phylo_kmer_db::omega() const noexcept
{
    return 1.0f;
}

std::string_view phylo_kmer_db::tree() const noexcept
{
    return _tree;
}

impl::search_result::search_result() noexcept
{}

impl::search_result::search_result(
    impl::search_result::const_iterator begin,
    impl::search_result::const_iterator end) noexcept
    : _begin{ begin }, _end{ end }
{}

impl::search_result::const_iterator impl::search_result::begin() const noexcept
{
    return _begin;
}

impl::search_result::const_iterator impl::search_result::end() const noexcept
{
    return _end;
}
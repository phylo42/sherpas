#include "core/phylo_kmer_db.h"

using namespace core;

phylo_kmer_db::phylo_kmer_db(size_t kmer_size) noexcept
    : _kmer_size{ kmer_size }
{}

void phylo_kmer_db::put(key_type key, inner_key_type branch, value_type score)
{
    if (auto it = _map.find(key); it != _map.end())
    {
        if (auto inner_it = it->second.find(branch); inner_it != it->second.end())
        {
            if (inner_it->second < score)
            {
                _map[key][branch] = score;
            }
        }
        else
        {
            _map[key][branch] = score;
        }
    }
    else
    {
        _map[key][branch] = score;
    }
}

#include <iostream>

std::optional<impl::search_result> phylo_kmer_db::search(key_type key) const
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

phylo_kmer_db::const_iterator phylo_kmer_db::begin() const
{
    return std::begin(_map);
}

phylo_kmer_db::const_iterator phylo_kmer_db::end() const
{
    return std::end(_map);
}

size_t phylo_kmer_db::size() const
{
    return _map.size();
}

size_t phylo_kmer_db::kmer_size() const
{
    return _kmer_size;
}

impl::search_result::search_result() noexcept
{}

impl::search_result::search_result(
    impl::search_result::const_iterator begin,
    impl::search_result::const_iterator end) noexcept
    : _begin{ begin }, _end{ end }
{}

impl::search_result::const_iterator impl::search_result::begin() const
{
    return _begin;
}

impl::search_result::const_iterator impl::search_result::end() const
{
    return _end;
}
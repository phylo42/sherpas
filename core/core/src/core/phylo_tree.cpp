#include <core/phylo_tree.h>
#include <core/phylo_kmer.h>
#include <algorithm>

using namespace core;
using namespace core::impl;
using std::vector;
using std::string;
using std::move;
using std::begin, std::end;

phylo_node::phylo_node()
{
    _clean();
}

phylo_node::~phylo_node() noexcept
{
    for (auto child : _children)
    {
        delete child;
    }
}

bool phylo_node::operator==(const phylo_node& rhs) const noexcept
{
    return (_preorder_id == rhs._preorder_id) && (_label == rhs._label);
}

bool phylo_node::operator!=(const phylo_node& rhs) const noexcept
{
    return !operator==(rhs);
}

std::string phylo_node::get_label() const noexcept
{
    return _label;
}

phylo_node* phylo_node::get_parent() const noexcept
{
    return _parent;
}

phylo_node::id_type phylo_node::get_preorder_id() const noexcept
{
    return _preorder_id;
}

phylo_node::id_type phylo_node::get_postorder_id() const noexcept
{
    return _postorder_id;
}

phylo_node::branch_length_type phylo_node::get_branch_length() const noexcept
{
    return _branch_length;
}

std::vector<phylo_node*> phylo_node::get_children() const
{
    return _children;
}

void phylo_node::_clean()
{
    _preorder_id = -1;
    _postorder_id = -1;
    _label = "";
    _branch_length = 0.0;
    _children.clear();
    _parent = nullptr;
}

void phylo_node::_add_children(phylo_node* node)
{
    _children.push_back(node);
}

phylo_node* impl::get_leftmost_leaf(phylo_node* root) noexcept
{
    while (root && !root->get_children().empty())
    {
        root = root->get_children()[0];
    }
    return root;
}

postorder_tree_iterator::postorder_tree_iterator() noexcept
    : postorder_tree_iterator{ nullptr }
{}

postorder_tree_iterator::postorder_tree_iterator(phylo_node* node) noexcept
    : _current{ node }, _postorder_id{ 0 }
{}

postorder_tree_iterator& postorder_tree_iterator::operator=(const postorder_tree_iterator& rhs) noexcept
{
    if (*this != rhs)
    {
        _current = rhs._current;
    }
    return *this;
}

postorder_tree_iterator::operator pointer() const noexcept
{
    return _current;
}

bool postorder_tree_iterator::operator==(const postorder_tree_iterator& rhs) const noexcept
{
    return _current == rhs._current;
}

bool postorder_tree_iterator::operator!=(const postorder_tree_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

postorder_tree_iterator& postorder_tree_iterator::operator++()
{
    /// Go upside down if necessary. We need to know the index of current node in the parent->children
    phylo_node* temp = _current->get_parent();
    int idx = _id_in_parent(_current);
    while (idx == -1 && temp)
    {
        temp = _current->get_parent();
        idx = _id_in_parent(_current);
    }

    /// the end of the tree
    if (temp == nullptr)
    {
        _current = nullptr;
    }
        /// visit the next sibling
    else if ((size_t) idx + 1 < temp->get_children().size())
    {
        _current = temp->get_children()[idx + 1];
        _current = get_leftmost_leaf(_current);
    }
        /// visit the parent
    else
    {
        _current = temp;
    }
    return *this;
}

postorder_tree_iterator::reference postorder_tree_iterator::operator*() const
{
    return *_current;
}

postorder_tree_iterator::pointer postorder_tree_iterator::operator->() const noexcept
{
    return _current;
}

phylo_node::id_type postorder_tree_iterator::_id_in_parent(const phylo_node* node) const
{
    if (node->get_parent() != nullptr)
    {
        /// WARNING:
        /// Here we perform a linear search to look for an index
        /// in parent's children list. This is definitely not the best way
        /// to do this. It is okay for small trees though.
        const auto& children = node->get_parent()->get_children();
        const auto it = std::find(begin(children), end(children), node);
        if (it != end(children))
        {
            return distance(begin(children), it);
        }
    }

    /// TODO: reimplement this with std::optional
    return -1;
}

phylo_tree::phylo_tree(phylo_node* root, size_t node_count)
    : _root{ root } , _node_count{ node_count }
{
    auto it = postorder_tree_iterator{ core::impl::get_leftmost_leaf(root) };
    const auto end = postorder_tree_iterator{ nullptr };
    phylo_node::id_type postorder_id = 0;
    for (; it != end; ++it)
    {
        _preorder_id_node_mapping[it->get_preorder_id()] = it;

        it->_postorder_id = postorder_id;
        _postorder_id_node_mapping[postorder_id] = it;
        ++postorder_id;
    }
}

phylo_tree::~phylo_tree() noexcept
{
    delete _root;
}

phylo_tree::const_iterator core::phylo_tree::begin() const noexcept
{
    return postorder_tree_iterator{ core::impl::get_leftmost_leaf(_root) };
}

phylo_tree::const_iterator core::phylo_tree::end() const noexcept
{
    return postorder_tree_iterator(nullptr);
}

size_t phylo_tree::get_node_count() const noexcept
{
    return _node_count;
}

phylo_tree::value_pointer phylo_tree::get_root() const noexcept
{
    return _root;
}

std::optional<phylo_node*> phylo_tree::get_by_preorder_id(phylo_node::id_type preorder_id) const noexcept
{
    if (const auto it = _preorder_id_node_mapping.at(preorder_id); it)
    {
        return { it };
    }
    else
    {
        return { std::nullopt };
    }
}

std::optional<phylo_node*> phylo_tree::get_by_postorder_id(phylo_node::id_type preorder_id) const noexcept
{
    if (const auto it = _postorder_id_node_mapping.at(preorder_id); it)
    {
        return { it };
    }
    else
    {
        return { std::nullopt };
    }
}
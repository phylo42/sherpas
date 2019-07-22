#ifndef RAPPAS_CPP_PHYLO_TREE_H
#define RAPPAS_CPP_PHYLO_TREE_H

#include <string>
#include <vector>
#include <unordered_map>

namespace core
{
    class phylo_tree;
}

namespace rappas
{
    namespace io
    {
        class newick_parser;
        core::phylo_tree load_newick(const std::string& file_name);
    }
}

namespace core
{
    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A node of a phylogenetic tree.
    class phylo_node
    {
        friend rappas::io::newick_parser;
        friend core::phylo_tree;

    public:
        /// Member types

        /// \brief Node post-/pre-order id type
        using id_type = int;

        /// \brief Branch length type
        using branch_length_type = double;

        phylo_node();
        phylo_node(const phylo_node& other) = delete;
        phylo_node& operator=(const phylo_node&) = delete;
        ~phylo_node() noexcept;

        /// WARNING: this operator only checks for the id and label fields
        bool operator==(const phylo_node& rhs) const noexcept;
        bool operator!=(const phylo_node& rhs) const noexcept;

        std::string get_label() const noexcept;
        phylo_node* get_parent() const noexcept;
        id_type get_preorder_id() const noexcept;
        id_type get_postorder_id() const noexcept;
        branch_length_type get_branch_length() const noexcept;

        std::vector<phylo_node*> get_children() const;

    private:
        /// Clean node and fill with the default values. Used in the default constructor
        void _clean();

        void _add_children(phylo_node* node);

    private:
        id_type _preorder_id;
        id_type _postorder_id;
        std::string _label;
        branch_length_type _branch_length;

        std::vector<phylo_node*> _children;
        phylo_node* _parent;
    };

    namespace impl
    {
        /// \brief Finds the leftmost leaf of a subtree.
        /// \details Used to start a depth-first search
        phylo_node* get_leftmost_leaf(phylo_node* root) noexcept;

        //////////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief A forward access non-const iterator for phylo_node objects. Performs a postorder depth-first
        /// search among a subtree of an input phylo_node.
        class postorder_tree_iterator
        {
        public:
            using iterator_category = std::forward_iterator_tag;
            using reference = phylo_node&;
            using pointer = phylo_node*;

        public:
            postorder_tree_iterator() noexcept;
            explicit postorder_tree_iterator(phylo_node* node) noexcept;
            postorder_tree_iterator(const postorder_tree_iterator& other) = default;
            postorder_tree_iterator& operator=(const postorder_tree_iterator& rhs) noexcept;
            postorder_tree_iterator& operator=(postorder_tree_iterator&&) = delete;
            ~postorder_tree_iterator() noexcept = default;

            /// \brief Converts an iterator to the const phylo_node* it points to.
            operator pointer() const noexcept;

            bool operator==(const postorder_tree_iterator& rhs) const noexcept;
            bool operator!=(const postorder_tree_iterator& rhs) const noexcept;

            /// \brief The increment operator. Contains the logic of depth-first search
            postorder_tree_iterator& operator++();

            /// Access
            reference operator*() const;
            pointer operator->() const noexcept;

        private:
            /// \brief Finds an index of this node in the parent's array of children
            phylo_node::id_type _id_in_parent(const phylo_node* node) const;

        private:
            phylo_node* _current;
            phylo_node::id_type _postorder_id;
        };
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A phylogenetic tree class
    /// \defails phylo_tree is only constructable by the rappas::io::load_newick function.
    /// Non-copyable. Phylo-nodes are not modifiable.
    /// \sa core::phylo_node, rappas::io::load_newick
    class phylo_tree
    {
    public:
        /// Member types
        /// \brief A const iterator type. Performs a post-order depth-first search
        using const_iterator = impl::postorder_tree_iterator;
        using value_pointer = phylo_node*;

        /// Ctors, dtor and operator=
        phylo_tree(value_pointer root, size_t node_count);
        phylo_tree(phylo_tree&&) noexcept = default;
        phylo_tree(const phylo_tree&) = delete;
        phylo_tree& operator=(const phylo_tree&) = delete;
        phylo_tree& operator=(phylo_tree&&) = delete;
        ~phylo_tree() noexcept;


        /// Iterators
        /// \brief Returns an iterator to the beginning
        const_iterator begin() const noexcept;
        /// \brief Returns an iterator to the beginning
        const_iterator end() const noexcept;


        /// Access
        /// \brief Returns the number of nodes in a tree.
        size_t get_node_count() const noexcept;
        /// \brief Returns a pointer to the root node.
        value_pointer get_root() const noexcept;
        /// \brief Returns a pointer to a phylo_node with a given preorder_id, if presented.
        /// \details This operation does not require any traversal and implemented in O(1).
        /// \sa get_by_postorder_id
        std::optional<core::phylo_node*> get_by_preorder_id(phylo_node::id_type preorder_id) const noexcept;
        /// \brief Returns a pointer to a phylo_node with a given postorder_id, if presented.
        /// \details This operation does not require any traversal and implemented in O(1).
        /// \sa get_by_preorder_id
        std::optional<core::phylo_node*> get_by_postorder_id(phylo_node::id_type preorder_id) const noexcept;
    private:
        /// \brief A root node.
        value_pointer _root;

        /// \brief A total number of nodes in a tree.
        /// \details WARNING: This class does not double check if the tree actually have this
        /// number of nodes, it is just passed as an argument in the constructor.
        size_t _node_count;

        /// \brief A map for phylo_node_preorder_id-> phylo_node for fast access
        std::unordered_map<phylo_node::id_type, core::phylo_node*> _preorder_id_node_mapping;
        /// \brief A map for phylo_node_postorder_id-> phylo_node for fast access
        std::unordered_map<phylo_node::id_type, core::phylo_node*> _postorder_id_node_mapping;
    };
}

#endif
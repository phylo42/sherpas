#include <core/phylo_tree.h>
#include <core/phylo_kmer.h>
#include <io/file_io.h>
#include <iostream>
#include <stack>
#include <string_view>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/tokenizer.hpp>

using namespace core;
using std::vector, std::stack;
using std::string;
using std::move;
using std::cout, std::endl;
using std::begin, std::end;
using std::string_view;

phylo_node::phylo_node()
{
    _clean();
}

phylo_node::phylo_node(int id, const std::string& label, float branch_length,
	const std::vector<phylo_node*>& children, phylo_node* parent)
	: _id(id)
	, _label(label)
	, _branch_length(branch_length)
	, _children(children)
	, _parent(parent)
{
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
    return (_id == rhs._id) && (_label == rhs._label);
}

bool phylo_node::operator!=(const phylo_node& rhs) const noexcept
{
    return !operator==(rhs);
}

int phylo_node::get_id() const
{
    return _id;
}

std::string phylo_node::get_label() const
{
    return _label;
}

phylo_node* phylo_node::get_parent() const
{
    return _parent;
}

float phylo_node::get_branch_length() const
{
    return _branch_length;
}

std::vector<phylo_node*> phylo_node::get_children() const
{
    return _children;
}

void phylo_node::_clean()
{
	_id = -1;
	_label = "";
	_branch_length = 0;
	_children.clear();
	_parent = nullptr;
}

void phylo_node::_add_children(phylo_node* node)
{
    _children.push_back(node);
}

core::phylo_tree::phylo_tree(core::phylo_node* root, size_t node_count) noexcept
    : _root{ root }, _node_count{ node_count }
{}

core::phylo_tree::~phylo_tree() noexcept
{
    delete _root;
}

phylo_node* core::get_leftmost_child(phylo_node* root)
{
    while (!root->get_children().empty())
    {
        root = root->get_children()[0];
    }
    return root;
}

core::phylo_tree::const_iterator core::phylo_tree::begin() const
{
    return phylo_tree_iterator<true>{ core::get_leftmost_child(_root) };
}

core::phylo_tree::const_iterator core::phylo_tree::end() const
{
    return phylo_tree_iterator<true>(nullptr);
}

size_t core::phylo_tree::get_node_count() const
{
    return _node_count;
}

namespace core
{
///
/// \brief A class for parsing .newick-formatted files.
/// \details This class parses phylogenetic trees in the newick format. It designed to support
///  a buffered reading from disk.
///
    class newick_parser
    {
    public:
        newick_parser();
        newick_parser(const newick_parser&) = delete;
        newick_parser(newick_parser&&) = delete;
        ~newick_parser() = default;

        /// Parse an input buffer. This function can be called more than once,
        /// during the buffered reading from disk.
        /// \param data A string variable containing the current buffer data to parse.
        void parse(const string_view& data);

        phylo_node* get_root() const;
        size_t get_node_count() const;

    private:
        /// Parse next symbol of input data.
        /// \param ch A character to parse
        void _parse_character(char ch);

        /// \brief Handles a left parenthesis in input data.
        /// \details A left parenthesis indicates that a new node with children should be created.
        /// We will create it later though, during the _handle_right_parenthesis call
        /// because the phylo_node class has no default constructor for some design reasons.
        /// \sa phylo_node::phylo_node, _handle_right_parenthesis
        void _handle_left_parenthesis();

        /// \details The list of children for "current" parent node is over.
        /// The next symbols are referred to the parent node
        void _handle_right_parenthesis();

        /// \details Node delimiter, we create a node from the text content we collected so far
        void _handle_comma();

        /// \details End of file, we take last node as root
        void _handle_semicolon();

        /// \details Keep reading the current node description
        /// \param ch A character to parse
        void _handle_text(char ch);

        void _start_node();
        phylo_node* _finish_node();
        void _parse_node_text();

    private:
        stack<phylo_node*> _node_stack;
        phylo_node* _root;
        int _node_index;
        string _node_text;

        bool _parsing_node;
        bool _end_of_file;
    };
}

newick_parser::newick_parser()
    : _root(nullptr)
    , _node_index(-1)
    , _parsing_node(false)
    , _end_of_file(false)
{
}

void newick_parser::parse(const string_view& data)
{
    for (char c : data)
    {
        _parse_character(c);

        if (_end_of_file)
        {
            break;
        }
    }
}

phylo_node* newick_parser::get_root() const
{
    return _root;
}

size_t newick_parser::get_node_count() const
{
    return (size_t) _node_index + 1;
}

void newick_parser::_parse_character(char ch)
{
    switch (ch)
    {
        case '(':
            _handle_left_parenthesis();
            break;
        case ')':
            _handle_right_parenthesis();
            break;
        case ',':
            _handle_comma();
            break;
        case ';':
            _handle_semicolon();
            break;
        default:
            _handle_text(ch);
            break;
    }
}

void newick_parser::_handle_left_parenthesis()
{
    _start_node();
    _parsing_node = false;
}

void newick_parser::_handle_right_parenthesis()
{
    _finish_node();
    _parsing_node = true;
}

void newick_parser::_handle_comma()
{
    _finish_node();
}

void newick_parser::_handle_semicolon()
{
    _root = _finish_node();
    _end_of_file = true;
}

void newick_parser::_handle_text(char ch)
{
    /// A node can start from a parenthesis (ex. "(A, B)C") or from a label (ex. "A").
    /// The second case is equal to "()A". In this case we need to create a node as soon
    /// as we read the first symbol
    if (!_parsing_node)
    {
        _start_node();
        _parsing_node = true;
    }

    /// Keep reading the node description
    _node_text.push_back(ch);
}

void newick_parser::_start_node()
{
    ++_node_index;
    phylo_node* parent = _node_stack.empty() ? nullptr : _node_stack.top();
    _node_stack.push(new phylo_node());
    _node_stack.top()->_id = _node_index;
    _node_stack.top()->_parent = parent;
}

phylo_node* newick_parser::_finish_node()
{
    _parse_node_text();

    /// Add the node to its parent's list
    phylo_node* current_node = _node_stack.top();
    _node_stack.pop();
    if (current_node->_parent != nullptr)
    {
        current_node->_parent->_add_children(current_node);
    }

    _parsing_node = false;
    return current_node;
}

void newick_parser::_parse_node_text()
{
    phylo_node* current_node = _node_stack.top();

    // the content can be like "node_label:branch_length", ":branch_length", "node_label" or just ""
    if (!_node_text.empty())
    {
        using tokenizer = boost::tokenizer<boost::char_separator<char>>;
        tokenizer tokens(_node_text, boost::char_separator<char>(":"));
        auto it = begin(tokens);

        // if node label presented
        if (!boost::starts_with(_node_text, ":"))
        {
            current_node->_label = *(it++);
        }

        if (it != end(tokens))
        {
            current_node->_branch_length = std::stof(*it);
        }
    }

    // the current node is over
    _node_text.clear();
}

core::phylo_tree core::load_newick(const string& file_name)
{
    cout << "Loading newick: " + file_name << endl;

    newick_parser parser;
    rappas::io::buffered_reader reader(file_name);
    if (reader.good())
    {
        while (!reader.empty())
        {
            auto chunk = reader.read_next_chunk();
            parser.parse(chunk);
        }
    }
    else
    {
        throw std::runtime_error("Cannot open file: " + file_name);
    }

    cout << "Loaded a tree of " << parser.get_node_count() << " nodes." << endl << endl;
    return { parser.get_root(), parser.get_node_count() };
}

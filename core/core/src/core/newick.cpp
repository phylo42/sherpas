#include <iostream>
#include <sstream>
#include <stack>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/tokenizer.hpp>
#include <core/newick.h>
#include <core/phylo_tree.h>
#include <core/phylo_kmer.h>
#include <utils/io/file_io.h>
#include <iomanip>

using std::string, std::string_view;
using core::phylo_node;

namespace rappas
{
    namespace io
    {
        /// \brief A class for parsing .newick-formatted files.
        /// \details This class parses phylogenetic trees in the newick format. It designed to support
        ///  a buffered reading from disk.
        class newick_parser
        {
        public:
            newick_parser();
            newick_parser(const newick_parser&) = delete;
            newick_parser(newick_parser&&) = delete;
            ~newick_parser() = default;

            /// \brief Parses an input buffer. This function can be called more than once,
            /// during the buffered reading from disk.
            /// \param data A string variable containing the current buffer data to parse.
            void parse(std::string_view data);

            core::phylo_node* get_root() const;
            size_t get_node_count() const;

        private:
            /// \brief Parses next symbol of input data.
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
            core::phylo_node* _finish_node();
            void _parse_node_text();

            std::stack<core::phylo_node*> _node_stack;
            core::phylo_node* _root;
            int _node_index;
            std::string _node_text;

            bool _parsing_node;
            bool _end_of_file;
        };
    }
}

using rappas::io::newick_parser;

newick_parser::newick_parser()
    : _root(nullptr)
    , _node_index(-1)
    , _parsing_node(false)
    , _end_of_file(false)
{}

void newick_parser::parse(string_view data)
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
    _node_stack.top()->_preorder_id = _node_index;
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

core::phylo_tree rappas::io::load_newick(const string& file_name)
{
    std::cout << "Loading newick: " + file_name << std::endl;

    /// Load a tree from file
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

    /// Assign post-order ids to the phylo_node's
    auto tree = core::phylo_tree{ parser.get_root(), parser.get_node_count() };
    std::cout << "Loaded a tree of " << parser.get_node_count() << " nodes.\n\n" << std::flush;
    return tree;
}

core::phylo_tree rappas::io::parse_newick(std::string_view newick_string)
{
    newick_parser parser;
    parser.parse(newick_string);
    return core::phylo_tree{ parser.get_root(), parser.get_node_count() };
}

std::ostream& operator<<(std::ostream& out, const core::phylo_node& node)
{
    const auto num_children = node.get_children().size();
    if (num_children > 0)
    {
        out << "(";
        size_t i = 0;
        for (; i < num_children - 1; ++i)
        {
            out << *node.get_children()[i] << ",";
        }
        out << *node.get_children()[num_children - 1] << ")";
    }

    if (!node.get_label().empty())
    {
        out << node.get_label();
    }
    out << ":" << std::setprecision(10) << node.get_branch_length();
    out << "{" << node.get_postorder_id() << "}";
    return out;
}

std::string rappas::io::to_newick(const core::phylo_tree& tree)
{
    std::ostringstream stream;
    stream << *tree.get_root() << ";";
    return stream.str();
}
#include <core/phylo_tree.h>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

void write_test_tree(const std::string& filename)
{
    std::ofstream out(filename);
    out << "(A:0.1,B:0.2,((C:0.1,D:0.2)Y:0.1,(E:0.1,F:0.2)Z:0.2)X:0.2)W:0.0;";
}

std::ostream& operator<<(std::ostream& out, const core::phylo_tree& tree)
{
    std::cout << "Nodes: " << tree.get_node_count() << '\n';
    for (const auto& node : tree)
    {
        out << node.get_label() << ": " << node.get_branch_length() << '\n';
    }
    return out;
}

int main()
{
    const auto filename = boost::filesystem::unique_path().string();
    write_test_tree(filename);
    const auto tree = core::load_newick(filename);
    std::cout << tree << std::endl;
}
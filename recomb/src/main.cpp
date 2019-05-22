//============================================================================
// Name        : main.cpp
// Author      : gs
// Version     :
// Copyright   :
// Description : @_y
//============================================================================

#define SEQ_TYPE_DNA
#define USE_SKA_FLAT_HASH_MAP
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <time.h>
//#include "k-mer.h"
//#include "arcs.h"
#include "query.h"
#include "write.h"
//#include "heap.h"
//#include "read.h"
//#include <boost/filesystem.hpp>
//#include "rapidjson/reader.h"
#include <io/fasta.h>
#include <core/phylo_tree.h>
#include <core/serialization.h>

using namespace std;


int main(int argc, char** argv) {
	std::ios::sync_with_stdio(false);
    	if (argc != 4)
    	{
        	std::cout << "Usage:\n\t" << argv[0] << " DATABASE_FILE TREE_FILE QUERY_FILE "  << std::endl;
        	return 1;
    	}
	cout << "Let's go ! @_y" << endl;
	clock_t t=clock();
	const auto db = core::load(std::string{ argv[1]});
	const auto tree = core::load_newick(std::string{ argv[2]});
    	const auto sequences = rappas::io::read_fasta(std::string{ argv[3]});
	int k=db.kmer_size();
	int db_size=db.size();
	int s=sequences.size();
	int tree_size=tree.get_node_count();
	cout << "db and infos loaded in " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
	std::vector<core::phylo_kmer_db::key_type> codes(0);
	std::vector<Arc> branches=getArcs(tree_size);
	std::vector<Arc*> read(0);
	Htree H(read);
	std::vector<std::vector<Arc*>> windows(0);
	for(int i=0; i<s; i++)
	{
		cout << "query: " << sequences[i].header() << endl;
		t=clock();
		codes=encode_string_views(sequences[i].sequence(),k);
		readQuery(codes, db, &branches, &H);
		H.getTop(10);
		H.print(10);
		cout << "time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
		clearBranches(&branches);
		H.Hempty();
		cout << endl;
		t=clock();
		slidingWindow(codes, 500, 5, db, &branches, &windows);
		printChange(windows,1);
		cout << "time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
		SciPlot(i+1,sequences[i].sequence(), branches, windows);
		clearBranches(&branches);
		windInit(&windows);
		cout << "-------" << endl;
	}
	std::vector<int> groups(tree_size+1);
	for(int i=1; i<tree_size+1; i++)
	{
		groups[i]=(i-1)/9;
	}
	std::vector<double> bg(9, core::score_threshold(db.kmer_size()));
	core::phylo_kmer_db db2=GroupDb(db, groups, 9);
	search(db,0);
	search(db2,0);
	int lab=1;
	for (const auto& node : tree)
   	{
        	cout << to_string(lab) << "-" << node.get_label() << ": " << node.get_branch_length() << '\n';
		lab++;
    	}
	cout << "The End" << endl;
	return 0;
}

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
#include "groups.h"
//#include "heap.h"
//#include "read.h"
//#include <boost/filesystem.hpp>
//#include "rapidjson/reader.h"
#include <io/fasta.h>
//#include <core/phylo_tree.h>
#include <core/serialization.h>

using namespace std;


int main(int argc, char** argv) {
	std::ios::sync_with_stdio(false);
    	if (argc < 4)
    	{
        	std::cout << "Usage:\n\t" << argv[0] << " DATABASE_FILE TREE_FILE QUERY_FILE (GROUPS_FILE)"  << std::endl;
        	return 1;
    	}
	cout << "Let's go ! @_y" << endl;
	clock_t t=clock();
	core::phylo_kmer_db db = core::load(std::string{ argv[1]});
	core::phylo_tree tree = core::load_newick(std::string{ argv[2]});
    	const auto sequences = rappas::io::read_fasta(std::string{ argv[3]});
	int k=db.kmer_size();
	int db_size=db.size();
	int s=sequences.size();
	int tree_size=tree.get_node_count();
	std::vector<core::phylo_kmer_db::key_type> codes(0);
	std::vector<Arc> branches=getArcs(tree_size);
	std::vector<Arc*> read(0);
	Htree H(read);
	std::vector<std::vector<Arc*>> windows(0);
	std::vector<int> leafMap;
	std::vector<std::string> ref(0);
	//if(argc==5)
	//{
		recordGroups(argv[4], &leafMap, &ref);
		std::vector<int> groups=makeGroups(tree,leafMap);
		core::phylo_kmer_db db2=GroupDb(db, groups);
	//}
	cout << "db and infos loaded in " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
	for(int i=0; i<s; i++)
	{
		cout << "query: " << sequences[i].header() << endl;
		t=clock();
		codes=encode_string_views(sequences[i].sequence(),k);
		readQuery(codes, db2, &branches, &H);
		H.getTop(-1);
		H.print(-1, ref);
		cout << "placement time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
		clearBranches(&branches);
		H.Hempty();
		cout << endl;
		t=clock();
		slidingWindow(codes, 500, 5, db2, &branches, &windows);
		//printChange(windows,1, ref);
		SciPlot(i+1,sequences[i].header(), ref, branches, windows);
		cout << "windows time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
		clearBranches(&branches);
		windInit(&windows);
		cout << "-------" << endl;
	}
	//getId("/home/guillaume/Documents/LIRMM/Samples/HIV/ref/alignment-ref.fasta","/home/guillaume/Documents/LIRMM/Samples/HIV/ref/");
	cout << "The End" << endl;
	return 0;
}

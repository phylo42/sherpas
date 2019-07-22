//============================================================================
// Name        : main.cpp
// Author      : gs
// Version     :
// Copyright   :
// Description : @_y
//============================================================================


#include <iostream>
#include <fstream>
#include <unordered_map>
#include <time.h>
//#include "k-mer.h"
//#include "arcs.h"
#include "query.h"
#include "write.h"
#include "groups.h"
#include "input.h"
//#include "heap.h"
//#include "read.h"
//#include <boost/filesystem.hpp>
//#include "rapidjson/reader.h"
#include <utils/io/fasta.h>
//#include <core/phylo_tree.h>
#include <core/serialization.h>

using namespace std;


int main(int argc, char** argv) {
	std::ios::sync_with_stdio(false);
    	if (argc < 3)
    	{
        	std::cout << "Usage:\n\t" << argv[0] << " DATABASE_FILE TREE_FILE QUERY_FILE"  << std::endl;
        	return 1;
    	}
	cout << "Let's go ! @_y" << endl;
	clock_t t=clock();
	core::phylo_kmer_db db = core::load(std::string{ argv[1]});
	core::phylo_tree tree = rappas::io::load_newick(std::string{ argv[2]});
    	std::vector<rappas::io::fasta> sequences = rappas::io::read_fasta(std::string{ argv[3]});
	int k=db.kmer_size();
	int db_size=db.size();
	int s=sequences.size();
	int tree_size=tree.get_node_count();
	int ws=500;
	int shift=(ws+k-1)/2;
	std::vector<core::phylo_kmer_db::key_type> codes(0);
	std::vector<Arc> branches=getArcs(tree_size);
	std::vector<Arc*> read(0);
	Htree H(read);
	std::vector<std::vector<Arc*>> windows(0);
	std::vector<std::string> ref_db2(0);
	std::vector<std::string> ref_arc(0);
	std::vector<int> group_id(0);
	getArcRef(tree, &ref_arc);
	cout << "-------" << endl;
	getDb2Ref(tree, &ref_db2, &group_id);
	core::phylo_kmer_db db2=GroupDb(db, group_id);
	cout << "db and infos loaded in " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
	t=clock();
	for(int i=0; i<s; i++)
	{
		cout << "query: " << sequences[i].header() << endl;
		//t=clock();
		codes=encode_string_views(sequences[i].sequence(),k);
		/*readQuery(codes, db, &branches, &H);
		H.getTop(5);
		H.print(5, ref_arc);
		//cout << "placement time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
		clearBranches(&branches);
		H.Hempty();
		cout << endl;*/
		//t=clock();
		slidingWindow(codes, ws, 5, db, &branches, &windows);
		printChange(windows, shift, ref_arc);
		//cout << "windows time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
		//SciPlot(i+1,(sequences[i].header()), ref_arc, branches, windows, shift);
		Csv(i+1,(sequences[i].header()).data(), ref_arc, windows, shift);
		clearBranches(&branches);
		windInit(&windows);
		cout << "-------" << endl;
	}
	cout << "overall time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
	/*std::vector<rappas::io::fasta> test= rappas::io::read_fasta("/home/guillaume/Documents/LIRMM/Samples/HIV/consensus-gapfree.fasta");
	std::vector<rappas::io::fasta> nq;
	std::vector<std::string> summary (0);
	random_q(&test, &nq, 15, &summary, 1, 1, 2);
		writeFasta(&nq,"/home/guillaume/Documents/LIRMM/Samples/HIV/compendium/random-recomb.fasta");
	writeInfo(&summary, "/home/guillaume/Documents/LIRMM/Samples/HIV/compendium/infos-random.txt");*/
	cout << "The End" << endl;
	return 0;
}

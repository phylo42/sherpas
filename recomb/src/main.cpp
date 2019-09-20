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
#include <core/newick.h>

using namespace std;


int main(int argc, char** argv) {
	std::ios::sync_with_stdio(false);
    	if (argc != 3)
    	{
		cout << argc << endl;
        	std::cout << "Usage:\n\t" << argv[0] << " DATABASE_FILE ; QUERY_FILE"  << std::endl;
        	return 1;
    	}
	cout << "Let's go ! @_y" << endl;
	clock_t t=clock();
	core::phylo_kmer_db db_rap = core::load(std::string{ argv[1]});
	size_t k=db_rap.kmer_size();
	core::phylo_kmer_db db_max {k, db_rap.omega(), std::string{db_rap.tree()} };
	core::phylo_kmer_db db_small {k, db_rap.omega(), std::string{db_rap.tree()} };
	core::phylo_tree tree = rappas::io::parse_newick(db_rap.tree());
	std::vector<rappas::io::fasta> sequences = rappas::io::read_fasta(std::string{ argv[2]});
	int s=sequences.size();
	int tree_size=tree.get_node_count();
	int ws=500;
	int top=5;
	char dbtype='A';
	//cout << "Window size" << endl;
	//cin >> ws;
	//cout << "Scores to record" << endl;
	//cin >> top;
	int shift=(ws+k-1)/2;
	std::vector<core::phylo_kmer_db::key_type> codes(0);
	std::vector<Arc> branches=getArcs(tree_size);
	std::vector<Arc*> read(0);
	Htree H(read);
	std::vector<std::vector<Arc*>> windows(0);
	std::vector<std::string> ref(0);
	std::vector<int> group_id(0);
	core::phylo_kmer_db* db=NULL;
	if(dbtype=='A')
	{
		cout << "Using RAPPAS database" << endl;
		db=&db_rap;
		getArcRef(tree, &ref);
	}
	if(dbtype=='B')
	{
		cout << "Using max-based database" << endl;
		db=&db_max; 
		getDb2Ref(tree, &ref, &group_id);
		GroupDb(db_rap, db, group_id);
	}
	//onlyRoot(*db, &db_small, &ref);
	rmTop(*db, &db_small, &ref);
	db=&db_small;
	cout << "db and infos loaded in " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
	t=clock();
	for(int i=0; i<s; i++)
	{
		cout << ">" << sequences[i].header() << endl;
		//t=clock();
		codes=encode_string_views(sequences[i].sequence(),k);
		/*readQuery(codes, *db, &branches, &H);
		H.getTop(top);
		H.print(top, ref);
		//cout << "placement time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
		clearBranches(&branches);
		H.Hempty();
		cout << endl;*/
		//t=clock();
		//slidingWindow(codes, ws, top, *db, &branches, &windows);
		//printChange(windows, shift, ref);
		//cout << "windows time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
		//Csv(i+1,(sequences[i].header()).data(), ref, windows, shift);
		clearBranches(&branches);
		windInit(&windows);
		//cout << "-------" << endl;
	}
	cout << "overall time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
	/*std::vector<rappas::io::fasta> test= rappas::io::read_fasta("/home/guillaume/Documents/LIRMM/Samples/HIV/autumn-toyset/pol/HIV1_ALL_2017_pol_DNA.fasta");
	//std::vector<rappas::io::fasta> comp= rappas::io::read_fasta("/home/guillaume/Documents/LIRMM/Samples/HIV/compendium/compendium-norec-align.fasta");
	std::vector<rappas::io::fasta> nq;
	std::vector<std::string> summary (0);
	random_q(&test, &nq, 25, &summary, 2, 2);
	std::vector<rappas::io::fasta> nq_gf=gapRm(&nq);
	writeFasta(&nq_gf,"/home/guillaume/Documents/LIRMM/Samples/HIV/compendium/random-recomb.fasta");
	writeInfo(&summary, "/home/guillaume/Documents/LIRMM/Samples/HIV/compendium/infos-random.txt");*/
	//tikzDoc("/home/guillaume/Documents/LIRMM/Samples/HIV/autumn-toyset/pol/infos-random.txt", "/home/guillaume/Documents/LIRMM/Samples/HIV/autumn-toyset/pol/resA-k10-notop.txt", "/home/guillaume/Documents/LIRMM/Samples/HIV/autumn-toyset/pol/res-jpHMM.txt");
	//compRes("/home/guillaume/Documents/LIRMM/Samples/HIV/autumn-toyset/pol/infos-random.txt", "/home/guillaume/Documents/LIRMM/Samples/HIV/autumn-toyset/pol/res-jpHMM.txt");
	compMosaic("/home/guillaume/Documents/LIRMM/Samples/HIV/autumn-toyset/pol/infos-random.txt", "/home/guillaume/Documents/LIRMM/Samples/HIV/autumn-toyset/pol/res-scueal.txt");
	//codes=encode_string_views("AGGTGGGCCTTGA", 8);
	//cout << codes.size() << endl;
	cout << "The End" << endl;
	return 0;
}

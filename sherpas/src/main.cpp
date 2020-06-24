//============================================================================
// Name        : main.cpp
// Author      : gs
// Version     :
// Copyright   :
// Description : @_y
//============================================================================


#include <iostream>
#include <time.h>
#include "query.h"
#include "output.h"
#include <utils/io/fasta.h>
#include <xpas/serialization.h>
#include <xpas/newick.h>

using namespace std;

int main(int argc, char** argv) {
	std::ios::sync_with_stdio(false);
	int dflag=0;
	int qflag=0;
	int gflag=0;
	int oflag=0;
	int mflag=0;
	int tflag=0;
	int wflag=0;
	int cflag=0;
	int lflag=0;
	int kflag=0;
	int arg=0;
	std::string dbadd="";
	std::string qadd="";
	std::string gadd="";
	std::string oadd="";
	int ws=300;
	double theta=100;
	char dbtype='F';
	cout << "Let's go ! @_y" << endl;
	clock_t t=clock();
	while (-1 != (arg = getopt(argc, argv, "d:q:g:o:m:t:w:clk")))
	{
		switch(arg)
		{
			//mandatory parameters
			case 'd':
				// path to database
				dflag=1;
				dbadd=optarg;
				break;
			case 'q':
				// path to query file
				qadd=optarg;
				qflag=1;
				break;
			case 'g':
				// path to group-assignment file --TO DO--
				gadd=optarg;
				gflag=1;
				break;
			case 'o':
				// path to output directory
				oadd=optarg;
				oflag=1;
				break;
			//optional parameters (default values set)
			case 'm':
				// sherpas method (R/F)
				dbtype=optarg[0];
				if(dbtype=='R' && tflag==0)
				{
					theta=0.99;
				}
				mflag=1;
				break;
			case 't':
				// threshold for post-control
				theta=stod(optarg);
				tflag=1;
				break;
			case 'w':
				// window size
				ws=stoi(optarg);
				wflag=1;
				break;
			//activable extra options
			case 'c':
				// circular if activated
				cflag=1;
				break;
			case 'l':
				// single-line output if activated
				lflag=1;
				break;
			case 'k':
				// no N/A-removal if activated
				kflag=1;
				break;
			default:
        			abort ();
		}
	}
	if (dflag==0)
	{
    		cout << "please provide a path to database (-d)" << endl;
		return 1;
	}
	if (qflag==0)
	{
    		cout << "please provide a path to queries sequences (-q)" << endl;
		return 1;
	}
	if (gflag==0)
	{
    		cout << "please provide a path to a group-assignment file (-g)" << endl;
		return 1;
	}
	if (oflag==0)
	{
    		cout << "please provide an output directory (-o)" << endl;
		return 1;
	}
	if (dbtype != 'F' && dbtype != 'R')
	{
    		cout << "please choose a valid method (F or R)" << endl;
		return 1;
	}
	if (ws<100)
	{
    		cout << "please choose a window size greater than 100" << endl;
		return 1;
	}
	if (theta<0)
	{
    		cout << "please choose a nonnegative threshold" << endl;
		return 1;
	}
	if (theta>1 && dbtype=='R')
	{
    		cout << "please note that the threshold must be lower than 1 when using method R" << endl;
		return 1;
	}
	cout << "DB: " << dbadd << endl;
	cout << "queries: " << qadd << endl;
	cout << "groups: " << gadd << endl;
	cout << "output: " << oadd << endl;
	cout << "method: " << dbtype;
	if(mflag ==0)
	{
		cout << " (default)";
	}
	cout << endl;
	cout << "threshold: " << theta;
	if(tflag ==0)
	{
		cout << " (default)";
	}
	cout << endl;
	cout << "window: " << ws;
	if(wflag ==0)
	{
		cout << " (default)";
	}
	cout << endl;
    xpas::phylo_kmer_db db_rap = xpas::load(dbadd);

	size_t k=db_rap.kmer_size();
    xpas::phylo_tree tree = xpas::io::parse_newick(db_rap.tree());

	int tree_size=tree.get_node_count();
    xpas::phylo_kmer_db db_small {k, db_rap.omega(), std::string{db_rap.tree()} };
	int wi=100;
	int top=2;
	int shift=(wi+k-1)/2;
	std::string qfile=fileName(qadd);
	if(cflag ==1)
	{
		cout << "circular genome" << endl;
		make_circu(qadd, oadd+qfile+"-circ"+to_string(ws)+".fasta", (ws+k-1)/2);
		wi=ws;
		shift=0;
		qadd=oadd+qfile+"-circ"+to_string(ws)+".fasta";
	}

	std::vector<xpas::io::fasta> sequences = xpas::io::read_fasta(qadd);

	// CODE REVIEW: see query.h
	//sequences=gapRm(&sequences);

	int s=sequences.size();
	std::vector<std::vector<xpas::phylo_kmer_db::key_type>> codes(0);
	std::vector<Arc> branches=getArcs(tree_size);
	std::vector<Arc*> read(0);
	Htree H(read);
	std::vector<std::vector<Arc*>> windows(0);
	std::vector<std::string> ref(0);
	std::vector<double> rat(0);
	std::vector<int> group_id(0);
    xpas::phylo_kmer_db* db=&db_rap;
	getArcRef(tree, &ref, gadd);
	if(dbtype=='F')
	{
		cout <<"Using full database" << endl;
	}
	if(dbtype=='R')
	{
		cout <<"Using reduced database" << endl;
		onlyRoot(*db, &db_small, &ref);
		db=&db_small;
	}
	std::vector<std::string> all=listGroups(ref);
	all.push_back("all");
	cout << "db and infos loaded in " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
	t=clock();
	std::string outfile=oadd+"res-"+qfile+".txt";
	remove(&(outfile[0]));
	ofstream writef(outfile, ios::app);
	printHead(qfile, dbtype, theta, ws, cflag, kflag, &writef);
	for(int i=0; i<s; i++)
	{
		cout << s-i << " - " <<  sequences[i].header() << endl;
		writef << ">" << sequences[i].header() << endl;
		codes=encode_ambiguous_string(sequences[i].sequence(),k);
		slidingVarWindow(codes, wi, ws, top, *db, branches, windows, rat, dbtype);
		mergeNA(printChange(windows, shift, ref, theta, rat, dbtype), cflag, lflag, kflag, &writef);

		// CODE REVIEW: see the declaration of clearBranches
		clearBranches(branches);

		windInit(windows);
		rat.clear();
	}
	writef.close();
	cout << "overall time " << float(clock()-t)/CLOCKS_PER_SEC<< " sec." <<endl;
	cout << "The End !" << endl;
	return 0;
}

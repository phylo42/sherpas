/*
 * query.cpp
 *
 *  Created on: 23 oct. 2018
 *      Author: scholz
 */

#include<iostream>
#include<fstream>
#include<vector>
#include <cstring>
#include<cmath>
#include "query.h"

using namespace std;
//the main algorithm


// CODE REVIEW: See the header file
/*
std::vector<xpas::io::fasta> gapRm(std::vector<xpas::io::fasta> *sequences)
{
	//gap-remover; useful to remove gaps (now that was unexpected)
	int s=(*sequences).size();

    // CODE REVIEW: would be better to change the input vector instead.
    std::vector<xpas::io::fasta> res;

	for(int i=0; i<s; i++)
	{
	    // CODE REVIEW: Nice try. I think there was an easier way though.
		(res).emplace_back(move((((*sequences)[i]).header()).data()), move(xpas::io::clean_sequence((((*sequences)[i]).sequence()).data())));
	}
	return res;
}*/

// CODE REVIEW: see header file
void make_circu(const std::string& res, const std::string& rep, int p)
{
	// adds p bases of the end to the beginning of the sequence, and p bases of the beginning to the end. p should be (ws+k-1)/2 for queries.
	// creates a new file. To be used until circularity is considered in the core for db building and queries reading (if that happens).
	ofstream write(rep);
	std::vector<xpas::io::fasta> seq= xpas::io::read_fasta(res);
	std::string prefix = "";
	std::string suffix = "";
	for(int i=0; i<seq.size(); i++)
	{
		prefix = (seq[i].sequence()).substr(0, p);
		suffix = (seq[i].sequence()).substr((seq[i].sequence()).length()-p);
		write << ">" << seq[i].header() << endl;
		write << suffix << seq[i].sequence() << prefix << endl;
	}
}

void windInit(std::vector<std::vector<Arc*>> &wind)
{
    // CODE REVIEW:
    // If you just want to empty the vector, there is an easier way:
    //   wind.clear();

	// empties windows vector
	const size_t s = wind.size();
	for(size_t i=0; i<s; i++)
	{
	    // CODE REVIEW: Are you sure you don't want to delete the pointers
	    // that are inside? I don't quite understand where do they come from and
	    // where they should be deallocated
		wind.pop_back();
	}
}

// tests whether a character is a nucleotide (1) or a gap/uncertain state (0).
int nucl(char c)
{
    // CODE REVIEW: What is the alphabet and how to interpret characters has to be all in one place,
    // and at the moment it is in xpas::seq_traits. Consider using xpas::key_to_code instead
    if (xpas::seq_traits::key_to_code(c))
    {
        return 1;
    }
    return 0;
    /*
    int res=0;
    if(c=='A'||c=='C'||c=='T'||c=='G')
    {
        res=1;
    }
    return res;*/
}

// CODE REVIEW: instead of passing by value, consider passing by a reference to const
// CODE REVIEW: type alias
std::vector<sherpas::key_codes> fixCcodes(const std::vector<sherpas::key_codes>& ccodes, std::string_view long_read, size_t kmer_size)
{
    // assigns empty code string to kmers with more than one ambiguity
    std::vector<sherpas::key_codes> res(0);
    std::vector<xpas::phylo_kmer_db::key_type> nope (0);
    nope.push_back(-1);

    int pos=0;

    // CODE REVIEW: size_t
    size_t check=0;
    for(size_t i=0; i<kmer_size; i++)
    {
        // CODE REVIEW: see query.h
        check=check+1-nucl(long_read[i]);
    }

    if(check>1)
    {
        res.push_back(nope);
    }
    else
    {
        res.push_back(ccodes[pos]);
        pos++;
    }

    for(int i=kmer_size; i<long_read.size(); i++)
    {
        check=check+nucl(long_read[i-kmer_size])-nucl(long_read[i]);
        if(check>1)
        {
            res.push_back(nope);
        }
        else
        {
            res.push_back(ccodes[pos]);
            pos++;
        }
    }
    return res;
}

std::vector<sherpas::key_codes> encode_ambiguous_string(std::string_view long_read, const size_t kmer_size)
{
	//reads the query as a list of kmers, then translated into a list of multi-codes. kmers with one ambiguity get codes for each possibility, kmers with more than one ambiguity are skipped.
	size_t position = 0;
	std::vector<sherpas::key_codes> res(0);

	// CODE REVIEW: I can't see any usage of this. We push back a value, but the vector is not use
	std::vector<xpas::phylo_kmer_db::key_type> nope (0);
	nope.push_back(-1);


	for (const auto& [kmer, codes] : xpas::to_kmers<xpas::one_ambiguity_policy>(long_read, kmer_size))
	{
	    res.push_back(codes);
		++position;
	}

	// CODE REVIEW: Possible bug: long_read.size() is size_t, which is unsigned, so as kmer_size.
	// If long_read + 1 < kmer_size, the right hand side overflows.
	// Example: long_read = 8, kmer_size = 10. res.size() == 0, but
	// long_read.size()-kmer_size+1 is a crazily big number.
	// Then we take the branch despite that there are no kmers.
	if (res.size() != long_read.size()-kmer_size+1)
	{
		res=fixCcodes(res, long_read, kmer_size);
	}
	return res;
}

// CODE REVIEW: passing by reference
void readQuery(const std::vector<std::vector<xpas::phylo_kmer_db::key_type>>& codes,
               const xpas::phylo_kmer_db& db,
               std::vector<Arc>& branches,
               Htree& H)
{
	//RAPPAS algorithm on a given query
	//result stored in heap *H
	std::vector<Arc*> res(0);
	size_t k=db.kmer_size();
	double thr = log10(xpas::score_threshold(db.omega(),k));
	int q=codes.size();
	std::vector<xpas::phylo_kmer_db::key_type> ka;
	std::vector<double> Lamb(branches.size());
	int w=0;

	// CODE REVIEW: unused variable
	//double add=0;

	for (int i=0; i<q; i++)
   	{
	    ka=codes[i];
		if (ka.size() == 1)
		{
			if (auto entries = db.search(ka[0]); entries)
            {
                for (const auto& [branch, score] : *entries)
                {
                    branches[branch].updateScore(score,1);
                    if (branches[branch].getKmers() == 1)
                    {
                        res.push_back(&branches[branch]);
                    }
                }
            }
		}
		else
		{
			w=ka.size();
			Lamb.clear();
			for (int p=0; p<w; p++)
			{
				if (auto entries = db.search(ka[p]); entries)
				{
					for (const auto& [branch, score] : *entries)
					{
						Lamb[branch]+=(pow(10,score)-pow(10,thr))/w;
					}
				}
			}

			for (int j=0; j<Lamb.size(); j++)
			{
				if (Lamb[j] !=0)
				{
					branches[j].updateScore(Lamb[j],1);
					if (branches[j].getKmers() == 1)
					{
						res.push_back(&branches[j]);
					}
				}
			}
		}
	}
	int s=(res).size();
	double up;
	for(int i=0; i<s; i++)
	{
		up = thr * (q - res[i]->getKmers());
		res[i]->updateScore(up,0);
		H.push(res[i]);
	}
}

void addKmer(const std::vector<xpas::phylo_kmer_db::key_type>& keys, Htree& H,
             const xpas::phylo_kmer_db& db,
             std::vector<Arc>& branches, double thr, int sw, int k, int move)
{
    // CODE REVIEW:
    // There are benefits of passing by reference here.
    // instead of passing vectors by value:
    //   - save time not copying data
    // instead of passing by pointer:
    //   - we don't have to check if pointer is null
    //   - avoid parenthesis hell. Also things like "(*ptr)." can written like "ptr->".
    //
    //   Expressions like
    //      "(*((*H).h(i))).getKmers()>0)"
    //   are less clear than
    //      "H.h(i)->getKmers() > 0".
    //
    //   By using references we make syntax more clear.
    //   With some refactoring we can make it more readable:
    //
    //   Arc* current_arc = H.h(i);
    //   bool has_kmers = current_arc->getKmers() > 0;
    //   if (has_kmers) { ... }
    //
    //   Compare this three lines of code to "(*((*H).h(i))).getKmers()>0)"

	//updates scores when a kmer is added.
	int w=keys.size();
	double add=0;
	if(w==1)
	{
		if (auto entries = db.search(keys[0]); entries)
		{
			for (const auto& [branch, score] : *entries)
			{
				branches[branch].updateScore(score-thr,1);

				if((branches[branch]).getKmers() == 1)
				{
					branches[branch].updateScore(((sw-1)*thr),0);
					H.push(&(branches[branch]));
				}
				else
				{
					H.heapUp((branches[branch].ArcCheck()));
				}
			}
		}
	}
	else if(w>1)
	{
		std::vector<double> Lamb(branches.size());
		for(int i=0; i<w; i++)
		{
			if (auto entries = db.search(keys[i]) ; entries)
			{
				for (const auto& [branch, score] : *entries)
				{
					Lamb[branch]+=(pow(10,score)-pow(10,thr))/w;
				}
			}
		}
		for(int j=0; j<Lamb.size(); j++)
		{
			if(Lamb[j] !=0)
			{
				add=log10(Lamb[j]+pow(10,thr))-thr;
				branches[j].updateScore(add,1);
				if((branches[j]).getKmers()==1)
				{
					branches[j].updateScore(((sw-1)*thr),0);
					H.push(&(branches[j]));
				}
				else
				{
					H.heapUp(branches[j].ArcCheck());
				}
			}
		}
	}

	if (move==0)
	{
		for(int i=0; i < H.size(); i++)
		{
			if(H.h(i)->getKmers() > 0)
			{
				H.h(i)->updateScore(thr,0);
			}
		}
	}
}

void rmKmer(const std::vector<xpas::phylo_kmer_db::key_type>& keys,
            Htree& H, const xpas::phylo_kmer_db& db,
            std::vector<Arc>& branches, double thr, int move)
{
	//updates scores when a kmer is removed.
	int w=keys.size();
	double rem=0;

	if (w==1)
	{
		if (auto entries = db.search(keys[0]) ; entries)
		{
			for (const auto& [branch, score] : *entries)
			{
				branches[branch].updateScore(thr-score,-1);

				if(branches[branch].getKmers()==0)
				{
					H.pop(branches[branch].ArcCheck());
				}
				else
				{
					H.heapDown(branches[branch].ArcCheck());
				}
			}
		}
	}
	else if (w>1)
	{
		std::vector<double> Lamb(branches.size());
		for(int i=0; i<w; i++)
		{
			if (auto entries = db.search(keys[i]); entries)
			{
				for (const auto& [branch, score] : *entries)
				{
					Lamb[branch]+=(pow(10,score)-pow(10,thr))/w;
				}
			}
		}

		for(int j=0; j<Lamb.size(); j++)
		{
			if(Lamb[j] !=0)
			{
				rem=thr-log10(Lamb[j]+pow(10,thr));
				branches[j].updateScore(rem,-1);
				if(branches[j].getKmers() == 0)
				{
					H.pop(branches[j].ArcCheck());
				}
				else
				{
					H.heapDown(branches[j].ArcCheck());
				}
			}
		}
	}

	if(move==0)
	{
		for(int i=0; i<H.size(); i++)
		{
			if(H.h(i)->getKmers() > 0)
			{
				H.h(i)->updateScore(-thr,0);
			}
		}
	}
}

void slidingVarWindow(const std::vector<std::vector<xpas::phylo_kmer_db::key_type>>& codes,
                      int wi, int sw, int m,
                      const xpas::phylo_kmer_db& db,
                      std::vector<Arc>& branches,
                      std::vector<std::vector<Arc*>>& res,
                      std::vector<double>& rat, char met)
{
	//the main algorithm.
	//basically a sliding window except the window size increases (from wi to sw) at the beginning of the query and decreases (from sw to wi) at the end.
	//equivalent to a sliding window if wi=sw (case of circular queries).
	//as the window size increases by 2 at each step, wi and ws should have same parity (I can't guarantee what happens otherwise, but this may induce a shift in positions somewhere).
	Arc *neutral=new Arc(-1,0);
	std::vector<Arc*> rest(0);
	Htree H(rest);
	std::vector<Arc*> temp(m);
	size_t k=db.kmer_size();
	double thr = log10(xpas::score_threshold(db.omega(),k));
	int q=codes.size();
	std::vector<xpas::phylo_kmer_db::key_type> ka;
	std::vector<xpas::phylo_kmer_db::key_type> kr;
	if(sw>q || wi>sw)
	{
		cout << "Dimension problems: check windows size vs k-mer size and/or windows size vs query length" << endl;
		cout << wi << "<" << sw << "<" << q << endl;
	}
	else
	{
		std::vector<std::vector<xpas::phylo_kmer_db::key_type>> firstw(codes.begin(),codes.begin()+wi-1);
		readQuery(firstw, db, branches, H);
		H.getTop(m);

		if(met=='R')
        {
            rat.push_back(H.lRatio(k));
        }
        else
        {
            rat.push_back(H.dRatio(k));
        }

        // CODE REVIEW: this is not a pessimization to get red of this variable.
        // The compiler is smart enough to inline the call of H.size() everywhere.
		// int s=H.size();

		// CODE REVIEW: the code of this loop copied 4 times in this function.
		// Consider making it a function
		for(int i=0; i<m; i++)
		{
		    // CODE REVIEW: use H.size()
			if(H.size() > i)
			{
			    // CODE REVIEW: What is the point of copying Arc data here?
			    // We create a copy of the arc, and most importantly we **allocate one arc in heap memory**.
			    // This is a quite expensive operation. Depending on the hardware and your luck, that may be
			    // orders of magnitude slower than:
			    //
			    //   Arc best = *H.h(i);
			    //

			    // I would question that it is needed to use
			    //   std::vector<Arc*> temp;
			    // instead of
			    //   std::vector<Arc> temp;
			    // (the same applies to the returning type).
			    // Consider switching to std::vector<Arc>, that will handle new allocations more wisely.
			    // The heap memory for std::vector is not allocated element by element like here, but
			    // in big chunks, which would make the total number of memory allocations less.

				Arc *best = new Arc(*(H.h(i)));
				temp[i]=best;
			}
			else
			{
				temp[i]=neutral;
			}
		}

        // CODE REVIEW: consider reserving memory for res in the beginning of this function.
        // If I understand it correctly, the result size is about m + (sw - wi - 1) + (q - sw) + m + (q - wi - 1 - q + sw)
        // (the number of res.push_back calls)
        // https://en.cppreference.com/w/cpp/container/vector/reserve
		res.push_back(temp);

		for(int j=wi; j<sw-1; j+=2)
		{
			ka=codes[j];
			kr=codes[j+1];
			addKmer(ka, H, db, branches, thr, j+1, k, 0);
			addKmer(kr, H, db, branches, thr, j+1, k, 0);

			// CODE REVIEW: In all three loops for different parts of the window, the code from here
			// to the end of the loop is the same. Correct me if I am wrong. If not, consider making it
			// a function
			H.getTop(m);

			if (met=='R')
			{
				rat.push_back(H.lRatio(k));
			}
			else
			{
				rat.push_back(H.dRatio(k));
			}

			// CODE REVIEW: see above
			//s=H.size();
			for(int i=0; i<m; i++)
			{
			    // CODE REVIEW: see above
				if(H.size()>i)
				{
				    // CODE REVIEW: see above
					Arc *best=new Arc(*(H.h(i)));
					temp[i]=best;
				}
				else
				{
					temp[i]=neutral;
				}
			}
			res.push_back(temp);
		}

		for(int j=sw; j<q; j++)
		{
			ka=codes[j];
			kr=codes[j-sw];

			addKmer(ka, H, db, branches, thr, j+1, k, 1);
			rmKmer(kr, H, db, branches, thr, 1);

            // CODE REVIEW: see the first loop
			H.getTop(m);
			if(met=='R')
			{
				rat.push_back(H.lRatio(k));
			}
			else
			{
				rat.push_back(H.dRatio(k));
			}

            // CODE REVIEW: see above
			//s=H.size();
			for(int i=0; i<m; i++)
			{
                // CODE REVIEW: see above
				if (H.size()>i)
				{
                    // CODE REVIEW: see above
					Arc * best = new Arc(*(H.h(i)));
					temp[i]=best;
				}
				else
				{
					temp[i]=neutral;
				}
			}
			res.push_back(temp);
		}

		for(int j=q-sw; j<q-wi-1; j+=2)
		{
			ka=codes[j];
			kr=codes[j+1];

			rmKmer(ka, H, db, branches, thr, 0);
			rmKmer(kr, H, db, branches, thr, 0);

            // CODE REVIEW: see the first loop
			H.getTop(m);

			if(met=='R')
			{
				rat.push_back(H.lRatio(k));
			}
			else
			{
				rat.push_back(H.dRatio(k));
			}

            // CODE REVIEW: see above
			//s=H.size();
			for(int i=0; i<m; i++)
			{
                // CODE REVIEW: see above
				if(H.size() > i)
				{
                    // CODE REVIEW: see above
					Arc *best=new Arc(*(H.h(i)));
					temp[i]=best;
				}
				else
				{
					temp[i]=neutral;
				}
			}
			res.push_back(temp);
		}
	}
}

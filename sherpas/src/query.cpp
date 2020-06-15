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

std::vector<rappas::io::fasta> gapRm(std::vector<rappas::io::fasta> *sequences)
{
	//gap-remover; useful to remove gaps (now that was unexpected)
	int s=(*sequences).size();
	std::vector<rappas::io::fasta> res;
	for(int i=0; i<s; i++)
	{
		(res).emplace_back(move((((*sequences)[i]).header()).data()), move(rappas::io::clean_sequence((((*sequences)[i]).sequence()).data())));
	}
	return res;
}

void make_circu(std::string res, std::string rep, int p)
{
	// adds p bases of the end to the beginning of the sequence, and p bases of the beginning to the end. p should be (ws+k-1)/2 for queries.
	// creates a new file. To be used until circularity is considered in the core for db building and queries reading (if that happens).
	ofstream write(rep);
	std::vector<rappas::io::fasta> seq= rappas::io::read_fasta(res);
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

void windInit(std::vector<std::vector<Arc*>> *wind)
{
	//empties windows vector
	int s=(*wind).size();
	for(int i=0; i<s; i++)
	{
		(*wind).pop_back();
	}
}

int nucl(char c)
//tests whether a character is a nucleotide (1) or a gap/uncertain state (0).
{
    int res=0;
    if(c=='A'||c=='C'||c=='T'||c=='G')
    {
        res=1;
    }
    return res;
}

std::vector<std::vector<core::phylo_kmer_db::key_type>> encode_ambiguous_string(std::string_view long_read, const size_t kmer_size)
{
	//reads the query as a list of kmers, then translated into a list of multi-codes. kmers with one ambiguity get codes for each possibility, kmers with more than one ambiguity are skipped.
	size_t position = 0;
	std::vector<std::vector<core::phylo_kmer_db::key_type>> res(0);
	std::vector<core::phylo_kmer_db::key_type> nope (0);
	nope.push_back(-1);
	for (const auto& [kmer, codes] : core::to_kmers<core::one_ambiguity_policy>(long_read, kmer_size))
	{
		{
			res.push_back(codes);
		}
		++position;
	}
	if(res.size() != long_read.size()-kmer_size+1)
	{
		res=fixCcodes(res, long_read, kmer_size);
	}
	return res;
}

std::vector<std::vector<core::phylo_kmer_db::key_type>> fixCcodes(std::vector<std::vector<core::phylo_kmer_db::key_type>> ccodes, std::string_view long_read, const size_t kmer_size)
{
	// assigns empty code string to kmers with more than one ambiguity
	std::vector<std::vector<core::phylo_kmer_db::key_type>> res(0);
	std::vector<core::phylo_kmer_db::key_type> nope (0);
	nope.push_back(-1);
	int check=0;
	int pos=0;
	for(int i=0; i<kmer_size; i++)
	{
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

void readQuery(std::vector<std::vector<core::phylo_kmer_db::key_type>> codes, const core::phylo_kmer_db& db, std::vector<Arc>* branches, Htree *H)
{
	//RAPPAS algorithm on a given query
	//result stored in heap *H
	std::vector<Arc*> res(0);
	size_t k=db.kmer_size();
	double thr = log10(core::score_threshold(db.omega(),k));
	int q=codes.size();
	std::vector<core::phylo_kmer_db::key_type> ka;
	std::vector<double> Lamb((*branches).size());
	int w=0;
	double add=0;
	for(int i=0; i<q; i++)
   	{
        	ka=codes[i];
		if(ka.size() == 1)
		{
			if (auto entries = db.search(ka[0]) ; entries)
    			{
        			for (const auto& [branch, score] : *entries)
        			{
					(*branches)[branch].updateScore(score,1);
					if(((*branches)[branch]).getKmers()==1)
					{
						(res).push_back(&((*branches)[branch]));
					}
        			}
			}
    		}
		else
		{
			w=ka.size();
			Lamb.clear();
			for(int p=0; p<w; p++)
			{
				if (auto entries = db.search(ka[p]) ; entries)
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
					(*branches)[j].updateScore(Lamb[j],1);
					if(((*branches)[j]).getKmers()==1)
					{
						(res).push_back(&((*branches)[j]));
					}
				}
			}
    		}
	}
	int s=(res).size();
	double up;
	for(int i=0; i<s; i++)
	{
		up=(thr)*(q-(*(res)[i]).getKmers());
		(*(res)[i]).updateScore(up,0);
		(*H).push(res[i]);
	}
}

void addKmer(std::vector<core::phylo_kmer_db::key_type> keys, Htree* H, const core::phylo_kmer_db& db, std::vector<Arc>* branches, double thr, int sw, int k, int move)
{
	//updates scores when a kmer is added.
	int w=keys.size();
	double add=0;
	if(w==1)
	{
		if (auto entries = db.search(keys[0]) ; entries)
		{
			for (const auto& [branch, score] : *entries)
			{
				(*branches)[branch].updateScore(score-thr,1);
				if(((*branches)[branch]).getKmers()==1)
				{
					(*branches)[branch].updateScore(((sw-1)*thr),0);
					(*H).push(&((*branches)[branch]));
				}
				else
				{
					(*H).heapUp(((*branches)[branch].ArcCheck()));
				}
			}
		}
	}
	else if(w>1)
	{
		std::vector<double> Lamb((*branches).size());
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
				(*branches)[j].updateScore(add,1);
				if(((*branches)[j]).getKmers()==1)
				{
					(*branches)[j].updateScore(((sw-1)*thr),0);
					(*H).push(&((*branches)[j]));
				}
				else
				{
					(*H).heapUp(((*branches)[j].ArcCheck()));
				}
			}
		}
	}
	if(move==0)
	{
		for(int i=0; i<(*H).size(); i++)
		{
			if((*((*H).h(i))).getKmers()>0)
			{
				(*((*H).h(i))).updateScore(thr,0);
			}
		}
	}
}

void rmKmer(std::vector<core::phylo_kmer_db::key_type> keys, Htree* H, const core::phylo_kmer_db& db, std::vector<Arc>* branches, double thr, int move)
{
	//updates scores when a kmer is removed.
	int w=keys.size();
	double rem=0;
	if(w==1)
	{
		if (auto entries = db.search(keys[0]) ; entries)
		{
			for (const auto& [branch, score] : *entries)
			{
				(*branches)[branch].updateScore(thr-score,-1);
				if(((*branches)[branch]).getKmers()==0)
				{
					(*H).pop(((*branches)[branch]).ArcCheck());
				}
				else
				{
					(*H).heapDown(((*branches)[branch].ArcCheck()));
				}
			}
		}
	}
	else if(w>1)
	{
		std::vector<double> Lamb((*branches).size());
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
				rem=thr-log10(Lamb[j]+pow(10,thr));
				(*branches)[j].updateScore(rem,-1);
				if(((*branches)[j]).getKmers()==0)
				{
					(*H).pop(((*branches)[j]).ArcCheck());
				}
				else
				{
					(*H).heapDown(((*branches)[j].ArcCheck()));
				}
			}
		}
	}
	if(move==0)
	{
		for(int i=0; i<(*H).size(); i++)
		{
			if((*((*H).h(i))).getKmers()>0)
			{
				(*((*H).h(i))).updateScore(-thr,0);
			}
		}
	}
}


void slidingVarWindow(std::vector<std::vector<core::phylo_kmer_db::key_type>> codes, int wi, int sw, int m, const core::phylo_kmer_db& db, std::vector<Arc>* branches, std::vector<std::vector<Arc*>> *res, std::vector<double> *rat, char met)
{
	//the main algorithm.
	//basically a sliding window except the window size increases (from wi to sw) at the beginning of the query and decreases (from sw to wi) at the end.
	//equivalent to a sliding window if wi=sw (case of circular queries).
	//as the window size increases by 2 at each step, wi and ws should have same parity (I can't guarantee what happens otherwise, but this may induce a shift in positions somewhere).
	//Arc *neutral=new Arc(-1,0);
	std::vector<Arc*> rest(0);
	Htree H(rest);
	std::vector<Arc*> temp(m);
	size_t k=db.kmer_size();
	double thr = log10(core::score_threshold(db.omega(),k));
	int q=codes.size();
	std::vector<core::phylo_kmer_db::key_type> ka;
	std::vector<core::phylo_kmer_db::key_type> kr;
	if(sw>q || wi>sw)
	{
		cout << "Dimension problems: check windows size vs k-mer size and/or windows size vs query length" << endl;
		cout << wi << "<" << sw << "<" << q << endl;
	}
	else
	{
		std::vector<std::vector<core::phylo_kmer_db::key_type>> firstw(codes.begin(),codes.begin()+wi-1);
		readQuery(firstw, db, branches, &H);
		H.getTop(m);
		if(met=='R')
			{
				(*rat).push_back(H.lRatio(k));
			}
			else
			{
				(*rat).push_back(H.dRatio(k));
			}
		int s=H.size();
		for(int i=0; i<m; i++)
		{
			if(s>i)
			{
				Arc *best=new Arc(*(H.h(i)));
				temp[i]=best;
			}
			else
			{
				Arc *neutral=new Arc(-1,0);
				temp[i]=neutral;
			}
		}
		(*res).push_back(temp);
		for(int j=wi; j<sw-1; j+=2)
		{
			ka=codes[j];
			kr=codes[j+1];
			addKmer(ka, &H, db, branches, thr, j+1, k, 0);
			addKmer(kr, &H, db, branches, thr, j+1, k, 0);
			H.getTop(m);
			if(met=='R')
			{
				(*rat).push_back(H.lRatio(k));
			}
			else
			{
				(*rat).push_back(H.dRatio(k));
			}
			s=H.size();
			for(int i=0; i<m; i++)
			{
				if(s>i)
				{
					Arc *best=new Arc(*(H.h(i)));
					temp[i]=best;
				}
				else
				{
					Arc *neutral=new Arc(-1,0);
					temp[i]=neutral;
				}
			}
			(*res).push_back(temp);
		}
		for(int j=sw; j<q; j++)
		{
			ka=codes[j];
			kr=codes[j-sw];
			addKmer(ka, &H, db, branches, thr, j+1, k, 1);
			rmKmer(kr, &H, db, branches, thr, 1);
			H.getTop(m);
			if(met=='R')
			{
				(*rat).push_back(H.lRatio(k));
			}
			else
			{
				(*rat).push_back(H.dRatio(k));
			}
			s=H.size();
			for(int i=0; i<m; i++)
			{
				if(s>i)
				{
					Arc *best=new Arc(*(H.h(i)));
					temp[i]=best;
				}
				else
				{
					Arc *neutral=new Arc(-1,0);
					temp[i]=neutral;
				}
			}
			(*res).push_back(temp);
		}
		for(int j=q-sw; j<q-wi-1; j+=2)
		{
			ka=codes[j];
			kr=codes[j+1];
			rmKmer(ka, &H, db, branches, thr, 0);
			rmKmer(kr, &H, db, branches, thr, 0);
			H.getTop(m);
			if(met=='R')
			{
				(*rat).push_back(H.lRatio(k));
			}
			else
			{
				(*rat).push_back(H.dRatio(k));
			}
			s=H.size();
			for(int i=0; i<m; i++)
			{
				if(s>i)
				{
					Arc *best=new Arc(*(H.h(i)));
					temp[i]=best;
				}
				else
				{
					Arc *neutral=new Arc(-1,0);
					temp[i]=neutral;
				}
			}
			(*res).push_back(temp);
		}
	}
}

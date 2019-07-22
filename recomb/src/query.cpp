/*
 * query.cpp
 *
 *  Created on: 23 oct. 2018
 *      Author: scholz
 */

#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include<cmath>
#include "query.h"

using namespace std;

std::vector<core::phylo_kmer_db::key_type> encode_string_views(std::string_view long_read, const size_t kmer_size)
{
	std::vector<core::phylo_kmer_db::key_type> res(0);
    	for (const auto& [kmer, code] : core::to_kmers(long_read, kmer_size))
    	{
		res.push_back(code);
		//std::cout << kmer << ": " << code << " " << std::endl;
    	}
	return res;
}

void search(const core::phylo_kmer_db& db, core::phylo_kmer_db::key_type key)
{
    if (auto entries = db.search(key); entries)
    {
        std::cout << "Found " << key << ":\n";
        for (const auto& [branch, score] : *entries)
        {
            std::cout << "\tbranch " << branch << ": " << score << '\n';
        }
    }
    else
    {
        std::cout << "Key " << key << " not found.\n";
    }
}

void readQuery(std::vector<core::phylo_kmer_db::key_type> codes, const core::phylo_kmer_db& db, std::vector<Arc>* branches, Htree *H)
{
	std::vector<Arc*> res(0);
	int k=db.kmer_size();
	double thr = core::score_threshold(k);
	int q=codes.size();
	for(int i=0; i<q; i++)
   	{
        	//cout << "adding " << codes[i] << endl;
		if (auto entries = db.search(codes[i]) ; entries)
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
		else
		{
			cout << "Nope " << codes[i] << endl;
		}
    	}
	int s=(res).size();
	double up;
	for(int i=0; i<s; i++)
	{
		up=/*log*/(thr)*(q-(*(res)[i]).getKmers());
		(*(res)[i]).updateScore(up,0);
		(*H).push(res[i]);
	}
}

void printScore(std::vector<Arc*> result, std::vector<std::string> ref)
{
	int s=result.size();
	for(int i=0; i<s; i++)
	{
		//cout << i+1 << ".";
		(*result[i]).printPlace(ref);
		cout << ": " << (*result[i]).getScore() << endl;
	}
}

void printChange(std::vector<std::vector<Arc*>> result, int shift, std::vector<std::string> ref)
{
	int s=result.size();
	int c=0;
	int aref=-1;
	for(int i=0; i<s; i++)
	{
		if(ref[(*result[i][0]).getPlace()] != ref[aref])
		{
			cout << i+shift << ",";
			//printScore(result[i], ref);
			aref=(*result[i][0]).getPlace();
			cout << ref[aref] << ",";
			c++;
		}
	}
	cout << endl;
	cout << "->" << c << " slices" << endl;
}

double lRatio(std::vector<Arc*> result, int i)
{
	int s=result.size();
	double q=0;
	for(int j=0; j<s;j++)
	{
		if(j!=i)
		{
			q+=exp((*result[j]).getScore());
		}
	}
	double ret=exp((*result[i]).getScore())/q;
	return ret;
}

void slidingWindow(std::vector<core::phylo_kmer_db::key_type> codes, int sw, int m, const core::phylo_kmer_db& db, std::vector<Arc>* branches, std::vector<std::vector<Arc*>> *res)
{
	Arc *neutral=new Arc(-1,0);
	std::vector<Arc*> rest(0);
	Htree H(rest);
	std::vector<Arc*> temp(m);
	int k=db.kmer_size();
	double thr = core::score_threshold(k);
	int q=codes.size();
	core::phylo_kmer_db::key_type ka;
	core::phylo_kmer_db::key_type kr;
	if(sw>q+k-1 || sw<k)
	{
		cout << "Dimension problems: check window size vs k-mer size and/or window size vs query length" << endl;
		cout << k << "<" << sw << "<" << q << endl;
	}
	else
	{
		std::vector<core::phylo_kmer_db::key_type> firstw(codes.begin(),codes.begin()+sw-k+1);
		readQuery(firstw, db, branches, &H);
		H.getTop(m);
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
				temp[i]=neutral;
			}
		}
		(*res).push_back(temp);
		for(int i=sw-k+1; i<q;i++)
		{
			ka=codes[i];
			kr=codes[i-sw+k-1];
			//cout << "adding " << ka << "; removing " << kr << endl;
			if (auto entries = db.search(ka) ; entries)
    			{
        			for (const auto& [branch, score] : *entries)
        			{
					(*branches)[branch].updateScore(score-thr,1);
					if(((*branches)[branch]).getKmers()==1)
					{
						(*branches)[branch].updateScore(((sw+1-k)*thr),0);
						H.push(&((*branches)[branch]));
					}
					else
					{
							H.heapUp(((*branches)[branch].ArcCheck()));
					}
        			}
    			}
			else
			{
				cout << "Nope " << ka << endl;
			}
			if (auto entries = db.search(kr) ; entries)
    			{
        			for (const auto& [branch, score] : *entries)
        			{
					(*branches)[branch].updateScore(thr-score,-1);
					if(((*branches)[branch]).getKmers()==0)
					{
						H.pop(((*branches)[branch]).ArcCheck());
					}
					else
					{
						H.heapDown(((*branches)[branch].ArcCheck()));
					}
        			}
    			}
			else
			{
				cout << "Nope " << kr << endl;
			}
			H.getTop(m);
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
					temp[i]=neutral;
				}
			}
			(*res).push_back(temp);
		}
	}
}

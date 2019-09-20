/*
 * groups.cpp
 *
 *  Created on: 05 jun. 2019
 *      Author: scholz
 */


#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include<cmath>
#include "input.h"

using namespace std;
//mainly, used to create new sequences from a set of sequences.
//two main constructions; extract subsequences and build artificial recombinants.

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


void incr(std::vector<int> *list, int a)
{
	//adds elt. a in an ordered vector *list; keeps *list ordered.
	int s=(*list).size();
	if(s==0 || a>= (*list)[s-1])
	{
		(*list).push_back(a);
	}
	else
	{
		int i=0;
		while((*list)[i]<a)
		{
			i++;
		}
		(*list).insert((*list).begin()+i,a);
	}
}

std::string getAccess(std::string leaf)
{
	//isolates access number from sequence name.
	int s=leaf.length();
	int dot=-1;
	for(int i=0; i<s; i++)
	{
		if(leaf[i]=='.')
		{
			dot=i;
		}
	}
	return leaf.substr(dot+1);
}

void writeFasta(std::vector<rappas::io::fasta> *sequences, std::string res)
{
	//writes a list of sequences in a .fasta file.
	//res is the folder in which the file "new-align.fasta" will appear.
	int s=(*sequences).size();
	ofstream write(res);
	for(int i=0; i<s; i++)
	{
		write << ">" << ((*sequences)[i].header())/*.substr(0,10)*/ << endl;
		write << (*sequences)[i].sequence();
		if(i<s-1)
		{
			write << endl;
		}
	}
}

std::vector<rappas::io::fasta> gapRm(std::vector<rappas::io::fasta> *sequences)
{
	//gap-remover; useful to create gap-free supplies from which to build recombinants
	int s=(*sequences).size();
	std::string seq="";
	std::string tmp="";
	std::string nm="";
	std::vector<rappas::io::fasta> res;
	for(int i=0; i<s; i++)
	{
		tmp="";
		seq=(((*sequences)[i]).sequence()).data();
		nm=(((*sequences)[i]).header()).data();
		for(int j=0; j<seq.length(); j++)
		{
			if(nucl(seq[j]))
			{
				tmp+=seq[j];
			}
		}
		(res).emplace_back(move(nm), move(tmp));
	}
	return res;
}

std::vector<int> splitGroups(std::vector<rappas::io::fasta> *sequences)
{
	std::string g="";
	int s=(*sequences).size();
	std::vector<int> res(0);
	for(int i=0; i<s; i++)
	{
		if(groupFromName((((*sequences)[i]).header()).data())!=g)
		{
			g=groupFromName((((*sequences)[i]).header()).data());
			res.push_back(i);
		}
	}
	res.push_back(s);
	return res;
}

void writeInfo(std::vector<string> *summary, std::string res)
{
	//writes a .txt file summarizing groups of breakpoints for a set of artificial recombinants generated with "random-q".
	//res is the folder in which the file "info-random.txt" will appear.
	int s=(*summary).size();
	ofstream write(res);
	for(int i=0; i<s; i++)
	{
		write << ">query_" << i+1 << endl;
		write << (*summary)[i];
		if(i<s-1)
		{
			write << endl;
		}
	}
}

void split_seq(std::vector<rappas::io::fasta> *sequences, std::vector<int> bp, int seq)
{
	//from a sequence seq and a list of breakpoints bp, creates one new sequence per interval.
	//replaces everything outside of the interval with gaps
	int prev=-1;
	int bits=bp.size()+1;
	std::string s= (((*sequences)[seq]).sequence()).data();
	int len=s.length();
	std::string cp;
	std::string name= (((*sequences)[seq]).header()).data();
	//name.pop_back();
	bp.push_back(len);
	for (int i=0; i< bits; i++)
	{
		if(bp[i] > prev && bp[i] < len+1)
		{
			cp=s;
			for(int j=0; j<len-1; j++)
			{
				if(j<(prev+1) || j>bp[i])
				{
					cp[j]='-';
				}
			}
		(*sequences).emplace_back(move(name+"_"+to_string(i+1)), move(cp));
		//cout << name << endl;
		//cout << move(name+"_"+to_string(i+1)) << endl;
		prev=bp[i];
		}
		else
		{
			cout << "breakpoints should be in increasing order and smaller than sequence size:" << endl;
			cout << prev << "; " << bp[i] << "; " << len << endl;
			i=bits;
		}
	}
	bp.pop_back();
}

void new_q(std::vector<rappas::io::fasta> *sequences, std::vector<rappas::io::fasta> *nq, std::vector<int> *bp, std::vector<int> *seq, int ref)
{
	//creates an artificial recombinant from breakpoints *bp and origin of segments *seq.
	int bits=(*bp).size()+1;
	if(bits==(*seq).size())
	{
		std::string res="";
		std::string s="";
		int prev=-1;
		int len=0;
		int diff=0;
		std::string name="query_"+to_string(ref);
		for (int i=0; i< bits; i++)
		{
			s= (((*sequences)[(*seq)[i]]).sequence()).data();
			len=s.length();
			if(len > prev)
			{
				if(i==bits-1)
				{
					(*bp).push_back(len-1);
				}
				if((*bp)[i] > len-1)
				{
					(*bp)[i]=len-1;
					cout << "sequence too short for 2nd breakpoint; breakpoint shifted left" << endl;
				}
				for(int j=0; j<len-1; j++)
				{
					if(!(j<(prev+1) || j>(*bp)[i]))
					{
						res.push_back(s[j]);
						if(nucl(s[j])==0)
						{
							diff++;
						}
					}
				}
				prev=(*bp)[i];
				(*bp)[i]=(*bp)[i]-diff+1;
				//diff=0;

			}
			else
			{
				(*seq)[i]=-1;
				cout << "sequence too short for 1st breakpoint; insertion ignored" << endl;
				if(i==bits-1)
				{
					(*bp).push_back(-1);
				}
			}
		}
		(*nq).emplace_back(move(name), move(res));
		/*if((*bp).size() == bits)
		{
			(*bp).pop_back();
		}*/
	}
	else
	{
		cout << "Number of bits must correspond to number of sequences:" << endl;
		cout << bits << "; " << (*seq).size() << endl;
	}
}

std::string record(std::vector<rappas::io::fasta> *sequences, std::vector<int> bp, std::vector<int> seq)
{
	//records sequence of groups and breakpoints in a string.
	int bits=bp.size();
	std::string res="";
	if(bits==seq.size())
	{
		//res.append(groupFromName((((*sequences)[seq[0]]).header()).data()));
		for(int i=0; i< bits; i++)
		{
			if(seq[i]> -1)
			{
				res.append(groupFromName((((*sequences)[seq[i]]).header()).data()));
				res.append(",");
				res.append(to_string(bp[i]));
				res.append(",");
			}
		}
		res.pop_back();
	}
	else
	{
		cout << "Number of bits must correspond to number of sequences:" << endl;
		cout << bits << "; " << seq.size() << endl;
	}
	return res;	
}

void random_q(std::vector<rappas::io::fasta> *sequences, std::vector<rappas::io::fasta> *nq, int n, std::vector<std::string> *summary, int b_max, int p_max)
{
	//generates n random artificial recombinants, records the information
	//b_min and b_max control min and max of breakpoints
	//p_max controls the max number of contributing sequences
	std::vector<int> bp(0);
	std::vector<int> seq(0);
	int s=(*sequences).size();
	int p=-1;
	int rec=0;
	int len=0;
	int g=-1;
	int g_prev=-1;
	int beg=0;
	std::vector<int> groups=splitGroups(sequences);
	for(int i=0; i<n; i++)
	{
		g=rand() %(groups.size()-1);
		seq.push_back(groups[g]+(rand() %(groups[g+1]-groups[g])));
		//seq.push_back(rand() %s);
		len=(((*sequences)[seq[0]]).sequence()).length();
		beg=0;
		while(! nucl((((*sequences)[seq[0]]).sequence())[beg]))
		{
			beg++;
		}
		for(int j=0; j<b_max; j++)
		{
			incr(&bp, beg+ rand() %(len+1-beg));
			if(j<p_max-1)
			{
				g_prev=g;
				while(g==g_prev)
				{
					g=rand() %(groups.size()-1);
				}
				p=groups[g]+(rand() %(groups[g+1]-groups[g]));
				//p=rand() %s;
			}
			else
			{
				p=seq[rand() %p_max];
				while(p==seq[j] && rec<10)
					{
						p=seq[rand() %p_max];
						rec++;
					}
				rec=0;
			}
			seq.push_back(p);
		}
		//cout << i+1 << endl;
		new_q(sequences, nq, &bp, &seq, i+1);
		(*summary).push_back(record(sequences, bp,seq));
		bp.clear();
		seq.clear();
	}
}

std::vector<std::vector<std::string>> readRes(std::string res)
{
	std::vector<std::vector<std::string>> read(0);
	std::string tmp="";
	std::vector<std::string> tab(0);
	ifstream data(res);
	if(data)
	{
		std::string line="";
		while(getline(data,line))
		{
			if(line[0] !='>')
			{
				for(int i=0; i<line.length(); i++)
				{
					if(line[i] !=',')
					{
						tmp+=line[i];				
					}
					else
					{
						tab.push_back(tmp);
						tmp="";
					}
				}
				tab.push_back(tmp);
				tmp="";
				if(tab.size() %2 == 0)
				{
					tab.insert(tab.begin(), "0");
				}
				read.push_back(tab);
				tab.clear();
			}
		}
		data.close();
	}
	else
	{
		cout << "No file there (or something like that)" << endl;
	}
	return read;
}

void compRes(std::string res1, std::string res2)
{
	std::vector<std::vector<std::string>> read1=readRes(res1);
	std::vector<std::vector<std::string>> read2=readRes(res2);
	int s=read2.size();
	int val=0;
	int t=-1;
	double total=0;
	double max=-1;
	int k1=-1;
	int k2=-1;
	double comp=-1;
	for(int i=0; i<s; i++)
	{
		t=(read2[i]).size();
		//max=stoi(read2[i][t-1])-stoi(read2[i][0])+1;
		max=0;
		k1=0;
		k2=0;
		comp=0;
		for(int j=stoi(read2[i][0]); j<stoi(read2[i][t-1])+1; j++)
		{
			if(j>stoi(read1[i][k1+2]))
			{
				k1+=2;
			}
			if(j>stoi(read2[i][k2+2]))
			{
				k2+=2;
			}
			//if(read2[i][k2+1] !="N/A")
			{
				max++;
				if(read1[i][k1+1]==read2[i][k2+1])
				{
					comp++;
				}
			}
		}
		//if(comp*100/max<98)
		{
			cout << i+1 << ": " << comp*100/max << endl;
		}
		if(max!=0)
		{
			val++;
			total+=comp*100/max;
		}
	}
	cout << "------------" << endl;
	cout << total/val << endl;
}

void compMosaic(std::string res1, std::string res2)
{
	std::vector<std::vector<std::string>> read1=readRes(res1);
	std::vector<std::vector<std::string>> read2=readRes(res2);
	int r=0;
	int s=0;
	int i=0;
	int j=0;
	double check=0;
	for(int k=0; k<read1.size(); k++)
	{
		r=(read1[k]).size();
		s=(read2[k]).size();
		i=1;
		j=1;
		check=0;
		if(r==s)
		{
			while(check==0 && i<s)
			{
				if(read1[k][i] !=read2[k][i])
				{
					check=nan("");
				}
				i+=2;
			}
		}
		if(r < s)
		{
			while(j<s && check > -1)
			{
				if(read1[k][i] == read2[k][j])
				{
					i+=2;
				}
				else if(i==1 || read1[k][i-2] != read2[k][j])
				{
					check++;
				}
				j+=2;
				if(j==s && i<r)
				{
					check=nan("");
				}
			}
		}
		if(s < r)
		{
			while(i<r && check < 1)
			{
				if(read1[k][i] == read2[k][j])
				{
					j+=2;
				}
				else if(j==1 || read1[k][i] != read2[k][j-2])
				{
					check--;
				}
				i+=2;
				if(i==r && j<s)
				{
					check=nan("");
				}
			}
		}
		if(check ==0)
		{
			int m=2;
			for(int l=2; l<r-1; l+=2)
			{
				while(read2[k][m-1]==read2[k][m+1])
				{
					m+=2;
				}
				cout << stoi(read2[k][m])-stoi(read1[k][l]) << ",";
				m+=2;
			}
		}
	}
	cout << endl;
}



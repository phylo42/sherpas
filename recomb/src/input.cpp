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

void incr(std::vector<int> *list, int a)
{
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

void getId(std::string align, std::string res) //only once; depends on the sequences name; current version works for names starting with group name followed with '.'
//res is the folder in which the file "leaf-summary.txt" will appear.
{
	const auto sequences = rappas::io::read_fasta(align);
	int i=0;
	int c=0;
	int s=sequences.size();
	std::string line="";
	ofstream write(res+"leaves-summary.txt");
	for(int j=0; j<s; j++)
	{
		line=sequences[j].header();
		while(line[i] !='.' || c<0)
		{
			if(line[i]=='.')
			{
				c++;
			}
			i++;
		}
		write << line.substr(0, 10) << " " << line.substr(0,i) << endl;
		i=0;
		c=0;
	}
}

void writeFasta(std::vector<rappas::io::fasta> *sequences, std::string res)
//res is the folder in which the file "new-align.fasta" will appear.
{
	int s=(*sequences).size();
	ofstream write(res);
	for(int i=0; i<s; i++)
	{
		write << ">" << ((*sequences)[i].header()).substr(0,10) << endl;
		write << (*sequences)[i].sequence();
		if(i<s-1)
		{
			write << endl;
		}
	}
}

void writeInfo(std::vector<string> *summary, std::string res)
//res is the folder in which the file "new-align.fasta" will appear.
{
	int s=(*summary).size();
	ofstream write(res);
	for(int i=0; i<s; i++)
	{
		write << "> query_" << i+1 << endl;
		write << (*summary)[i];
		if(i<s-1)
		{
			write << endl;
		}
	}
}

std::string cutsec(rappas::io::fasta seq, int cut, int side)
{
	std::string res="";
	/*if(cut>(seq.sequence()).size() || cut<0)
	{
		cout << "breakpoint out of bounds" << endl; 
	}
	else
	{
		if(side==1)
		{
			res=seq(0,cut);
		}
		if(side==2)
		{
			res=seq(cut);
		}
	}*/
	return res;
}

void split_seq(std::vector<rappas::io::fasta> *sequences, std::vector<int> bp, int seq)
{
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
	int bits=(*bp).size()+1;
	if(bits==(*seq).size())
	{
		std::string res="";
		std::string s="";
		int prev=-1;
		int len=0;
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
					}
				}
				prev=(*bp)[i];
			}
			else
			{
				(*seq)[i]=-1;
				cout << "sequence too short for 1st breakpoint; insertion ignored" << endl;
			}
		}
		(*nq).emplace_back(move(name), move(res));
		if((*bp).size() == bits)
		{
			(*bp).pop_back();
		}
	}
	else
	{
		cout << "Number of bits must correspond to number of sequences:" << endl;
		cout << bits << "; " << (*seq).size() << endl;
	}
}

std::string record(std::vector<rappas::io::fasta> *sequences, std::vector<int> bp, std::vector<int> seq)
{
	int bits=bp.size();
	std::string res="";
	if(bits==seq.size()-1)
	{
		res.append((((*sequences)[seq[0]]).header()));
		res.pop_back();
		for(int i=1; i< bits+1; i++)
		{
			if(seq[i]> -1)
			{
				res.append(",");
				res.append(to_string(bp[i-1]));
				res.append(",");
				res.append((((*sequences)[seq[i]]).header()));
				res.pop_back();
			}
		}
	}
	else
	{
		cout << "Number of bits must correspond to number of sequences:" << endl;
		cout << bits << "; " << seq.size() << endl;
	}
	return res;	
}

void random_q(std::vector<rappas::io::fasta> *sequences, std::vector<rappas::io::fasta> *nq, int n, std::vector<std::string> *summary, int b_min, int b_max, int p_max)
{
	std::vector<int> bp(0);
	std::vector<int> seq(0);
	int s=(*sequences).size();
	int len=-1;
	int b=-1;
	int p=-1;
	int rec=0;	
	for(int i=0; i<n; i++)
	{
		seq.push_back(rand() %s);
		len=(((*sequences)[seq[0]]).sequence()).length();
		if(b_min != b_max)
		{
			b=rand() %(b_max-b_min+1) +b_min;
		}
		else
		{
			b=b_min;
		}
		for(int j=0; j<b; j++)
		{
			incr(&bp, rand() %len);
			if(j<p_max-1)
			{
				p=rand() %s;
			}
			else
			{
				p=seq[rand() %p_max];
				while(p==seq[j] || rec<10)
					{
						p=seq[rand() %p_max];
						rec++;
					}
			rec=0;
			}
			seq.push_back(p);
		}
		//cout << i << endl;
		new_q(sequences, nq, &bp, &seq, i+1);
		(*summary).push_back(record(sequences, bp,seq));
		bp.clear();
		seq.clear();
	}
}


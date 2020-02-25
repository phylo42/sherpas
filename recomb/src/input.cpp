/*
 * groups.cpp
 *
 *  Created on: 05 jun. 2019
 *      Author: scholz
 */


#include<iostream>
#include<fstream>
#include <random>
#include<vector>
#include <string>
#include<cmath>
#include "input.h"

using namespace std;
//used for all kind of tests (simulations and measures). Not supposed to be part of the released version.

void incr(std::vector<int> *list, int a)
{
	//adds elt. a in an ordered vector *list by keeping it ordered.
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

std::string groupFromName(std::string leaf)
{
	//gets group from sequence name.
	//WARNING: depends on how the sequence are named: version 1 works for names starting with group name followed with '.' (HIV); version 2 works for names ending with "_[gr_name]" (HBV).
	//not needed for main algo (replaced by groupFromTab) but still needed in the simulating recombinants part (just because it is easier just to leave it here than to update that part).
	int i=0;
	std::string res="";
	/*while(leaf[i] !='.' && i<leaf.length()) // beginning of version 1 (HIV)
	{
		i++;
	}
	res=leaf.substr(0,i);*/ //end of version 1 (HIV)
	/*if(res[0]=='A' || (res[0]=='0' && res[1]=='1'))
	{
		res="A";
	}*/
	while(leaf[i] !='_') // beginning of version 2 (HBV)
	{
		i++;
	}
	res=leaf.substr(i+1); // end of version 2 (HBV)
	return res;
}

int bpDist(std::vector<int> bp, int max)
{
	//computes the minimal distance between two breakpoints, and between breakpoints and extremities (0 and max).
	//assumes the breakpoints are ordered and max greater than them all (otherwise I don't know what happens).
	//int res=bp[0]+max-bp.back(); //for circular
	int res=bp[0]; //for linear
	bp.push_back(max);//for linear
	for(int i=0; i<bp.size()-1; i++)
	{
		if(bp[i+1]-bp[i] < res || res<0)
		{
			res=bp[i+1]-bp[i];
		}
	}
	return res;
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
	int bits=(*bp).size();
	if(bits==(*seq).size())
	{
		std::string res="";
		std::string s="";
		int prev=0;
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
					(*bp).push_back(len);
				}
				if((*bp)[i] > len)
				{
					(*bp)[i]=len;
					cout << "sequence too short for 2nd breakpoint; breakpoint shifted left" << endl;
				}
				for(int j=0; j<len; j++)
				{
					if(j>prev-1 && j<(*bp)[i])
					{
						res.push_back(s[j]);
						if(s[j]=='-')
						{
							diff++;
						}
					}
				}
				prev=(*bp)[i];
				(*bp)[i]=(*bp)[i]-diff;
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
		while((*bp).size() > bits)
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
	//records sequence of groups and breakpoints in a string.
	int bits=bp.size();
	std::string res="0,";
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

void random_q(std::vector<rappas::io::fasta> *sequences, std::vector<rappas::io::fasta> *nq, int n, std::vector<std::string> *summary)
{
	//generates n random artificial recombinants, records the information
	//b_min and b_max control min and max of breakpoints
	//p_max controls the max number of contributing sequences
	//pre-queries must have same length (otherwise things might go wrong)
	std::vector<int> bp(0);
	std::vector<int> seq(0);
	int s=(*sequences).size();
	int p=-1;
	int rec=0;
	int len=0;
	int g=-1;
	int g_prev=-1;
	int b_max=0;
	int p_max=0;
	int i=0;
	std::geometric_distribution<> dist_bp(0.2);
	std::geometric_distribution<> dist_par(0.8);
	std::random_device rd;
	std::vector<int> groups=splitGroups(sequences);
	while(i<n)
	{
		//g=rand() %(groups.size()-1);
		//seq.push_back(groups[g]+(rand() %(groups[g+1]-groups[g])));
		p=rand() %s;
		while(p+1==88 || p+1==165 || p+1==177 || p+1==188 || p+1==391 || p+1==710 || p+1==711 || p+1==941 || p+1==1129 || p+1==1155 || p+1==1175 || p+1==1295 || p+1==1505 || p+1==1515 || p+1==1517 || p+1==1518 || p+1==1520 || p+1==1523 || p+1==1524 || p+1==1556 || p+1==1575 || p+1==1577 || p+1==1587 || p+1==1613 || p+1==1727 || p+1==1729 || p+1==1754)
		{
			p=rand() %s;
		}
		seq.push_back(p);
		len=(((*sequences)[seq[0]]).sequence()).length();
		b_max=dist_bp(rd)+1;
		p_max=dist_par(rd)+2;
		if(b_max>p_max-2)
		{
			cout << b_max << " bp and " << p_max << " parents." << endl;
			for(int j=0; j<b_max; j++) //b_max if linear ; b_max-1 if circular
			{
				incr(&bp, rand() %(len+1));
				if(j<p_max-1)
				{
					p=rand() %s;
					while(groupFromName((((*sequences)[seq[0]]).header()).data())==groupFromName((((*sequences)[p]).header()).data()) /*|| p+1==88 || p+1==165 || p+1==177 || p+1==188 || p+1==391 || p+1==710 || p+1==711 || p+1==941 || p+1==1129 || p+1==1155 || p+1==1175 || p+1==1295 || p+1==1505 || p+1==1515 || p+1==1517 || p+1==1518 || p+1==1520 || p+1==1523 || p+1==1524 || p+1==1556 || p+1==1575 || p+1==1577 || p+1==1587 || p+1==1613 || p+1==1727 || p+1==1729 || p+1==1754*/)
					{
						p=rand() %s;
					}
				}
				else
				{
					p=seq[rand() %p_max];
					while(p==seq[j])
						{
							p=seq[rand() %p_max];
							rec++;
						}
					rec=0;
				}
				seq.push_back(p);
			}
			//incr(&bp, rand() %(len+1));//only if circular
			//seq.push_back(seq[0]); //only if circular
			cout << bpDist(bp,(((*sequences)[p]).sequence()).length()) << " -> ";
			bp.push_back((((*sequences)[p]).sequence()).length());
			new_q(sequences, nq, &bp, &seq, i+1);
			cout << bpDist(bp,(((*sequences)[p]).sequence()).length()) << endl;
			if(bpDist(bp,(((*sequences)[p]).sequence()).length()) >100)
			{
				(*summary).push_back(record(sequences, bp,seq));			
				i++;
				cout << "query " << i << " built" << endl;
			}
			else
			{
				(*nq).pop_back();
			}
		}
		bp.clear();
		seq.clear();
	}
}

std::vector<std::vector<std::string>> readRes(std::string res)
{
	//reads result file as produced by sherpas
	std::vector<std::vector<std::string>> read(0);
	std::string tmp="";
	std::vector<std::string> tab(0);
	ifstream data(res);
	if(data)
	{
		std::string line="";
		int check=0;
		while(getline(data,line))
		{
			if(check==0)
			{
				check=(line[0]=='>');
			}
			if(check==1)
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
					/*if(tab.size() %2 == 0)
					{
						tab.insert(tab.begin(), "0");
					}*/
					read.push_back(tab);
					tab.clear();
				}
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

void compRes(std::string res1, std::string res2, std::string group)
{
	//computes sensitivity and precision
	std::vector<std::vector<std::string>> read1=readRes(res1);
	std::vector<std::vector<std::string>> read2=readRes(res2);
	int s=read2.size();
	double valS=0;
	double totalS=0;
	double maxS=0;
	double compS=0;
	double valP=0;
	double totalP=0;
	double maxP=0;
	double compP=0;
	int k1=-1;
	int k2=-1;
	
	for(int i=0; i<s; i++)
	{
		maxS=0;
		compS=0;
		maxP=0;
		compP=0;
		k1=0;
		k2=0;
		for(int j=stoi(read2[i][0])+1; j<stoi((read2[i]).back())+1; j++)
		{
			if(j>stoi(read1[i][k1+2]))
			{
				k1+=2;
			}
			if(j>stoi(read2[i][k2+2]))
			{
				k2+=2;
			}
			if(read1[i][k1+1] == group || group=="all") //sensitivity
			{
				maxS++;
				if(read1[i][k1+1]==read2[i][k2+1])
				{
					compS++;
				}
			}
			if(read2[i][k2+1] == group || (group=="all" && read2[i][k2+1] !="N/A")) //precision
			{
				maxP++;
				if(read1[i][k1+1]==read2[i][k2+1])
				{
					compP++;
				}
			}
		}
		if(maxS !=0 /*&& (!isGroup(read1[i], "A") && !isGroup(read1[i], "A1"))*/)
		{
			valS+=maxS;
			totalS+=compS*100;
		}
		if(maxP !=0 /*&& (!isGroup(read1[i], "A") && !isGroup(read1[i], "A1"))*/)
		{
			valP+=maxP;
			totalP+=compP*100;
		}
	}
	cout << totalS/valS << " & " << totalP/valP << " & ";// << endl;
}

void compMosaic(std::string res1, std::string res2)
{
	//compares simulated mosaic with true mosaic.
	std::vector<std::vector<std::string>> read1=readRes(res1);
	std::vector<std::vector<std::string>> read2=readRes(res2);
	int r=0;
	int s=0;
	int i=0;
	int j=0;
	int l=0;
	double check=0;
	std::vector<double> count (4);
	int p=0;
	int e=0;
	for(int k=0; k<read1.size(); k++)
	{
		l=3;
		//s=(read2[k]).size();
		while(l<(read1[k]).size())
		{
			if(read1[k][l]==read1[k][l-2])
			{
				(read1[k]).erase((read1[k]).begin()+l-1);
				(read1[k]).erase((read1[k]).begin()+l-1);
			}
			else
			{
				l+=2;
			}			
		}
		/*while(stoi(read1[k][2])<stoi(read2[k][0]))
		{
			(read1[k]).erase((read1[k]).begin());
			(read1[k]).erase((read1[k]).begin());
		}
		while(stoi(read1[k][(read1[k]).size()-3])>stoi((read2[k]).back()))
		{
			(read1[k]).pop_back();
			(read1[k]).pop_back();
		}*/
		l=1;
		while(l<(read2[k]).size())
		{
			if(read2[k][l]=="N/A") //ratio-control version (l from 1 to size)
			{
				(read2[k]).erase((read2[k]).begin()+l-1);
				(read2[k]).erase((read2[k]).begin()+l-1);
				if(l>1 && l<read2[k].size() && read2[k][l]==read2[k][l-2])
				{
					(read2[k]).erase((read2[k]).begin()+l-1);
					(read2[k]).erase((read2[k]).begin()+l-1);
				}
			}
			else
			{
				l+=2;
			}			
		}
		s=(read2[k]).size();
		r=(read1[k]).size();
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
		/*if(check ==0)
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
		}*/
		//if(k+1<2232 || k+1 > 3877)
		//if(k+1==37 || k+1==201 || k+1==326 || k+1==511 || k+1==663 || k+1==1081 || k+1==1393 || k+1==1661 || k+1==2222 || k+1==3000 || k+1==3404 || k+1==3737 || k+1==4712 || k+1==5210 || k+1==5992 || k+1==6181 || k+1==6752 || k+1==7113 || k+1==7880 || k+1==8024 || k+1==8066 || k+1==8614 || k+1==9001 || k+1==9216 || k+1==9970)
		{
			if(check==0)
			{
				count[0]++;
				//cout << k+1 << ",";
			}
			if(check>0)
			{
				count[1]++;
				//cout << k+1 << ",";
			}
			if(check<0)
			{
				count[2]++;
				//cout << k+1 << ",";
			}
			if(isnan(check))
			{
				count[3]++;
				//cout << k+1 << ",";
			}
		}
	}
	//cout << endl << "------------" << endl;
	cout << /*thr << " & " <<*/ count[0]*100/read1.size() << " & " << count[1]*100/read1.size() << " & " << count[2]*100/read1.size() << " & " << count[3]*100/read1.size() << " \\\\" << endl;
}

void compMosaicCirc(std::string res1, std::string res2)
{
	//same as above for circular genomes
	std::vector<std::vector<std::string>> read1=readRes(res1);
	std::vector<std::vector<std::string>> read2=readRes(res2);
	int r=0;
	int s=0;
	int i=0;
	int j=0;
	int l=0;
	double check=0;
	std::vector<double> count (4);
	int p=0;
	int e=0;
	for(int k=0; k<read1.size(); k++)
	{
		l=3;
		//s=(read2[k]).size();
		while(l<(read1[k]).size())
		{
			if(read1[k][l]==read1[k][l-2])
			{
				(read1[k]).erase((read1[k]).begin()+l-1);
				(read1[k]).erase((read1[k]).begin()+l-1);
			}
			else
			{
				l+=2;
			}			
		}
		l=1;
		while(l<(read2[k]).size())
		{
			if(read2[k][l]=="N/A") //ratio-control version (l from 1 to size)
			{
				(read2[k]).erase((read2[k]).begin()+l-1);
				(read2[k]).erase((read2[k]).begin()+l-1);
				if(l>1 && l<read2[k].size() && read2[k][l]==read2[k][l-2])
				{
					(read2[k]).erase((read2[k]).begin()+l-1);
					(read2[k]).erase((read2[k]).begin()+l-1);
				}
			}
			else
			{
				l+=2;
			}			
		}
		if((read1[k]).size() > 3 && read1[k][1]==read1[k][(read1[k]).size()-2]) //circular genome only
		{
			(read1[k]).pop_back();
			(read1[k]).pop_back();
		}
		if((read2[k]).size() > 3 && read2[k][1]==read2[k][(read2[k]).size()-2]) //circular genome only
		{
			(read2[k]).pop_back();
			(read2[k]).pop_back();
		}
		s=(read2[k]).size();
		r=(read1[k]).size();
		check=nan("");
		l=1;
		if(r==s)
		{
			while(isnan(check) && l<s)
			{
				i=1;
				j=1;
				check=0;
				if(l>1)
				{
					(read2[k]).push_back(read2[k][1]);
					(read2[k]).push_back(read2[k][0]);
					(read2[k]).erase((read2[k]).begin());
					(read2[k]).erase((read2[k]).begin());
				}
				l+=2;		
				while(check==0 && i<s)
				{
					if(read1[k][i] !=read2[k][i])
					{
						check=nan("");
					}
					i+=2;
				}
			}
		}
		if(r < s)
		{
			while(isnan(check) && l<s)
			{
				i=1;
				j=1;
				check=0;
				if(l>1)
				{
					(read2[k]).push_back(read2[k][1]);
					(read2[k]).push_back(read2[k][0]);
					(read2[k]).erase((read2[k]).begin());
					(read2[k]).erase((read2[k]).begin());
				}
				l+=2;
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
		}
		if(s < r)
		{
			while(isnan(check) && l<s)
			{
				i=1;
				j=1;
				check=0;
				if(l>1)
				{
					(read2[k]).push_back(read2[k][1]);
					(read2[k]).push_back(read2[k][0]);
					(read2[k]).erase((read2[k]).begin());
					(read2[k]).erase((read2[k]).begin());
				}
				l+=2;
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
		}
		//if(k+1<2232 || k+1 > 3877)
		{
			if(check==0)
			{
				count[0]++;
				//cout << k+1 << ",";
			}
			if(check>0)
			{
				count[1]++;
				//cout << k+1 << ",";
			}
			if(check<0)
			{
				count[2]++;
				//cout << k+1 << ",";
			}
			if(isnan(check))
			{
				count[3]++;
				//cout << k+1 << ",";
			}
		}
	}
	//cout << endl << "------------" << endl;
	cout << /*thr << " & " <<*/ count[0]*100/read1.size() << " & " << count[1]*100/read1.size() << " & " << count[2]*100/read1.size() << " & " << count[3]*100/read1.size() << " \\\\" << endl;
}

void compComp(std::string res1, std::string res2, int shift)
{
	//like compMosaic except the order and multiplicity of the groups are not taken into account
	std::vector<std::vector<std::string>> read1=readRes(res1);
	std::vector<std::vector<std::string>> read2=readRes(res2);
	std::vector<std::string> compo1 (0);
	std::vector<std::string> compo2 (0);
	std::vector<int> count (4);
	int s=0;
	int r=0;
	double check=0;
	for(int k=0; k<read1.size(); k++)
	{
		if(shift>0)
		{
			while(stoi(read1[k][2])<shift)
			{
				(read1[k]).erase((read1[k]).begin());
				(read1[k]).erase((read1[k]).begin());
			}
			while(stoi(read1[k][(read1[k]).size()-3])>stoi((read1[k]).back())-shift)
			{
				(read1[k]).pop_back();
				(read1[k]).pop_back();
			}
		}
		compo1.clear();
		compo2.clear();
		for(int i=1; i<(read1[k]).size(); i+=2)
		{
			if(isGroup(compo1, read1[k][i])==0)
			{
				compo1.push_back(read1[k][i]);
				//cout << read1[k][i] << ",";
			}
		}
		//cout << " // ";
		for(int i=0; i<(read2[k]).size(); i++)
		{
			//if(((read2[k][i]).length()<3 || read2[k][i]=="CPZ" || read2[k][i]=="01_AE") && isGroup(compo2, read2[k][i])==0) // comet not considering crf
			//if(isGroup(compo2, read2[k][i])==0) //comet considering crf
			if(isGroup(compo2, read2[k][i])==0 && (i %2)==1 && read2[k][i] !="N/A") //Sh/Sc/jpHMM version
			{
				compo2.push_back(read2[k][i]);
				//cout << read2[k][i] << ",";
			}
		}
		r=compo1.size();
		s=compo2.size();
		check=s-r;
		if(r>=s)
		{
			for(int i=0; i<s; i++)
			{
				if(isGroup(compo1, compo2[i])==0)
				{
					check=nan("");
				}
			}
		}
		else
		{
			for(int i=0; i<r; i++)
			{
				if(isGroup(compo2, compo1[i])==0)
				{
					check=nan("");
				}
			}
		}
		if(check==0)
		{
			count[0]++;
			//cout << k+1 << ",";
		}
		if(check>0)
		{
			count[1]++;
		}
		if(check<0)
		{
			count[2]++;
		}
		if(isnan(check))
		{
			count[3]++;
		}
		//cout << " -> " << check << endl;
	}
	cout << count[0]*100/read1.size() << " & " << count[1]*100/read1.size() << " & " << count[2]*100/read1.size() << " & " << count[3]*100/read1.size() << " \\\\" << endl;
}

void interpret(std::vector<std::string> *res, std::string seq, int k)
{
	//translates coordinates within gapfree sequences into coordinates within the alignment.
	// Careful about the way '?' are treated (gaps in the current version)
	//for fixres use only
	int p=0;
	int add=0;
	int j=1;
	int i=0;
	while(i<seq.length())
	{
		while(i>=stoi((*res)[p])+add && p+2<(*res).size())
		{
			(*res)[p]=to_string(stoi((*res)[p])+add);
			p+=2;
		}
		if(seq[i]=='-'  || seq[i]=='?')
		{
			add++;
		}
		/*if(i<k/2 && nucl(seq[i])==0 && seq[i] !='-'  && seq[i] != '?' )
		{
			i=i-k/2;
		}
		if(i<seq.length()-k/2)
		{
			if(nucl(seq[i+k/2])==0 && seq[i+k/2] !='-'  && seq[i+k/2] != '?')
			{
				add+=k;
				if(i<k/2)
				{
					add=add-k/2+i+1;
				}
				i+=k/2+1;
				j=0;
				while(j<k && i+j<seq.length())
				{
					if(nucl(seq[i+j])==0 && seq[i+j] !='-'  && seq[i+j] != '?' )
					{
						add+=j+1;
						i+=j+1;
						j=0;
					}
					else
					{
						j++;
					}
				}
				if(i+j>=seq.length())
				{
					add+=j+1-k;
				}
				i+=k/2-1;
			}
		}*/
		i++;
	}
	while(p<(*res).size())
	{
		(*res)[p]=to_string(stoi((*res)[p])+add);
		p+=2;
	}
}

void fixRes(std::string res, std::string rep, std::string seq, int k)
{
	//translates breakpoints coordinates within sequence into breakpoints coordinates with regard to the alignment seq. All sequences of the alignment must have been processed and must appear in the same order.
	//useful when bp have not been recorded wrt the gapfree sequence (scueal dataset).
	//keep for release purpose ?
	std::vector<std::vector<std::string>> read=readRes(res);
	std::vector<rappas::io::fasta> sequences= rappas::io::read_fasta(seq);
	std::vector<std::string> Rtmp (0);
	std::string Stmp="";
	ofstream writef(rep+"-fixed.txt");
	//ofstream writep(rep+"-pruned.txt");
	if(read.size() == sequences.size())
	{
		for(int i=0; i<read.size(); i++)
		{
			writef << ">" << (sequences[i]).header() << endl;
			//writep << ">" << (sequences[i]).header() << endl;
			Rtmp=read[i];
			Stmp=((sequences[i]).sequence()).data();
			if((Stmp).length() != stoi(Rtmp.back()))
			{
				interpret(&Rtmp, Stmp, k);
			}
			writef << Rtmp[0] << ",";
			//writep << stoi(Rtmp[0])+shift << ",";
			for(int j=1; j<Rtmp.size()-1; j++)
			{
				writef << Rtmp[j] << ",";
				//writep << Rtmp[j] << ",";
			}
			writef << Rtmp.back() << endl;
			//writep << stoi(Rtmp.back())-shift << endl;
		}
	}
	else
	{
		cout << "files for res and sequences disagree in size" << endl;
		cout << read.size() << "; " << sequences.size() << endl;
	}
}

void statsNA(std::string res)
{
	//shows how many N/A segments a result contains, and the average length of such segments
	std::vector<std::vector<std::string>> read=readRes(res);
	int l=0;
	double length=0;
	double seg=0;
	for(int i=0; i<read.size(); i++)
	{
		l=1;
		while(l<(read[i]).size())
		{
			if(read[i][l]=="N/A")
			{
				seg++;
				length+=stoi(read[i][l+1])-stoi(read[i][l-1]);
			}
			l+=2;
		}
	}
	cout << seg << " & " << length/seg << endl;
}

std::string splitM(std::string res, int p)
{
	//gets p^th group of a result of the form g_1*...*g_n
	//useless?
	int k=0;
	int l=0;
	std::string r="";
	while(l+2<p || k<res.length())
	{
		if(res[k]=='*')
		{
			l++;
		}
		if(l+1==p && res[k]!='*')
		{
			r+=res[k];
		}
		k++;
	}
	return r;
}

int isGroup(std::vector<std::string> res, std::string gr)
{
	int r=0;
	int i=0;
	while(i<res.size() && r==0)
	{
		if(res[i]==gr /*|| (res[i]=="A" && gr[0]=='A') || (res[i]=="N" && gr=="CPZ")*/)
		{
			r=1;
		}
		i++;
	}
	return r;
}

void read_jstyle(std::string res, std::string rep)
{
	//translates jpHMM output format into ours
	//useless once we synchronize our output format with the jpHMM one, but current one easier for tests (all measure fonctions are written wrt our format)
	ofstream write(rep);
	ifstream data(res);
	if(data)
	{
		int i=0;
		std::string line="";
		std::string tmp="";
		std::vector<std::string> rec (0);
		int s=0;
		//int check=0;
		while(getline(data,line))
		{
			if(line[0]=='>')
			{
				i=0;
				s++;
				while(line[i] != '(' && i<line.length())
				{
					i++;
				}
				write << line.substr(0,i-1) << endl;
			}
			else if(line.length() !=0 && s>0)
			{
				for(int j=0; j<line.length(); j++)
				{
					if(line[j]=='\t')
					{
						rec.push_back(tmp);
						tmp="";
					}
					else
					{
						tmp+=line[j];
					}
					if(j+1==line.length())
					{
						rec.push_back(tmp);
						tmp="";
					}
				}
				//write << line << "\t" << line.length() << endl;
			}
			if(line.length()==0 && s>0)
			{
				rec[0]="0";
				for(int k=0; k<rec.size()-2; k+=3)
				{
					if(rec[k+2].length() > 5)
					{
						rec[k+2]="N/A";
					}
					write << rec[k] << "," << rec[k+2] << ",";
				}
				write << rec[rec.size()-2] << endl;
				rec.clear();
			}
		}
	}
	else
	{
		cout << "No file there (or something like that)" << endl;
	}
}

void read_comet(std::string res, std::string rep)
{
	//careful: comet skips colons in sequence names (and maybe other symbols as well).
	ofstream write(rep);
	ifstream data(res);
	if(data)
	{
		int l=0;
		std::string line="";
		int sp=0;
		while(getline(data,line))
		{
			if(l>0)
			{
				write << ">";
				sp=0;
				for(int i=0;sp<3; i++)
				{
					if(line[i]=='\t')
					{
						sp++;
					}
					if((sp==0 || sp==2))
					{
						if(line[i]=='\t')
						{
							write << endl;
						}
						else
						{
							write << line[i];
						}
					}			
				}
				write << endl;
			}
			l++;
		}
		
	}
	else
	{
		cout << "No file there (or something like that)" << endl;
	}
}

void order_res(std::string res, std::string rep, std::string seq)
{
	//when the results are not ordered the way the sequences given are ordered (eg by comet), rearranges the results wrt the sequence order.
	ofstream write(rep);
	ifstream data(res);
	std::vector<rappas::io::fasta> sequences= rappas::io::read_fasta(seq);
	if(data)
	{
		std::vector<std::string> tab (0);
		std::string line="";
		int j=0;
		while(getline(data,line))
		{
			tab.push_back(line);
		}
		for(int i=0;i<sequences.size(); i++)
		{
			j=0;
			while(j<tab.size() && (tab[j]).substr(1) !=(sequences[i]).header())
			{
				j+=2;
			}
			//cout << (sequences[i]).header() << " is " << tab[j-2] << endl;
			write << tab[j] << endl << tab[j+1] << endl;
		}
	}
	else
	{
		cout << "No file there (or something like that)" << endl;
	}
}

void prune(std::string res, std::string seq, int shift)
{
	//basically removes the first and last "shift" bases in a result file.
	//tends to become useless as variable window size replaces fixed window size
	std::vector<rappas::io::fasta> test= rappas::io::read_fasta(seq);
	std::vector<std::vector<std::string>> tmp=readRes(res);
	int j=0;
	for(int i=0; i<tmp.size();i++)
	{
		cout << ">" << (test[i]).header() << endl;
		cout << shift << ",";
		j=0;
		while(stoi(tmp[i][j])<shift)
		{
			j+=2;
		}
		while(stoi(tmp[i][j])<stoi((tmp[i]).back())-shift)
		{
			cout << tmp[i][j-1] << "," << tmp[i][j] << ",";
			j+=2;
		}
		cout << tmp[i][j-1] << "," << stoi((tmp[i]).back())+1-shift << endl;
	}
}

std::vector<int> readRead(std::string name)
{
	//stores relevant infos from the name of a nanosim output in a table of ints
	//0 -> number of the query
	//1 -> start of alignment
	//2 -> 1 if F, 0 if R
	//3 -> first junk
	//4 -> nbr of sites noised
	//5 -> last junk
	std::vector<int> res (0);
	std::string tmp="";
	int i=0;
	int pos=0;
	while(name[i] !='-')
	{
		i++;
	}
	i++;
	while(i<name.size()+1)
	{
		if(name[i]=='_' || i==name.size())
		{
			if(pos <2 || pos>4)
			{
				res.push_back(stoi(tmp));
			}
			if(pos==4)
			{
				res.push_back(tmp=="F");
			}
			pos++;
			tmp="";
		}
		else
		{
			tmp+=name[i];
		}
		i++;
	}
	return res;
}

std::vector<int> pickN(std::vector<rappas::io::fasta> *list, int max)
{
	//picks one noisy read for each query from 1 to max.
	std::vector<int> res;
	std::vector<int> infos (6);
	int check=0;
	int j=0;
	for(int i=0; i<max; i++)
	{
		check=0;
		j=0;
		while(check==0 && j<(*list).size())
		{
			infos=readRead((((*list)[j]).header()).data());
			if(infos[0]== i+1 && infos[2]==1)
			{
				check=1;
			}
			else
			{
				j++;
			}
		}
		if(check==1)
		{
			res.push_back(j);
			//cout << "* " << i+1 << endl;
		}
		else
		{
			res.push_back(-1);
			cout << "* " << i+1 << endl;
		}
	}	
	return res;
}

std::vector<std::string> lineErr(std::string line)
{
	std::vector<std::string> res (0);
	int i=0;
	std::string tmp="";
	while(i<line.length())
	{
		if(line[i]=='	')
		{
			res.push_back(tmp);
			tmp="";
		}
		else
		{
			tmp+=line[i];
		}
		i++;
	}
	return res;
}

void readErr(std::string res, std::string rep)
{
	ofstream write(rep);
	ifstream data(res);
	if(data)
	{
		int i=0;
		int place=0;
		std::string line="";
		std::string nm="";
		std::vector<std::string> info (4);
		while(getline(data,line))
		{
			info=lineErr(line);
			if(info[0] !=nm && info[0] !="Seq_name")
			{
				cout << info[0] << endl;
				write << endl << ">" << info[0] << endl;
				nm=info[0];
			}
			if(info[2]=="ins")
			{
				write << info[1] << "," << info[3] << ",";
			}
			if(info[2]=="del")
			{
				write << info[1] << ",-" << info[3] << ",";
			}
		}
	}
	else
	{
		cout << "No file there (or something like that)" << endl;
	}
}

void sortErr(std::string res, std::string rep, std::vector<rappas::io::fasta>* queries)
{
	ofstream write(rep);
	std::string nm="";
	std::string line="";
	int check=0;
	int seq=0;
	int pos=0;
	for(int i=0; i<(*queries).size(); i++)
	{
		check=0;
		nm=((*queries)[i].header()).data();
		cout << i+1 << ": " << nm << endl;
		seq=readRead(nm)[0];
		pos=readRead(nm)[1];
		ifstream full(res);
		while(getline(full, line) && check<2)
		{
			if(check==1)
			{
				write << line << endl;
				check++;
			}
			if(line[0]=='>')
			{
				if(readRead(line)[0]==seq && readRead(line)[1]==pos)
				{
					write << ">" << nm << endl;
					check++;
				}
			}
		}
		full.close();
	}
}

void infoReads(std::string clean, std::string err, std::string res, std::vector<rappas::io::fasta>* reads)
{
	std::vector<std::vector<std::string>> info=readRes(clean);
	std::vector<std::vector<std::string>> move=readRes(err);
	std::vector<string> nm= readNm(err);
	ofstream write(res);
	int init=0;
	int bp=0;
	int i=0;
	int j=0;
	int shift=0;
	std::string seq="";
	std::vector<int> read (6);
	int finit=0;
	for(int k=0; k<nm.size(); k++)
	{
		cout << "query " << k+1 << endl;
		i=2;
		j=(move[k]).size()-3;
		bp=0;
		shift=0;
		read=readRead(nm[k]);
		init=read[1]-read[3];
		seq=(((*reads)[k]).sequence()).data();
		finit=seq.size();
		while(stoi(info[k][i])<init)
		{
			i+=2;
		}
		write << nm[k] << endl;
		write << "0,";
		while(i<(info[k]).size())
		{
			bp=stoi(info[k][i])-init;
			while(stoi(move[k][j])<bp && j<(move[k]).size()-2)
			{
				shift+=stoi(move[k][j+1]);
				j=j-2;
				if(j<0)
				{
					j=(move[k]).size()-2;
				}
			}
			if(bp+shift<finit && i<(info[k]).size()-2)
			{
				write << info[k][i-1] << "," << bp+shift << ",";
			}
			else
			{
				write << info[k][i-1] << "," << finit << endl;
				i=(info[k]).size()-2;
			}
			i+=2;
		}
	}
}

std::string round(std::string val, int arr)
{
	std::string res="";
	int j=0;
	int i=0;
	while(val[j] !='.' && j<val.length()-1)
	{
		res+=val[j];
		j++;
	}
	while(i<arr+1 && j+i<val.length()-1)
	{
		res+=val[j+i];
		i++;
	}
	res+=' ';
	return res;
}

void shrinkTable(std::string tab, std::string res, std::vector<std::string> kp, int arr, int col)
{
	//don't ask...
	ofstream write(res);
	ifstream data(tab);
	if(data)
	{
		std::string line="";
		int check=0;
		std::string tmp="";
		std::vector<std::string> ltmp (0);
		int j=0;
		while(getline(data,line))
		{
			ltmp.clear();
			if(j<46 || line[0]=='\\' || line.length()==0 || line[0]=='P')
			{
				write << line << endl;
			}
			else
			{
				for(int i=0; i<line.length(); i++)
				{
					if(line[i]!='&')
					{
						tmp+=line[i];
					}
					else
					{
						ltmp.push_back(tmp);
						tmp="";
					}
				}
				ltmp.push_back(tmp);
				tmp="";
				if(ltmp[0]==kp[check])
				{
					check++;
					write << ltmp[0];
					for(int k=1; k<col+1; k++)
					{
						write << "&" << round(ltmp[k], arr);
					}
					write << "\\\\" << endl;
				}
				//write << "line here" << endl;
			}
			j++;
		}
	}
	else
	{
		cout << "No file there (or something like that)" << endl;
	}
}

/*
 * output.cpp
 *
 *  Created on: 11 feb. 2020
 *      Author: scholz
 */


#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include<cmath>
#include "output.h"

using namespace std;
//deals with output and writes files.

std::vector<std::string> readNm(std::string res)
{
	//gets list of names from a result file
	std::vector<std::string> read(0);
	ifstream data(res);
	if(data)
	{
		std::string line="";
		while(getline(data,line))
		{
			if(line[0] =='>')
			{
				read.push_back(line);
			}
		}
		data.close();
	}
	else
	{
		cout << "No file found at " << res << endl;
	}
	return read;
}

std::string fileName(std::string add)
{
	//Extracts file name from its address.
	std::string res=add;
	int p=add.length()-1;
	int e=0;
	while(add[p] != '/' && p>0)
	{
		if(add[p]=='.')
		{
			e=p;
		}
		p--;
	}
	if(p>0)
	{
		res=add.substr(p+1, e-p-1);
	}
	return res;
}

void printHead(std::string qfile, char dbtype, double theta, int ws, int cflag, int kflag, std::ofstream* writef)
{
	//prints header of a output file, comtaining informations on the parameters used.
	(*writef) << "# SHERPAS output for queries \"" << qfile << "\"" << endl;
	(*writef) << "#" << endl;
	(*writef) << "# method: " << dbtype << "; threshold: " << theta << "; window: " << ws << endl;
	if(cflag==1)
	{
		(*writef) << "#" << endl << "# circular genome" << endl;
	}
	if(kflag==1)
	{
		(*writef) << "#" << endl << "# no post-treatment of N/A-segments" << endl;
	}
	(*writef) << "#" << endl << endl; 
}

double lRatio(std::vector<Arc*> result, int i, int k)
{
	//likelihood ratio
	int s=result.size();
	double q=0;
	for(int j=0; j<s;j++)
	{
		q+=(pow(10,((*result[j]).getScore()-(*result[i]).getScore())/k))*((*result[j]).getScore() !=0);
	}
	return 1/q;
}

std::vector<std::string> printChange(std::vector<std::vector<Arc*>> result, int shift, std::vector<std::string> ref, double thr, int k, char m)
{
	//from a list of result vectors, the segments with breakpoints and values between (does not print anything, despite the name).
	std::vector<std::string> res(0);
	int s=result.size();
	int rs=ref.size();
	int aref=rs;
	res.push_back("0");
	ref.push_back("N/A");
	double rat=0;
	for(int i=0; i<s; i++)
	{
		if(m=='F')
		{
			rat=pow(10,((*result[i][0]).getScore()-(*result[i][1]).getScore())/k);
		}
		else
		{
			rat=lRatio(result[i],0,k);
		}
		if(ref[(*result[i][0]).getPlace()] != ref[(*result[i][1]).getPlace()] && (i==0 || aref<rs) && rat<thr)
		{
			if(i>0)
			{
				res.push_back(to_string(i+shift));
			}
			aref=rs;
			res.push_back("N/A");
		}
		else if(ref[(*result[i][0]).getPlace()] != ref[aref] && (ref[(*result[i][0]).getPlace()] == ref[(*result[i][1]).getPlace()] || rat>=thr))
		{
			if(i>0)
			{
				res.push_back(to_string(i+shift));
			}
			aref=((*result[i][0]).getPlace());
			if(isTop(ref[aref]))
			{
				res.push_back("N/A");
			}
			else
			{
				res.push_back(ref[aref]);
			}
		}
	}
	res.push_back(to_string(s+2*shift));
	return res;
}

void mergeNA(std::vector<std::vector<std::string>> read, std::string seq, int circ, int lin, int keep, std::ofstream* write)
{
	//Remove "N/A" section when both adjacent sections are the same.
	std::vector<string> nm= readNm(seq);
	int l=0;
	for(int i=0; i<read.size(); i++)
	{
		(*write) << nm[i] << endl;
		if(keep==0)
		{
			l=3;
			while(l<(read[i]).size()-2)
			{
				if(read[i][l]=="N/A" && read[i][l-2]==read[i][l+2])
				{
					for(int rec=0; rec<4; rec++)
					{
						(read[i]).erase((read[i]).begin()+l-1);
					}
				}
				else if(read[i][l]=="N/A" && read[i][l+2]=="N/A")
				{
					for(int rec=0; rec<2; rec++)
					{
						(read[i]).erase((read[i]).begin()+l);
					}
				}
				else
				{
					l+=2;
				}
			}
			if(circ==1)
			{
				if(read[i][1]=="N/A" && read[i][(read[i]).size()-2]=="N/A" && read[i]	[3]==read[i][(read[i]).size()-4])
				{
					(read[i]).erase((read[i]).end()-2);
					(read[i]).erase((read[i]).end()-2);
					(read[i]).erase((read[i]).begin()+1);
					(read[i]).erase((read[i]).begin()+1);
				}
			}
		}
		if(lin==1)
		{
			for(int j=0; j<(read[i]).size()-1; j++)
			{
				(*write) << read[i][j] << ",";
			}
			(*write) << (read[i]).back() << endl;
		}
		else
		{
			(*write) << "1\t";
			for(int j=1; j<(read[i]).size()-2; j+=2)
			{
				(*write) << stoi(read[i][j+1])-1 << "\t";
				(*write) << read[i][j] << endl << read[i][j+1] << "\t";
			}
			(*write) << (read[i]).back()  << "\t" << read[i][(read[i]).size()-2] << endl << endl;
		}
	}
}

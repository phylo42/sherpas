/*
 * output.cpp
 *
 *  Created on: 11 feb. 2020
 *      Author: scholz
 */


#include<iostream>
#include<vector>
#include <string>
#include <boost/filesystem.hpp>
#include "output.h"

using namespace std;
namespace fs = boost::filesystem;

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

void outdir(std::string add)
{
	int i=add.length()-1;
	while(add[i] != '/')
	{
		i--;
	}
	cout << add.substr(0,i+1) << endl;
	if(fs::exists(add.substr(0,i+1))==0)
	{
		fs::create_directory(add.substr(0,i+1));
	}
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

std::vector<std::string> printChange(std::vector<std::vector<Arc*>> result, int shift, std::vector<std::string> ref, double thr, std::vector<double> rat, char m)
{
	//from a list of result vectors, the segments with breakpoints and values between (does not print anything, despite the name).
	std::vector<std::string> res(0);
	int s=result.size();
	int rs=ref.size();
	int aref=rs;
	res.push_back("0");
	ref.push_back("N/A");
	for(int i=0; i<s; i++)
	{
		if((ref[(*result[i][0]).getPlace()] != ref[(*result[i][1]).getPlace()] || m=='R') && (i==0 || aref<rs) && rat[i]<thr)
		{
			if(i>0)
			{
				res.push_back(to_string(i+shift));
			}
			aref=rs;
			res.push_back("N/A");
		}
		else if((ref[(*result[i][0]).getPlace()] != ref[aref]) && ((ref[(*result[i][0]).getPlace()] == ref[(*result[i][1]).getPlace()] && m=='F') || rat[i]>=thr))
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
		delete result[i][0];
		delete result[i][1];
	}
	res.push_back(to_string(s+2*shift));
	return res;
}

void mergeNA(std::vector<std::string> read, int circ, int lin, int keep, std::ofstream* write)
{
	//Remove "N/A" section when both adjacent sections are the same.
	int l=0;
	if(keep==0)
	{
		l=3;
		while(l<(read).size()-2)
		{
			if(read[l]=="N/A" && read[l-2]==read[l+2])
			{
				for(int rec=0; rec<4; rec++)
				{
					(read).erase((read).begin()+l-1);
				}
			}
			else if(read[l]=="N/A" && read[l+2]=="N/A")
			{
				for(int rec=0; rec<2; rec++)
				{
					(read).erase((read).begin()+l);
				}
			}
			else
			{
				l+=2;
			}
		}
			if(circ==1)
		{
			if(read[1]=="N/A" && read[(read).size()-2]=="N/A" && read[3]==read[(read).size()-4])
			{
				(read).erase((read).end()-2);
				(read).erase((read).end()-2);
				(read).erase((read).begin()+1);
				(read).erase((read).begin()+1);
			}
		}
	}
	if(lin==1)
	{
		for(int j=0; j<(read).size()-1; j++)
		{
			(*write) << read[j] << ",";
		}
		(*write) << (read).back() << endl;
	}
	else
	{
		(*write) << "1\t";
		for(int j=1; j<(read).size()-2; j+=2)
		{
			(*write) << stoi(read[j+1])-1 << "\t";
			(*write) << read[j] << endl << read[j+1] << "\t";
		}
		(*write) << (read).back()  << "\t" << read[(read).size()-2] << endl << endl;
	}
}

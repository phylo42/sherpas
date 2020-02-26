/*
 * groups.cpp
 *
 *  Created on: 22 may 2019
 *      Author: scholz
 */


#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include<cmath>
#include "groups.h"

using namespace std;
//everything needed to put the notion of groups at the heart of the output.
//does not affect the way the main algorithm works.

std::vector<std::string> readGrFile(std::string csv)
{
	//reads the .csv file with groups information
	std::vector<std::string> tab (0);
	ifstream data(csv);
	if(data)
	{
		std::string line="";
		int check=0;
		int p=0;
		while(getline(data,line))
		{
			p=0;
			while(line[p] !=',')
			{
					p++;
			}
			tab.push_back(line.substr(0,p));
			tab.push_back(line.substr(p+1));
		}
		data.close();
	}
	else
	{
		cout << "No file found at " << csv << endl;
	}
	return tab;
}

std::string groupFromTab(std::string nm, std::vector<std::string> tab)
{
	//uses the table read from the .csv file to get the group of reference sequence.
	std::string res="";
	int j=0;
	int check=0;
	while(j<tab.size() && check==0)
	{
		if(tab[j]==nm)
		{
			check=1;
			res=tab[j+1];
		}
		j+=2;
	}
	if(check==0)
	{
		cout << "Leaf \"" << nm << "\" not found in the reference tree." << endl;
		res="N/A";
	}
	return res;
}

void getArcRef(core::phylo_tree& tree, std::vector<std::string>* ref, std::string gfile)
{
	//maps each branch to the corresponding group in *ref
	if((*ref).size()>0)
	{
		cout << "Warning: vectors to be filled-up are nonempty" << endl;
	}
	int i=-1;
	std::string group="";
	std::vector<std::string> tab=readGrFile(gfile);
	for (const auto& node : tree)
   	{
		i=node.get_postorder_id();
		if((node.get_children()).size()==0)
		{
			group=groupFromTab(node.get_label(), tab);
		}
		else
		{
			if((*ref)[(*node.get_children()[0]).get_postorder_id()] == (*ref)[(*node.get_children()[1]).get_postorder_id()])
			{
				group=(*ref)[(*node.get_children()[1]).get_postorder_id()];
			}
			else
			{
				group=(*ref)[(*node.get_children()[0]).get_postorder_id()]+ "*" + (*ref)[(*node.get_children()[1]).get_postorder_id()];
			}
		}
		(*ref).push_back(group);
    	}
}

std::vector<std::string> listGroups(std::vector<std::string> ref)
{
	std::vector<std::string> res (0);
	int s=0;
	int j=0;
	for(int i=0; i<ref.size(); i++)
	{
		if(ref[i][1] != '*' && ref[i][2] != '*' && ref[i][3] != '*')
		{
			j=0;
			while(j<s)
			{
				if(res[j] != ref[i])
				{
					j++;
				}
				else
				{
					j=s+1;
				}
			}
			if(j==s)
			{
				res.push_back(ref[i]);
				s++;
			}
		}
	}
	return res;
}

int isTop(std::string gr)
{
	//check from its label if the branch of a tree is a "top branch" (ie not in a clade).
	int res=0;
	int i=0;
	while(i<gr.length() && res==0)
	{
		if(gr[i]=='*')
		{
			res=1;
		}
		i++;
	}
	return res;
}

void onlyRoot(const core::phylo_kmer_db& db, core::phylo_kmer_db *db2, std::vector<std::string>* ref)
{
	//shrinks a database by keeping only the score of the root of each group
	double thr=core::score_threshold(db.omega(), db.kmer_size());
	for(const auto& [key, entries] : db)
	{
		for (const auto& [branch, score] : entries)
		{
			if((*ref)[branch+1] != (*ref)[branch] && !isTop((*ref)[branch]))
			{
				(*db2).insert(key, {branch, score});
			}
		}
	}
}


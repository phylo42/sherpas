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

int max(std::vector<int> list)
{
	int res=-1;
	for(int i=0;i<list.size(); i++)
	{
		if(list[i]>res)
		{
			res=list[i];
		}
	}
	return res;
}

void getId(std::string align, std::string res) //only once; depends on the sequences name; current version works for names starting with "Ref.[group_id]."
{
	ifstream read(align);
	int i=0;
	int c=0;
	if(read)
	{
		ofstream write(res+"leaves-summary.txt");
		std::string line;
		while(getline(read,line))
		{
			if(line[0]=='>')
			{
				while(line[i] !='.' || c<1)
				{
					if(line[i]=='.')
					{
						c++;
					}
					i++;
				}
				write << line.substr(1) << " " << line.substr(5,i-5) << endl;
				i=0;
				c=0;
			}
		}
	}
	else
	{
		cout << "No file there (or something like that)" << endl;
	}
}

void recordGroups(std::string doc,std::vector<int>* leafMap, std::vector<std::string>* ref)
{
	if((*leafMap).size() > 0 || (*ref).size() > 0)
	{
		cout << "careful; references nonempty" << endl;
	}
	ifstream read(doc);
	if(read)
	{
		//(*leafMap).push_back(0);
		(*ref).push_back("out");
		std::string line;
		std::string group;
		int i=0;
		int g=0;
		while(getline(read,line))
		{
			while((line[i] !=' ' || line[i+1]==' ') && i<line.length())
			{
				i++;
			}
			group=line.substr(i+1);
			//cout << group << endl;
			if(group != (*ref).back())
			{
				(*ref).push_back(group);
				g++;
			}
			(*leafMap).push_back(g);
			i=0;
		}
	}
	else
	{
		cout << "No file there (or something like that)" << endl;
	}
}

std::vector<int> makeGroups(core::phylo_tree& tree, std::vector<int> leafMap)
{
	int tree_size=tree.get_node_count();
	std::vector<int> groups(tree_size,-1);
	int l=0;
	int i=0;
	int m=max(leafMap)+1;
	for (const auto& node : tree)
   	{
		i=node.get_id();
        	if(node.get_label()!="")
		{
			groups[i]=leafMap[l];
			l++;
		}
		else
		{
			if(groups[(*node.get_children()[0]).get_id()] == groups[(*node.get_children()[1]).get_id()])
			{
				groups[i]=groups[(*node.get_children()[1]).get_id()];
			}
			else
			{
				groups[i]=m;
				//m++; //active -> top arcs grouped individually ; commented out -> top arcs grouped all together
			}
		}
    	}
	return groups;
}

core::phylo_kmer_db GroupDb(const core::phylo_kmer_db& db, std::vector<int> groups)
{
	core::phylo_kmer_db db2 {groups.size()};
	double thr=core::score_threshold(db.kmer_size());
	int g=max(groups)+1;
	std::vector<double> bg(g, thr);
	for(int key=0; key<db.size(); key++)
		if (auto entries = db.search(key); entries)
		{
			for (const auto& [branch, score] : *entries)
			{
				if(score>bg[groups[branch]])
				{
					bg[groups[branch]]=score;
				}
			}
			for(int i=0; i<g; i++)
			{
				if(bg[i] > thr)
					{
						db2.put(key, i, bg[i]);
						bg[i]=thr;
					}
			}
		}
		else
		{
		std::cout << "Key " << key << " not found.\n";
		}
	return db2;
}




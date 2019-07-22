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

std::string groupFromName(std::string leaf)
{
	int i=0;
	while(leaf[i] !='.' && i<leaf.length())
	{
		i++;
	}
	return leaf.substr(0,i);
}

void getArcRef(core::phylo_tree& tree, std::vector<std::string>* ref)
{
	if((*ref).size()>0)
	{
		cout << "Warning: vectors to be filled-up are nonempty" << endl;
	}
	int i=-1;
	std::string group="";
	for (const auto& node : tree)
   	{
		i=node.get_postorder_id();
		if((node.get_children()).size()==0)
		{
			group=groupFromName(node.get_label());
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
		cout << i << ": " << (*ref)[i] << endl;
    	}
}

void getDb2Ref(core::phylo_tree& tree, std::vector<std::string>* ref, std::vector<int>* group_id)
{
	if((*ref).size()>0 ||(*group_id).size()>0 )
	{
		cout << "Warning: vectors to be filled-up are nonempty" << endl;
	}
	int i=-1;
	std::string group="";
	int l=-1;
	std::string gtemp="";
	for (const auto& node : tree)
   	{
		i=node.get_postorder_id();
        	if((node.get_children()).size()==0)
		{
			if(groupFromName(node.get_label()) != gtemp)
			{
				l++;
				gtemp=groupFromName(node.get_label());
				(*ref).push_back(gtemp);
				//cout << l << ": " << (*ref)[l] << endl;
			}
		}
		else
		{
			if((*group_id)[(*node.get_children()[0]).get_postorder_id()] != (*group_id)[(*node.get_children()[1]).get_postorder_id()])
			{
				l++;
				(*ref).push_back((*ref)[(*group_id)[(*node.get_children()[0]).get_postorder_id()]]+ "*" + (*ref)[(*group_id)[(*node.get_children()[1]).get_postorder_id()]]);
				//cout << l << ": " << (*ref)[l] << endl;
			}
		}
		(*group_id).push_back(l);
		//cout << i << ": " << (*group_id)[i] << endl;
    	}
}

core::phylo_kmer_db GroupDb(const core::phylo_kmer_db& db, std::vector<int> groups)
{
	size_t k=db.kmer_size();
	core::phylo_kmer_db db2 {k, std::string{db.tree()} };
	double thr=core::score_threshold(k);
	int g=max(groups)+1;
	//cout << g << endl;
	std::vector<double> bg(g, thr);
	for(const auto& [key, entries] : db)
	{
		//cout << key << endl;
		for (const auto& [branch, score] : entries)
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
					db2.insert(key, {i, bg[i]});
					bg[i]=thr;
				}
		}
	}
	return db2;
}




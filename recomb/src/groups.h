/*
 * groups.h
 *
 *  Created on: 22 may 2019
 *      Author: scholz
 */

#ifndef GROUPS_H_
#define GROUPS_H_
#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include <core/phylo_kmer_db.h>
#include <core/phylo_tree.h>
//#include "query.h"

int max(std::vector<int> list);

void recordGroups(std::string doc,std::vector<int>* leafMap, std::vector<std::string>* ref);

std::vector<int> makeGroups(core::phylo_tree& tree, std::vector<int> leafMap);

core::phylo_kmer_db GroupDb(const core::phylo_kmer_db& db, std::vector<int> groups);


#endif /* GROUPS_H_ */

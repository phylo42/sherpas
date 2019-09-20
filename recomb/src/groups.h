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

std::string groupFromName(std::string leaf);

void getArcRef(core::phylo_tree& tree, std::vector<std::string>* ref);

void getDb2Ref(core::phylo_tree& tree, std::vector<std::string>* ref, std::vector<int>* group_id);

void GroupDb(const core::phylo_kmer_db& db, core::phylo_kmer_db *db2, std::vector<int> groups);

void rmTop(const core::phylo_kmer_db& db, core::phylo_kmer_db *db2, std::vector<std::string>* ref);

void onlyRoot(const core::phylo_kmer_db& db, core::phylo_kmer_db *db2, std::vector<std::string>* ref);

#endif /* GROUPS_H_ */

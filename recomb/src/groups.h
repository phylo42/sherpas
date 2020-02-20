/*
 * groups.h
 *
 *  Created on: 22 may 2019
 *      Author: scholz
 */

#ifndef GROUPS_H_
#define GROUPS_H_
#include <iostream>
#include <fstream>
#include<vector>
#include <string>
#include <core/phylo_kmer_db.h>
#include <core/phylo_tree.h>

std::vector<std::string> readGrFile(std::string csv);

std::string groupFromTab(std::string nm, std::vector<std::string> tab);

void getArcRef(core::phylo_tree& tree, std::vector<std::string>* ref, std::string gfile);

std::vector<std::string> listGroups(std::vector<std::string> ref);

int isTop(std::string gr);

void onlyRoot(const core::phylo_kmer_db& db, core::phylo_kmer_db *db2, std::vector<std::string>* ref);

#endif /* GROUPS_H_ */

/*
 * write.h
 *
 *  Created on: 6 nov. 2018
 *      Author: scholz
 */

#ifndef WRITE_H_
#define WRITE_H_
#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include "query.h"

int findArc(int a, std::string t);

double lengthArc(int pos, std::string t);

std::string subtree(int pos, std::string t);

double distNode(int pos, std::string t);

double sumPaths(std::string t);

int nbLeaves(std::string t);

double lengthPendant(int pos, std::string t);

void resInit(std::vector<Arc*> *res);

void windInit(std::vector<std::vector<Arc*>> *wind);

//void jplace(std::vector<std::string> q, std::string t, std::vector<std::vector<Kmer*>> list, std::vector<Arc>* branches, std::vector<std::string> infos, std::vector<Arc*> *res);

void SciPlot(int n, std::string_view query, std::vector<Arc> branches, std::vector<std::vector<Arc*>> res);


#endif /* WRITE_H_ */

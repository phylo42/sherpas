/*
 * query.h
 *
 *  Created on: 23 oct. 2018
 *      Author: scholz
 */

#ifndef QUERY_H_
#define QUERY_H_
#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include "arcs.h"
#include "heap.h"
#include <core/phylo_kmer_db.h>
#include <core/kmer_iterator.h>

std::vector<core::phylo_kmer_db::key_type> encode_string_views(std::string_view long_read, const size_t kmer_size);

void search(const core::phylo_kmer_db& db, core::phylo_kmer_db::key_type key);

void readQuery(std::vector<core::phylo_kmer_db::key_type> codes, const core::phylo_kmer_db& db, std::vector<Arc>* branches, Htree *H);

void printScore(std::vector<Arc*> result, std::vector<std::string> ref);

void printChange(std::vector<std::vector<Arc*>> result, int shift, std::vector<std::string> ref);

double lRatio(std::vector<Arc*> result, int i);

void slidingWindow(std::vector<core::phylo_kmer_db::key_type> codes, int sw, int m, const core::phylo_kmer_db& db, std::vector<Arc>* branches, std::vector<std::vector<Arc*>> *res);

#endif /* QUERY_H_ */

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
#include "heap.h"
#include <core/phylo_kmer_db.h>
#include <core/kmer_iterator.h>
#include <utils/io/fasta.h>

int nucl(char c);

std::vector<rappas::io::fasta> gapRm(std::vector<rappas::io::fasta> *sequences);

void make_circu(std::string res, std::string rep, int p);

void windInit(std::vector<std::vector<Arc*>> *wind);

std::vector<std::vector<core::phylo_kmer_db::key_type>> encode_ambiguous_string(std::string_view long_read, const size_t kmer_size);

std::vector<std::vector<core::phylo_kmer_db::key_type>> fixCcodes(std::vector<std::vector<core::phylo_kmer_db::key_type>> ccodes, std::string_view long_read, const size_t kmer_size);

void readQuery(std::vector<std::vector<core::phylo_kmer_db::key_type>> codes, const core::phylo_kmer_db& db, std::vector<Arc>* branches, Htree *H);

void addKmer(std::vector<core::phylo_kmer_db::key_type> keys, Htree* H, const core::phylo_kmer_db& db, std::vector<Arc>* branches, double thr, int sw, int k, int move);

void rmKmer(std::vector<core::phylo_kmer_db::key_type> keys, Htree* H, const core::phylo_kmer_db& db, std::vector<Arc>* branches, double thr, int move);

void slidingVarWindow(std::vector<std::vector<core::phylo_kmer_db::key_type>> codes, int wi, int sw, int m, const core::phylo_kmer_db& db, std::vector<Arc>* branches, std::vector<std::vector<Arc*>> *res);

#endif /* QUERY_H_ */

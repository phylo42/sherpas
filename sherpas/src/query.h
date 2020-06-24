/*
 * query.h
 *
 *  Created on: 23 oct. 2018
 *      Author: scholz
 */

#ifndef QUERY_H_
#define QUERY_H_
#include<iostream>
#include <fstream>
#include <vector>
#include <string>
#include "heap.h"
#include <xpas/phylo_kmer_db.h>
#include <xpas/kmer_iterator.h>
#include <utils/io/fasta.h>

// CODE REVIEW: no need to declare here. This function is local to query.cpp
//int nucl(char c);

// CODE REVIEW: passing by pointer is C-way, consider using a reference. The most important difference is
// that it can not be null.
// CODE REVIEW: Actually xpas::io::read_fasta already removes the gaps. I am quite sure I told you when I've
// changed this. Anyway, this method is not needed anymore.
//std::vector<xpas::io::fasta> gapRm(std::vector<xpas::io::fasta> *sequences);

// CODE REVIEW: consider passing by reference
void make_circu(const std::string& res, const std::string& rep, int p);

// CODE REVIEW: consider passing by reference
void windInit(std::vector<std::vector<Arc*>>& wind);


// CODE REVIEW: I wanna make your life simpler.
// What type is key_codes depends on which version of xpas::to_kmers you use though...
// Here I assume that we always use xpas::to_kmers<xpas::one_ambiguity_policy>.
// That can be changed in future.
namespace sherpas
{
    using key_codes = std::vector<xpas::phylo_kmer_db::key_type>;
}

// CODE REVIEW: const size_t here does not do anything
// CODE REVIEW: type alias
std::vector<sherpas::key_codes> encode_ambiguous_string(std::string_view long_read, size_t kmer_size);


// CODE REVIEW: no need to declare here
//std::vector<sherpas::key_codes> fixCcodes(const std::vector<sherpas::key_codes>& ccodes,
//                                          std::string_view long_read,
//                                          size_t kmer_size);

// CODE REVIEW: no need to declare here
//void readQuery(std::vector<std::vector<xpas::phylo_kmer_db::key_type>> codes, const xpas::phylo_kmer_db& db, std::vector<Arc>* branches, Htree *H);

// CODE REVIEW: consider using a reference to const instead of passing by value
// CODE REVIEW: consider using a reference instead of a pointer
void addKmer(const std::vector<xpas::phylo_kmer_db::key_type>& keys, Htree& H,
             const xpas::phylo_kmer_db& db, std::vector<Arc>& branches,
             double thr, int sw, int k, int move);

// CODE REVIEW: consider passing by reference
void rmKmer(const std::vector<xpas::phylo_kmer_db::key_type>& keys, Htree& H,
            const xpas::phylo_kmer_db& db,
            std::vector<Arc>& branches, double thr, int move);

// CODE REVIEW: consider passing by reference
void slidingVarWindow(const std::vector<std::vector<xpas::phylo_kmer_db::key_type>>& codes, int wi, int sw, int m,
                      const xpas::phylo_kmer_db& db,
                      std::vector<Arc>& branches,
                      std::vector<std::vector<Arc*>>& res,
                      std::vector<double>& rat, char met);

#endif /* QUERY_H_ */

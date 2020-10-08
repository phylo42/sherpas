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


void make_circu(const std::string& res, const std::string& rep, int p);

namespace sherpas
{
    using key_codes = std::vector<xpas::phylo_kmer_db::key_type>;
}

std::vector<sherpas::key_codes> encode_ambiguous_string(std::string_view long_read, size_t kmer_size);


void addKmer(const std::vector<xpas::phylo_kmer_db::key_type>& keys, Htree& H,
             const xpas::phylo_kmer_db& db, std::vector<Arc>& branches,
             double thr, int sw, int k, int move);

void rmKmer(const std::vector<xpas::phylo_kmer_db::key_type>& keys, Htree& H,
            const xpas::phylo_kmer_db& db,
            std::vector<Arc>& branches, double thr, int move);

void windOut(Htree& H, size_t k, int m, char met, std::vector<std::vector<Arc*>>& res, std::vector<double>& rat, Arc* neutral);

void slidingVarWindow(const std::vector<std::vector<xpas::phylo_kmer_db::key_type>>& codes, int wi, int sw, int m,
                      const xpas::phylo_kmer_db& db,
                      std::vector<Arc>& branches,
                      std::vector<std::vector<Arc*>>& res,
                      std::vector<double>& rat, char met);

#endif /* QUERY_H_ */

/*
 * output.h
 *
 *  Created on: 11 feb. 2020
 *      Author: scholz
 */

#ifndef WRITE_H_
#define WRITE_H_
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "arcs.h"
#include "groups.h"
#include <xpas/phylo_kmer_db.h>
#include <xpas/kmer_iterator.h>
#include <utils/io/fasta.h>

std::vector<std::string> readNm(std::string res);

std::string fileName(std::string add);

void printHead(std::string qfile, char dbtype, double theta, int ws, int cflag, int kflag, std::ofstream* writef);

std::vector<std::string> printChange(std::vector<std::vector<Arc*>> result, int shift, std::vector<std::string> ref, double thr, std::vector<double> rat, char m);

void mergeNA(std::vector<std::string> read, int circ, int lin, int keep, std::ofstream* write);

#endif /* OUTPUT_H_ */

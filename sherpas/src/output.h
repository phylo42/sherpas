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
#include "xpas/fasta.h"

std::vector<std::string> readNm(std::string res);

// Creates a directory if it does not exists
void createDirectory(const std::string& filename);

// Check if input file exists. If not, throws an exception
void checkFileExists(const std::string& filename);

// Extracts the file name from its full name
std::string getFileName(const std::string& filename);

void printHead(std::string qfile, char dbtype, double theta, int ws, int cflag, int kflag, std::ofstream* writef);

std::vector<std::string> printChange(std::vector<std::vector<Arc*>> result, int shift, std::vector<std::string> ref, double thr, std::vector<double> rat, char m);

void mergeNA(std::vector<std::string> read, int circ, int lin, int keep, std::ofstream* write, std::string seq_header);

#endif /* OUTPUT_H_ */

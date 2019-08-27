/*
 * input.h
 *
 *  Created on: 05 jun. 2019
 *      Author: scholz
 */

#ifndef INPUT_H_
#define INPUT_H_
#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include <core/phylo_tree.h>
#include <utils/io/fasta.h>
#include "groups.h"

int nucl(char c);

void incr(std::vector<int> *list, int a);

void writeFasta(std::vector<rappas::io::fasta> *sequences, std::string res);

std::vector<rappas::io::fasta> gapRm(std::vector<rappas::io::fasta> *sequences);

std::vector<int> splitGroups(std::vector<rappas::io::fasta> *sequences);

void writeInfo(std::vector<std::string> *summary, std::string res);

void split_seq(std::vector<rappas::io::fasta> *sequences, std::vector<int> bp, int seq);

void new_q(std::vector<rappas::io::fasta> *sequences, std::vector<rappas::io::fasta> *nq, std::vector<int> *bp, std::vector<int> *seq, int ref);

std::string record(std::vector<rappas::io::fasta> *sequences, std::vector<int> bp, std::vector<int> seq);

void random_q(std::vector<rappas::io::fasta> *sequences, std::vector<rappas::io::fasta> *nq, int n, std::vector<std::string> *summary, int b_max, int p_max);


#endif /* INPUT_H_ */

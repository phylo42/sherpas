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
#include "query.h"
#include "output.h"

void incr(std::vector<int> *list, int a);

int bpDist(std::vector<int> bp, int max);

std::string getAccess(std::string leaf);

void writeFasta(std::vector<rappas::io::fasta> *sequences, std::string res);

std::vector<int> splitGroups(std::vector<rappas::io::fasta> *sequences);

void writeInfo(std::vector<std::string> *summary, std::string res);

void split_seq(std::vector<rappas::io::fasta> *sequences, std::vector<int> bp, int seq);

void new_q(std::vector<rappas::io::fasta> *sequences, std::vector<rappas::io::fasta> *nq, std::vector<int> *bp, std::vector<int> *seq, int ref);

std::string record(std::vector<rappas::io::fasta> *sequences, std::vector<int> bp, std::vector<int> seq);

void random_q(std::vector<rappas::io::fasta> *sequences, std::vector<rappas::io::fasta> *nq, int n, std::vector<std::string> *summary);

std::vector<std::vector<std::string>> readRes(std::string res);

void compRes(std::string res1, std::string res2, std::string group);

void compMosaic(std::string res1, std::string res2);

void compMosaicCirc(std::string res1, std::string res2);

void compComp(std::string res1, std::string res2, int shift);

void interpret(std::vector<std::string> *res, std::string seq, int k);

void fixRes(std::string res, std::string rep, std::string seq, int k);

void statsNA(std::string res);

std::string splitM(std::string res, int p);

int isGroup(std::vector<std::string> res, std::string gr);

void read_jstyle(std::string res, std::string rep);

void read_comet(std::string res, std::string rep);

void order_res(std::string res, std::string rep, std::string seq);

void prune(std::string res, std::string seq, int shift);

std::vector<int> readRead(std::string name);

std::vector<int> pickN(std::vector<rappas::io::fasta> *list, int max);

std::vector<std::string> lineErr(std::string line);

void readErr(std::string res, std::string rep);

void sortErr(std::string res, std::string rep, std::vector<rappas::io::fasta>* queries);

void infoReads(std::string clean, std::string err, std::string res, std::vector<rappas::io::fasta>* reads);

#endif /* INPUT_H_ */

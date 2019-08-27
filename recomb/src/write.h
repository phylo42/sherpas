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


void resInit(std::vector<Arc*> *res);

void windInit(std::vector<std::vector<Arc*>> *wind);

void SciPlot(int n, std::string_view query, std::vector<std::string> ref, std::vector<Arc> branches, std::vector<std::vector<Arc*>> res, int shift);

void Csv(int n, std::string query, std::vector<std::string> ref, std::vector<std::vector<Arc*>> res, int shift);

int stars(std::string gr);

std::string color(std::string gr);

std::string div(std::string bp);

void tikzLine(std::string rec, int y);

void tikzDoc(std::string real, std::string resA, std::string resB);


#endif /* WRITE_H_ */

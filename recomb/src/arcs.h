/*
 * arcs.h
 *
 *  Created on: 23 oct. 2018
 *      Author: scholz
 */

#ifndef ARCS_H_
#define ARCS_H_
#include<iostream>
#include<fstream>
#include<vector>
#include <string>

class Arc
{
public:
	Arc(int a, double p);
	void reinit();
	void updateScore(double p, int i);
	void updateCheck(int i);
	int ArcCheck();
	int getPlace();
	void printPlace();
	double getScore();
	int getKmers();
	int compareArc(int a);
	int BinSearch(std::vector<Arc*> list, int s);



private:
	int m_place;
	double m_score;
	int m_kmers;
	int m_check;
};

std::vector<Arc> getArcs(int s);

void clearBranches (std::vector<Arc>* b);


#endif /* ARCS_H_ */

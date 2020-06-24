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
	void printPlace(std::vector<std::string> ref);
	double getScore();
	int getKmers();
	int compareArc(int a);

private:
	int m_place;
	double m_score;
	int m_kmers;
	int m_check;
};

std::vector<Arc> getArcs(int s);

// CODE REVIEW: Prefer references to pointers to pass an argument
void clearBranches(std::vector<Arc>& b);


#endif /* ARCS_H_ */

/*
 * arcs.cpp
 *
 *  Created on: 23 oct. 2018
 *      Author: scholz
 */


#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include "arcs.h"

using namespace std;

Arc::Arc(int a, double p)
{
	m_place=a;
	m_score=p;
	m_check=0;
	if(p==0)
	{
		m_kmers=0;
	}
	else
	{
		m_kmers=1;
	}
}

void Arc::reinit()
{
	m_score=0;
	m_kmers=0;
	m_check=-1;
}

void Arc::updateScore(double p, int i)
{
	m_score+=p;
	m_kmers+=i;
}

void Arc::updateCheck(int i)
{
	m_check=i;
}

int Arc::ArcCheck()
{
	return m_check;
}

int Arc::getPlace()
{
	return m_place;
}

void Arc::printPlace(std::vector<std::string> ref)
{
	int s=ref.size();
	if(s>0 && m_place <s)
	{
		cout << ref[m_place] << " ";
	}
	else
	{
		cout << m_place << " ";
	}
	//cout << endl;
}

double Arc::getScore()
{
	return m_score;
}

int Arc::getKmers()
{
	return m_kmers;
}

int Arc::compareArc(int a)
{
	int r=0;
	if (m_place==a)
	{
		r=1;
	}
	return r;
}



std::vector<Arc> getArcs(int s)
{
	std::vector<Arc> res;
	for(int i=0;i<s+1;i++)
	{
		Arc a=Arc(i,0);
		res.push_back(a);
	}
	return res;
}

void clearBranches (std::vector<Arc>* b)
{
	int s=(*b).size();
	for(int i=0; i<s; i++)
	{
		(*b)[i].reinit();
	}
}

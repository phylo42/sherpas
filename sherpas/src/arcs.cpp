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
//basic operation to record at each step of the main algorithm, key information on arcs.

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

	// CODE REVIEW: Consider reserving memory for vectors of known size.
	// That makes push_back, emplace_back operations faster in practice.
	// See: https://en.cppreference.com/w/cpp/container/vector/reserve
	res.reserve(s + 1);

	for(int i=0;i<s+1;i++)
	{
	    // CODE REVIEW: push_back makes an unnecessary copy of object a in the worst case or a move in the best case.
	    // Both can be avoided with emplace_back.
	    // See: https://en.cppreference.com/w/cpp/container/vector/emplace_back
		//Arc a=Arc(i,0);
		//res.push_back(a);
		res.emplace_back(i, 0);
	}
	return res;
}

void clearBranches(std::vector<Arc>& b)
{
    // CODE REVIEW: vector::size actually returns size_t.
    // CODE REVIEW: const
	const size_t s = b.size();
	for(size_t i=0; i<s; i++)
	{
		b[i].reinit();
	}
}

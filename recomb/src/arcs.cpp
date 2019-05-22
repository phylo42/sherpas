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

void Arc::printPlace()
{
	cout << m_place << " ";
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

int Arc::BinSearch(std::vector<Arc*> list, int s)
{
	int p=list.size();
	//cout << p << "; " << s << endl;
	if(p>0)
	{
		if(s<1)
		{
			s=list.size();
		}
		if(s<p)
		{
			p=s;
		}
		if(m_score<(*list[p-1]).getScore())
		{
			p++;
		}
		else
		{
			int sup=p;
			int inf=0;
			p=p/2;
			int toobig=1;
			int toosmall=1;
			int ls=list.size();
			while((toobig+toosmall)>0)
			{
				toobig=(m_score>(*list[p]).getScore());
				if(toobig==1 && p==0)
				{
					p=-1;
					toobig=0;
				}
				if(p+1<ls)
				{
					toosmall=(m_score<(*list[p+1]).getScore());
				}
				else
				{
					toosmall=0;
				}
				if(toobig==1)
				{
					sup=p;
					p=(p+inf)/2;
				}
				if(toosmall==1)
				{
					inf=p;
					p=(p+sup)/2;
				}
			}
			p++;
		}
	}
	return p;
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

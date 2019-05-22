/*
 * heap.h
 *
 *  Created on: 19 févr. 2019
 *      Author: guillaume
 */

#ifndef HEAP_H_
#define HEAP_H_
#include<iostream>
#include<fstream>
#include<vector>
#include <string>
#include "arcs.h"

//https://www.techiedelight.com/min-heap-max-heap-implementation-c/

class Htree
{
public:
	Htree(std::vector<Arc*> res);
	void Hempty();
	Arc* h(int i);
	int hPlace(int i);
	double hScore(int i);
	int parent(int i);
	int left(int i);
	int right(int i);
	void print(int i);
	int size();
	void swap(int i, int j);
	void push(Arc* a);
	void pop(int i);
	void heapDown(int i);
	void heapUp(int i);
	void getTop(int m); //pb: currently emptying H...


private:
	std::vector<Arc*> m_res;
};


#endif /* HEAP_H_ */

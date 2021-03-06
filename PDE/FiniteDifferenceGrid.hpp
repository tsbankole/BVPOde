#pragma once
#include"Node.h"
#include<vector>
#include<cassert>
#include<iostream>

class FiniteDifferenceGrid {
public:
	friend class BvpOde;
private:
	double deltax;
	std::vector<double> mNodes;
public:
	FiniteDifferenceGrid(int numNodes, double xMin, double xMax) {
		deltax = (xMax - xMin) / ( numNodes -1 );
		for (int i = 0; i < numNodes; i++) {
			mNodes.push_back(xMin + i * (deltax));
		}
	}
};
#pragma once
#include"Node.h"
#include<vector>

class FiniteDifferenceGrid {
public:
	friend class BvpOde;
private:
	double deltax;
	std::vector<double> mNodes;
public:
	FiniteDifferenceGrid(int numNodes, double xMin, double xMax) {
		deltax = (xMax - xMin) / numNodes;
		for (int i = 0; i < numNodes; i++) {
			mNodes.push_back(xMin + i * (deltax));
		}
	}
};
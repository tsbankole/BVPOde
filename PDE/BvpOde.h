#pragma once
#include"LinearSystem.h"
#include<iostream>
#include <string>
#include "Matrix2D.h"
#include "Vector.h"
#include "SecondOrderOde.h"
#include "BoundaryConditions.hpp"
#include "FiniteDifferenceGrid.hpp"

class BvpOde {
private:
	BvpOde( const BvpOde& OtherBvpOde){}

	int mNumNodes;

	SecondOrderOde* mpOde;

	BoundaryConditions* mpBconds;

	Vector* mpSolVec;
	
	Vector* mpRhsVec;

	Matrix2D* mpLhsMat;

	LinearSystem* mpLinearSystem;

	FiniteDifferenceGrid* mpFDGrid;

	std::string mFilename;

	void PopulateMatrix();
	void PopulateVector();
	void ApplyBoundaryConditions();

public:
	BvpOde(SecondOrderOde* pOde, BoundaryConditions* pBcs, int numNodes);
	~BvpOde();

	void SetFilename(const std::string& name) {
		mFilename = name;
		/*std::cout << mpLhsMat->mData[0][0] << " " << mpLhsMat->mData[0][1] << "\n";
		int i = 238;
		std::cout << mpLhsMat->mData[i][i-1] << " " << mpLhsMat->mData[i][i] << "  " << mpLhsMat->mData[i][i+1] <<  "+++" << mpRhsVec->mData[i] <<  " xval " << mpFDGrid->mNodes[i] <<  "\n";*/
	}

	void Solve();
	void WriteSolutionFile();


};
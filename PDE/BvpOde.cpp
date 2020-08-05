// See E:\Sci_comp\scientific computing_Python_CCpp\SciCompCC++\X\creating in constructor for removing the address type instantiation of pointers in constructors
#include "BvpOde.h"
#include <cassert>
#include <iostream>
#include <fstream>

BvpOde::BvpOde(SecondOrderOde* pOde, BoundaryConditions* pBcs, int numNodes) {
	assert(numNodes > 2);
	assert((pBcs->mLhsBcIsDirichlet ^ pBcs->mLhsBcIsNeumann) && (pBcs->mRhsBcIsDirichlet ^ pBcs->mRhsBcIsNeumann));
	mNumNodes = numNodes;
	mpOde = pOde;
	mpBconds = pBcs;
	mpSolVec = new Vector(numNodes);
	//mpSolVec = &solvec;
	mpRhsVec = new Vector(numNodes);
	//mpRhsVec = &rhsvec;
	mpLhsMat = new Matrix2D(numNodes, numNodes);
	//mpLhsMat = &lhsmat;
	
	//mpLinearSystem = &ls;
	mpFDGrid = new FiniteDifferenceGrid(numNodes, mpOde->mXmin, mpOde->mXmax);
	PopulateMatrix();
	PopulateVector();
	ApplyBoundaryConditions();
	mpLinearSystem = new LinearSystem(*mpLhsMat, *mpRhsVec);
	mFilename = "bvpodefilename.dat";
	
}

BvpOde::~BvpOde() {
	delete mpSolVec;
	delete mpRhsVec;
	delete mpLhsMat;
	delete mpLinearSystem;
	delete mpFDGrid;
}

void BvpOde::PopulateMatrix() {
	double alph, beta, gamma, dx, Uxx, Ux, U;
	dx = mpFDGrid->deltax;
	Uxx = mpOde->mCoeffOfUxx;
	Ux = mpOde->mCoeffOfUx;
	U = mpOde->mCoeffOfU;
	alph = (Uxx/dx/dx + Ux/2/dx );
	beta = (-2 * Uxx/dx/dx + U);
	gamma = (Uxx/dx/dx - Ux/2/dx);

	// central difference for all rows except first and last, which are to be populated by Applying boundary conditions.
	for (int i = 1; i < mNumNodes - 1; i++) {
		mpLhsMat->mData[i][i + 1] = alph;
		mpLhsMat->mData[i][i - 1] = gamma;
		mpLhsMat->mData[i][i] = beta;
	}
}

void BvpOde::PopulateVector() {
	for (int i = 1; i < mNumNodes - 1; i++) {
		mpRhsVec->mData[i] = mpOde->mpRhsFunc(mpFDGrid->mNodes[i]);
	}
}

void BvpOde::ApplyBoundaryConditions() {
	// fill in the vector first and last elements according to Bconditions
	mpRhsVec->mData[0] = mpBconds->mLhsBcValue;
	mpRhsVec->mData[mNumNodes - 1] = mpBconds->mRhsBcValue;


   // fill in the first and last row entries for the LHS matrix
	if (mpBconds->mLhsBcIsDirichlet) {
		mpLhsMat->mData[0][0] = 1.0;
	}
	else if (mpBconds->mLhsBcIsNeumann) {
		mpLhsMat->mData[0][0] = -1.0 / mpFDGrid->deltax;
		mpLhsMat->mData[0][1] = 1.0 / mpFDGrid->deltax;
	}
	if (mpBconds->mRhsBcIsDirichlet) {
		mpLhsMat->mData[mNumNodes-1][mNumNodes-1] = 1.0;
	}
	else if (mpBconds->mRhsBcIsNeumann) {
		mpLhsMat->mData[mNumNodes - 1][mNumNodes - 2] = -1.0 / mpFDGrid->deltax;
		mpLhsMat->mData[mNumNodes - 1][mNumNodes - 1] = 1.0 / mpFDGrid->deltax;
	}
}

void BvpOde::WriteSolutionFile() {
	assert(mFilename != "");
	std::ofstream write_output(mFilename);
	assert(write_output.is_open());
	for (int i = 0; i < mNumNodes; i++) {
		write_output << mpSolVec->mData[i] << "\n";
	}
	write_output.close();
}

void BvpOde::Solve() {
	*mpSolVec = mpLinearSystem->Solve();
	WriteSolutionFile();
}
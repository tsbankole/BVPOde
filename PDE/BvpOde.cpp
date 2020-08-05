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
	Vector solvec = Vector(numNodes);
	mpSolVec = &solvec;
	Vector rhsvec = Vector(numNodes);
	mpRhsVec = &rhsvec;
	Matrix2D lhsmat = Matrix2D(numNodes, numNodes);
	mpLhsMat = &lhsmat;
	LinearSystem ls = LinearSystem(lhsmat, rhsvec);
	mpLinearSystem = &ls;
	mpFDGrid = &FiniteDifferenceGrid(numNodes, mpOde->mXmin, mpOde->mXmax);
	PopulateMatrix();
	PopulateVector();
	ApplyBoundaryConditions();
}

void BvpOde::PopulateMatrix() {
	double alph, beta, gamma, dx, Uxx, Ux, U;
	dx = mpFDGrid->deltax;
	Uxx = mpOde->mCoeffOfUxx;
	Ux = mpOde->mCoeffOfUx;
	U = mpOde->mCoeffOfU;
	alph = (Uxx / dx / dx + Ux/ dx );
	beta = (-2 * Uxx / dx/ dx + U);
	gamma = (Uxx / dx / dx - Ux / dx);

	// central difference for all rows except first and last, which are to be populated by Applying boundary conditions.
	for (int i = 1; i < mNumNodes - 1; i++) {
		mpLhsMat->mData[i][i + 1] = alph;
		mpLhsMat->mData[i][i - 1] = beta;
		mpLhsMat->mData[i][i] = gamma;
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
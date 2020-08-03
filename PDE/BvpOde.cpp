#include "BvpOde.h"
#include <cassert>

BvpOde::BvpOde(SecondOrderOde* pOde, BoundaryConditions* pBcs, int numNodes) {
	assert(numNodes > 2);
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
	mpFDGrid = &FiniteDifferenceGrid(numNodes, mpOde->mXmin, mpOde->mXmax);
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

	for (int i = 1; i < mNumNodes - 1; i++) {
		mpLhsMat->mData[i][i + 1] = alph;
		mpLhsMat->mData[i][i - 1] = beta;
		mpLhsMat->mData[i][i] = gamma;
	}
}

void BvpOde::PopulateVector() {
	mpRhsVec->mData[0] = mpBconds->mLhsBcValue;
	mpRhsVec->mData[mNumNodes - 1] = mpBconds->mRhsBcValue;
	for (int i = 1; i < mNumNodes - 1; i++) {
		mpRhsVec->mData[i] = mpOde->mpRhsFunc(mpFDGrid->mNodes[i]);
	}

}
#include "BoundaryConditions.hpp"
// set Dirichlet Boundary conditions Left and Right
void BoundaryConditions::SetLhsDirichletBc(double lhsValue) {
	mLhsBcIsDirichlet = true;
	mLhsBcIsNeumann = false;
	mLhsBcValue = lhsValue;
}

void BoundaryConditions::SetRhsDirichletBc(double rhsValue) {
	mRhsBcIsDirichlet = true;
	mRhsBcIsNeumann = false;
	mRhsBcValue = rhsValue;
}

// set Neumann Boundary conditions Left and Right
void BoundaryConditions::SetLhsNeumannBc(double lhsValue) {
	mLhsBcIsDirichlet = false;
	mLhsBcIsNeumann = true;
	mLhsBcValue = lhsValue;
}

void BoundaryConditions::SetRhsNeumannBc(double rhsValue) {
	mRhsBcIsDirichlet = false;
	mRhsBcIsNeumann = true;
	mRhsBcValue = rhsValue;
}
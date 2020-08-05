#include "BoundaryConditions.hpp"
// set Dirichlet Boundary conditions Left and Right
void BoundaryConditions::SetLhsDirichletBc(double lhsValue) {
	mLhsBcIsDirichlet = true;
	mLhsBcValue = lhsValue;
}

void BoundaryConditions::SetRhsDirichletBc(double rhsValue) {
	mRhsBcIsDirichlet = true;
	mRhsBcValue = rhsValue;
}

// set Neumann Boundary conditions Left and Right
void BoundaryConditions::SetLhsNeumannBc(double lhsValue) {
	mLhsBcIsNeumann = true;
	mLhsBcValue = lhsValue;
}

void BoundaryConditions::SetRhsNeumannBc(double rhsValue) {
	mRhsBcIsNeumann = true;
	mRhsBcValue = rhsValue;
}
#pragma once
class BoundaryConditions
{
public:
	friend class BvpOde;

private:
	bool mLhsBcIsDirichlet;
	bool mRhsBcIsDirichlet;
	bool mLhsBcIsNeumann;
	bool mRhsBcIsNeumann;
	double mLhsBcValue;
	double mRhsBcValue;

public:
	BoundaryConditions() {
		mLhsBcIsDirichlet = false;
		mLhsBcIsDirichlet = false;
		mRhsBcIsDirichlet = false;
		mRhsBcIsNeumann = false;
	};
	void SetLhsDirichletBc(double lhsValue);
	void SetRhsDirichletBc(double rhsValue);
	void SetLhsNeumannBc(double lhsDerivValue);
	void SetRhsNeumannBc(double rhsDerivValue);
};
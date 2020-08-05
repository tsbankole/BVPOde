#pragma once
#include<cassert>
class SecondOrderOde
{
	friend class BvpOde;
private:
	double mCoeffOfUxx;
	double mCoeffOfUx;
	double mCoeffOfU;

	double(*mpRhsFunc)(double x);

	double mXmin;
	double mXmax;

public:
	SecondOrderOde(double coeffUxx, double coeffUx, double coeffU,
		double(*rhs)(double), double xMinimum, double xMaximum)
	{
		assert(xMaximum > xMinimum);
		mCoeffOfUxx = coeffUxx;
		mCoeffOfUx = coeffUx;
		mCoeffOfU = coeffU;
		mpRhsFunc = rhs;
		mXmax = xMaximum;
		mXmin = xMinimum;
	}
};
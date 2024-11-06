#include "FEMConfig.h"
#include <cassert>
#include <cmath>

double FEMConfig::mu(int type)
{
	switch (type)
	{
	case 1: return 1;
	default: assert(false);
	}
}
double FEMConfig::gamma(int type)
{
	switch (type)
	{
	case 1: return 1;
	default: assert(false);
	}
}

Point FEMConfig::F(double x, double y, int type)
{
	double xRes, yRes;
	Point result;
	switch (type)
	{
	case 1:
		xRes = 2 * cos(y);
		yRes = 2 * sin(x);
		break;
	default: assert(false);
	}

	result.x = xRes;
	result.y = yRes;

	return result;
}

Point FEMConfig::S1(double x, double y, int type)
{
	double xRes, yRes;
	Point result;
	switch (type)
	{
	case 1:
		xRes = cos(y);
		yRes = sin(x);
		break;
	default: assert(false);
	}

	result.x = xRes;
	result.y = yRes;

	return result;
}
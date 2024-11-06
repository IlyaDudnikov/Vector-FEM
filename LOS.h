#pragma once
#include "operationsSolver.h"
#include <iostream>

void LOS_sq(int N, int maxit, type e, vector<double>& x, vector<int> ia, vector<int> ja,
	vector<double> di, vector<double> al, vector<double> au, vector<double> vec); // x - начальное приближение
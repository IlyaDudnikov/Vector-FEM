#pragma once
//#include <stdio.h>
//#include <math.h>
//#include <conio.h>
//#include <stdlib.h>
//#include <fstream>
//#include <iostream>
//#include <iomanip>
#include <vector>

using namespace std;

typedef double type;

void sum_vector(vector<double> a, vector<double> b, vector<double>& res, int N);

double scal_prod(vector<double> a, vector<double> b, int N);  // скалярное произведение

void mult_vector(vector<double> a, vector<double> b, vector<double>& res, int N);

void mult_coef(vector<double> vec, double k, vector<double>& res, int N);

void mult_matrix(vector<int> ia, vector<int> ja, vector<double> di,
	vector<double> al, vector<double> au, vector<double> vec, vector<double>& res, int N);

void LU_sq(int N, vector<double> al, vector<double> di, vector<double> au, vector<int> ia, vector<int> ja,
	vector<double>& al_sq, vector<double>& di_sq, vector<double>& au_sq);

void mult_L_obr(vector<double> aa, vector<double> b, vector<double>& y, int N, vector<int> ia, vector<int> ja, vector<double> di_sq);

void mult_U_obr(vector<double> aa, vector<double> b, vector<double>& y, int N, vector<int> ia, vector<int> ja, vector<double> di_sq);

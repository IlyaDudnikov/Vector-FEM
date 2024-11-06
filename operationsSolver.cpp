#include "operationsSolver.h"

void sum_vector(vector<double> a, vector<double> b, vector<double>& res, int N)
{
	for (int i = 0; i < N; i++)
		res[i] = a[i] + b[i];
}

double scal_prod(vector<double> a, vector<double> b, int N)  // скал€рное произведение
{
	double sum = 0;
	for (int i = 0; i < N; i++)
		sum += a[i] * b[i];
	return sum;
}

void mult_vector(vector<double> a, vector<double> b, vector<double>& res, int N)
{
	for (int i = 0; i < N; i++)
		res[i] = a[i] * b[i];
}

void mult_coef(vector<double> vec, double k, vector<double>& res, int N)
{
	for (int i = 0; i < N; i++)
		res[i] = k * vec[i];
}

void mult_matrix(vector<int> ia, vector<int> ja, vector<double> di,
	vector<double> al, vector<double> au, vector<double> vec, vector<double>& res, int N)
{
	for (int i = 0; i < N; i++)
	{
		int i0 = ia[i];
		int i1 = ia[i + 1];
		res[i] = di[i] * vec[i];
		for (int k = i0; k < i1; k++)
		{
			int j = ja[k];
			res[i] += al[k] * vec[j];
			res[j] += au[k] * vec[i];
		}
	}
}

void LU_sq(int N, vector<double> al, vector<double> di, vector<double> au, vector<int> ia, vector<int> ja,
	vector<double>& al_sq, vector<double>& di_sq, vector<double>& au_sq)
{
	//копирование-инициализаци€
	for (int i = 0; i < N; i++)
		di_sq[i] = di[i];

	for (int i = 0; i < ia[N]; i++)
	{
		al_sq[i] = al[i];
		au_sq[i] = au[i];
	}

	for (int i = 0; i < N; i++)
	{
		type sd = 0; //переменные суммировани€

		int i0 = ia[i];
		int i1 = ia[i + 1];

		for (int k = i0; k < i1; k++) {
			int j = ja[k];
			type sl = 0, su = 0;
			int j0 = ia[j];
			int j1 = ia[j + 1];
			int ki = i0;
			int kj = j0;

			for (; ki < k && kj < j1;) {
				int jl = ja[ki];
				int ju = ja[kj];
				if (jl == ju) {
					sl += au_sq[kj] * al_sq[ki];
					su += al_sq[kj] * au_sq[ki];
					ki++; kj++;
				}
				else if (jl < ju) ki++;
				else kj++;
			}
			au_sq[k] = (au_sq[k] - su) / di_sq[j];
			al_sq[k] = (al_sq[k] - sl) / di_sq[j];
			sd += au_sq[k] * al_sq[k];
		}

		di_sq[i] = sqrt(di_sq[i] - sd);
	}
}

void mult_L_obr(vector<double> aa, vector<double> b, vector<double>& y, int N, vector<int> ia, vector<int> ja, vector<double> di_sq)
{
	for (int i = 0; i < N; i++) {
		type s = 0; //переменные суммировани€

		int i0 = ia[i];//индекс 1го элемента в iтой строке
		int i1 = ia[i + 1];

		for (int k = i0; k < i1; k++) {
			int j = ja[k];
			s += y[j] * aa[k];
		}
		y[i] = (b[i] - s) / di_sq[i];
	}
}

void mult_U_obr(vector<double> aa, vector<double> b, vector<double>& y, int N, vector<int> ia, vector<int> ja, vector<double> di_sq)
{
	for (int i = 0; i < N; i++)
		y[i] = b[i];
	for (int i = N - 1; i >= 0; i--) {
		int i0 = ia[i];//индекс 1го элемента в iтой строке
		int i1 = ia[i + 1];

		y[i] /= di_sq[i];

		for (int k = i1 - 1; k >= i0; k--) {
			int j = ja[k];
			y[j] -= y[i] * aa[k];
		}
	}
}

#include "LOS.h"

void LOS_sq(int N, int maxit, type e, vector<double>& x, vector<int> ia, vector<int> ja,
	vector<double> di, vector<double> al, vector<double> au, vector<double> vec) // x - начальное приближение
{
	type a, b, skal1, skal2;
	int k;
	vector<double> r(N), p(N), z(N), al_sq(ia[N]), au_sq(ia[N]), di_sq(N), res1(N), res2(N), LAUr(N);
	LU_sq(N, al, di, au, ia, ja, al_sq, di_sq, au_sq);
	mult_matrix(ia, ja, di, al, au, x, res1, N);
	for (int i = 0; i < N; i++) {
		res2[i] = vec[i] - res1[i];
	}
	mult_L_obr(al_sq, res2, r, N, ia, ja, di_sq);
	mult_U_obr(au_sq, r, z, N, ia, ja, di_sq);
	mult_matrix(ia, ja, di, al, au, z, res1, N);
	mult_L_obr(al_sq, res1, p, N, ia, ja, di_sq);

	type nev = scal_prod(r, r, N);

	for (k = 0; k < maxit && nev > e; k++)
	{
		skal1 = scal_prod(p, r, N);
		skal2 = scal_prod(p, p, N);
		a = skal1 / skal2;
		mult_coef(z, a, res1, N);
		sum_vector(x, res1, x, N);
		mult_coef(p, -a, res1, N);
		sum_vector(r, res1, r, N);
		mult_U_obr(au_sq, r, LAUr, N, ia, ja, di_sq);
		mult_matrix(ia, ja, di, al, au, LAUr, res2, N);
		mult_L_obr(al_sq, res2, LAUr, N, ia, ja, di_sq);
		skal1 = scal_prod(p, LAUr, N);
		b = -skal1 / skal2;
		mult_U_obr(au_sq, r, res2, N, ia, ja, di_sq);
		mult_coef(z, b, res1, N);
		sum_vector(res2, res1, z, N);
		mult_coef(p, b, res1, N);
		sum_vector(LAUr, res1, p, N);
		nev = scal_prod(r, r, N);
	}
	/*std::cout << "iteration: " << k << "\n";
	std::cout << "squared norm of the residual: " << nev << "\n";*/
}

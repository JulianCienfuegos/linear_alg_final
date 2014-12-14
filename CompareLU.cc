/* 
 * This code surveys C++ linear system solvers which use the LU
 * Factorization. We will consider:
 * 
 * BLAS
 * An LU factorization and solver which I wrote.
 * Armadillo
 * Eigen
 * 
 * Compile with: g++ -I /usr/local/include/eigen3 CompareLU.cc -larmadillo -llapack -lblas -o CompareLU
 */

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <armadillo>
#include <Eigen/Dense>

int N = 1024;
const int num_iters = 1;


extern "C" void dgetrf_(int* dim1, int* dim2, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );

using namespace std;
using namespace arma;
using Eigen::MatrixXd;
using Eigen::PartialPivLU;

void lufac(double **, double **);
void fSub(double **, double *, double *);
void bSub(double **, double *, double *);
void MVMult(double **, double *, double *);

int main()
{
	/********************** ALLOCATE MEMORY ***************************/
	/*                   ARRAYS FOR MY SOLVER                         */
	double ** A = new double * [N];
	double ** P = new double * [N];
	double ** L = new double * [N];
	double ** U = new double * [N];
	double * b = new double [N];
	double * b_hat = new double [N];
	double * y = new double [N];
	double *x = new double [N];

	for(int i = 0; i < N; i++)
	{
		A[i] = new double [N];
		P[i] = new double [N];
		L[i] = new double [N];
		U[i] = new double [N];
	}
	
	srand((signed)time(NULL));
	
	for (int i = 0; i < N; i++)
	{
	  for (int j = 0; j < N; j++)
	  {
		  A[i][j] = rand() % 100 + 1;
		  L[i][j] = 0;
		  U[i][j] = 0;
		  if(i == j) // Make the identity matrix
			P[i][j] = 1;
		  else
			P[i][j] = 0;
	  }
	}

	for (int i = 0; i < N; i++)
	{
		b[i] = rand() % 100 + 1;
	}
	/*                    BLAS SOLVER ARRAYS                          */
	double * A_blas = new double [N*N];
	double * b_blas = new double [N];
	int ipiv[N];
	for (int i = 0; i < N*N; i++)
		A_blas[i] = (double)(rand()%100);
	for (int i = 0; i < N; i++)
		b_blas[i] = (double)(rand()%100);
	/*                  ARMADILLO MATRICES                           */
	mat A_arma = randu<mat>(N, N);
	mat P_arma, L_arma, U_arma;
	mat b_arma = zeros<mat>(N,1);
	mat x_arma = zeros<mat>(N,1);
	
	/*                    EIGEN MATRICES                             */
	MatrixXd A_eigen = MatrixXd::Random(N,N);
	MatrixXd x_eigen = MatrixXd::Zero(N,1);
	MatrixXd b_eigen = MatrixXd::Random(N,1);
	/*                    OTHER VARIABLES                            */
	char trans = 'N';
	int dim = N;
	int nrhs = 1;
	int LDA = dim;
	int LDB = dim;
	int info;
	double run_times[4];
	clock_t start, stop;

	/************************** MY SOLVER *****************************/
int yes = 0;
if(yes == 1)
{
	printf("Running My Solver...\n");
	start = clock();
	for(int iter = 0; iter < num_iters; iter++)
	{
		lufac(A,P);
	
		/* Get U */
		for (int i = 0; i < N; i++)
		{
			for(int j = i; j < N; j++)
			{
				U[i][j] = A[i][j];
			}
		}
	
		/* Get L */
		for(int i = 0; i < N; i++)
		{
			L[i][i] = 1;
			for(int j = i-1; j >=0; j--)
			{
				L[i][j] = A[i][j];
			}
		}
			
		/*Now we multiply b by P and solve PAx = Pb => LUx = Pb*/
		MVMult(P, b, b_hat);
	
		/* Solve Ly = b_hat */
		fSub(L, y, b_hat);
	
		/* Solve Ux =  y */
		bSub(U, x, y);
		
		/* Now we can do what we want with the x! */
	}
	stop = clock();
	run_times[0] = (stop - start)/(num_iters*CLOCKS_PER_SEC);
}	
	/* We can clean up a bit */
	for (int i = 0; i < N; i++)
	{
		delete[] A[i];
		delete[] P[i];
		delete[] L[i];
		delete[] U[i];
	}
	delete[] A;
	delete[] P;
	delete[] L;
	delete[] U;
	delete[] b;
	delete[] b_hat;
	delete[] x;
	delete[] y;
	
	/*********************** BLAS SOLVER ******************************/
	printf("Running BLAS...\n ");
	start = clock();
	for(int iter = 0; iter < num_iters; iter++)
	{
		for (int i = 0; i < N; i++)
			b_blas[i] = (double)(rand()%100);
		dgetrf_(&dim, &dim, A_blas, &LDA, ipiv, &info);
		dgetrs_(&trans, &dim, &nrhs, A_blas, &LDA, ipiv, b_blas, &LDB, &info);
	}
	stop = clock();
	clock_t start1 = clock();
	for(int i = 0; i < num_iters; i++)
	{
		for (int i = 0; i < N; i++)
			b_blas[i] = (double)(rand()%100);
	}
	clock_t stop1 = clock();
	run_times[1] = (double)(stop - start)/(num_iters*CLOCKS_PER_SEC);
	run_times[1] -= (double)(stop1 - start)/(num_iters*CLOCKS_PER_SEC);
	delete[] A_blas;
	delete[] b_blas;
	
	/************************ ARMADILLO *******************************/
	printf("Running Armadillo...\n");
	start = clock();
	for(int iter = 0; iter < num_iters; iter++)
	{
		lu(L_arma, U_arma, P_arma, A_arma);
		x_arma = solve(trimatu(U_arma), solve(trimatl(L_arma), P_arma*b_arma));	
	}
	stop = clock();
	run_times[2] = (double)(stop - start)/(num_iters*CLOCKS_PER_SEC);
	
	/************************** EIGEN *********************************/
	printf("Running Eigen...\n");
	/*
	start = clock();
	for(int iter = 0; iter < num_iters; iter++)
	{
		x_eigen = A_eigen.lu().solve(b_eigen);
	}
	stop = clock();
	run_times[3] = (double)(stop - start)/(num_iters*CLOCKS_PER_SEC);
	*/
	/********************* PRINT RESULTS ******************************/
	printf("My Solver: %+.0e\nBLAS: %+.0e \nArmadillo: %+.0e\nEigen: %+.0e\n", \
		run_times[0], run_times[1], run_times[2], run_times[3]);
}

void lufac(double ** A, double ** P)
{
	for(int i = 0; i < N-1; i++) // loop over all columns except the last.
	{
		
		double max = abs(A[i][i]);
		int maxIdx = i; //row i has the max right now
		for(int j = i+1; j < N; j++) //loop over rows
		{
			if(abs(A[j][i]) > max){maxIdx = j;} // Find max entry
		}
		if(maxIdx != i) // Perform a swap if the max isn't on the diagonal.
		{
			for(int col = 0; col < N; col++) //swap rows of A and P
			{
				double temp = A[i][col];
				A[i][col] = A[maxIdx][col];
				A[maxIdx][col] = temp;
				double temp2 = P[i][col];
				P[i][col] = P[maxIdx][col];
				P[maxIdx][col] = temp2;
			}
		}
		
		for (int j = i+1; j < N; j++)
		{
			A[j][i] = A[j][i]/A[i][i];
		}
		for(int j = i+1; j < N; j++)
		{
			for (int k = i+1; k < N; k++)
			{
				A[j][k] = A[j][k] - A[j][i]*A[i][k];
			}
		}
	}
}

void bSub(double ** U, double * x, double * b)
{
	for(int i = N-1; i >=0 ; i--)
	{
		for(int j = i+1; j < N; j++)
		{
			x[i] -= x[j]*U[i][j]; 
		}
		x[i] = (x[i] + b[i])/U[i][i];
	}

}

void fSub(double ** L, double * x, double * b)
{
	for(int i = 0; i < N; i++)
	{
		for (int j = 1; j < i; j++)
		{
			b[i] -= L[i][j]*x[j];
		}
		x[i] = b[i]/L[i][i];
	}
}

void MVMult(double ** M, double * x, double * b)
{
	for(int i = 0; i < N; i++)
	{
		b[i] = 0; /* This could be done elsehwere*/
	}
	/*Could this be optimized?*/
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			b[i] += M[i][j]*x[j];
		}
	}
}

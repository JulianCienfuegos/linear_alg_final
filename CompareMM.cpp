/*
 * This code surveys matrix multiplication in C++
 * 
 * 
 * BLAS
 * Tiled MM on one core
 * Tiled MM on multiple cores
 * Armadillo MM
 * Eigen MM
 * 
 * are what we will experiment with.
 * 
 * Compilation Instructions: 
 * g++ -I /usr/local/include/eigen3 CompareMM.cpp -larmadillo -lblas -o CompareMM
 * 
 * Where -I /usr/local/include/eigen3 is the path to your Eigen directory.
 * 
 */  

#include <stdio.h>
#include <ctime>
#include <armadillo>
#include <omp.h>
#include <Eigen/Dense>

#define bSize 64
#define NITER 1
int N = 1024;

/*
 * The function dgemm() performs one multiplication of the form 
 * C = alpha*A*B + beta*C. Below we will set alpha to 1 and beta to zero.
 */

extern "C" void dgemm_(char *transa, char *transb, int *m, int *n, int* k, 
			 double *alpha, double *A, int *lda, double *B, int *ldb,
			 double *beta,  double *C, int *ldc ); 


void BMM(double **, double **, double **);
void omp_BMM(double **, double **, double **);

using namespace std;
using namespace arma;
using Eigen::MatrixXd;

int main(int argc, char * argv[])
{
/************************ MAKE MATRICES **************************/
/*                         C++ 2D ARRAY                          */
double ** A_block = new double * [N];
double ** B_block = new double * [N];
double ** C_block = new double * [N];

for(int i = 0; i < N; i++)
{
	A_block[i] = new double [N];
	B_block[i] = new double [N];
	C_block[i] = new double [N];
}
srand(time(NULL));
for(int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	{
		A_block[i][j] = (double)(rand()%100);
		B_block[i][j] = (double)(rand()%100);
	}	
	
/*                        C++ 1D ARRAY                           */
double * A_blas = new double [N*N];
double * B_blas = new double [N*N];
double * C_blas = new double [N*N];
for (int i = 0; i < N*N; i++)
{
	A_blas[i] = (double)(rand()%100);
	B_blas[i] = (double)(rand()%100);
}

/*                         ARMADILLO MATRICES                    */
mat A = randu<mat>(N, N);
mat B = randu<mat>(N, N);
mat C = zeros(N,N);

/*                          EIGEN MATRICES                       */

MatrixXd A_eigen = MatrixXd::Random(N,N);
MatrixXd B_eigen = MatrixXd::Random(N,N);
MatrixXd C_eigen = MatrixXd::Zero(N,N);

/*                         OTHER VARIABLES                       */
double run_times[5];
clock_t start, stop;
char trans = 'N';
double alpha = 1;
double beta = 0;
/*                           BEGIN COMPARE!                      */

printf("Running BMM...\n");
start = clock();
for(int i = 0; i < NITER; i++)
	BMM(A_block, B_block, C_block);
stop = clock();
run_times[0] = (double)(stop - start)/(NITER*CLOCKS_PER_SEC);

printf("Running omp_BMM...\n");
start = clock();
for(int i = 0; i < NITER; i++)
	omp_BMM(A_block, B_block, C_block);
stop = clock();
run_times[1] = (double)(stop - start)/(NITER*CLOCKS_PER_SEC);

printf("Running Armadillo MM...\n");
start = clock();
for(int i = 0; i < NITER; i++)
	C = A*B;
stop = clock();
run_times[2] = (double)(stop - start)/(NITER*CLOCKS_PER_SEC);

printf("Running BLAS...\n");
start = clock();
for(int i = 0; i < NITER; i++)
	dgemm_(&trans, &trans, &N, &N, &N, &alpha, A_blas, &N, B_blas,\
		&N, &beta, C_blas, &N);
stop = clock();
run_times[3] = (double)(stop - start)/(NITER*CLOCKS_PER_SEC);

printf("Running Eigen...\n");
start = clock();
for(int i = 0; i < NITER; i++)
	C_eigen = A_eigen*B_eigen;
stop = clock();
run_times[4] = (double)(stop - start)/(NITER*CLOCKS_PER_SEC);
	
printf("BMM: %f\nomp_BMM: %f \nArmadillo: %f\nBLAS: %f\nEigen: %f\n", \
	run_times[0], run_times[1], run_times[2], run_times[3], run_times[4]);

}

void BMM(double ** A, double ** B, double ** C)
{
/* Block Matrix Multiply */
for(int i = 0; i < N; i += bSize)
	for(int k = 0; k < N; k += bSize)
		for(int j = 0; j < N; j += bSize)
			for(int it = i; it < i+bSize; it++)
				for(int kt = k; kt < k+bSize; kt++)
					for(int jt = j; jt < j+bSize; jt++)
					{
						C[it][jt] += A[it][kt]*B[kt][jt];
					}
}


void omp_BMM(double ** A, double ** B, double ** C)
{
int i, j, k, it, jt, kt;
/* openMP version of BMM */
#pragma omp parallel shared(A,B,C), private(i,j,k, it, jt, kt)
{
#pragma omp for schedule(static)
for(int i = 0; i < N; i += bSize)
	for(int k = 0; k < N; k += bSize)
		for(int j = 0; j < N; j += bSize)
			for(int it = i; it < i+bSize; it++)
				for(int kt = k; kt < k+bSize; kt++)
					for(int jt = j; jt < j+bSize; jt++)
						C[it][jt] += A[it][kt]*B[kt][jt];
}
}

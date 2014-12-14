// LAPACK test code
// compile with: g++ main.cpp -llapack -lblas -o testprog

#include <iostream>
#include <vector>

using namespace std;

extern "C" void dgetrf_(int* dim1, int* dim2, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );

int main()
{
  char trans = 'T';
  int dim = 3;
  int nrhs = 1;
  int LDA = dim;
  int LDB = dim;
  int info;

  double a[dim*dim];
  double b[dim];
  int ipiv[dim];
  	for(int i = 0; i < dim*dim; i++)
	{
		std::cout << a[i] << " ";
		if ((i+1)%dim == 0)
			std::cout << std::endl;
	}
	for(int i = 0; i < dim; i++)
	{
		std::cout << b[i] << " ";
	}

  dgetrf_(&dim, &dim, a, &LDA, ipiv, &info);
  dgetrs_(&trans, &dim, &nrhs, a, &LDA, ipiv, b, &LDB, &info);

	for(int i = 0; i < dim*dim; i++)
	{
		std::cout << a[i] << " ";
		if ((i+1)%dim == 0)
			std::cout << std::endl;
	}
	for(int i = 0; i < dim; i++)
	{
		std::cout << b[i] << " ";
	}

  return(0);
}

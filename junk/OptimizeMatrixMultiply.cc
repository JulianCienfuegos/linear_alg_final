#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#define N 500

using namespace std;

void ijk(int A[N][N], int B[N][N], int C[N][N])
{
	// Should be second best
	for(int i=0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

void ikj(int A[N][N], int B[N][N], int C[N][N])
{
	// Should be the best
	for(int i=0; i < N; i++)
	{
		for (int k = 0; k < N; k++)
		{
			for (int j = 0; j < N; j++)
			{
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}


void jki(int A[N][N], int B[N][N], int C[N][N])
{
	// Should be the worst
	for(int j=0; j < N; j++)
	{
		for (int k = 0; k < N; k++)
		{
			for (int i = 0; i < N; i++)
			{
				C[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

void MakeArray(int** M)
{
	M = new int* [N];
	for(int i = 0; i < N; i++)
	{
		M[i] = new int [N];
	}
}

int main()
{
/*
  int A[N][N];
  int B[N][N];
  int C[N][N];
  int D[N][N];
  int E[N][N];
*/
  
  int **A = new int *[N];
  int **B = new int *[N];
  int **C = new int *[N];	
  int **D = new int *[N];
  int **E = new int *[N];
  MakeArray(A);
  MakeArray(B);
  MakeArray(C);
  MakeArray(D);
  MakeArray(E);
  
  srand((unsigned)time(NULL));

  for (int i = 0; i < N; i++)
  {
      for (int j = 0; j < N; j++)
      {
		  A[i][j] = rand() % 10;
		  B[i][j] = rand() % 10;
		  C[i][j] = 0;
		  D[i][j] = 0;
		  E[i][j] = 0;
      }
  }
  clock_t begin = clock();
  ijk(**A, **B, **C);
  clock_t end = clock();
  cout << "(2nd?)Elapsed time for ijk is " << (double)(end - begin)/CLOCKS_PER_SEC << endl;
  
  begin = clock();
  ikj(**A, **B, **D);
  end = clock();
  cout << "(1st?)Elapsed time for ikj is " << (double)(end - begin)/CLOCKS_PER_SEC << endl;
  
  begin = clock();
  jki(**A, **B, **E);
  end = clock();
  cout << "(3rd?)Elapsed time for jki is " << (double)(end - begin)/CLOCKS_PER_SEC << endl;
  
}

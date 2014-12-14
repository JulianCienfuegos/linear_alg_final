#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>

#define N  1024
#define num_iters  2
#define bSize  32

using namespace std;

void blockMult(float ** A, float ** B, float ** C)
{
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

void lufac(float ** A, float ** P)
{
	for(int i = 0; i < N-1; i++) // loop over all columns except the last.
	{
		
		float max = abs(A[i][i]);
		int maxIdx = i; //row i has the max right now
		for(int j = i+1; j < N; j++) //loop over rows
		{
			if(abs(A[j][i]) > max){maxIdx = j;} // Find max entry
		}
		if(maxIdx != i) // Perform a swap if the max isn't on the diagonal.
		{
			for(int col = 0; col < N; col++) //swap rows of A and P
			{
				float temp = A[i][col];
				A[i][col] = A[maxIdx][col];
				A[maxIdx][col] = temp;
				float temp2 = P[i][col];
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

void MVMult(float ** M, float * x, float * b)
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

void bSub(float ** U, float * x, float * b)
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

void fSub(float ** L, float * x, float * b)
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

int main()
{
	/* Allocate memory and set up arrays*/
	float ** A = new float * [N];
	float ** P = new float * [N];
	float ** L = new float * [N];
	float ** U = new float * [N];
	float ** A_const = new float * [N];
	float ** C = new float * [N];
	float * b = new float [N];

	for(int i = 0; i < N; i++)
	{
		A[i] = new float [N];
		P[i] = new float [N];
		L[i] = new float [N];
		U[i] = new float [N];
		A_const[i] = new float [N];
		C[i] = new float [N];
	}
	
	srand((signed)time(NULL));
	
	for (int i = 0; i < N; i++)
	{
	  for (int j = 0; j < N; j++)
	  {
		  A[i][j] = rand() % 100 + 1;
		  L[i][j] = 0;
		  U[i][j] = 0;
		  A_const[i][j] = A[i][j];
		  C[i][j] = 0;
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
	
	/* Now we perform the LU factorization of A */
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
	
	/* We can clean up a bit */
	for (int i = 0; i < N; i++)
	{
		delete[] A[i];
		delete[] A_const[i];
		delete[] C[i];
	}
	delete[] A;
	delete[] A_const;
	delete[] C;

	/* Set up some things */
	float * b_hat = new float [N];
	float * y = new float [N];
	float *x = new float [N];
	for (int i = 0; i < N; i++)
	{
		b_hat[i] = 0;
		y[i] = 0;
		x[i] = 0;
	}

	/*Now we multiply b by P and solve PAx = Pb => LUx = Pb*/
	MVMult(P, b, b_hat);

	/* Solve Ly = b_hat */
	fSub(L, y, b_hat);

	/* Solve Ux =  y */
	bSub(U, x, y);
	
	/* Now we can do what we want with the x! */
	
	/* We can clean up a bit */
	for (int i = 0; i < N; i++)
	{
		delete[] P[i];
		delete[] L[i];
		delete[] U[i];
	}
	delete[] P;
	delete[] L;
	delete[] U;
	delete[] b;
	delete[] b_hat;
	delete[] x;
	delete[] y;
}

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

#define N 5
int main(int argc, char** argv)
{
	mat A = randu<mat>(N, N);
	mat B = randu<mat>(N, N);
	mat C = zeros<mat>(N, N);
	cout << C << endl;
	C = A*B;
	
	return 0;
	
}

#include <iostream>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using namespace std;
#define N 3
int main()
{

MatrixXi x = MatrixXi::Random(N,N);
MatrixXi y = MatrixXi::Random(N,N);
MatrixXi z = MatrixXi::Zero(N,N);
cout << z << endl;

z = x*y;

cout << z<<endl;
}

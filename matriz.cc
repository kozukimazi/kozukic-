#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <unsupported/Eigen/KroneckerProduct>

using Eigen::MatrixXd;
using namespace std;

Eigen:: MatrixXcd anticonmutator(Eigen::MatrixXcd,Eigen::MatrixXcd);
Eigen:: MatrixXcd JordanWigner(Eigen::MatrixXcd,Eigen::MatrixXcd, int,int);

int main()
{
  Eigen::MatrixXcd sigmax(2,2);
    sigmax(0,0) = 0;
    sigmax(0,1) = 1;
    sigmax(1,0) = 1;
    sigmax(1,1) = 0;
  Eigen::MatrixXcd sigmay(2,2);
    sigmay(0,0) = 0;
    sigmay(0,1) = complex<double>(0,-1.0);
    sigmay(1,0) = complex<double>(0,1.0);
    sigmay(1,1) = 0;
  Eigen::Matrix2d sigmaz(2,2);
  sigmaz(0,0) = 1;
  sigmaz(1,1) = -1;  
  sigmaz(0,1) = 0;
  sigmaz(1,0) = 0;  
  
  cout<< anticonmutator(sigmaz,sigmaz) << endl;
  Eigen::Matrix3d A;
    A << 1, 0,0,
         0, 0,0,
         0, 0,1;

  Eigen::Matrix2d B;
    B << 0, 0,
         0, 1;

  Eigen::MatrixXcd K = Eigen::kroneckerProduct(A, B).eval();

  std::cout << "Kronecker Product:\n" << K << std::endl;
  return 0;

}

Eigen::MatrixXcd anticonmutator(Eigen::MatrixXcd A, Eigen::MatrixXcd B){
  return (A*B) + (B*A);
}

Eigen::MatrixXcd JordanWigner(Eigen::MatrixXcd sigmax, Eigen::MatrixXcd sigmaz, int n,int m){
  
  return sigmax;
}



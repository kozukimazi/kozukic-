#include "Qtermo.h"


Qtermo::Qtermo(const MatrixXcd &H,int n){
    Hamiltonian = H;
    Hilbert = n;
}

MatrixXcd Qtermo::anticonmutator(const MatrixXcd &A, const MatrixXcd & B){
    return A*B + B*A;
}


VectorXcd Qtermo::reshape(MatrixXcd &A){
   VectorXcd vecA = Map<VectorXcd>(A.data(), A.size());
  return vecA;
}

MatrixXcd Qtermo::reshapeback(VectorXcd &A, int n){
    MatrixXcd matA = Map<MatrixXcd>(A.data(),n,n);
  return matA;
}

MatrixXcd Qtermo::isolatedevol(const MatrixXcd &H){
  int d = H.rows();
  MatrixXcd Id = MatrixXcd::Identity(d, d);
  MatrixXcd totji = (kroneckerProduct(Id,H).eval() - kroneckerProduct(H.transpose(),Id).eval());

  return complex<double>(0,-1.0)*totji;
}


MatrixXcd Qtermo::evol(const MatrixXcd &rho0, double t,int n){
   MatrixXcd liouville = t*isolatedevol(Hamiltonian);
  //cout<<"la matriz:"<<endl;
  //cout<<liouville.exp()<<endl;
  MatrixXcd m0 = rho0;
  VectorXcd init0 = reshape(m0);
  //cout<<init0<<endl;
  VectorXcd initf = liouville.exp() * init0;
  

  return reshapeback(initf,n);
  };





#include "jordanwigner.h"

void jordanwigner:: setdim( int hil){
    hilbert = hil;
}

Eigen::MatrixXcd jordanwigner::tensor(Eigen::Matrix2d M, int N){
    Eigen::MatrixXcd Fin = M;
    for (int i = 1; i<N; i++){ 
       Fin = Eigen::kroneckerProduct(M, Fin).eval();
    }
    return Fin;
}


Eigen::MatrixXcd jordanwigner::opcreator(int l){
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
    Eigen::MatrixXcd sigmam(2,2);
    sigmam(0,0) = 0;
    sigmam(0,1) = 1;
    sigmam(1,0) = 0;
    sigmam(1,1) = 0;  
  Eigen::Matrix2d sigmaz(2,2);
  sigmaz(0,0) = 1;
  sigmaz(1,1) = -1;  
  sigmaz(0,1) = 0;
  sigmaz(1,0) = 0;
  Eigen::Matrix2d I2 = Eigen::Matrix2d::Identity();
  Eigen::MatrixXcd Iden = tensor(I2, hilbert-l);
  Eigen::MatrixXcd tot;
  if(l > 1){
//here we defined the tensorial product
  if (l!=hilbert){
  Eigen::MatrixXcd sigmazs = tensor(sigmaz, l-1);
  Eigen::MatrixXcd aux = Eigen::kroneckerProduct(sigmazs,sigmam).eval();
   tot = Eigen::kroneckerProduct(aux,Iden).eval();}
  else{
    Eigen::MatrixXcd sigmazs = tensor(sigmaz, l-1);
    Eigen::MatrixXcd aux = Eigen::kroneckerProduct(sigmazs,sigmam).eval();
    tot = aux;
  } 
  }

  else{
    if (l!=hilbert){      
    tot = Eigen::kroneckerProduct(sigmam,Iden).eval();}
    else{
      tot = sigmam;
    }
  }

  return tot;
  
};
#include <iostream>
#include "jordanwigner.h"

Eigen::MatrixXcd anticonmutator(const Eigen::MatrixXcd& , const Eigen::MatrixXcd&  );
Eigen::VectorXcd reshape(Eigen::MatrixXcd& );
Eigen::MatrixXcd reshapeback( Eigen::VectorXcd& );

int main()
{
  jordanwigner tot;
  tot.setdim(2);  
  Eigen::Matrix2d sigmaz(2,2);
  sigmaz(0,0) = 1;
  sigmaz(1,1) = -1;  
  sigmaz(0,1) = 0;
  sigmaz(1,0) = 0;
  Eigen::MatrixXcd Df1 = tot.opcreator(1);
  Eigen::MatrixXcd Df2 = tot.opcreator(2);
  //Eigen::MatrixXcd Df3 = tot.opcreator(3);
  Eigen::VectorXcd vec1 = reshape(Df1);
  //cout<< anticonmutator(Df2.adjoint(),Df2)<< endl;
  cout<< Df1 <<endl;
  cout<< vec1<<endl;
  return 0;
}

Eigen::MatrixXcd anticonmutator(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd&  B){
  return A*B + B*A;
}


Eigen::VectorXcd reshape(Eigen::MatrixXcd& A){
  Eigen::VectorXcd vecA = Eigen::Map<Eigen::VectorXcd>(A.data(), A.size());
  return vecA;
}




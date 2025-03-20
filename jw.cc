#include <iostream>
#include "jordanwigner.h"
#include <fstream>
#include <iomanip>



//all big matrices for references
Eigen::MatrixXcd anticonmutator(const Eigen::MatrixXcd& , const Eigen::MatrixXcd&  );
Eigen::VectorXcd reshape(Eigen::MatrixXcd& );
Eigen::MatrixXcd reshapeback( Eigen::VectorXcd&,int);
Eigen::MatrixXcd isolatedevol(const Eigen::MatrixXcd& );
Eigen::MatrixXcd evol(const Eigen::MatrixXcd&, const Eigen::MatrixXcd&, double,int);


int main()
{
  jordanwigner tot;
  double eps = 0.4;
  //double t = 0.1;
  tot.setdim(1);  
  //Eigen::Matrix4d rho0;  // 4x4 matrix of doubles

    // Initialize matrix elements
    //rho0 << 0.25, 0, 0, 0,
      //      0, 0.25, 0, 0,
        //    0, 0, 0.25, 0,
          //  0, 0, 0, 0.25;
  Eigen::Matrix2d rho0;

    rho0 << 0.5, 0.2, 
            0.2, 0.5;



  Eigen::MatrixXcd Df1 = tot.opcreator(1);
  //Eigen::MatrixXcd Df2 = tot.opcreator(2);
  Eigen::MatrixXcd H = eps*(Df1.adjoint()*Df1 + Df1.adjoint()*Df1);
  Eigen::VectorXcd vec1 = reshape(Df1);

  //Now here we make the file
  //cout<<H<<endl;
  //cout<<"rand"<<endl;
  ofstream outFile("datos1.dat");

  outFile<<fixed<<setprecision(6);
  
  double t0 = 0;
  double dt = 0.2;

  for(int i =0; i <4;i++){
    
    //cout<<t0<<endl;

    Eigen::MatrixXcd Mat = evol(rho0,H,t0,2);
    cout<<"la matriz es:"<<endl;
    cout<<Mat<<endl;
    
    outFile<<t0<< "\t"<<Mat(1,1).real()<<endl;
    t0 = t0+dt;
  }

   outFile.close();   

  return 0;
}

Eigen::MatrixXcd anticonmutator(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd&  B){
  return A*B + B*A;
}


Eigen::VectorXcd reshape(Eigen::MatrixXcd& A){
  Eigen::VectorXcd vecA = Eigen::Map<Eigen::VectorXcd>(A.data(), A.size());
  return vecA;
}

Eigen::MatrixXcd reshapeback(Eigen::VectorXcd& A,int n){
  Eigen::MatrixXcd matA = Eigen::Map<Eigen::MatrixXcd>(A.data(),n,n);
  return matA;
}


Eigen::MatrixXcd isolatedevol(const Eigen::MatrixXcd & H){
  int d = H.rows();
  Eigen::MatrixXcd Id = Eigen::MatrixXcd::Identity(d, d);
  Eigen::MatrixXcd totji = (Eigen::kroneckerProduct(Id,H).eval() - Eigen::kroneckerProduct(H.transpose(),Id).eval());

  return complex<double>(0,-1.0)*totji;
}

Eigen::MatrixXcd evol(const Eigen::MatrixXcd &rho0,const Eigen::MatrixXcd & H, double t,int n){
  
  Eigen::MatrixXcd liouville = t*isolatedevol(H);
  cout<<"la matriz:"<<endl;
  //cout<<liouville.exp()<<endl;
  Eigen::MatrixXcd m0 = rho0;
  Eigen::VectorXcd init0 = reshape(m0);
  //cout<<init0<<endl;
  Eigen::VectorXcd initf = liouville.exp() * init0;
  

  return reshapeback(initf,n);

}

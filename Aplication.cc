#include <iostream>
#include "jordanwigner.h"
#include "Qtermo.h"
#include <fstream>
#include <iomanip>


int main()
{
  jordanwigner tot;
  double eps = 0.4;
  //double t = 0.1;
  tot.setdim(2);  
  Eigen::Matrix4d rho0;  // 4x4 matrix of doubles

    // Initialize matrix elements
    rho0 << 0.25, 0, 0, 0,
            0, 0.25, 0.1, 0,
            0, 0.1, 0.25, 0,
            0, 0, 0, 0.25;


  MatrixXcd Df1 = tot.opcreator(1);
  MatrixXcd Df2 = tot.opcreator(2);
  MatrixXcd H = 0.1*Df1.adjoint()*Df1 +eps*(Df1.adjoint()*Df2 + Df2.adjoint()*Df1);
  
  Qtermo prod(H,H.rows());

  //Now here we make the file
  //cout<<H<<endl;
  //cout<<"rand"<<endl;
  ofstream outFile("datosap.dat");

  outFile<<fixed<<setprecision(6);
  
  double t0 = 0;
  double dt = 0.2;

  for(int i =0; i <100;i++){
    
    //cout<<t0<<endl;

    MatrixXcd Mat = prod.evol(rho0,t0,4);
    //cout<<"la matriz es:"<<endl;
    //cout<<Mat<<endl;
    
    outFile<<t0<< "\t"<<Mat(1,2).real()<<endl;
    t0 = t0+dt;
  }

   outFile.close();   

  return 0;
}

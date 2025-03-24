#ifndef QTERMO_H
#define QTERMO_H

#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions> 

using namespace std;
using namespace Eigen;

class Qtermo{
    private:
    MatrixXcd Hamiltonian;
    int Hilbert;
    public:
    //constructor declaration 
    Qtermo(const MatrixXcd &,int);
    MatrixXcd anticonmutator(const MatrixXcd& ,const MatrixXcd&);
    VectorXcd reshape(MatrixXcd& );
    MatrixXcd reshapeback(VectorXcd&,int);
    MatrixXcd isolatedevol(const MatrixXcd& );  
    MatrixXcd evol(const MatrixXcd&, double,int);

    
};

#endif
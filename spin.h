#ifndef SPIN_H
#define SPIN_H

#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions> 

using namespace Eigen;
using namespace std;

class spin{

    private:
    int numspin;
    double S;
    int d;//dimension (2s+1)
    MatrixXcd Sx(d,d);
    MatrixXcd Sy(d,d);
    MatrixXcd Sz(d,d);
    // Private helper functions
    void createSz();
    void createup();
    void createdown();
    void createSxSy();

    public:
    spin(int,double,int);    
    MatrixXcd tensor(MatrixXcdd, int);
    MatrixXcd Szf(int);
    MatrixXcd Syf(int);
    MatrixXcd Sxf(int);


};




#endif
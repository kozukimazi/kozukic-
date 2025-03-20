#ifndef JORDANWIGNER_H
#define JORDANWIGNER_H

#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions> 

using Eigen::MatrixXd;
using namespace std;


class jordanwigner{

    private:
    int hilbert;
    public:
    void setdim(int);
    Eigen::MatrixXcd tensor(Eigen::Matrix2d, int);
    Eigen::MatrixXcd opcreator(int);
    

};

#endif

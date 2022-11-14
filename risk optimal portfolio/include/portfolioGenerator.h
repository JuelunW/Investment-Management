#pragma once
#include "std_lib_facilities.h"

class portfolioGenerator {
private:
    //the number of stocks
    const int n = 7;
    //weights contribution
    vector<double> weights;
    //shares contribution
    vector<double> shares;
    //stock price
    vector<vector<double>> price;
    //daily value of portfolio
    vector<double> dailyVoP;
    //inFile
    void loadFile();
    //random weight
    void setRandomWeight();
    //set shares
    void setShares();
    //get portfolio daily returns
    void setDailyVoP();
    //
    //yearly log returns
    vector<vector<double>> dailyReturn();
    //mean yearly log returns
    vector<double> meanReturn();
    //yearly covariance matrix
    vector<vector<double>> cov();
    //
    vector<vector<double>> getCofactor(vector<vector<double>>& A, int p, int q, int nn) const;
    double determinant(vector<vector<double>>& A, int nn);
    vector<vector<double>> adjoint(vector<vector<double>>& A);
    vector<vector<double>> inverse(vector<vector<double>>& A);
    //
    //multiply two matrix
    vector<vector<double>> matrixM(vector<vector<double>> A, vector<vector<double>> B);
    //matrix * vector
    vector<double> matrixM(vector<vector<double>> A, vector<double> B);
    //vector * matrix
    vector<double> matrixM(vector<double> A, vector<vector<double>> B);
    //vector * vector
    vector<vector<double>> matrixM(vector<double> A, vector<double> B);
    double matrixMn(vector<double> A, vector<double> B);
    //vector * number
    vector<double> matrixM(vector<double> A, double B);
    //subtraction of two matrix
    vector<vector<double>> matrixS(vector<vector<double>> A, vector<vector<double>> B);
    //Addition of two vector
    vector<double> matrixP(vector<double> A, vector<double> B);

public:
    //initial capital
    double cap = 1000000;
    //random weight ctor
    portfolioGenerator();
    //optimal
    explicit portfolioGenerator(int a);
    //other weight ctor
    explicit portfolioGenerator(vector<double>& w);
    //get daily value of portfolio
    vector<double> getDailyVoP() const;
    //get weights
    void getWeights() const;
};

ostream& operator<< (ostream& os, const portfolioGenerator& p);
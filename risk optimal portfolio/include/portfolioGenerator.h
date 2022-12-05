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
    vector<vector<double>> matrixM(vector<vector<double>> matrixA, vector<vector<double>> matrixB);
    //matrix * vector
    vector<double> matrixM(vector<vector<double>> matrixA, vector<double> verticalB);
    //vector * matrix
    vector<double> matrixM(vector<double> horizontalA, vector<vector<double>> matrixB);
    //vector * vector
    vector<vector<double>> matrixM(vector<double> verticalA, vector<double> horizontalB);
    double matrixMn(vector<double> horizontalA, vector<double> verticalB);
    //vector * number
    vector<double> matrixM(vector<double> verticalA, double B);
    //subtraction of two matrix
    vector<vector<double>> matrixS(vector<vector<double>> matrixA, vector<vector<double>> matrixB);
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
    //get draw-down
    vector<double> getDrawDown() const;
    //get maximum draw-down
    void getMDD() const;
};

ostream& operator<< (ostream& os, const portfolioGenerator& p);
#include "std_lib_facilities.h"
#include "portfolioGenerator.h"

int main()
{
    //test case 1
    vector<double> w1;
    w1.assign(7, 1/7.0);
    portfolioGenerator p1(w1);
    cout << p1;

    //test case 2
    vector<double> w2;
    w2.assign(7, 0);
    w2[0] = 1;
    portfolioGenerator p2(w2);
    cout << p2;


    //random weight
    portfolioGenerator randomPortfolio;
    cout << randomPortfolio;

    //My strategy is maximum sharp ratio
    portfolioGenerator optimalPortfolio(1);
    cout << "The weight of risk optimal portfolio is:" << endl;
    optimalPortfolio.getWeights();
    cout << optimalPortfolio;
    //outFile the daily values of optimal portfolio
    vector<double> out = optimalPortfolio.getDailyVoP();
    ofstream outFile;
    outFile.open("../../output/optimalPortfolio.csv");

    outFile << "daily portfolio value" << endl;
    for(auto i:out){
        outFile << i << endl;
    }

    outFile.close();
    optimalPortfolio.getMDD();
    return 0;
}
#include "std_lib_facilities.h"
#include "portfolioGenerator.h"

int main()
{
    //test case 1
    vector<double> w1;
    w1.assign(7, 1/7.0);
    portfolioGenerator p1(w1);

    vector<double> tmp1 = p1.getDailyVoP();
    cout << "The profit of p1 is\n$" << tmp1.back() - tmp1.front() << endl;
    cout << (tmp1.back()/tmp1.front() - 1)*100 << "%" << endl;
    cout << endl;


    //test case 2
    vector<double> w2;
    w2.assign(7, 0);
    w2[0] = 1;
    portfolioGenerator p2(w2);

    vector<double> tmp2 = p2.getDailyVoP();
    cout << "The profit of p2 is\n$" << tmp2.back() - tmp2.front() << endl;
    cout << (tmp2.back()/tmp2.front() - 1)*100 << "%" << endl;
    cout << endl;


    //random weight
    portfolioGenerator randomPortfolio;
    //cout << p1;
    vector<double> tmp = randomPortfolio.getDailyVoP();
    cout << "The profit of random weight is\n$" << tmp.back() - tmp.front() << endl;
    cout << (tmp.back()/tmp.front() - 1)*100 << "%" << endl;
    cout << endl;

    //My strategy is maximum sharp ratio
    portfolioGenerator optimalPortfolio(1);
    cout << "The weight of risk optimal portfolio is:" << endl;
    optimalPortfolio.getWeights();

    //outFile the daily values of optimal portfolio
    vector<double> out = optimalPortfolio.getDailyVoP();

    ofstream outFile;
    outFile.open("../../output/optimalPortfolio.csv");
    for(auto i:out){
        outFile << i << endl;
    }
    outFile.close();

    return 0;
}
#include "portfolioGenerator.h"
#include "std_lib_facilities.h"
portfolioGenerator::portfolioGenerator(){
    loadFile();
    setRandomWeight();
    setShares();
    setDailyVoP();
}
portfolioGenerator::portfolioGenerator(int a){
    loadFile();

    //setOptimalWeight();
    vector<double> mu = meanReturn();
    vector<vector<double>> Sigma = cov();

    vector<double> l; //should be vertical
    l.assign(n, 1);
    vector<vector<double>> I;//diagonal
    for(int i = 0; i < n; i++){
        vector<double> tmp;
        for(int j = 0; j < n; j++){
            if(j == i) tmp.push_back(1);
            else tmp.push_back(0);
        }
        I.push_back(tmp);
    }


    vector<double> W_0 = matrixM(matrixM(inverse(Sigma), l),
                                 (1/matrixMn(matrixM(l, inverse(Sigma)), l)));

    vector<vector<double>> B = matrixM(inverse(Sigma), matrixS(I, matrixM(l, W_0)));
    vector<double> W_1 = matrixM(B, mu);

    vector<vector<double>> ds;
    double compare = 0;
    int aa;
    for(int A = 1; A <= 100; A++){
        vector<double> tmp;
        tmp.push_back(A);

        vector<double> W_A = matrixP(W_0, matrixM(W_1, (1/A)));
        for(auto i:W_A) tmp.push_back(i);

        double mu_A = matrixMn(W_A, mu);
        double sig_A = sqrt(matrixMn(W_A, matrixM(Sigma, W_A)));
        tmp.push_back(mu_A);
        tmp.push_back(sig_A);

        double SR = mu_A/sig_A;//assume risk-free rate = 0%
        tmp.push_back(SR);
        if(SR > compare){
            compare = SR;
            aa = A - 1;
        }

        ds.push_back(tmp);
    }

    for(int i = 1; i < n + 1; i++){
        weights.push_back(ds[aa][i]);
    }
    //
    //
    setShares();
    setDailyVoP();
}
portfolioGenerator::portfolioGenerator(vector<double>& w){
    loadFile();
    weights.clear();
    for(auto i:w) weights.push_back(i);
    setShares();
    setDailyVoP();
}

void portfolioGenerator::loadFile(){
    ifstream inFile;
    inFile.open("../../input/Book1.csv");
    double dailyPrice;

    inFile >> dailyPrice;
    inFile.clear();
    inFile.ignore(200, '\n');

    while (!inFile.eof()) {
        vector<double> tmp;
        for(int i = 0; i < n; i++){
            inFile >> dailyPrice;
            tmp.push_back(dailyPrice);
        }
        price.push_back(tmp);
    }
}

void portfolioGenerator::setRandomWeight(){
    random_device seed;
    uniform_real_distribution<double> uniform(0, 1);

    weights.clear();
    for(int i = 0; i < n; i++) weights.push_back(uniform(seed));
    double sum = 0;
    for(int i = 0; i < n; i++) sum += weights[i];
    for(int i = 0; i < n; i++) weights[i] /= sum;
}

void portfolioGenerator::setShares(){
    shares.clear();
    for(int i = 0; i < n; i++) shares.push_back(cap*weights[i]/price[0][i]);
}

void portfolioGenerator::setDailyVoP(){
    dailyVoP.clear();
    for(int row = 0; row < price.size() - 1; row++){
        double tmp = 0;
        for(int col = 0; col < n; col++) tmp += shares[col]*price[row][col];
        dailyVoP.push_back(tmp);
    }
}

vector<double> portfolioGenerator::getDailyVoP() const {
    return dailyVoP;
}

void portfolioGenerator::getWeights() const{
    for(auto i:weights) cout << i << endl;
}

ostream& operator<<(ostream& os, const portfolioGenerator& p)
{
    vector<double> tmp = p.getDailyVoP();
    os << "The profit of random weight is\n$" << tmp.back() - tmp.front() << endl;
    os << (tmp.back()/tmp.front() - 1)*100 << "%" << endl;
    os << endl;
    return os;
}

//optimal portfolio
vector<vector<double>> portfolioGenerator::dailyReturn() {
    vector<vector<double>> dailyReturns;
    for(int row = 0; row < price.size() - 1; row++){
        vector<double> tmp;
        for(int col = 0; col < n; col++) tmp.push_back(log(price[row + 1][col]/price[row][col]));
        dailyReturns.push_back(tmp);
    }
    return dailyReturns;
}

vector<double> portfolioGenerator::meanReturn() {
    vector<vector<double>> ret = dailyReturn();
    vector<double> result;
    for(int i = 0; i < n; i++){
        double sum = 0;
        for(int j = 0; j < ret.size(); j++){
            sum += ret[j][i];
        }
        result.push_back(252*sum/ret.size());
    }
    return result;
}

vector<vector<double>> portfolioGenerator::cov(){
    vector<vector<double>> ret = dailyReturn();
    vector<double> meanReturns = meanReturn();
    vector<vector<double>> result;

    for(int row = 0; row < n; row++){
        vector<double> tmp;
        for(int col = 0; col < n; col++){
            double tt = 0;
            for(int i = 0; i < ret.size(); i++){
                tt += (meanReturns[row] - ret[i][row])*(meanReturns[col] - ret[i][col])
                        /ret.size();
            }
            tmp.push_back(tt);
        }
        result.push_back(tmp);
    }
    return result;
}

vector<vector<double>> portfolioGenerator::getCofactor(vector<vector<double>>& A, int p, int q, int nn) const {
    vector<vector<double>> result(n, vector<double>(n));
    int i = 0, j = 0;
    // Looping for each element of the matrix
    for (int row = 0; row < nn; row++) {
        for (int col = 0; col < nn; col++) {
            //  Copying into temporary matrix only those
            //  element which are not in given row and
            //  column
            if (row != p && col != q) {
                result[i][j++] = A[row][col];

                // Row is filled, so increase row index and
                // reset col index
                if (j == nn - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
    return result;
}

double portfolioGenerator::determinant(vector<vector<double>>& A, int nn) {
    double D = 0; // Initialize result

    //  Base case : if matrix contains single element
    if (nn == 1) return A[0][0];

    double sign = 1; // To store sign multiplier

    // Iterate for each element of first row

    for (int f = 0; f < nn; f++) {
        // Getting Cofactor of A[0][f]
        vector<vector<double>> temp = getCofactor(A, 0, f, nn);
        D += sign * A[0][f] * determinant(temp, nn - 1);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

vector<vector<double>> portfolioGenerator::adjoint(vector<vector<double>>& A){
    vector<vector<double>> adj(n, vector<double> (n));

    // temp is used to store cofactors of A[][]
    int sign = 1;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Get cofactor of A[i][j]
            vector<vector<double>> temp = getCofactor(A, i, j, n);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign) * (determinant(temp, n - 1));
        }
    }
    return adj;
}

vector<vector<double>> portfolioGenerator::inverse(vector<vector<double>>& A){
    // Find determinant of A[][]
    double det = determinant(A, n);
    if (det == 0) {
        cout << "Singular matrix, can't find its inverse";
        exit(1);
    }

    // Find adjoint
    vector<vector<double>> adj = adjoint(A);

    vector<vector<double>> inverse(n, vector<double> (n));
    // Find Inverse using formula "inverse(A) =
    // adj(A)/det(A)"
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inverse[i][j] = adj[i][j] / det;

    return inverse;
}

vector<vector<double>> portfolioGenerator::matrixM(vector<vector<double>> A, vector<vector<double>> B){
    vector<vector<double>> result;
    vector<double> tmp;

    for (int row = 0; row < A.size(); row++){
        for (int col = 0; col < B[0].size(); col++){
            double tt = 0;
            for(int i = 0; i < B.size(); i++) tt += A[row][i] * B[i][col];
            tmp.push_back(tt);
        }
        result.push_back(tmp);
        tmp.clear();
    }

    return result;
}
vector<double> portfolioGenerator::matrixM(vector<vector<double>> A, vector<double> B){
    vector<double> result;
    for (int row = 0; row < A.size(); row++){
        double tt = 0;
        for (int col = 0; col < B.size(); col++){
            tt += A[row][col] * B[col];
        }
        result.push_back(tt);
    }
    return result;
}
vector<double> portfolioGenerator::matrixM(vector<double> A, vector<vector<double>> B){
    vector<double> result;
    for (int col = 0; col < A.size(); col++){
        double tt = 0;
        for (int row = 0; row < B.size(); row++){
            tt += A[row]* B[row][col];
        }
        result.push_back(tt);
    }
    return result;
}
vector<vector<double>> portfolioGenerator::matrixM(vector<double> A, vector<double> B){
    vector<vector<double>> result;
    for (int row = 0; row < A.size(); row++){
        vector<double> tt;
        for (int col = 0; col < B.size(); col++) tt.push_back(A[row]* B[col]);
        result.push_back(tt);
    }
    return result;
}
double portfolioGenerator::matrixMn(vector<double> A, vector<double> B){
    double tt = 0;
    for (int i = 0; i < A.size(); i++) tt += A[i] * B[i];
    return tt;
}
vector<double> portfolioGenerator::matrixM(vector<double> A, double B){
    vector<double> result;

    for (int i = 0; i < A.size(); i++){
        double tt = A[i] * B;
        result.push_back(tt);
    }

    return result;
}

vector<vector<double>> portfolioGenerator::matrixS(vector<vector<double>> A, vector<vector<double>> B){
    vector<vector<double>> result;
    vector<double> tmp;

    for (int row = 0; row < A.size(); row++){
        for (int col = 0; col < B[0].size(); col++) tmp.push_back(A[row][col] - B[row][col]);
        result.push_back(tmp);
        tmp.clear();
    }
    return result;
}
vector<double> portfolioGenerator::matrixP(vector<double> A, vector<double> B){
    vector<double> tmp;
    for (int i = 0; i < A.size(); i++) tmp.push_back(A[i] + B[i]);
    return tmp;
}
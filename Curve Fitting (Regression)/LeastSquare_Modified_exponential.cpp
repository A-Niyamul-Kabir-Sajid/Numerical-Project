// Modified Exponential Regression: T = a + b * e^(t/4)
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

// Function to compute regression coefficients and predict value
//T=a+b*e^(t/4)=y;
//Y=a+b*f(t)    
//Y=A + bX
//Y=y, X=e^(t/4), A=a, b=b
pair<double,double> modifiedExponentialRegression(const vector<double>& t, const vector<double>& T, double predictT) {
    int n = t.size();
    if(n < 2) {
        cout << "At least 2 data points are required.\n";
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }

    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;

    for(int i = 0; i < n; i++) {
        double X = exp(t[i]/4.0); // Transform t to X
        Sx += X;
        Sy += T[i];
        Sxx += X * X;
        Sxy += X * T[i];
    }

    double denom = n*Sxx - Sx*Sx;
    if(fabs(denom) < 1e-12) {
        cout << "Cannot perform regression (all transformed X values may be identical).\n";
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }

    double b = (n*Sxy - Sx*Sy)/denom;
    double a = (Sy - b*Sx)/n;

    double predictedY = a + b * exp(predictT/4.0);

    cout << fixed << setprecision(6);
    cout << "Regression equation: T = " << a << " + " << b << " * e^(t/4)\n";
    cout << "Predicted value at t = " << predictT << " is T = " << predictedY << "\n";

    return {a, b};
}

int main() {
    freopen("input/LeastSquare_Modified_exponential.txt", "r", stdin);
    freopen("output/LeastSquare_Modified_exponential_output.txt", "w", stdout);

    int n;
    //cout << "Enter number of data points: ";
    if(!(cin >> n) || n < 2) return 0;

    vector<double> t(n), T(n);
    //cout << "Enter data points (t T) one per line:\n";
    for(int i = 0; i < n; i++) {
        cin >> t[i] >> T[i];
    }

    double predictT;
    //cout << "Enter the value of t to predict T: ";
    if(!(cin >> predictT)) return 0;

    modifiedExponentialRegression(t, T, predictT);
    
    return 0;
}
/*
Input:
5
0 2
1 2.7
2 3.7
3 5.0
4 7.4
5

Output:
Regression equation: T = -1.294980 + 3.110722 * e^(t/4)
Predicted value at t = 5.000000 is T = 9.562507
*/
// Simple Linear Regression (y = a + b*x)
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

// Function to compute linear regression coefficients a and b
pair<double,double> LinearRegression(const vector<double>& x, const vector<double>& y, double predictX) {
    int n = x.size();
    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;

    for (int i = 0; i < n; i++) {
        Sx += x[i];
        Sy += y[i];
        Sxx += x[i] * x[i];
        Sxy += x[i] * y[i];
    }

    double denom = n * Sxx - Sx * Sx;
    if (fabs(denom) < 1e-12) {  // Prevent division by zero
        cout << "Cannot perform linear regression (all x values might be the same).\n";
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }

    double b = (n * Sxy - Sx * Sy) / denom;
    double a = (Sy - b * Sx) / n;

    cout << fixed << setprecision(6);
    cout << "Regression Line: y = " << a << " + " << b << "*x\n";
    cout << "Predicted value at x = " << predictX << " is y = " << (a + b * predictX) << "\n";

    return {a, b};
}

int main() {
    freopen("input/LeastSquare_Linear.txt", "r", stdin);
    freopen("output/LeastSquare_Linear_output.txt", "w", stdout);

    int n;
    //cout << "Enter number of data points: ";
    if (!(cin >> n) || n < 2) {
        cout << "Invalid input. At least 2 data points are required.\n";
        return 0;
    }

    vector<double> x(n), y(n);
    //cout << "Enter the data points (x y) one per line:\n";
    for (int i = 0; i < n; i++) {
        cin >> x[i] >> y[i];
    }

    double predictX;
    //cout << "Enter the value of x to predict y: ";
    if (!(cin >> predictX)) {
        cout << "Invalid input for prediction.\n";
        return 0;
    }

    LinearRegression(x, y, predictX);
    
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
Regression Line: y = 1.540000 + 1.310000*x
Predicted value at x = 5.000000 is y = 8.090000
*/
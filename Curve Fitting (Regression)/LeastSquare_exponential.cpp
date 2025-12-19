#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

//P = P0 * e^(kt);
//ln P = ln P0 + kt
//Y = A + bX
//Y = ln P, X = t, A = ln P0, b = k
pair<double, double> exponential(const vector<double>& x, const vector<double>& y, double predictX) {
    int n = x.size();
    if (n < 2) {
        cout << "At least 2 data points are required.\n";
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }

    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;

    for (int i = 0; i < n; i++) {
        if (y[i] <= 0) {
            cout << "Error: All Y values must be positive for logarithm.\n";
            return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
        }
        double ly = log(y[i]);
        Sx += x[i];
        Sy += ly;
        Sxx += x[i] * x[i];
        Sxy += x[i] * ly;
    }

    double denom = n * Sxx - Sx * Sx;
    if (fabs(denom) < 1e-12) {
        cout << "Cannot perform regression (denominator too small).\n";
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }

    double b = (n * Sxy - Sx * Sy) / denom;
    double A = (Sy - b * Sx) / n;
    double a = exp(A);

    double predictedY = a * exp(b * predictX);

    cout << fixed << setprecision(6);
    cout << "Regression equation: P = " << a << " * e^(" << b << " * t)\n";
    cout << "Predicted value at t = " << predictX << " is P = " << predictedY << endl;

    return {a, b};
}

int main() {
    freopen("input/LeastSquare_exponential.txt", "r", stdin);
    freopen("output/LeastSquare_exponential_output.txt", "w", stdout);

    int n;
    //cout << "Enter number of data points: ";
    if (!(cin >> n) || n < 2) return 0;

    vector<double> x(n), y(n);
    //cout << "Enter data points (t P) one per line:\n";
    for (int i = 0; i < n; i++) {
        cin >> x[i] >> y[i];
    }

    double predictX;
    //cout << "Enter the value of t to predict P: ";
    if (!(cin >> predictX)) return 0;

    exponential(x, y, predictX);
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
Regression equation: P = 1.963168 * e^(0.323285 * t)
Predicted value at t = 5.000000 is P = 9.884674*/
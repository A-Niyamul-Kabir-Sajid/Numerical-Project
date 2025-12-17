// Newton Divided Difference Interpolation Method (Unequally Spaced)
// Includes error estimation

#include <bits/stdc++.h>
using namespace std;

// Divided Difference Interpolation
double dividedDifference(const vector<double>& x, const vector<double>& y, double valX, double& error) {
    int n = x.size();
    vector<vector<double>> table(n, vector<double>(n, 0.0));

    // Fill first column
    for(int i=0; i<n; i++) table[i][0] = y[i];

    // Build divided difference table
    for(int j=1; j<n; j++)
        for(int i=0; i<n-j; i++)
            table[i][j] = (table[i+1][j-1] - table[i][j-1]) / (x[i+j] - x[i]);

    // Evaluate
    double result = table[0][0];
    double xTerm = 1.0;
    for(int i=1; i<n; i++) {
        xTerm *= (valX - x[i-1]);
        result += table[0][i] * xTerm;
    }

    // Error estimation: use last term of table if exists
    error = fabs(table[0][n-1] * xTerm);

    return result;
}

int main() {
    int n;
    cout << "Enter number of data points: ";
    cin >> n;
    vector<double> x(n), y(n);
    cout << "Enter data points (x y):\n";
    for(int i=0; i<n; i++) cin >> x[i] >> y[i];

    double valX;
    cout << "Enter value to interpolate: ";
    cin >> valX;

    double error;
    double valY = dividedDifference(x, y, valX, error);
    cout << "Interpolated value at x = " << valX << " is y = " << valY << endl;
    cout << "Estimated interpolation error: " << error << endl;

    return 0;
}

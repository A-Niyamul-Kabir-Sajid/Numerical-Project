// Newton Forward Interpolation Method (Equally Spaced)
// Includes error estimation

#include <bits/stdc++.h>
using namespace std;

// Print vector differences (optional for debugging)
void printDifferences(const vector<double>& v, int order) {
    cout << "Order " << order << " differences: ";
    for (auto val : v) cout << val << " ";
    cout << endl;
}

// Newton Forward Interpolation
double newtonForward(const vector<double>& x, const vector<double>& y, double valX, double& error) {
    int n = x.size();
    double h = x[1] - x[0];
    double u = (valX - x[0]) / h;

    vector<double> delY = y;
    double result = y[0];
    double uTerm = 1.0;
    double factorial = 1.0;

    // Build differences
    for (int i = 1; i < n; i++) {
        uTerm *= (u - (i - 1));
        factorial *= i;

        vector<double> newDelY;
        for (int j = 0; j < delY.size() - 1; j++)
            newDelY.push_back(delY[j+1] - delY[j]);

        delY = newDelY;
        if(delY.empty()) break;

        result += (uTerm * delY[0]) / factorial;
    }

    // Error estimation: use next term if exists
    if (!delY.empty()) error = fabs((uTerm * delY[0]) / factorial);
    else error = 0.0;

    return result;
}

int main() {
    int n;
    //cout << "Enter number of equally spaced data points: ";
    cin >> n;
    vector<double> x(n), y(n);
    //cout << "Enter data points (x y):\n";
    for(int i=0; i<n; i++) cin >> x[i] >> y[i];

    double valX;
    //cout << "Enter value to interpolate: ";
    cin >> valX;

    double error;
    double valY = newtonForward(x, y, valX, error);
    cout << "Interpolated value at x = " << valX << " is y = " << valY << endl;
    cout << "Estimated interpolation error: " << error << endl;

    return 0;
}
/*
Input:
4
1 1
2 4
3 9
4 16
2.5
Output:
Interpolated value at x = 2.5 is y = 6.25
Estimated interpolation error: 0
*/
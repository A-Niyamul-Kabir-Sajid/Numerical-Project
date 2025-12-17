// Polynomial Least Squares Regression with Gaussian Elimination
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

int main() {
    int n;
    cout << "Enter number of data points: ";
    cin >> n;
    if (n <= 1) {
        cout << "At least two data points are required.\n";
        return 0;
    }

    vector<double> x(n), y(n);
    cout << "Enter the data points (x y) one per line:\n";
    for (int i = 0; i < n; i++) {
        cin >> x[i] >> y[i];
    }

    int m;
    cout << "Enter degree of polynomial (m < n): ";
    cin >> m;
    if (m >= n) {
        cout << "Polynomial degree must be less than number of data points.\n";
        return 0;
    }

    double TheNewX;
    cout << "Enter value of x to predict y: ";
    cin >> TheNewX;

    // Construct normal equations A * coeffs = B
    int size = m + 1;
    vector<vector<double>> A(size, vector<double>(size, 0));
    vector<double> B(size, 0);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < n; k++) {
                A[i][j] += pow(x[k], i + j);
            }
        }
        for (int k = 0; k < n; k++) {
            B[i] += y[k] * pow(x[k], i);
        }
    }

    // Gaussian elimination with partial pivoting
    for (int i = 0; i < size; i++) {
        // Partial pivoting
        int maxRow = i;
        for (int k = i + 1; k < size; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) maxRow = k;
        }
        if (fabs(A[maxRow][i]) < 1e-12) {
            cout << "Matrix is singular or nearly singular. Cannot solve.\n";
            return 0;
        }
        if (maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(B[i], B[maxRow]);
        }

        // Elimination
        for (int j = i + 1; j < size; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < size; k++) {
                A[j][k] -= factor * A[i][k];
            }
            B[j] -= factor * B[i];
        }
    }

    // Back substitution
    vector<double> coeffs(size, 0);
    for (int i = size - 1; i >= 0; i--) {
        coeffs[i] = B[i];
        for (int j = i + 1; j < size; j++) {
            coeffs[i] -= A[i][j] * coeffs[j];
        }
        coeffs[i] /= A[i][i];
    }

    // Output coefficients
    cout << "\nPolynomial coefficients:\n";
    for (int i = 0; i < size; i++) {
        cout << "Coefficient of x^" << i << " = " << fixed << setprecision(6) << coeffs[i] << endl;
    }

    // Predict value at TheNewX
    double TheNewY = 0.0;
    for (int i = 0; i < size; i++) {
        TheNewY += coeffs[i] * pow(TheNewX, i);
    }

    cout << "\nPredicted value at x = " << TheNewX << " is y = " << fixed << setprecision(6) << TheNewY << endl;

    return 0;
}

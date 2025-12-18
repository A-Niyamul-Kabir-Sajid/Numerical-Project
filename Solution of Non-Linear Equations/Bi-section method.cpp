#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

// Constants for numerical stability
const double EPS = 1e-9;      // Precision for root finding
const double STEP = 0.5;     // Step size to scan for sign changes
const double RANGE = 100.0;  // Range to scan (-100 to 100)
const int MAX_ITER = 100;    // Safety cap for bisection iterations

// Evaluate polynomial using Horner's Method: O(n)
double evaluate(const vector<double>& poly, double x) {
    double res = 0;
    for (double coeff : poly)
        res = res * x + coeff;
    return res;
}

// Function to print the polynomial nicely
void printPolynomial(const vector<double>& poly) {
    int n = poly.size() - 1;
    cout << "f(x) = ";
    for (int i = 0; i <= n; i++) {
        if (poly[i] == 0) continue;
        if (poly[i] > 0 && i != 0) cout << " + ";
        if (poly[i] < 0) cout << " - ";
        
        double val = abs(poly[i]);
        if (val != 1 || i == n) cout << val;
        
        int power = n - i;
        if (power > 1) cout << "x^" << power;
        else if (power == 1) cout << "x";
    }
    cout << endl;
}

int main() {
    int degree;
    //cout << "Enter the degree of the polynomial: ";
    if (!(cin >> degree)) return 1;

    vector<double> coeffs(degree + 1);
    //cout << "Enter coefficients (from highest power to constant):" << endl;
    for (int i = 0; i <= degree; i++) {
        cin >> coeffs[i];
    }

    printPolynomial(coeffs);

    vector<double> roots;

    // 1. Scanning Phase: Look for intervals [a, b] where f(a) and f(b) have different signs
    for (double a = -RANGE; a < RANGE; a += STEP) {
        double b = a + STEP;
        double fa = evaluate(coeffs, a);
        double fb = evaluate(coeffs, b);

        // Check if 'a' itself is a root to avoid missing it
        if (abs(fa) < EPS) {
            roots.push_back(a);
            continue;
        }

        // Bisection condition: signs must be opposite
        // Using (fa > 0) != (fb > 0) is safer than fa * fb < 0 to prevent underflow
        if ((fa > 0) != (fb > 0)) {
            double low = a;
            double high = b;

            for (int i = 0; i < MAX_ITER; i++) {
                double mid = low + (high - low) / 2.0;
                double fmid = evaluate(coeffs, mid);

                if (abs(fmid) < EPS || (high - low) / 2.0 < EPS) {
                    roots.push_back(mid);
                    break;
                }

                // Check which sub-interval contains the sign change
                if ((evaluate(coeffs, low) > 0) != (fmid > 0)) {
                    high = mid;
                } else {
                    low = mid;
                }
            }
        }
    }

    // 2. Cleaning Phase: Sort and remove duplicates caused by overlapping intervals
    if (roots.empty()) {
        cout << "No real roots found in the range [-" << RANGE << ", " << RANGE << "]." << endl;
    } else {
        sort(roots.begin(), roots.end());
        
        // Remove duplicates where roots are closer than a threshold
        auto it = unique(roots.begin(), roots.end(), [](double x, double y) {
            return abs(x - y) < 1e-5; 
        });
        roots.erase(it, roots.end());

        cout << "Real roots found:" << endl;int i=1;
        for (double r : roots) {
            // Clean up near-zero results like -0.000000
            if (abs(r) < 1e-10) r = 0.0;
            cout << "x"<<i<<" = " << fixed << setprecision(6) << r << endl;
            i++;
        }
    }

    return 0;
}


/*
no root
input 
2
1 0 1   
output

f(x) = x^2 + 1
No real roots found in the range [-100, 100].

*/

/*
root found
input
2
1 -3 2
output
f(x) = x^2 - 3x + 2
Real roots found:
x1 = 1.000000
x2 = 2.000000
*/




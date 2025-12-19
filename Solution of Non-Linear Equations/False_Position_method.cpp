// Full implementation of the False Position Method (Regula Falsi) for Polynomial Roots
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

// Numerical Constants
const double EPS = 1e-9;      // Tolerance for root precision
const int MAX_ITER = 100;    // Maximum iterations per bracket
const double STEP = 0.5;     // Scanning step to find sign changes
const double RANGE = 1000.0; // The x-axis range to search (-1000 to 1000)

// Function to evaluate polynomial using Horner's Method
double evaluate(const vector<double>& cof, double x) {
    double res = 0;
    for (double c : cof)
        res = res * x + c;
    return res;
}

int main() {
    freopen("input/False_Position_method_root_found.txt", "r", stdin);
    freopen("output/False_Position_method_root_found.txt", "w", stdout);

    int n;
    //cout << "Enter the degree of the polynomial: ";
    if (!(cin >> n)) return 0;

    // Polynomial coefficients from highest power to constant
    vector<double> cof(n + 1);
   // cout << "Enter the " << n + 1 << " coefficients: ";
    for (int i = 0; i <= n; i++) cin >> cof[i];

    vector<double> found_roots;

    // 1. SCANNING PHASE
    // We look for intervals [a, b] where the function crosses the x-axis
    for (double a = -RANGE; a < RANGE; a += STEP) {
        double b = a + STEP;
        double fa = evaluate(cof, a);
        double fb = evaluate(cof, b);

        // Check if the start of the interval is already a root
        if (abs(fa) < EPS) {
            found_roots.push_back(a);
            continue;
        }

        // If signs are different, a root exists in [a, b]
        if ((fa > 0) != (fb > 0)) {
            double x0 = a, x1 = b;
            
            // 2. FALSE POSITION ITERATION
            for (int i = 0; i < MAX_ITER; i++) {
                double f0 = evaluate(cof, x0);
                double f1 = evaluate(cof, x1);

                // Prevent division by zero if the function is very flat
                if (abs(f1 - f0) < 1e-15) break;

                // Regula Falsi Formula
                double xr = x1 - f1 * (x0 - x1) / (f0 - f1);
                double fr = evaluate(cof, xr);

                // Check for convergence
                if (abs(fr) < EPS || abs(x1 - x0) < EPS) {
                    found_roots.push_back(xr);
                    break;
                }

                // Update the bracket
                if ((f0 > 0) != (fr > 0)) x1 = xr;
                else x0 = xr;
            }
        }
    }

    // 3. CLEANING PHASE
    if (found_roots.empty()) {
        cout << "\nNo real roots found in the range." << endl;
    } else {
        // Sort the roots
        sort(found_roots.begin(), found_roots.end());
        
        // Remove duplicate detections caused by adjacent scanning intervals
        auto it = unique(found_roots.begin(), found_roots.end(), [](double r1, double r2) {
            return abs(r1 - r2) < 1e-5; 
        });
        found_roots.erase(it, found_roots.end());

        cout << "\nUnique real roots found:" << endl;int i = 1;
        for (double r : found_roots) {
            
            // Clean up visual noise like -0.000000
            if (abs(r) < 1e-10) r = 0.0;
            cout << "x"<<i<<" = " << fixed << setprecision(6) << r << endl;
            i++;
        }
    }
    
    return 0;
}

/*
No root
input:
2
1 0 1
output:
No real roots found in the range.

*/

/*
With root
input:
2
1 -3 2
output:
Unique real roots found:
x1 = 1.000000
x2 = 2.000000
*/
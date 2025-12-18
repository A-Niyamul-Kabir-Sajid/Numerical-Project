// Full implementation of the Newton-Raphson Method for Polynomial Roots
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-8;    // Convergence tolerance
const int MAX_ITER = 100;    // Max iterations per root search
const double STEP = 0.5;     // Step size for scanning initial guesses

// Print polynomial nicely
void printfunc(const vector<int>& cof) {
    int n = cof.size();
    bool first = true;
    for(int i = 0; i < n; i++) {
        if(cof[i] == 0) continue;
        
        if(!first && cof[i] > 0) cout << '+';
        if(cof[i] == -1 && n-i-1 != 0) cout << '-';
        else if(cof[i] != 1 || n-i-1 == 0) cout << cof[i];

        int power = n-i-1;
        if(power > 1) cout << "x^" << power;
        else if(power == 1) cout << "x";
        first = false;
    }
    if(first) cout << "0";
    cout << endl;
}

// Evaluate polynomial at x using Horner's method
double func(const vector<int>& cof, double x) {
    double res = 0;
    for(int i = 0; i < cof.size(); i++)
        res = res * x + cof[i];
    return res;
}

// Evaluate derivative of polynomial at x using Horner's method
double deriv(const vector<int>& cof, double x) {
    double res = 0;
    int n = cof.size();
    for(int i = 0; i < n - 1; i++)
        res = res * x + (double)cof[i] * (n - i - 1);
    return res;
}

int main() {
    int n;
    if (!(cin >> n)) return 0;

    vector<int> cof(n + 1);
    for(int i = 0; i <= n; i++)
        cin >> cof[i];

    printfunc(cof);

    vector<double> roots;
    double start = -1000, end = 1000;

    // Scan the range to find initial guesses for Newton-Raphson
    for(double i = start; i <= end; i += STEP) {
        double x0 = i;
        int iter = 0;

        while(iter < MAX_ITER) {
            double f = func(cof, x0);
            double df = deriv(cof, x0);

            // Avoid division by zero
            if(fabs(df) < 1e-12) break; 

            double x1 = x0 - f / df;

            // Check if the sequence has converged
            if(fabs(x1 - x0) < EPS) {
                // Verify that x1 is actually a root (Newton can converge to non-roots if it diverges)
                if(fabs(func(cof, x1)) < 1e-6) {
                    roots.push_back(x1);
                }
                break;
            }

            x0 = x1;
            iter++;
        }
    }

    if(roots.empty()) {
        cout << "No real roots found.\n";
    } else {
        // Sort roots to prepare for deduplication
        sort(roots.begin(), roots.end());
        
        // Remove duplicate roots found from different initial guesses
        auto it = unique(roots.begin(), roots.end(), [](double a, double b) {
            return fabs(a - b) < 1e-4; // Threshold for considering roots "the same"
        });
        roots.erase(it, roots.end());

        cout << "Roots found:\n";
        for(double r : roots) {
            int i=1;
            // Fix sign of zero if result is -0.000000
            if(fabs(r) < 1e-10) r = 0.0;
            cout <<"x"<<i<<" = "<< fixed << setprecision(6) << r <<endl;
            i++;
        }
        cout << endl;
    }

    return 0;
}

/*
No root
input:
2
1 0 1
output:
1x^2+1
No real roots found.

*/

/*
Root
input:
3
1 -6 11 -6

output:
x^3-6x^2+11x-6
Roots found:
1.000000 2.000000 3.000000 


*/
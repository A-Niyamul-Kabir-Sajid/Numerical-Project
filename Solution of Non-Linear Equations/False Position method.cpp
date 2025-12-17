#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-8;     // convergence tolerance
const int MAX_ITER = 100;     // max iterations per root
const double STEP = 0.5;      // interval step for scanning

// Print polynomial nicely
void printfunc(const vector<int>& cof) {
    int n = cof.size();
    for(int i = 0; i < n; i++) {
        if(cof[i] > 0 && i != 0) cout << '+';
        if(cof[i] != 0) {
            cout << cof[i];
            if(n-i-1 > 1) cout << "x^" << n-i-1;
            else if(n-i-1 == 1) cout << "x";
        }
    }
    cout << endl;
}

// Evaluate polynomial at x using Horner's method
double func(const vector<int>& cof, double x) {
    double res = 0;
    for(int i = 0; i < cof.size(); i++)
        res = res*x + cof[i];
    return res;
}

int main() {
    int n;
    cin >> n;

    vector<int> cof(n+1);
    for(int i = 0; i <= n; i++)
        cin >> cof[i];

    printfunc(cof);

    vector<double> roots;
    double start = -1000, end = 1000;

    // Scan intervals to detect sign changes (possible roots)
    for(double a = start; a < end; a += STEP) {
        double b = a + STEP;
        double fa = func(cof, a);
        double fb = func(cof, b);

        if(fa*fb > 0) continue; // no root in interval

        double x0 = a, x1 = b;
        int iter = 0;
        set<long long> seen;

        while(iter < MAX_ITER) {
            double f0 = func(cof, x0);
            double f1 = func(cof, x1);

            if(fabs(f1 - f0) < EPS) break; // avoid division by zero

            double xr = x1 - f1*(x0 - x1)/(f0 - f1); // false position formula

            if(fabs(func(cof, xr)) < EPS || seen.count(round(xr*1e6))) {
                roots.push_back(xr);
                break;
            }

            // Update interval
            if(f0*func(cof, xr) < 0) x1 = xr;
            else x0 = xr;

            seen.insert(round(xr*1e6));
            iter++;
        }
    }

    if(roots.empty()) {
        cout << "No real roots found.\n";
    } else {
        sort(roots.begin(), roots.end());
        cout << "Roots found:\n";
        for(double r : roots)
            cout << fixed << setprecision(6) << r << " ";
        cout << endl;
    }

    return 0;
}
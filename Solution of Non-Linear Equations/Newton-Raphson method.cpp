#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-8;   // convergence tolerance
const int MAX_ITER = 100;   // max iterations per root
const double STEP = 0.5;    // interval step for scanning

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

// Evaluate derivative of polynomial at x
double deriv(const vector<int>& cof, double x) {
    double res = 0;
    int n = cof.size();
    for(int i = 0; i < n-1; i++)
        res = res*x + cof[i]*(n-i-1);
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

    for(double i = start; i < end; i += STEP) {
        double x0 = i;
        int iter = 0;
        set<long long> seen;

        while(iter < MAX_ITER) {
            double f = func(cof, x0);
            double df = deriv(cof, x0);

            if(fabs(df) < EPS) break; // derivative too small, avoid division by zero

            double x1 = x0 - f/df;

            if(fabs(x1 - x0) < EPS || seen.count(round(x1*1e6))) {
                roots.push_back(x1);
                break;
            }

            seen.insert(round(x1*1e6));
            x0 = x1;
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
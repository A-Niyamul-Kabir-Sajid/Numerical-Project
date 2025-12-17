// Bisection Method for Finding Polynomial Roots

#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-8;     // Convergence tolerance
const int MAX_ITER = 100;     // Maximum iterations per root
const double STEP = 0.5;      // Interval step for scanning possible roots
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

    // Scan intervals for sign changes
    for(double a = start; a < end; a += STEP) {
        double b = a + STEP;
        double fa = func(cof, a);
        double fb = func(cof, b);

        if(fa*fb > 0) continue; // no root in this interval

        double left = a, right = b;
        int iter = 0;
        set<long long> seen; // Avoid duplicate roots

        // Bisection iterations
        while(iter < MAX_ITER && fabs(right - left) > EPS) {
            double mid = (left + right)/2.0;
            double fmid = func(cof, mid);

            if(fabs(fmid) < EPS || seen.count(round(mid*1e6))) {
                roots.push_back(mid);
                break;
            }

            if(fa * fmid < 0) right = mid;
            else {
                left = mid;
                fa = fmid; // update fa for next iteration
            }

            seen.insert(round(mid*1e6));
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
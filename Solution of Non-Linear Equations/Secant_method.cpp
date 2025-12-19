#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-8;   // convergence tolerance
const int MAX_ITER = 100;   // max iterations per root
const double STEP = 0.5;    // interval step for root scanning

// Print polynomial nicely
void printfunc(vector<int> cof) {
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

// Evaluate polynomial at x
double func(const vector<int>& cof, double x) {
    double ret = 0;
    int n = cof.size();
    for(int i = 0; i < n; i++)
        ret = ret*x + cof[i];  // Horner's method
    return ret;
}

int main() {
    freopen("input/Secant_method_root_found.txt", "r", stdin);
    freopen("output/Secant_method_root_found.txt", "w", stdout);

    int n;
    cin >> n;

    vector<int> cof(n+1);
    for(int i = 0; i <= n; i++)
        cin >> cof[i];

    printfunc(cof);

    vector<double> roots;
    double start = -1000;  // scanning range
    double end = 1000;

    // scan intervals for sign changes
    for(double i = start; i < end; i += STEP) {
        double f1 = func(cof, i);
        double f2 = func(cof, i + STEP);
        if(f1*f2 <= 0) {  // possible root
            double x0 = i, x1 = i+STEP;
            int iter = 0;
            set<long long> seen;  // to avoid duplicates

            while(iter < MAX_ITER) {
                double f0 = func(cof, x0);
                double f1val = func(cof, x1);

                if(fabs(f1val - f0) < EPS) break;  // avoid division by zero

                double xr = x1 - f1val*(x0 - x1)/(f0 - f1val);

                if(fabs(xr - x1) < EPS || seen.count(round(xr*1e6))) {
                    roots.push_back(xr);
                    break;
                }

                seen.insert(round(xr*1e6));
                x0 = x1;
                x1 = xr;
                iter++;
            }
        }
    }

    if(roots.empty()) {
        cout << "No real roots found.\n";
    } else {
        cout << "Roots found:\n";
        sort(roots.begin(), roots.end());
        for(double r : roots)
            cout << fixed << setprecision(6) << r << " ";
        cout << endl;
    }
    
    return 0;
}


/*
No Root
input:
2
1 0 1
output:
1x^2+1
No real roots found.
*/

/*
root
input:
2
1 0 -4
output:
1x^2-4
Roots found:
-2.000000 -2.000000 2.000000 2.000000 
*/
#include <bits/stdc++.h>
using namespace std;

// Polynomial evaluation
double f(vector<double>& a, double x)
{
    double s = 0;
    for(int i = 0; i < a.size(); i++)
        s += a[i] * pow(x, i);
    return s;
}
int main()
{
    freopen("input/Second_order_derivative.txt", "r", stdin);
    freopen("output/Second_order_derivative_output.txt", "w", stdout);
    
    int degree;
    //cout << "Enter degree of polynomial: ";
    cin >> degree;
    vector<double> a(degree + 1);
    //cout << "Enter coefficients a_n to a_0:\n";
    for(int i = degree; i >= 0; i--)
        cin >> a[i];
    double p;
    //cout << "Enter value of p: ";
    cin >> p;
    double h = 1.0;
    double fpp = (f(a, p + h) - 2*f(a, p) + f(a, p - h)) / (h * h);
    cout << "f\"(p) = " << fpp << endl;
    
    return 0;
}
/*
Input:
3
1 0 -1 1
2
Output:
f"(p) = 12
*/
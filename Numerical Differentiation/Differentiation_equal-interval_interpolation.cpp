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
    freopen("input/Differentiation_equal-interval_interpolation.txt", "r", stdin);
    freopen("output/Differentiation_equal-interval_interpolation_output.txt", "w", stdout);
    
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
    vector<double> y(degree + 1);
    for(int i = 0; i <= degree; i++)
        y[i] = f(a, i);
    // Forward differences
    vector<vector<double>> d(degree + 1, vector<double>(degree + 1));
    for(int i = 0; i <= degree; i++)
        d[i][0] = y[i];
    for(int j = 1; j <= degree; j++)
        for(int i = 0; i <= degree - j; i++)
            d[i][j] = d[i+1][j-1] - d[i][j-1];
    double u = (p - 0) / h;
    double fp = d[0][1];
    if(degree >= 2) fp += ((2*u - 1)/2.0) * d[0][2];
    if(degree >= 3) fp += ((3*u*u - 6*u + 2)/6.0) * d[0][3];
    fp /= h;
    cout << "f'(p) = " << fp << endl;
    
    return 0;
}
/*
Input:
2
1 3 2
1
Output:
f'(p)=5
*/

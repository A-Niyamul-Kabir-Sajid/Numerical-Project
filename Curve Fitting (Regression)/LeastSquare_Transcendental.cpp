#include<iostream>
#include<vector>
#include<cmath>
#include<iomanip>
#include<limits>

using namespace std;


//y=ax^b
//lny=ln a + b ln x
//Y=A + bX
//Y=ln y, X=ln x, A=ln a
pair<double,double> transcendental(const vector<double>& x,const vector<double>& y,double predictX){
    int n= x.size();
    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    for(int i = 0; i < n; i++) {
        if(x[i] <= 0 || y[i] <= 0) return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
        double lx = log(x[i]);
        double ly = log(y[i]);
        Sx  += lx;
        Sy  += ly;
        Sxx += lx * lx;
        Sxy += lx * ly;
    }
    double denom = n*Sxx - Sx*Sx;
    if(fabs(denom) < 1e-12) return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    double b= (n*Sxy-Sx*Sy)/denom;
    double A=(Sy-b*Sx)/n;
    double a=exp(A);
    cout << "Predicted value at x = " << predictX << " is y = ";
    cout << fixed << setprecision(6) << (a*pow(predictX,b)) << "\n";

    return {a,b};
}



int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if(!(cin>>n)) return 0;
    vector<double> x(n), y(n);
    for(int i = 0; i < n; i++) {
        cin>>x[i]>>y[i];
    }
    double predictX;
    if(!(cin>>predictX)) return 0;

    
    transcendental(x,y,predictX);
   
}
/*
Input:
4
1 2.7
2 3.7
3 5.0
4 7.4
5

Output:
Predicted value at x = 5 is y = 9.110219
*/
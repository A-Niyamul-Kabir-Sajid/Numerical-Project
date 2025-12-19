#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
double f(double x){
    return x*x;
}
int main(){
    freopen("input/2_Simpsons_three-eighths_rule.txt", "r", stdin);
    freopen("output/2_Simpsons_three-eighths_rule_output.txt", "w", stdout);
    
    double a,b;
    int n;
    //cout<<"Enter lower limit a: ";
    cin>>a;
    //cout<<"Enter upper limit b: ";
    cin>>b;
    //cout<<"Enter number of intervals (must be multiple of 3): ";
    cin>>n;
    double h=(b-a)/n;
    double sum=f(a)+f(b);
    for(int i=1;i<n;i++){
        double x=a+i*h;
        if(i%3==0) sum+=2*f(x);
        else sum+=3*f(x);
    }
    double result=(3*h/8)*sum;
    cout<<"\nSimpson 3/8 Result = "<<result<<endl;
}
/*
Input:
0 2 4
Output:
Simpson 3/8 Result = 2.29688

*/
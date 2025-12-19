#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
double f(double x){
    return x*x;
}
int main(){
    freopen("input/1_Simpsons_One-third_rule.txt", "r", stdin);
    freopen("output/1_Simpsons_One-third_rule_output.txt", "w", stdout);
    
    double a,b;
    int n;
    //cout<<"Enter lower limit a: ";
    cin>>a;
    //cout<<"Enter upper limit b: ";
    cin>>b;
    //cout<<"Enter number of intervals (must be even): ";
    cin>>n;
    double h=(b-a)/n;
    double sum=f(a)+f(b);
    for(int i=1;i<n;i++){
        double x=a+i*h;
        if(i%2==0) sum+=2*f(x);
        else sum+=4*f(x);
    }
    double result=(h/3)*sum;
    cout<<"\nSimpson 1/3 Result = "<<result<<endl;
}
/*
Input:
0 2 4
Output:
Simpson 1/3 Result = 2.66667
*/

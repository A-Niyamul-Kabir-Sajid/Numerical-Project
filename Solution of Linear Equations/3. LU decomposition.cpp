#include <bits/stdc++.h>
using namespace std;
int main(){
    int n;
    cout<<"Enter order of matrix: ";
    cin>>n;
    vector<vector<double>>A(n,vector<double>(n));
    vector<vector<double>>L(n,vector<double>(n,0));
    vector<vector<double>>U(n,vector<double>(n,0));
    vector<double>B(n),Y(n),X(n);
    const double EPS=1e-12;
    cout<<"Enter coefficient matrix A:\n";
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            cin>>A[i][j];
    cout<<"Enter constant matrix B:\n";
    for(int i=0;i<n;i++)
        cin>>B[i];
    for(int i=0;i<n;i++){
        for(int k=i;k<n;k++){
            double sum=0;
            for(int j=0;j<i;j++)
                sum+=L[i][j]*U[j][k];
            U[i][k]=A[i][k]-sum;
        }
        if(fabs(U[i][i])<EPS){
            cout<<"No unique solution exists (zero pivot found)\n";
            return 0;
        }
        for(int k=i;k<n;k++){
            if(i==k)
                L[i][i]=1;
            else{
                double sum=0;
                for(int j=0;j<i;j++)
                    sum+=L[k][j]*U[j][i];
                L[k][i]=(A[k][i]-sum)/U[i][i];
            }
        }
    }
    for(int i=0;i<n;i++){
        double sum=0;
        for(int j=0;j<i;j++)
            sum+=L[i][j]*Y[j];
        Y[i]=B[i]-sum;
    }
    for(int i=n-1;i>=0;i--){
        double sum=0;
        for(int j=i+1;j<n;j++)
            sum+=U[i][j]*X[j];
        if(fabs(U[i][i])<EPS){
            if(fabs(Y[i]-sum)<EPS)
                cout<<"Infinite solutions exist\n";
            else
                cout<<"No solution exists\n";
            return 0;
        }
        X[i]=(Y[i]-sum)/U[i][i];
    }
    cout<<"\nL Matrix:\n";
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            cout<<setw(10)<<L[i][j]<<" ";
        cout<<endl;
    }
    cout<<"\nU Matrix:\n";
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            cout<<setw(10)<<U[i][j]<<" ";
        cout<<endl;
    }
    cout<<"\nSolution Vector X:\n";
    for(int i=0;i<n;i++)
        cout<<"x"<<i+1<<" = "<<X[i]<<endl;
    return 0;
}

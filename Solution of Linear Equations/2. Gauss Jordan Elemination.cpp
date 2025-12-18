#include <bits/stdc++.h>
using namespace std;

void printMat(const vector<vector<double>>&mat){
    for(auto &r:mat){
        for(double v:r) cout<<v<<" ";
        cout<<"\n";
    }
    cout<<"\n";
}

int main(){
    int n;
    cout<<"Enter number of variables: ";
    cin>>n;

    vector<vector<double>>A(n,vector<double>(n+1));
    cout<<"Enter augmented matrix row wise:\n";
    for(int i=0;i<n;i++)
        for(int j=0;j<=n;j++)
            cin>>A[i][j];

    const double EPS=1e-12;

    for(int i=0;i<n;i++){
        if(fabs(A[i][i])<EPS){
            for(int r=i+1;r<n;r++){
                if(fabs(A[r][i])>EPS){
                    swap(A[i],A[r]);
                    break;
                }
            }
            if(fabs(A[i][i])<EPS){
                bool coeffZero=true;
                for(int k=0;k<n;k++)
                    if(fabs(A[i][k])>EPS) coeffZero=false;

                if(coeffZero && fabs(A[i][n])<EPS){
                    cout<<"Infinite solutions exist\n";
                    return 0;
                }
                if(coeffZero){
                    cout<<"No solution exists\n";
                    return 0;
                }
            }
        }

        double temp=A[i][i];
        for(int k=0;k<=n;k++)
            A[i][k]/=temp;

        for(int j=0;j<n;j++){
            if(j==i) continue;
            double factor=A[j][i];
            for(int k=0;k<=n;k++)
                A[j][k]-=factor*A[i][k];
        }
    }

    printMat(A);

    cout<<"The solution is:\n";
    for(int i=0;i<n;i++)
        cout<<"x"<<i+1<<" = "<<A[i][n]<<"\n";

    return 0;
}

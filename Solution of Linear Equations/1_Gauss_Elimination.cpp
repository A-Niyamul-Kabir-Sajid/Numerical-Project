#include <bits/stdc++.h>
using namespace std;

void printMat(const vector<vector<double>>&mat){
    for(auto &row:mat){
        for(double v:row) cout<<v<<" ";
        cout<<"\n";
    }
    cout<<"\n";
}
int main(){
    freopen("input/1_Gauss_Elimination_infinite_solutions.txt", "r", stdin);
    freopen("output/1_Gauss_Elimination_infinite_solutions_output.txt", "w", stdout);
    
    int n;
    //cout<<"Enter number of variables: ";
    cin>>n;

    vector<vector<double>>A(n,vector<double>(n+1));
    //cout<<"Enter augmented matrix row wise:\n";
    for(int i=0;i<n;i++)
        for(int j=0;j<=n;j++)
            cin>>A[i][j];

    const double EPS=1e-12;
    int row=0;

    for(int col=0;col<n && row<n;col++){
        int sel=row;
        for(int i=row;i<n;i++)
            if(fabs(A[i][col])>fabs(A[sel][col]))
                sel=i;

        if(fabs(A[sel][col])<EPS) continue;

        swap(A[sel],A[row]);

        for(int i=row+1;i<n;i++){
            double factor=A[i][col]/A[row][col];
            for(int j=col;j<=n;j++)
                A[i][j]-=factor*A[row][j];
        }
        row++;
    }

    printMat(A);

    int rankA=0,rankAug=0;
    for(int i=0;i<n;i++){
        bool nzA=false;
        for(int j=0;j<n;j++)
            if(fabs(A[i][j])>EPS) nzA=true;
        if(nzA) rankA++;
        if(nzA||fabs(A[i][n])>EPS) rankAug++;
    }

    if(rankA<rankAug){
        cout<<"infinite solution exists\n";
        return 0;
    }
    if(rankA<n){
        cout<<"Infinite solutions exist\n";
        return 0;
    }

    vector<double>x(n);
    for(int i=n-1;i>=0;i--){
        double sum=A[i][n];
        for(int j=i+1;j<n;j++)
            sum-=A[i][j]*x[j];
        x[i]=sum/A[i][i];
    }

    cout<<"The solution is:\n";
    for(int i=0;i<n;i++)
        cout<<"x"<<i+1<<"="<<x[i]<<"\n";
    
    return 0;
}
/*
infinite solution
input 
2
1 2 3
2 4 6  
output
2 4 6 
0 0 0 

Infinite solutions exist


*/

/*
root found
input
2
2 3 8
4 1 10
output
4 1 10 
0 2.5 3 

The solution is:
x1=2.2
x2=1.2

*/
/*
infinite solution
input
2
1 2 3
2 4 10
output
2 4 10 
0 0 -2 

infinite solution exists


*/
#include<bits/stdc++.h>
using namespace std;

const double EPS=1e-10;

int rankOfMatrix(vector<vector<double>> mat){
    int n=mat.size();
    int m=mat[0].size();
    int rank=0;
    for(int col=0,row=0;col<m&&row<n;col++){
        int sel=row;
        for(int i=row;i<n;i++)
            if(fabs(mat[i][col])>fabs(mat[sel][col]))
                sel=i;
        if(fabs(mat[sel][col])<EPS) continue;
        swap(mat[sel],mat[row]);
        for(int i=row+1;i<n;i++){
            double factor=mat[i][col]/mat[row][col];
            for(int j=col;j<m;j++)
                mat[i][j]-=factor*mat[row][j];
        }
        row++;
        rank++;
    }
    return rank;
}

double determinant(vector<vector<double>> mat){
    int n=mat.size();
    double det=1.0;
    for(int i=0;i<n;i++){
        int maxRow=i;
        for(int k=i+1;k<n;k++)
            if(fabs(mat[k][i])>fabs(mat[maxRow][i]))
                maxRow=k;
        if(fabs(mat[maxRow][i])<EPS)
            return 0.0;
        if(maxRow!=i){
            swap(mat[i],mat[maxRow]);
            det*=-1;
        }
        det*=mat[i][i];
        for(int k=i+1;k<n;k++){
            double factor=mat[k][i]/mat[i][i];
            for(int j=i;j<n;j++)
                mat[k][j]-=mat[i][j]*factor;
        }
    }
    return det;
}

vector<vector<double>> cofactor(vector<vector<double>>& mat){
    int n=mat.size();
    vector<vector<double>> cof(n,vector<double>(n));
    for(int r=0;r<n;r++){
        for(int c=0;c<n;c++){
            vector<vector<double>> sub;
            for(int i=0;i<n;i++){
                if(i==r) continue;
                vector<double> row;
                for(int j=0;j<n;j++){
                    if(j==c) continue;
                    row.push_back(mat[i][j]);
                }
                sub.push_back(row);
            }
            double minor=determinant(sub);
            cof[r][c]=((r+c)%2==0?1:-1)*minor;
        }
    }
    return cof;
}

vector<vector<double>> transpose(vector<vector<double>>& mat){
    int n=mat.size();
    vector<vector<double>> t(n,vector<double>(n));
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            t[j][i]=mat[i][j];
    return t;
}

vector<vector<double>> numeric_multiply(vector<vector<double>>& mat,double val){
    int n=mat.size();
    vector<vector<double>> res(n,vector<double>(n));
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            res[i][j]=mat[i][j]*val;
    return res;
}

void print(vector<vector<double>>& mat){
    int n=mat.size();
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            cout<<fixed<<setprecision(5)<<mat[i][j]<<" ";
        cout<<"\n";
    }
}

int main(){
    freopen("input/4_Matrix_Inversion_singular.txt", "r", stdin);
    freopen("output/4_Matrix_Inversion_singular_output.txt", "w", stdout);
    
    int n;
    cin>>n;

    vector<vector<double>> A(n,vector<double>(n));
    vector<double> B(n);

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            cin>>A[i][j];
        cin>>B[i];
    }

    double det=determinant(A);

    if(fabs(det)>EPS){
        cout<<"Single(unique) solution exists\n";
        cout<<"Matrix is singular\n";
        cout<<"Determinant:"<<fixed<<setprecision(5)<<det<<"\n";

        vector<vector<double>> cof=cofactor(A);
        vector<vector<double>> adj=transpose(cof);
        vector<vector<double>> inv=numeric_multiply(adj,1.0/det);

        cout<<"Inverse Matrix:\n";
        print(inv);

        vector<double> X(n,0.0);
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                X[i]+=inv[i][j]*B[j];

        cout<<"Solution:\n";
        for(int i=0;i<n;i++)
            cout<<"x"<<i+1<<"="<<fixed<<setprecision(5)<<X[i]<<"\n";
        return 0;
    }

    vector<vector<double>> Aug(n,vector<double>(n+1));
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            Aug[i][j]=A[i][j];
        Aug[i][n]=B[i];
    }

    int rankA=rankOfMatrix(A);
    int rankAug=rankOfMatrix(Aug);

    if(rankA!=rankAug)
        cout<<"No solution exists\n";
    else
        cout<<"Many solutions exists\n";
    
    return 0;
}

/*
infinite solution or no solution
input 
2
1 2 3
2 4 6
output
Many solutions exists

L Matrix:
 1  0 
 2  1 

U Matrix:
 1  2 
 0  0 

*/

/*
root found
input
2
2 1 8
1 3 10
output

Single(unique) solution exists
Matrix is singular
Determinant:5.00000
Inverse Matrix:
0.60000 -0.20000 
-0.20000 0.40000 
Solution:
x1=2.80000
x2=2.40000

*/

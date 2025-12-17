#include <bits/stdc++.h>
using namespace std;
void printMat(const vector<vector<double>>&mat)
{
    int n=mat.size();
    int m=mat[0].size();
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            cout<<mat[i][j]<<" ";
        }
        cout<<"\n";
    }
    cout<<endl;
}
int main()
{
    int n;
    cout<<"Enter number of variables: ";
    cin>>n;
    cout<<"Enter augmented matrix row wise:\n";
    vector<vector<double>>A(n,vector<double>(n+1));
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n+1;j++)
        {
            cin>>A[i][j];
        }
    }
    for(int i=0;i<n;i++)
    {
        if(fabs(A[i][i])<1e-12)
        {
            for(int r=i+1;r<n;r++)
            {
                if(fabs(A[r][i])>1e-12)
                {
                    swap(A[i],A[r]);
                    break;
                }
            }
            if(fabs(A[i][i])<1e-12)
            {
                bool coeffZero=true;
                for(int k=0;k<n;k++)
                {
                    if(fabs(A[i][k])>1e-12)
                    {
                        coeffZero=false;
                        break;
                    }
                }
                if(coeffZero&&fabs(A[i][n])<1e-12)
                {
                    cout << "Infinite solutions exist\n";
                }
                else if(coeffZero)
                {
                    cout<<"No solution exists\n";
                }
            }
        }
        double temp=A[i][i];
        for(int k=0;k<n+1;k++)
        {
            A[i][k]=A[i][k]/temp;
        }
        for(int j=i+1;j<n;j++)
        {
            double factor=A[j][i]/A[i][i];
            for(int k=0;k<n+1;k++)
            {
                A[j][k]=A[j][k]-factor*A[i][k];
            }
        }
        for(int j=i-1;j>=0;j--)
        {
            double factor=A[j][i]/A[i][i];
            for(int k=0;k<n+1;k++)
            {
                A[j][k]=A[j][k]-factor*A[i][k];
            }
        }
    }
    printMat(A);
    cout<<"The solution is:\n";
    for(int i=0;i<n;i++)
    {
        cout<<"x"<<i+1<<" = "<<A[i][n]<<"\n";
    }
}

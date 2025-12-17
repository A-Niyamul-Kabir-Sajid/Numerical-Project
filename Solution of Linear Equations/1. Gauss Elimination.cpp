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
        for(int j=0;j<=n;j++)
        {
            cin>>A[i][j];
        }
    }
    const double EPS=1e-12;
    for(int i=0;i<n-1;i++)
    {
        if(fabs(A[i][i])<EPS)
        {
            for(int r=i+1;r<n;r++)
            {
                if(fabs(A[r][i])>EPS)
                {
                    swap(A[i],A[r]);
                    break;
                }
            }
        }
        if(fabs(A[i][i])<EPS)
        {
            bool allZero=true;
            for(int k=0;k<n;k++)
            {
                if(fabs(A[i][k])>EPS)
                {
                    allZero=false;
                    break;
                }
            }
            if(allZero&&fabs(A[i][n])<EPS)
            {
                cout<<"Infinite solutions exist\n";
            }
            else
            {
                cout<<"No solution exists\n";
            }
            return 0;
        }
        for(int j=i+1;j<n;j++)
        {
            double factor=A[j][i]/A[i][i];
            for(int k=i;k<=n;k++)
            {
                A[j][k]-=factor*A[i][k];
            }
        }
    }
    printMat(A);
    vector<double>x(n);
    for(int i=n-1;i>=0;i--)
    {
        if(fabs(A[i][i])<EPS)
        {
            cout<<"No unique solution exists\n";
            return 0;
        }
        double sum=A[i][n];
        for(int j=i+1;j<n;j++)
        {
            sum-=A[i][j]*x[j];
        }
        x[i]=sum/A[i][i];
    }
    cout<<"The solution is:\n";
    for(int i=0;i<n;i++)
    {
        cout<<"x"<<i+1<<" = "<<x[i]<<"\n";
    }
    return 0;
}

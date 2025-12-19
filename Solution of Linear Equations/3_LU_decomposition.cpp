#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;

int matrixRank(vector<vector<double>> mat, int n, int m)
{
    int rank = 0;
    for (int row = 0, col = 0; row < n && col < m; col++)
    {
        int sel = row;
        for (int i = row; i < n; i++)
            if (fabs(mat[i][col]) > fabs(mat[sel][col]))
                sel = i;

        if (fabs(mat[sel][col]) < EPS) continue;

        swap(mat[sel], mat[row]);

        for (int i = 0; i < n; i++)
        {
            if (i != row)
            {
                double factor = mat[i][col] / mat[row][col];
                for (int j = col; j < m; j++)
                    mat[i][j] -= factor * mat[row][j];
            }
        }
        row++;
        rank++;
    }
    return rank;
}

int main()
{
    freopen("input/3_LU_decomposition_infinite_solutions.txt", "r", stdin);
    freopen("output/3_LU_decomposition_infinite_solutions_output.txt", "w", stdout);
    
    int n;
    //cout << "Enter number of variables: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n), Y(n), X(n, 0);
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    //cout << "Enter augmented matrix row wise:\n";
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cin >> A[i][j];
        cin >> B[i];
    }

    for(int i = 0; i < n; i++)
    {
        for(int k = i; k < n; k++)
        {
            double sum = 0;
            for(int j = 0; j < i; j++)
                sum += L[i][j] * U[j][k];
            U[i][k] = A[i][k] - sum;
        }

        for(int k = i; k < n; k++)
        {
            if(i == k)
                L[i][i] = 1;
            else
            {
                double sum = 0;
                for(int j = 0; j < i; j++)
                    sum += L[k][j] * U[j][i];
                if(fabs(U[i][i]) < EPS)
                    L[k][i] = 0;
                else
                    L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }

    double det = 1;
    for(int i = 0; i < n; i++)
        det *= U[i][i];

    if(fabs(det) > EPS)
    {
        for(int i = 0; i < n; i++)
        {
            double sum = 0;
            for(int j = 0; j < i; j++)
                sum += L[i][j] * Y[j];
            Y[i] = B[i] - sum;
        }

        for(int i = n - 1; i >= 0; i--)
        {
            double sum = 0;
            for(int j = i + 1; j < n; j++)
                sum += U[i][j] * X[j];
            X[i] = (Y[i] - sum) / U[i][i];
        }

        cout << "\ninfinite solution exists.\n";
        cout << "Solution Vector X:\n";
        for(int i = 0; i < n; i++)
            cout << "x" << i+1 << " = " << fixed << setprecision(6) << X[i] << endl;
    }
    else
    {
        vector<vector<double>> aug(n, vector<double>(n+1));
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
                aug[i][j] = A[i][j];
            aug[i][n] = B[i];
        }

        int rankA = matrixRank(A, n, n);
        int rankAug = matrixRank(aug, n, n+1);

        cout << "\nRank(A) = " << rankA << ", Rank(A|B) = " << rankAug << endl;

        if(rankA == rankAug && rankA < n)
            cout << "Infinite solutions exist.\n";
        else
            cout << "infinite solution exists.\n";
    }

    cout << "\nL Matrix:\n";
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout << setw(2) << L[i][j] << " ";
        cout << endl;
    }

    cout << "\nU Matrix:\n";
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            cout << setw(2) << U[i][j] << " ";
        cout << endl;
    }
    
    return 0;
}
/*
infinite solution
input 
2
1 2 3
2 4 6
output

Rank(A) = 1, Rank(A|B) = 1
Infinite solutions exist.

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
2 3 8
4 1 10
output


infinite solution exists.
Solution Vector X:
x1 = 2.200000
x2 = 1.200000

L Matrix:
1.000000 0.000000 
2.000000 1.000000 

U Matrix:
2.000000 3.000000 
0.000000 -5.000000 



*/
/*
infinite solution
input
2
1 2 3
2 4 10
output

Rank(A) = 1, Rank(A|B) = 2
infinite solution exists.

L Matrix:
 1  0 
 2  1 

U Matrix:
 1  2 
 0  0 


*/
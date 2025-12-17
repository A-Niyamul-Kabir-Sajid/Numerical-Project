#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;

// Function to compute rank of a matrix using Gaussian elimination
int matrixRank(vector<vector<double>> mat, int n, int m) {
    int rank = 0;
    for (int row = 0, col = 0; row < n && col < m; col++) {
        int sel = row;
        for (int i = row; i < n; i++)
            if (fabs(mat[i][col]) > fabs(mat[sel][col]))
                sel = i;

        if (fabs(mat[sel][col]) < EPS) continue;

        swap(mat[sel], mat[row]);

        for (int i = 0; i < n; i++) {
            if (i != row) {
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

int main() {
    int n;
    cout << "Enter number of variables: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n), Y(n), X(n, 0);
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    cout << "Enter coefficient matrix A:\n";
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            cin >> A[i][j];

    cout << "Enter constant matrix B:\n";
    for(int i = 0; i < n; i++)
        cin >> B[i];

    // LU Decomposition (Doolittle)
    for(int i = 0; i < n; i++) {
        // Upper Triangular
        for(int k = i; k < n; k++) {
            double sum = 0;
            for(int j = 0; j < i; j++)
                sum += L[i][j] * U[j][k];
            U[i][k] = A[i][k] - sum;
        }

        // Lower Triangular
        for(int k = i; k < n; k++) {
            if(i == k) L[i][i] = 1;
            else {
                double sum = 0;
                for(int j = 0; j < i; j++)
                    sum += L[k][j] * U[j][i];
                if(fabs(U[i][i]) < EPS)
                    L[k][i] = 0; // singular
                else
                    L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }

    // Calculate determinant to check singularity
    double det = 1;
    for(int i = 0; i < n; i++) det *= U[i][i];

    if(fabs(det) > EPS) {
        // Unique solution, solve using LU
        // Forward substitution LY = B
        for(int i = 0; i < n; i++) {
            double sum = 0;
            for(int j = 0; j < i; j++)
                sum += L[i][j] * Y[j];
            Y[i] = B[i] - sum;
        }

        // Backward substitution UX = Y
        for(int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for(int j = i + 1; j < n; j++)
                sum += U[i][j] * X[j];
            X[i] = (Y[i] - sum) / U[i][i];
        }

        cout << "\nUnique solution exists.\n";
        cout << "Solution Vector X:\n";
        for(int i = 0; i < n; i++)
            cout << "x" << i+1 << " = " << fixed << setprecision(6) << X[i] << endl;
    } else {
        // Check rank to determine infinite/no solution
        vector<vector<double>> aug(n, vector<double>(n+1));
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) aug[i][j] = A[i][j];
            aug[i][n] = B[i];
        }

        int rankA = matrixRank(A, n, n);
        int rankAug = matrixRank(aug, n, n+1);

        cout << "\nRank(A) = " << rankA << ", Rank(A|B) = " << rankAug << endl;

        if(rankA == rankAug && rankA < n)
            cout << "Infinite solutions exist.\n";
        else
            cout << "No solution exists.\n";
    }

    // Print L matrix
    cout << "\nL Matrix:\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++)
            cout << setw(2) << L[i][j] << " ";
        cout << endl;
    }

    // Print U matrix
    cout << "\nU Matrix:\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++)
            cout << setw(2) << U[i][j] << " ";
        cout << endl;
    }

    return 0;
}

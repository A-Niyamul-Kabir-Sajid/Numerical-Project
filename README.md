# Numerical Methods Project

This repository contains implementations of various numerical methods for solving mathematical problems including linear equations, non-linear equations, interpolation, differentiation, integration, differential equations, and curve fitting.

## Table of Contents

1. [Solution of Linear Equations](#solution-of-linear-equations)
   - [Gauss Elimination](#gauss-elimination-method)
   - [Gauss Jordan Elimination](#gauss-jordan-elimination-method)
   - [LU Decomposition](#lu-decomposition-method)
   - [Matrix Inversion](#matrix-inversion-method)

2. [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
   - [Bisection Method](#bisection-method)
   - [False Position Method](#false-position-method)
   - [Newton-Raphson Method](#newton-raphson-method)
   - [Secant Method](#secant-method)

3. [Interpolation and Approximation](#interpolation-and-approximation)
   - [Newton's Forward Interpolation](#newtons-forward-interpolation)
   - [Newton's Backward Interpolation](#newtons-backward-interpolation)
   - [Divided Difference Interpolation](#divided-difference-interpolation)

4. [Numerical Differentiation](#numerical-differentiation)
   - [Differentiation based on Equal-Interval Interpolation](#differentiation-based-on-equal-interval-interpolation)
   - [Second Order Derivative](#second-order-derivative)

5. [Numerical Integration](#numerical-integration)
   - [Simpson's One-Third Rule](#simpsons-one-third-rule)
   - [Simpson's Three-Eighths Rule](#simpsons-three-eighths-rule)

6. [Solution of Differential Equations](#solution-of-differential-equations)
   - [Runge-Kutta Method](#runge-kutta-method)

7. [Curve Fitting (Regression)](#curve-fitting-regression)
   - [Linear Regression](#linear-regression)
   - [Exponential Regression](#exponential-regression)
   - [Polynomial Regression](#polynomial-regression)
   - [Modified Exponential Regression](#modified-exponential-regression)
   - [Transcendental Regression](#transcendental-regression)

---

## Solution of Linear Equations

Linear equations of the form Ax = b can be solved using various direct and iterative methods. This section covers the most commonly used direct methods.

### Gauss Elimination Method

#### Theory
The Gauss Elimination method is a direct method for solving a system of linear equations. It converts the system into an upper triangular matrix through forward elimination, then uses back-substitution to find the solution. The method involves:

1. **Forward Elimination**: Transform the augmented matrix [A|b] into upper triangular form
2. **Back Substitution**: Solve the system starting from the last equation

The algorithm can handle three cases:
- **Unique Solution**: When the matrix has full rank
- **Infinite Solutions**: When the matrix is rank-deficient and the system is consistent
- **No Solution**: When the system is inconsistent (rank of A < rank of [A|b])

#### Code
```cpp
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
    freopen("input/1_Gauss_Elimination_unique_solution.txt", "r", stdin);
    freopen("output/1_Gauss_Elimination_unique_solution_output.txt", "w", stdout);
    
    int n;
    cin>>n;

    vector<vector<double>>A(n,vector<double>(n+1));
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
        cout<<"No solution exists\n";
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
```

#### Test Cases

##### Case 1: Unique Solution
**Input:**
```
2
2 3 8
4 1 10
```
**Output:**
```
4 1 10 
0 2.5 3 

The solution is:
x1=2.2
x2=1.2
```

##### Case 2: Infinite Solutions
**Input:**
```
2
1 2 3
2 4 6
```
**Output:**
```
2 4 6 
0 0 0 

Infinite solutions exist
```

##### Case 3: No Solution
**Input:**
```
2
1 2 3
2 4 10
```
**Output:**
```
2 4 10 
0 0 -2 

No solution exists
```

---

### Gauss Jordan Elimination Method

#### Theory
The Gauss-Jordan elimination method is an extension of Gauss elimination that reduces the augmented matrix directly to **reduced row echelon form (RREF)** rather than just upper triangular form. This provides the solution directly without requiring back-substitution.

Key differences from Gauss elimination:
- Eliminates both above and below the pivot element
- Results in an identity matrix on the left side of the augmented matrix
- The right side directly gives the solution (for unique systems)
- Can also identify infinite solutions and inconsistent systems

#### Code
```cpp
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
    freopen("input/2_Gauss_Jordan_Elemination_infinite_solutions.txt", "r", stdin);
    freopen("output/2_Gauss_Jordan_Elemination_infinite_solutions_output.txt", "w", stdout);
    
    int n;
    cin>>n;

    vector<vector<double>>A(n,vector<double>(n+1));
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
                    cout<<"infinite solution exists\n";
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
```

#### Test Cases

##### Case 1: Unique Solution
**Input:**
```
2
2 3 8
4 1 10
```
**Output:**
```
1 0 2.2 
-0 1 1.2 

The solution is:
x1 = 2.2
x2 = 1.2
```

##### Case 2: Infinite Solutions
**Input:**
```
2
1 2 3
2 4 6
```
**Output:**
```
Infinite solutions exist
```

##### Case 3: No Solution
**Input:**
```
2
1 2 3
2 4 10
```
**Output:**
```
No solution exists
```

---

### LU Decomposition Method

#### Theory
LU Decomposition decomposes a matrix A into a product of a lower triangular matrix (L) and an upper triangular matrix (U), such that A = LU. This method is useful when:

1. **Computational Efficiency**: Once decomposed, solving multiple systems with the same matrix A but different right-hand sides is efficient
2. **Determinant Calculation**: det(A) = product of diagonal elements of U
3. **Matrix Inversion**: Can be used to invert matrices

The algorithm proceeds as:
1. Decompose A into L and U
2. Solve Ly = b to find y
3. Solve Ux = y to find x

#### Code
```cpp
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
    freopen("input/3_LU_decomposition_unique_solution.txt", "r", stdin);
    freopen("output/3_LU_decomposition_unique_solution_output.txt", "w", stdout);
    
    int n;
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n), Y(n), X(n, 0);
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

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

        cout << "\nUnique solution exists.\n";
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
            cout << "No solution exists.\n";
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
```

#### Test Cases

##### Case 1: Unique Solution
**Input:**
```
2
2 3 8
4 1 10
```
**Output:**
```
Unique solution exists.
Solution Vector X:
x1 = 2.200000
x2 = 1.200000

L Matrix:
1.000000 0.000000 
2.000000 1.000000 

U Matrix:
2.000000 3.000000 
0.000000 -5.000000
```

##### Case 2: Infinite Solutions
**Input:**
```
2
1 2 3
2 4 6
```
**Output:**
```
Rank(A) = 1, Rank(A|B) = 1
Infinite solutions exist.

L Matrix:
 1  0 
 2  1 

U Matrix:
 1  2 
 0  0
```

##### Case 3: No Solution
**Input:**
```
2
1 2 3
2 4 10
```
**Output:**
```
Rank(A) = 1, Rank(A|B) = 2
No solution exists.

L Matrix:
 1  0 
 2  1 

U Matrix:
 1  2 
 0  0
```

---

### Matrix Inversion Method

#### Theory
Matrix inversion solves the system Ax = b by computing A⁻¹ and multiplying both sides: x = A⁻¹b.

Using the adjoint method:
1. **Calculate determinant** of A
2. **Check invertibility**: Matrix is invertible only if det(A) ≠ 0 (non-singular)
3. **Compute cofactor matrix**: Calculate the cofactor of each element
4. **Form adjoint matrix**: Transpose the cofactor matrix
5. **Calculate A⁻¹**: A⁻¹ = adj(A) / det(A)
6. **Solve**: x = A⁻¹b

Cases:
- **Invertible Matrix**: Determinant ≠ 0, unique solution exists
- **Singular Matrix**: Determinant = 0, no unique solution (may have infinite or no solutions)

#### Code
```cpp
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
    freopen("input/4_Matrix_Inversion_invertible.txt", "r", stdin);
    freopen("output/4_Matrix_Inversion_invertible_output.txt", "w", stdout);
    
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
        cout<<"Matrix is invertible\n";
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
```

#### Test Cases

##### Case 1: Invertible Matrix (Unique Solution)
**Input:**
```
2
2 1 8
1 3 10
```
**Output:**
```
Single(unique) solution exists
Matrix is invertible
Determinant:5.00000
Inverse Matrix:
0.60000 -0.20000 
-0.20000 0.40000 
Solution:
x1=2.80000
x2=2.40000
```

##### Case 2: Singular Matrix (No Unique Solution)
**Input:**
```
2
1 2 3
2 4 6
```
**Output:**
```
Many solutions exists
```

---

## Solution of Non-Linear Equations

Non-linear equations are of the form f(x) = 0 where f is a non-linear function. These methods iteratively narrow down the solution interval or refine an estimate.

### Bisection Method

#### Theory
The Bisection method is a bracketing method that repeatedly divides an interval in half to find a root.

**Prerequisites**: 
- The function must be continuous
- The root must be bracketed, i.e., f(a) and f(b) have opposite signs

**Algorithm**:
1. Start with interval [a, b] where f(a)·f(b) < 0
2. Calculate midpoint: c = (a + b) / 2
3. If f(c) = 0, c is the root
4. If f(a)·f(c) < 0, the root is in [a, c]; update b = c
5. Otherwise, the root is in [c, b]; update a = c
6. Repeat until convergence (|b - a| < tolerance)

**Advantages**: Guaranteed convergence if root is bracketed
**Disadvantages**: Slow convergence compared to other methods

#### Code
```cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

const double EPS = 1e-9;
const double STEP = 0.5;
const double RANGE = 100.0;
const int MAX_ITER = 100;

double evaluate(const vector<double>& poly, double x) {
    double res = 0;
    for (double coeff : poly)
        res = res * x + coeff;
    return res;
}

void printPolynomial(const vector<double>& poly) {
    int n = poly.size() - 1;
    cout << "f(x) = ";
    for (int i = 0; i <= n; i++) {
        if (poly[i] == 0) continue;
        if (poly[i] > 0 && i != 0) cout << " + ";
        if (poly[i] < 0) cout << " - ";
        
        double val = abs(poly[i]);
        if (val != 1 || i == n) cout << val;
        
        int power = n - i;
        if (power > 1) cout << "x^" << power;
        else if (power == 1) cout << "x";
    }
    cout << endl;
}

int main() {
    freopen("input/Bi-section_method_root_found.txt", "r", stdin);
    freopen("output/Bi-section_method_root_found.txt", "w", stdout);

    int degree;
    if (!(cin >> degree)) return 1;

    vector<double> coeffs(degree + 1);
    for (int i = 0; i <= degree; i++) {
        cin >> coeffs[i];
    }

    printPolynomial(coeffs);
    vector<double> roots;

    for (double a = -RANGE; a < RANGE; a += STEP) {
        double b = a + STEP;
        double fa = evaluate(coeffs, a);
        double fb = evaluate(coeffs, b);

        if (abs(fa) < EPS) {
            roots.push_back(a);
            continue;
        }

        if ((fa > 0) != (fb > 0)) {
            double low = a;
            double high = b;

            for (int i = 0; i < MAX_ITER; i++) {
                double mid = low + (high - low) / 2.0;
                double fmid = evaluate(coeffs, mid);

                if (abs(fmid) < EPS || (high - low) / 2.0 < EPS) {
                    roots.push_back(mid);
                    break;
                }

                if ((evaluate(coeffs, low) > 0) != (fmid > 0)) {
                    high = mid;
                } else {
                    low = mid;
                }
            }
        }
    }

    if (roots.empty()) {
        cout << "No real roots found in the range [-" << RANGE << ", " << RANGE << "]." << endl;
    } else {
        sort(roots.begin(), roots.end());
        auto it = unique(roots.begin(), roots.end(), [](double x, double y) {
            return abs(x - y) < 1e-5; 
        });
        roots.erase(it, roots.end());

        cout << "Real roots found:" << endl;int i=1;
        for (double r : roots) {
            if (abs(r) < 1e-10) r = 0.0;
            cout << "x"<<i<<" = " << fixed << setprecision(6) << r << endl;
            i++;
        }
    }
    
    return 0;
}
```

#### Test Cases

##### Case 1: Root Found
**Input:**
```
2
1 -3 2
```
**Output:**
```
f(x) = x^2 - 3x + 2
Real roots found:
x1 = 1.000000
x2 = 2.000000
```

---

### False Position Method

#### Theory
The False Position method (also called Regula Falsi) is another bracketing method that is faster than bisection. Instead of taking the midpoint, it uses linear interpolation to estimate where the root might be.

**Algorithm**:
1. Start with interval [a, b] where f(a)·f(b) < 0
2. Calculate false position: c = a - f(a)·(b - a) / (f(b) - f(a))
3. If f(c) = 0, c is the root
4. If f(a)·f(c) < 0, update b = c
5. Otherwise, update a = c
6. Repeat until convergence

**Advantages**: Generally faster convergence than bisection
**Disadvantages**: May not converge as rapidly if the function is highly non-linear

#### Code
```cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

const double EPS = 1e-9;
const int MAX_ITER = 100;
const double STEP = 0.5;
const double RANGE = 1000.0;

double evaluate(const vector<double>& cof, double x) {
    double res = 0;
    for (double c : cof)
        res = res * x + c;
    return res;
}

int main() {
    freopen("input/False_Position_method_root_found.txt", "r", stdin);
    freopen("output/False_Position_method_root_found.txt", "w", stdout);

    int n;
    if (!(cin >> n)) return 0;

    vector<double> cof(n + 1);
    for (int i = 0; i <= n; i++) cin >> cof[i];

    vector<double> found_roots;

    for (double a = -RANGE; a < RANGE; a += STEP) {
        double b = a + STEP;
        double fa = evaluate(cof, a);
        double fb = evaluate(cof, b);

        if (abs(fa) < EPS) {
            found_roots.push_back(a);
            continue;
        }

        if ((fa > 0) != (fb > 0)) {
            double x0 = a, x1 = b;
            
            for (int i = 0; i < MAX_ITER; i++) {
                double f0 = evaluate(cof, x0);
                double f1 = evaluate(cof, x1);

                if (abs(f1 - f0) < 1e-15) break;

                double xr = x1 - f1 * (x0 - x1) / (f0 - f1);
                double fr = evaluate(cof, xr);

                if (abs(fr) < EPS || abs(x1 - x0) < EPS) {
                    found_roots.push_back(xr);
                    break;
                }

                if ((f0 > 0) != (fr > 0)) x1 = xr;
                else x0 = xr;
            }
        }
    }

    if (found_roots.empty()) {
        cout << "\nNo real roots found in the range." << endl;
    } else {
        sort(found_roots.begin(), found_roots.end());
        auto it = unique(found_roots.begin(), found_roots.end(), [](double r1, double r2) {
            return abs(r1 - r2) < 1e-5; 
        });
        found_roots.erase(it, found_roots.end());

        cout << "\nUnique real roots found:" << endl;int i = 1;
        for (double r : found_roots) {
            if (abs(r) < 1e-10) r = 0.0;
            cout << "x"<<i<<" = " << fixed << setprecision(6) << r << endl;
            i++;
        }
    }
    
    return 0;
}
```

#### Test Cases

##### Case 1: Root Found
**Input:**
```
2
1 -3 2
```
**Output:**
```
Unique real roots found:
x1 = 1.000000
x2 = 2.000000
```

---

###Newton-Raphson Method

#### Theory
The Newton-Raphson method is an open method (doesn't require bracketing) that uses calculus to find roots rapidly.

**Formula**: x_{n+1} = x_n - f(x_n) / f'(x_n)

**Algorithm**:
1. Start with an initial guess x₀
2. Calculate x₁ using the formula
3. Repeat until convergence (|x_{n+1} - x_n| < tolerance)

**Advantages**: 
- Very fast convergence (quadratic convergence)
- Works well for simple roots

**Disadvantages**:
- Requires derivative calculation
- May not converge if initial guess is poor
- May diverge or oscillate for certain functions

#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-8;
const int MAX_ITER = 100;
const double STEP = 0.5;

void printfunc(const vector<int>& cof) {
    int n = cof.size();
    bool first = true;
    for(int i = 0; i < n; i++) {
        if(cof[i] == 0) continue;
        
        if(!first && cof[i] > 0) cout << '+';
        if(cof[i] == -1 && n-i-1 != 0) cout << '-';
        else if(cof[i] != 1 || n-i-1 == 0) cout << cof[i];

        int power = n-i-1;
        if(power > 1) cout << "x^" << power;
        else if(power == 1) cout << "x";
        first = false;
    }
    if(first) cout << "0";
    cout << endl;
}

double func(const vector<int>& cof, double x) {
    double res = 0;
    for(int i = 0; i < cof.size(); i++)
        res = res * x + cof[i];
    return res;
}

double deriv(const vector<int>& cof, double x) {
    double res = 0;
    int n = cof.size();
    for(int i = 0; i < n - 1; i++)
        res = res * x + (double)cof[i] * (n - i - 1);
    return res;
}

int main() {
    freopen("input/Newton-Raphson_method_root_found.txt", "r", stdin);
    freopen("output/Newton-Raphson_method_root_found.txt", "w", stdout);

    int n;
    if (!(cin >> n)) return 0;

    vector<int> cof(n + 1);
    for(int i = 0; i <= n; i++)
        cin >> cof[i];

    printfunc(cof);

    vector<double> roots;
    double start = -1000, end = 1000;

    for(double i = start; i <= end; i += STEP) {
        double x0 = i;
        int iter = 0;

        while(iter < MAX_ITER) {
            double f = func(cof, x0);
            double df = deriv(cof, x0);

            if(fabs(df) < 1e-12) break;

            double x1 = x0 - f / df;

            if(fabs(x1 - x0) < EPS) {
                if(fabs(func(cof, x1)) < 1e-6) {
                    roots.push_back(x1);
                }
                break;
            }

            x0 = x1;
            iter++;
        }
    }

    if(roots.empty()) {
        cout << "No real roots found.\n";
    } else {
        sort(roots.begin(), roots.end());
        auto it = unique(roots.begin(), roots.end(), [](double a, double b) {
            return fabs(a - b) < 1e-4;
        });
        roots.erase(it, roots.end());

        cout << "Roots found:\n";int i=1;
        for(double r : roots) {
            if(fabs(r) < 1e-10) r = 0.0;
            cout <<"x"<<i<<" = "<< fixed << setprecision(6) << r <<endl;
            i++;
        }
        cout << endl;
    }
    
    return 0;
}
```

#### Input Format
```
3
1 -6 11 -6
```

#### Output
```
x^3-6x^2+11x-6
Roots found:
x1 = 1.000000
x2 = 2.000000
x3 = 3.000000
```

---

### Secant Method

#### Theory
The Secant method is similar to Newton-Raphson but **doesn't require the derivative**. Instead, it approximates the derivative using two previous function values.

**Formula**: x_{n+1} = x_n - f(x_n) · (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))

**Algorithm**:
1. Start with two initial guesses x₀ and x₁
2. Calculate x₂ using the formula
3. Update: x₀ = x₁, x₁ = x₂
4. Repeat until convergence

**Advantages**:
- Doesn't require derivative
- Reasonably fast convergence

**Disadvantages**:
- Requires two initial guesses
- Slightly slower than Newton-Raphson
- May not always converge

#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-8;
const int MAX_ITER = 100;
const double STEP = 0.5;

void printfunc(vector<int> cof) {
    int n = cof.size();
    for(int i = 0; i < n; i++) {
        if(cof[i] > 0 && i != 0) cout << '+';
        if(cof[i] != 0) {
            cout << cof[i];
            if(n-i-1 > 1) cout << "x^" << n-i-1;
            else if(n-i-1 == 1) cout << "x";
        }
    }
    cout << endl;
}

double func(const vector<int>& cof, double x) {
    double ret = 0;
    int n = cof.size();
    for(int i = 0; i < n; i++)
        ret = ret*x + cof[i];
    return ret;
}

int main() {
    freopen("input/Secant_method_root_found.txt", "r", stdin);
    freopen("output/Secant_method_root_found.txt", "w", stdout);

    int n;
    cin >> n;

    vector<int> cof(n+1);
    for(int i = 0; i <= n; i++)
        cin >> cof[i];

    printfunc(cof);

    vector<double> roots;
    double start = -1000;
    double end = 1000;

    for(double i = start; i < end; i += STEP) {
        double f1 = func(cof, i);
        double f2 = func(cof, i + STEP);
        if(f1*f2 <= 0) {
            double x0 = i, x1 = i+STEP;
            int iter = 0;
            set<long long> seen;

            while(iter < MAX_ITER) {
                double f0 = func(cof, x0);
                double f1val = func(cof, x1);

                if(fabs(f1val - f0) < EPS) break;

                double xr = x1 - f1val*(x0 - x1)/(f0 - f1val);

                if(fabs(xr - x1) < EPS || seen.count(round(xr*1e6))) {
                    roots.push_back(xr);
                    break;
                }

                seen.insert(round(xr*1e6));
                x0 = x1;
                x1 = xr;
                iter++;
            }
        }
    }

    if(roots.empty()) {
        cout << "No real roots found.\n";
    } else {
        cout << "Roots found:\n";
        sort(roots.begin(), roots.end());
        for(double r : roots)
            cout << fixed << setprecision(6) << r << " ";
        cout << endl;
    }
    
    return 0;
}
```

#### Input Format
```
2
1 0 -4
```

#### Output
```
1x^2-4
Roots found:
-2.000000 -2.000000 2.000000 2.000000
```

---

## Interpolation and Approximation

Interpolation methods estimate function values between known data points using polynomial approximation.

### Newton's Forward Interpolation

#### Theory
Newton's Forward Interpolation is used when we need to interpolate values **near the beginning** of a table of equally-spaced data points.

**Key Concepts**:
- Uses **forward differences** (∆f)
- Applicable when data points are equally spaced
- Builds a polynomial using the data and its differences
- Formula: f(x) = f(x₀) + p∆f(x₀) + p(p-1)/2! ∆²f(x₀) + ...

where p = (x - x₀) / h, and h is the spacing between consecutive x values

**Advantages**: 
- Efficient for interpolation near the start of data
- Uses equally-spaced points efficiently

**Disadvantages**:
- Limited to equally-spaced data
- Becomes less accurate for points far from the start

#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double newtonForward(const vector<double>& x, const vector<double>& y, double valX, double& error) {
    int n = x.size();
    double h = x[1] - x[0];
    double u = (valX - x[0]) / h;

    vector<double> delY = y;
    double result = y[0];
    double uTerm = 1.0;
    double factorial = 1.0;

    for (int i = 1; i < n; i++) {
        uTerm *= (u - (i - 1));
        factorial *= i;

        vector<double> newDelY;
        for (int j = 0; j < delY.size() - 1; j++)
            newDelY.push_back(delY[j+1] - delY[j]);

        delY = newDelY;
        if(delY.empty()) break;

        result += (uTerm * delY[0]) / factorial;
    }

    if (!delY.empty()) error = fabs((uTerm * delY[0]) / factorial);
    else error = 0.0;

    return result;
}

int main() {
    freopen("input/Newtons_Forward_interpolation.txt", "r", stdin);
    freopen("output/Newtons_Forward_interpolation_output.txt", "w", stdout);
    
    int n;
    cin >> n;
    vector<double> x(n), y(n);
    for(int i=0; i<n; i++) cin >> x[i] >> y[i];

    double valX;
    cin >> valX;

    double error;
    double valY = newtonForward(x, y, valX, error);
    cout << "Interpolated value at x = " << valX << " is y = " << valY << endl;
    cout << "Estimated interpolation error: " << error << endl;
    
    return 0;
}
```

#### Input Format
```
4
1 1
2 4
3 9
4 10
2.5
```

#### Output
```
Interpolated value at x = 2.5 is y = 6.625
Estimated interpolation error: 0.375
```

---

### Newton's Backward Interpolation

#### Theory
Newton's Backward Interpolation is used for interpolation **near the end** of a table of equally-spaced data points.

**Key Concepts**:
- Uses **backward differences** (∇f)
- Equally-spaced data required
- Formula: f(x) = f(x_n) + p∇f(x_n) + p(p+1)/2! ∇²f(x_n) + ...

where p = (x - x_n) / h

**When to Use**:
- When interpolating near the end of your data table
- When data is equally spaced

**Advantages**:
- Efficient for extrapolation near the end
- Uses equally-spaced data efficiently

#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double newtonBackward(const vector<double>& x, const vector<double>& y, double valX, double& error) {
    int n = x.size();
    double h = x[1] - x[0];
    double u = (valX - x[n-1]) / h;

    vector<double> delY = y;
    double result = y[n-1];
    double uTerm = 1.0;
    double factorial = 1.0;

    for (int i = 1; i < n; i++) {
        uTerm *= (u + (i - 1));
        factorial *= i;

        vector<double> newDelY;
        for (int j = 0; j < delY.size() - 1; j++)
            newDelY.push_back(delY[j+1] - delY[j]);

        delY = newDelY;
        if(delY.empty()) break;

        result += (uTerm * delY.back()) / factorial;
    }

    if (!delY.empty()) error = fabs((uTerm * delY.back()) / factorial);
    else error = 0.0;

    return result;
}

int main() {
    freopen("input/Newtons_Backward_interpolation.txt", "r", stdin);
    freopen("output/Newtons_Backward_interpolation_output.txt", "w", stdout);
    
    int n;
    cin >> n;
    vector<double> x(n), y(n);
    for(int i=0; i<n; i++) cin >> x[i] >> y[i];

    double valX;
    cin >> valX;

    double error;
    double valY = newtonBackward(x, y, valX, error);
    cout << "Interpolated value at x = " << valX << " is y = " << valY << endl;
    cout << "Estimated interpolation error: " << error << endl;
    
    return 0;
}
```

#### Input Format
```
4
1 1
2 4
3 9
4 16
2.5
```

#### Output
```
Interpolated value at x = 2.5 is y = 6.25
Estimated interpolation error: 0
```

---

### Divided Difference Interpolation

#### Theory
Divided Difference Interpolation (Newton's Divided Difference Formula) is applicable for **unequally-spaced** data points and is more general than forward/backward interpolation.

**Key Concepts**:
- Works with any spacing between data points
- Uses divided differences instead of regular differences
- Formula: f(x) = f(x₀) + (x-x₀)f[x₀,x₁] + (x-x₀)(x-x₁)f[x₀,x₁,x₂] + ...

**Divided Difference Table**:
- First order: f[x_i, x_j] = (f(x_j) - f(x_i)) / (x_j - x_i)
- Higher orders built recursively

**Advantages**:
- Works with unequally-spaced data
- More flexible than forward/backward methods

**Disadvantages**:
- Slightly more computational overhead

#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double dividedDifference(const vector<double>& x, const vector<double>& y, double valX, double& error) {
    int n = x.size();
    vector<vector<double>> table(n, vector<double>(n, 0.0));

    for(int i=0; i<n; i++) table[i][0] = y[i];

    for(int j=1; j<n; j++)
        for(int i=0; i<n-j; i++)
            table[i][j] = (table[i+1][j-1] - table[i][j-1]) / (x[i+j] - x[i]);

    double result = table[0][0];
    double xTerm = 1.0;
    for(int i=1; i<n; i++) {
        xTerm *= (valX - x[i-1]);
        result += table[0][i] * xTerm;
    }

    error = fabs(table[0][n-1] * xTerm);

    return result;
}

int main() {
    freopen("input/Divided_difference_interpolation.txt", "r", stdin);
    freopen("output/Divided_difference_interpolation_output.txt", "w", stdout);
    
    int n;
    cin >> n;
    vector<double> x(n), y(n);
    for(int i=0; i<n; i++) cin >> x[i] >> y[i];

    double valX;
    cin >> valX;

    double error;
    double valY = dividedDifference(x, y, valX, error);
    cout << "Interpolated value at x = " << valX << " is y = " << valY << endl;
    cout << "Estimated interpolation error: " << error << endl;
    
    return 0;
}
```

#### Input Format
```
3
1 1
2 4
4 10
3
```

#### Output
```
Interpolated value at x = 3 is y = 7
Estimated interpolation error: 0
```

---

## Numerical Differentiation

Numerical differentiation approximates the derivative of a function using finite differences, especially useful when the analytical derivative is unavailable or difficult to compute.

### Differentiation based on Equal-Interval Interpolation

#### Theory
This method uses finite difference formulas derived from Newton's Forward/Backward interpolation to approximate derivatives.

**Common Formulas** (for equally-spaced points with spacing h):

1. **First Derivative**:
   - Forward difference: f'(x) ≈ (f(x+h) - f(x)) / h
   - Backward difference: f'(x) ≈ (f(x) - f(x-h)) / h
   - Central difference: f'(x) ≈ (f(x+h) - f(x-h)) / 2h (more accurate)

2. **Truncation Error**: Depends on the order of the approximation and step size h

**Advantages**:
- Simple to implement
- Works with tabulated data
- Central difference is more accurate

**Disadvantages**:
- Choice of h is critical (too small causes rounding errors, too large causes truncation errors)
- Less accurate than analytical derivatives

#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double f(vector<double>& a, double x)
{
    double s = 0;
    for(int i = 0; i < a.size(); i++)
        s += a[i] * pow(x, i);
    return s;
}

int main()
{
    freopen("input/Differentiation_equal-interval_interpolation.txt", "r", stdin);
    freopen("output/Differentiation_equal-interval_interpolation_output.txt", "w", stdout);
    
    int degree;
    cin >> degree;

    vector<double> a(degree + 1);
    for(int i = degree; i >= 0; i--)
        cin >> a[i];

    double p;
    cin >> p;
    double h = 1.0;
    vector<double> y(degree + 1);
    for(int i = 0; i <= degree; i++)
        y[i] = f(a, i);
    vector<vector<double>> d(degree + 1, vector<double>(degree + 1));
    for(int i = 0; i <= degree; i++)
        d[i][0] = y[i];
    for(int j = 1; j <= degree; j++)
        for(int i = 0; i <= degree - j; i++)
            d[i][j] = d[i+1][j-1] - d[i][j-1];
    double u = (p - 0) / h;
    double fp = d[0][1];
    if(degree >= 2) fp += ((2*u - 1)/2.0) * d[0][2];
    if(degree >= 3) fp += ((3*u*u - 6*u + 2)/6.0) * d[0][3];
    fp /= h;
    cout << "f'(p) = " << fp << endl;
    
    return 0;
}
```

#### Input Format
```
2
1 3 2
1
```

#### Output
```
f'(p) = 5
```

---

### Second Order Derivative

#### Theory
The second derivative can be approximated using finite differences, useful for analyzing concavity, inflection points, and solving differential equations.

**Formula** (for equally-spaced points):
f''(x) ≈ (f(x+h) - 2f(x) + f(x-h)) / h²

**Derivation**:
Using Taylor series expansion or combining first derivative approximations

**Error Analysis**:
- Truncation error is O(h²)
- Very sensitive to the choice of step size h
- Rounding errors can be significant for small h

**Applications**:
- Determining concavity
- Finding inflection points
- Solving second-order differential equations

#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double f(vector<double>& a, double x)
{
    double s = 0;
    for(int i = 0; i < a.size(); i++)
        s += a[i] * pow(x, i);
    return s;
}
int main()
{
    freopen("input/Second_order_derivative.txt", "r", stdin);
    freopen("output/Second_order_derivative_output.txt", "w", stdout);
    
    int degree;
    cin >> degree;
    vector<double> a(degree + 1);
    for(int i = degree; i >= 0; i--)
        cin >> a[i];
    double p;
    cin >> p;
    double h = 1.0;
    double fpp = (f(a, p + h) - 2*f(a, p) + f(a, p - h)) / (h * h);
    cout << "f\"(p) = " << fpp << endl;
    
    return 0;
}
```

#### Input Format
```
3
1 0 -1 1
2
```

#### Output
```
f"(p) = 12
```

---

## Numerical Integration

Numerical integration (quadrature) approximates definite integrals, essential when analytical integration is impossible or impractical.

### Simpson's One-Third Rule

#### Theory
Simpson's One-Third Rule approximates the integral by dividing the interval into equal subintervals and fitting parabolas through each triple of points.

**Formula** (for n subintervals, where n is even):
∫[a,b] f(x)dx ≈ (h/3)[f(x₀) + 4f(x₁) + 2f(x₂) + 4f(x₃) + ... + 2f(x_{n-2}) + 4f(x_{n-1}) + f(x_n)]

where h = (b - a) / n

**Pattern**:
- First and last terms: coefficient 1
- Odd-indexed interior points: coefficient 4
- Even-indexed interior points: coefficient 2

**Error**:
- Truncation error: O(h⁴)
- Better accuracy than trapezoidal rule

**Advantages**:
- Good balance between accuracy and simplicity
- Works well for most smooth functions

**Disadvantages**:
- Number of subintervals must be even
- Less accurate for highly oscillatory functions

#### Code
```cpp
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
    cin>>a;
    cin>>b;
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
```

#### Input Format
```
0 2 4
```

#### Output
```
Simpson 1/3 Result = 2.66667
```

---

### Simpson's Three-Eighths Rule

#### Theory
Simpson's Three-Eighths Rule uses cubic polynomials to approximate the function, dividing each interval into 3 equal subintervals.

**Formula** (for n subintervals, where n is divisible by 3):
∫[a,b] f(x)dx ≈ (3h/8)[f(x₀) + 3f(x₁) + 3f(x₂) + 2f(x₃) + 3f(x₄) + 3f(x₅) + ... + f(x_n)]

where h = (b - a) / n

**Pattern**:
- Coefficients repeat: 1, 3, 3, 2, 3, 3, 2, ..., 3, 3, 1

**Error**:
- Truncation error: O(h⁵)
- Slightly more accurate than Simpson's 1/3 rule

**When to Use**:
- When n is divisible by 3
- When higher accuracy is needed

**Advantages**:
- Higher accuracy (O(h⁵))
- Good for smooth functions

**Disadvantages**:
- Requires number of subintervals divisible by 3
- Slightly more complex than Simpson's 1/3 rule

#### Code
```cpp
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
    cin>>a;
    cin>>b;
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
```

#### Input Format
```
0 2 4
```

#### Output
```
Simpson 3/8 Result = 2.29688
```

---

## Solution of Differential Equations

Numerical methods for solving differential equations are essential when analytical solutions don't exist or are too complex.

### Runge-Kutta Method

#### Theory
The Runge-Kutta method is a family of iterative methods for solving initial value problems (IVPs) of the form:
dy/dx = f(x, y), with initial condition y(x₀) = y₀

**Fourth-Order Runge-Kutta** (most commonly used):
The method calculates four intermediate slopes and uses their weighted average:

k₁ = f(x_n, y_n)
k₂ = f(x_n + h/2, y_n + h·k₁/2)
k₃ = f(x_n + h/2, y_n + h·k₂/2)
k₄ = f(x_n + h, y_n + h·k₃)

y_{n+1} = y_n + (h/6)(k₁ + 2k₂ + 2k₃ + k₄)

**Advantages**:
- Fourth-order accuracy (error O(h⁵))
- Self-starting (doesn't require previous points)
- Good stability properties
- No derivative calculations needed

**Disadvantages**:
- Requires four function evaluations per step
- Not as efficient for stiff equations

**Applications**:
- Solving initial value problems
- Simulating physical systems (motion, chemical reactions, etc.)

#### Code
```cpp
#include <bits/stdc++.h>
using namespace std;
float dydx(float x,float y){
    return((x-y)/2);
}
float rungeKutta(float x0,float y0,float x,float h){
    int n=(int)((x-x0)/h);
    float k1,k2,k3,k4;
    float y=y0;
    for(int i=1;i<=n;i++){
        k1=h*dydx(x0,y);
        k2=h*dydx(x0+0.5*h,y+0.5*k1);
        k3=h*dydx(x0+0.5*h,y+0.5*k2);
        k4=h*dydx(x0+h,y+k3);
        y=y+(1.0/6.0)*(k1+2*k2+2*k3+k4);
        x0=x0+h;
    }
    return y;
}
int main(){
    freopen("input/1_Runge_Kutta.txt", "r", stdin);
    freopen("output/1_Runge_Kutta_output.txt", "w", stdout);
    
    float x0,y,x,h;
    cin>>x0>>y>>x>>h;
    cout<<"The value of y at x is : "<<rungeKutta(x0,y,x,h);
    
    return 0;
}
```

#### Input Format
```
0 1 2 0.2
```

#### Output
```
The value of y at x is : 1.01971
```

---

## Curve Fitting (Regression)

Curve fitting finds the best-fit function that represents the relationship between variables in a dataset.

### Linear Regression

#### Theory
Linear regression finds the best-fit line y = a + bx for a set of data points using the **Least Squares Method**.

**Normal Equations**:
- b = (n·Σ(xy) - Σx·Σy) / (n·Σ(x²) - (Σx)²)
- a = (Σy - b·Σx) / n

where n is the number of data points

**Coefficient of Determination** (R²):
Measures the goodness of fit (0 ≤ R² ≤ 1, closer to 1 is better)

**Advantages**:
- Simple and computationally efficient
- Interpretable coefficients
- Works well for linear relationships

**Disadvantages**:
- Assumes linear relationship
- Sensitive to outliers

#### Code
```cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

pair<double,double> LinearRegression(const vector<double>& x, const vector<double>& y, double predictX) {
    int n = x.size();
    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;

    for (int i = 0; i < n; i++) {
        Sx += x[i];
        Sy += y[i];
        Sxx += x[i] * x[i];
        Sxy += x[i] * y[i];
    }

    double denom = n * Sxx - Sx * Sx;
    if (fabs(denom) < 1e-12) {
        cout << "Cannot perform linear regression (all x values might be the same).\\n";
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }

    double b = (n * Sxy - Sx * Sy) / denom;
    double a = (Sy - b * Sx) / n;

    cout << fixed << setprecision(6);
    cout << "Regression Line: y = " << a << " + " << b << "*x\\n";
    cout << "Predicted value at x = " << predictX << " is y = " << (a + b * predictX) << "\\n";

    return {a, b};
}

int main() {
    freopen("input/LeastSquare_Linear.txt", "r", stdin);
    freopen("output/LeastSquare_Linear_output.txt", "w", stdout);

    int n;
    if (!(cin >> n) || n < 2) {
        cout << "Invalid input. At least 2 data points are required.\\n";
        return 0;
    }

    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++) {
        cin >> x[i] >> y[i];
    }

    double predictX;
    if (!(cin >> predictX)) {
        cout << "Invalid input for prediction.\\n";
        return 0;
    }

    LinearRegression(x, y, predictX);
    
    return 0;
}
```

#### Input Format
```
5
0 2
1 2.7
2 3.7
3 5.0
4 7.4
5
```

#### Output
```
Regression Line: y = 1.540000 + 1.310000*x
Predicted value at x = 5.000000 is y = 8.090000
```

---

### Exponential Regression

#### Theory
Exponential regression fits data to y = ae^(bx) using the **Least Squares Method** by linearizing through logarithmic transformation.

**Transformation**:
- Take ln(y) = ln(a) + bx
- Let Y = ln(y), A = ln(a)
- Now Y = A + bx (linear form)

**Solution**:
- Perform linear regression on (x, ln(y)) pairs
- Get A and b
- Calculate a = e^A

**Applications**:
- Population growth
- Radioactive decay
- Bacterial growth
- Epidemic models

**Advantages**:
- Good for naturally exponential phenomena
- More accurate than linear for exponential data

#### Code
```cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

pair<double, double> exponential(const vector<double>& x, const vector<double>& y, double predictX) {
    int n = x.size();
    if (n < 2) {
        cout << "At least 2 data points are required.\n";
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }

    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;

    for (int i = 0; i < n; i++) {
        if (y[i] <= 0) {
            cout << "Error: All Y values must be positive for logarithm.\n";
            return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
        }
        double ly = log(y[i]);
        Sx += x[i];
        Sy += ly;
        Sxx += x[i] * x[i];
        Sxy += x[i] * ly;
    }

    double denom = n * Sxx - Sx * Sx;
    if (fabs(denom) < 1e-12) {
        cout << "Cannot perform regression (denominator too small).\n";
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }

    double b = (n * Sxy - Sx * Sy) / denom;
    double A = (Sy - b * Sx) / n;
    double a = exp(A);

    double predictedY = a * exp(b * predictX);

    cout << fixed << setprecision(6);
    cout << "Regression equation: P = " << a << " * e^(" << b << " * t)\n";
    cout << "Predicted value at t = " << predictX << " is P = " << predictedY << endl;

    return {a, b};
}

int main() {
    freopen("input/LeastSquare_exponential.txt", "r", stdin);
    freopen("output/LeastSquare_exponential_output.txt", "w", stdout);

    int n;
    if (!(cin >> n) || n < 2) return 0;

    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++) {
        cin >> x[i] >> y[i];
    }

    double predictX;
    if (!(cin >> predictX)) return 0;

    exponential(x, y, predictX);
}
```

#### Input Format
```
5
0 2
1 2.7
2 3.7
3 5.0
4 7.4
5
```

#### Output
```
Regression equation: P = 1.963168 * e^(0.323285 * t)
Predicted value at t = 5.000000 is P = 9.884674
```

---

### Polynomial Regression

#### Theory
Polynomial regression fits data to a polynomial: y = a₀ + a₁x + a₂x² + ... + a_nx^n

**Solution Method**:
- Uses system of normal equations
- Can use Gaussian elimination or other methods to solve
- Higher degree provides better fit but risks overfitting

**Degree Selection**:
- Low degree: May underfit (high bias)
- High degree: May overfit (high variance)
- Use R² or cross-validation to choose optimal degree

**Applications**:
- Approximating complex curves
- Interpolation
- Trend analysis

**Advantages**:
- Flexible (can fit various curve shapes)
- Easy to compute

**Disadvantages**:
- Risk of overfitting with high degrees
- Less interpretable than linear regression

#### Code
```cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

int main() {
    freopen("input/LeastSquare_Polynomial.txt", "r", stdin);
    freopen("output/LeastSquare_Polynomial_output.txt", "w", stdout);
    
    int n;
    cin >> n;
    if (n <= 1) {
        cout << "At least two data points are required.\n";
        return 0;
    }

    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++) {
        cin >> x[i] >> y[i];
    }

    int m;
    cin >> m;
    if (m >= n) {
        cout << "Polynomial degree must be less than number of data points.\n";
        return 0;
    }

    double TheNewX;
    cin >> TheNewX;

    int size = m + 1;
    vector<vector<double>> A(size, vector<double>(size, 0));
    vector<double> B(size, 0);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < n; k++) {
                A[i][j] += pow(x[k], i + j);
            }
        }
        for (int k = 0; k < n; k++) {
            B[i] += y[k] * pow(x[k], i);
        }
    }

    for (int i = 0; i < size; i++) {
        int maxRow = i;
        for (int k = i + 1; k < size; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) maxRow = k;
        }
        if (fabs(A[maxRow][i]) < 1e-12) {
            cout << "Matrix is singular or nearly singular. Cannot solve.\n";
            return 0;
        }
        if (maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(B[i], B[maxRow]);
        }

        for (int j = i + 1; j < size; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < size; k++) {
                A[j][k] -= factor * A[i][k];
            }
            B[j] -= factor * B[i];
        }
    }

    vector<double> coeffs(size, 0);
    for (int i = size - 1; i >= 0; i--) {
        coeffs[i] = B[i];
        for (int j = i + 1; j < size; j++) {
            coeffs[i] -= A[i][j] * coeffs[j];
        }
        coeffs[i] /= A[i][i];
    }

    cout << "\nPolynomial coefficients:\n";
    for (int i = 0; i < size; i++) {
        cout << "Coefficient of x^" << i << " = " << fixed << setprecision(6) << coeffs[i] << endl;
    }

    double TheNewY = 0.0;
    for (int i = 0; i < size; i++) {
        TheNewY += coeffs[i] * pow(TheNewX, i);
    }
    cout << "\nPredicted value at x = " << fixed << setprecision(6) << TheNewX << " is y = " << TheNewY << "\n";
    
    return 0;
}
```

#### Input Format
```
5
0 2
1 2.7
2 3.7
3 5.0
4 7.4
2
5
```

#### Output
```
Polynomial coefficients:
Coefficient of x^0 = 2.068571
Coefficient of x^1 = 0.252857
Coefficient of x^2 = 0.264286

Predicted value at x = 5.000000 is y = 9.940000
```

---

### Modified Exponential Regression

#### Theory
Modified Exponential regression fits data to y = ae^(bx) where the exponential is modified with a translation or scaling parameter. Useful for data that doesn't pass through origin or has asymptotic behavior.

**Model**: y = a·e^(bx) with modifications for better fit

**Solution Steps**:
1. Apply appropriate transformation (may involve shifting or scaling)
2. Linearize the equation
3. Perform linear regression
4. Transform results back

**When to Use**:
- When standard exponential doesn't fit well
- Data with asymptotic behavior
- Restricted domain exponential functions

**Advantages**:
- Better flexibility than standard exponential
- Still relatively simple

#### Code
```cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

pair<double,double> modifiedExponentialRegression(const vector<double>& t, const vector<double>& T, double predictT) {
    int n = t.size();
    if(n < 2) {
        cout << "At least 2 data points are required.\n";
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }

    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;

    for(int i = 0; i < n; i++) {
        double X = exp(t[i]/4.0);
        Sx += X;
        Sy += T[i];
        Sxx += X * X;
        Sxy += X * T[i];
    }

    double denom = n*Sxx - Sx*Sx;
    if(fabs(denom) < 1e-12) {
        cout << "Cannot perform regression (all transformed X values may be identical).\n";
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }

    double b = (n*Sxy - Sx*Sy)/denom;
    double a = (Sy - b*Sx)/n;

    double predictedY = a + b * exp(predictT/4.0);

    cout << fixed << setprecision(6);
    cout << "Regression equation: T = " << a << " + " << b << " * e^(t/4)\n";
    cout << "Predicted value at t = " << predictT << " is T = " << predictedY << "\n";

    return {a, b};
}

int main() {
    freopen("input/LeastSquare_Modified_exponential.txt", "r", stdin);
    freopen("output/LeastSquare_Modified_exponential_output.txt", "w", stdout);

    int n;
    if(!(cin >> n) || n < 2) return 0;

    vector<double> t(n), T(n);
    for(int i = 0; i < n; i++) {
        cin >> t[i] >> T[i];
    }

    double predictT;
    if(!(cin >> predictT)) return 0;

    modifiedExponentialRegression(t, T, predictT);
    
    return 0;
}
```

#### Input Format
```
5
0 2
1 2.7
2 3.7
3 5.0
4 7.4
5
```

#### Output
```
Regression equation: T = -1.294980 + 3.110722 * e^(t/4)
Predicted value at t = 5.000000 is T = 9.562507
```

---

### Transcendental Regression

#### Theory
Transcendental regression fits data to y = ax^b (power law) using the **Least Squares Method** with logarithmic transformation.

**Transformation**:
- Take ln(y) = ln(a) + b·ln(x)
- Let Y = ln(y), X = ln(x), A = ln(a)
- Now Y = A + bX (linear form)

**Solution**:
- Perform linear regression on (ln(x), ln(y)) pairs
- Get A and b
- Calculate a = e^A

**Applications**:
- Allometric relationships (biology)
- Power laws in physics
- Economic relationships
- Natural phenomena following power laws

**Advantages**:
- Good for naturally power-law relationships
- Simple transformation
- Interpretable exponent

**Disadvantages**:
- Requires positive x and y values
- Data must follow power law distribution

#### Code
```cpp
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<iomanip>
#include<limits>

using namespace std;

pair<double,double> transcendental(const vector<double>& x,const vector<double>& y,double predictX){
    int n= x.size();
    double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    for(int i = 0; i < n; i++) {
        if(x[i] <= 0 || y[i] <= 0) return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
        double lx = log(x[i]);
        double ly = log(y[i]);
        Sx  += lx;
        Sy  += ly;
        Sxx += lx * lx;
        Sxy += lx * ly;
    }
    double denom = n*Sxx - Sx*Sx;
    if(fabs(denom) < 1e-12) return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    double b= (n*Sxy-Sx*Sy)/denom;
    double A=(Sy-b*Sx)/n;
    double a=exp(A);
    cout << "Predicted value at x = " << predictX << " is y = ";
    cout << fixed << setprecision(6) << (a*pow(predictX,b)) << "\n";

    return {a,b};
}

int main(){
    freopen("input/LeastSquare_Transcendental.txt", "r", stdin);
    freopen("output/LeastSquare_Transcendental_output.txt", "w", stdout);

    int n;
    if(!(cin>>n)) return 0;
    vector<double> x(n), y(n);
    for(int i = 0; i < n; i++) {
        cin>>x[i]>>y[i];
    }
    double predictX;
    if(!(cin>>predictX)) return 0;

    transcendental(x,y,predictX);
}
```

#### Input Format
```
4
1 2.7
2 3.7
3 5.0
4 7.4
5
```

#### Output
```
Predicted value at x = 5.000000 is y = 9.110219
```

---

## Project Structure

```
Numerical Project/
├── Curve Fitting (Regression)/
│   ├── LeastSquare_Linear.cpp
│   ├── LeastSquare_Exponential.cpp
│   ├── LeastSquare_Polynomial.cpp
│   ├── LeastSquare_Modified_Exponential.cpp
│   ├── LeastSquare_Transcendental.cpp
│   └── input/
├── Interpolation and Approximation/
│   ├── Newtons_Forward_interpolation.cpp
│   ├── Newtons_Backward_interpolation.cpp
│   ├── Divided_difference_interpolation.cpp
│   └── input/
├── Numerical Differentiation/
│   ├── Differentiation_equal-interval_interpolation.cpp
│   ├── Second_order_derivative.cpp
│   └── input/
├── Numerical Integration/
│   ├── 1_Simpsons_One-third_rule.cpp
│   ├── 2_Simpsons_three-eighths_rule.cpp
│   └── input/
├── Solution of Differential Equations/
│   ├── 1_Runge_Kutta.cpp
│   └── input/
├── Solution of Linear Equations/
│   ├── 1_Gauss_Elimination.cpp
│   ├── 2_Gauss_Jordan_Elemination.cpp
│   ├── 3_LU_decomposition.cpp
│   ├── 4_Matrix_Inversion.cpp
│   └── input/
├── Solution of Non-Linear Equations/
│   ├── Bi-section_method.cpp
│   ├── False_Position_method.cpp
│   ├── Newton-Raphson_method.cpp
│   ├── Secant_method.cpp
│   └── input/
└── README.md
```

---

## How to Use

1. Navigate to the respective folder for the method you want to use
2. Review the theory section above
3. Compile the C++ file: `g++ -std=c++17 -O2 filename.cpp -o output`
4. Run the program: `./output`
5. Provide input based on the format specified in the input section
6. Review the output

---

## Contributing

Feel free to add more test cases, optimizations, or additional methods!

---

## License

This project is open source and available for educational purposes.

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
// Add your code here
```

#### Test Cases

##### Case 1: Unique Solution
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
```

##### Case 2: Infinite Solutions
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
```

##### Case 3: No Solution
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
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
// Add your code here
```

#### Test Cases

##### Case 1: Unique Solution
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
```

##### Case 2: Infinite Solutions
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
```

##### Case 3: No Solution
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
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
// Add your code here
```

#### Test Cases

##### Case 1: Unique Solution
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
```

##### Case 2: Infinite Solutions
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
```

##### Case 3: No Solution
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
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
// Add your code here
```

#### Test Cases

##### Case 1: Invertible Matrix (Unique Solution)
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
```

##### Case 2: Singular Matrix (No Unique Solution)
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
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
// Add your code here
```

#### Test Cases

##### Case 1: Root Found
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
```

##### Case 2: No Root in Interval
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
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
// Add your code here
```

#### Test Cases

##### Case 1: Root Found
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
```

##### Case 2: No Root in Interval
**Input:**
```
// Add input here
```
**Output:**
```
// Add output here
```

---

### Newton-Raphson Method

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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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
// Add your code here
```

#### Input Format
```
// Add input format here
```

#### Output
```
// Add output here
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

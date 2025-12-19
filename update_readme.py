#!/usr/bin/env python3
"""
Script to automatically update README.md with code and test cases
"""

import os
import glob

# Define the project structure
project_root = r"c:\Users\Shahariar Emon\OneDrive\Desktop\Num\Numerical-Project"

# Map of section titles to their corresponding files
file_mappings = {
    "Gauss Elimination": {
        "cpp": "Solution of Linear Equations/1_Gauss_Elimination.cpp",
        "tests": [
            ("Unique Solution", "1_Gauss_Elimination_unique_solution"),
            ("Infinite Solutions", "1_Gauss_Elimination_infinite_solutions"),
            ("No Solution", "1_Gauss_Elimination_no_solution")
        ]
    },
    "Gauss Jordan": {
        "cpp": "Solution of Linear Equations/2_Gauss_Jordan_Elemination.cpp",
        "tests": [
            ("Unique Solution", "2_Gauss_Jordan_Elemination_unique_solution"),
            ("Infinite Solutions", "2_Gauss_Jordan_Elemination_infinite_solutions"),
            ("No Solution", "2_Gauss_Jordan_Elemination_no_solution")
        ]
    },
    "LU Decomposition": {
        "cpp": "Solution of Linear Equations/3_LU_decomposition.cpp",
        "tests": [
            ("Unique Solution", "3_LU_decomposition_unique_solution"),
            ("Infinite Solutions", "3_LU_decomposition_infinite_solutions"),
            ("No Solution", "3_LU_decomposition_no_solution")
        ]
    },
    "Matrix Inversion": {
        "cpp": "Solution of Linear Equations/4_Matrix_Inversion.cpp",
        "tests": [
            ("Invertible Matrix", "4_Matrix_Inversion_invertible"),
            ("Singular Matrix", "4_Matrix_Inversion_singular")
        ]
    },
    "Bisection": {
        "cpp": "Solution of Non-Linear Equations/Bi-section_method.cpp",
        "tests": [("Test Case", "Bi-section_method")]
    },
    "False Position": {
        "cpp": "Solution of Non-Linear Equations/False_Position_method.cpp",
        "tests": [("Test Case", "False_Position_method")]
    },
    "Newton-Raphson": {
        "cpp": "Solution of Non-Linear Equations/Newton-Raphson_method.cpp",
        "tests": [("Test Case", "Newton-Raphson_method")]
    },
    "Secant": {
        "cpp": "Solution of Non-Linear Equations/Secant_method.cpp",
        "tests": [("Test Case", "Secant_method")]
    },
    "Newton Forward": {
        "cpp": "Interpolation and Approximation/Newtons_Forward_interpolation.cpp",
        "tests": [("Test Case", "Newtons_Forward_interpolation")]
    },
    "Newton Backward": {
        "cpp": "Interpolation and Approximation/Newtons_Backward_interpolation.cpp",
        "tests": [("Test Case", "Newtons_Backward_interpolation")]
    },
    "Divided Difference": {
        "cpp": "Interpolation and Approximation/Divided_difference_interpolation.cpp",
        "tests": [("Test Case", "Divided_difference_interpolation")]
    },
    "Differentiation Equal Interval": {
        "cpp": "Numerical Differentiation/Differentiation_equal-interval_interpolation.cpp",
        "tests": [("Test Case", "Differentiation_equal-interval_interpolation")]
    },
    "Second Derivative": {
        "cpp": "Numerical Differentiation/Second_order_derivative.cpp",
        "tests": [("Test Case", "Second_order_derivative")]
    },
    "Simpson One-Third": {
        "cpp": "Numerical Integration/1_Simpsons_One-third_rule.cpp",
        "tests": [("Test Case", "1_Simpsons_One-third_rule")]
    },
    "Simpson Three-Eighths": {
        "cpp": "Numerical Integration/2_Simpsons_three-eighths_rule.cpp",
        "tests": [("Test Case", "2_Simpsons_three-eighths_rule")]
    },
    "Runge Kutta": {
        "cpp": "Solution of Differential Equations/1_Runge_Kutta.cpp",
        "tests": [("Test Case", "1_Runge_Kutta")]
    },
    "Linear Regression": {
        "cpp": "Curve Fitting (Regression)/LeastSquare_Linear.cpp",
        "tests": [("Test Case", "LeastSquare_Linear")]
    },
    "Exponential Regression": {
        "cpp": "Curve Fitting (Regression)/LeastSquare_exponential.cpp",
        "tests": [("Test Case", "LeastSquare_exponential")]
    },
    "Polynomial Regression": {
        "cpp": "Curve Fitting (Regression)/LeastSquare_Polynomial.cpp",
        "tests": [("Test Case", "LeastSquare_Polynomial")]
    },
    "Modified Exponential": {
        "cpp": "Curve Fitting (Regression)/LeastSquare_Modified_exponential.cpp",
        "tests": [("Test Case", "LeastSquare_Modified_exponential")]
    },
    "Transcendental": {
        "cpp": "Curve Fitting (Regression)/LeastSquare_Transcendental.cpp",
        "tests": [("Test Case", "LeastSquare_Transcendental")]
    }
}

def read_file(filepath):
    """Read file content"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            return f.read()
    except:
        return "// File not found"

def extract_code_section(cpp_content):
    """Extract just the code without comments at the end"""
    lines = cpp_content.split('\n')
    code_lines = []
    for line in lines:
        if line.strip().startswith('/*') and 'input' in line.lower():
            break
        code_lines.append(line)
    return '\n'.join(code_lines).rstrip()

print("README update helper")
print("=" * 50)

# Example: Print code and test cases for Gauss Elimination
for method_name, info in file_mappings.items():
    cpp_path = os.path.join(project_root, info["cpp"])
    
    print(f"\n### {method_name}")
    print(f"CPP File: {cpp_path}")
    
    if os.path.exists(cpp_path):
        code = read_file(cpp_path)
        clean_code = extract_code_section(code)
        print(f"Code length: {len(clean_code)} characters")
        
        # Get directory for input/output
        dir_path = os.path.dirname(cpp_path)
        
        for test_name, test_file in info["tests"]:
            input_path = os.path.join(dir_path, "input", f"{test_file}.txt")
            output_path = os.path.join(dir_path, "output", f"{test_file}_output.txt")
            
            if os.path.exists(input_path):
                input_content = read_file(input_path)
                print(f"  {test_name} - Input: {len(input_content)} chars")
            
            if os.path.exists(output_path):
                output_content = read_file(output_path)
                print(f"  {test_name} - Output: {len(output_content)} chars")
    else:
        print(f"  ERROR: File not found - {cpp_path}")

print("\n" + "=" * 50)
print("Due to the large size, manual insertion into README is recommended.")
print("Copy code blocks from the C++ files and paste into README.md")

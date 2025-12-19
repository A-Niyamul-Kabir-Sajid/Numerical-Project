# PowerShell script to collect all code, inputs, and outputs
# Run this to see all content that needs to go into README.md

$projectRoot = "c:\Users\Shahariar Emon\OneDrive\Desktop\Num\Numerical-Project"

Write-Host "========== GAUSS ELIMINATION ==========" -ForegroundColor Green
Write-Host "`n#### Code"
Get-Content "$projectRoot\Solution of Linear Equations\1_Gauss_Elimination.cpp" -Raw | Write-Host

Write-Host "`n##### Case 1: Unique Solution"
Write-Host "**Input:**"
Get-Content "$projectRoot\Solution of Linear Equations\input\1_Gauss_Elimination_unique_solution.txt" -Raw | Write-Host
Write-Host "**Output:**"
Get-Content "$projectRoot\Solution of Linear Equations\output\1_Gauss_Elimination_unique_solution_output.txt" -Raw | Write-Host

Write-Host "`n##### Case 2: Infinite Solutions"
Write-Host "**Input:**"
Get-Content "$projectRoot\Solution of Linear Equations\input\1_Gauss_Elimination_infinite_solutions.txt" -Raw | Write-Host
Write-Host "**Output:**"
Get-Content "$projectRoot\Solution of Linear Equations\output\1_Gauss_Elimination_infinite_solutions_output.txt" -Raw | Write-Host

Write-Host "`n##### Case 3: No Solution"
Write-Host "**Input:**"
Get-Content "$projectRoot\Solution of Linear Equations\input\1_Gauss_Elimination_no_solution.txt" -Raw | Write-Host
Write-Host "**Output:**"
Get-Content "$projectRoot\Solution of Linear Equations\output\1_Gauss_Elimination_no_solution_output.txt" -Raw | Write-Host

Write-Host "`n`n========== INSTRUCTIONS ==========" -ForegroundColor Cyan
Write-Host "1. Run this script and redirect output to a text file:"
Write-Host "   .\collect_all_codes.ps1 > all_content.txt"
Write-Host "`n2. Open all_content.txt and README.md side by side"
Write-Host "`n3. Copy-paste the code blocks into the appropriate sections"
Write-Host "`n4. Repeat for all other methods by adding similar blocks to this script"

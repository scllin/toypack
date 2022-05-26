# mrMed
mrMed (**MR-based mediation analysis**) is a method that performs causal medaiton analysis based on Genome-wide summary statistics using Mendelian Randomization (MR).

mrMed provides three main functions to perform MR-based mediation analysis and provide estimates of total effect (TE), direct effect (DE), indirect effect (IE), and mediation proportion (rho):
1. Diff_IVW
2. Prod_IVW
3. Prod_Median

### Reference
Causal mediation analysis: a MR approach
<https://...>

### 1. Install and load mrMed
To install and load the latest version from GitHub, run the following lines:
```r
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("scllin/mrMed")
```
Load mrMed 
```r
library(mrMed)
```

### 2. Example
```r
# Load the dataset with exposure:WHR, mediator:T2D, and outcome:CAD
data(WHR_T2D_CAD)

# Run mrMed with the default methods: Diff_IVW, Prod_IVW, Prod_Median
mrMed(WHR_T2D_CAD)

# Run mrMed with the the other methods
mrMed(WHR_T2D_CAD, method_list=c("Diff_IVW_0","Prod_IVW_0"))

# Run Diff-Median method (may require certain run time due to the bootstrap procedure)
mrMed(WHR_T2D_CAD, method_list=c("Diff_Median"))
```
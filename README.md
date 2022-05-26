# mrMed
mrMed (**MR-based mediation analysis**) is a method that performs causal medaiton analysis based on summary-data Mendelian Randomization (MR).

mrMed provides three main functions to perform MR-based mediation analysis and provide estimates of total effect (TE), direct effect (DE), indirect effect (IE), mediation proportion (rho):
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

# Run mrMed with the default three methods Diff_IVW, Prod_IVW, Prod_Median
mr_presso(WHR_T2D_CAD)

# Run MR-PRESSO on a multi-variable MR (MMR) model specifying several exposures
mr_presso(BetaOutcome = "Y_effect", BetaExposure = c("E1_effect", "E2_effect"), SdOutcome = "Y_se", SdExposure = c("E1_se", "E2_se"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)
```
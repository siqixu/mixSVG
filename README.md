# mixSVG: An Omnibus Test of Mixed Effects for Detecting Spatially Variable Genes with Spatial Transcriptomics Data
The "mixSVG" package implements an omnibus test to detect spatially variable genes with spatial transcriptomics data, which tests for fixed effect (mean-level) and randome effect (variance-level) of spatial heteogeneity of gene expression simultaneously.

## Setup
The "mixSVG" package depends on R (>= 3.5.0). Use the following command in R to install the "mixSVG" package:
```
library(devtools)
install_github("siqixu/mixSVG",ref="main") # install the "DeepMed" package
```
## Usage
```
mixSVG(count,
       coord,
       X = NULL,
       libsize_inc = TRUE,
       libsize = NULL,
       vtest_zero_prop = 0.995,
       ncore = 10,
       n_perm = 1000,
       sig = 0.05)
```
`count`: 
A q by n numeric matrix of gene expression counts for q genes and n spots.

`coord`: 
A n by 2 numeric matrix of two-dimensioanl spatial coordinates for n spots.

`X`: 
A n by p numeric matrix of p covariates for n spots.

`libsize_inc`: 
Whether to account for the library size. `libsize_inc=TRUE` by default.

`libsize`: 
A numeric vector of the library size for n spots. If `libsize_inc=TRUE`, then libsize will be the total expression counts of all genes on each spot by default. If `libsize_inc=FALSE`, then libsize=1 by default.

`vtest_zero_prop`: A numeric value bewteen 0 and 1. The mixed effects (fixed effect and random effect) will be tested when the proportion of zero expression counts across all spots is less than `vtest_zero_prop`. Otherwise, only the fixed effect will be tested.

`ncore`:
Number of cores used for parallel computation

`n_perm`:	
Number of permutation replicates for testing the random effect.

`sig`:	
Significance level for dectecting spatially variable genes based on the adjusted P-values of Benjamini-Hochberg method.

## Value
`results`:	A list of results for all genes, including:

`model0`: The estimation results under the null model, which contain the estimated coefficients of covariates (beta), the estimated variance of residual (tau), indicator of the estimation convergence (converge), numbder of iterations (iter);

`pval`: The combined P-value of mixSVG based on 13 kinds of transformations of spatial coordinates accounting for different spatial patterns;

`pval_pat`: A 13 by 3 matrix of P-values of mixSVG, where 13 rows represent 13 transformations of spatial coordinates, and three columns represent the test for mixed effect, fixed effect and random effect, respectively.

`pval_all`:	
A matrix of P-values for all genes. The first column contains the original P-values of mixSVG, and the second column contains the P-values adjusted by Benjamini-Hochberg method.

`pval_sig`:	
A matrix of P-values for the detected spatially variable genes based on the adjusted P-values of Benjamini-Hochberg method. The first column contains the original P-values of mixSVG, and the second column contains the P-values adjusted by Benjamini-Hochberg method.











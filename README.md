# mixSVG: An Omnibus Test of Mixed Effects for Detecting Spatially Variable Genes with Spatial Transcriptomics Data
The "mixSVG" package implements an omnibus test to detect spatially variable genes with spatial transcriptomics data, which tests for fixed effect (mean-level) and random effect (variance-level) of spatial heterogeneity of gene expression simultaneously.

## Setup
The "mixSVG" package depends on R (>= 3.5.0). Use the following command in R to install the "mixSVG" package:
```
library(devtools)
install_github("siqixu/mixSVG",ref="main") # install the "mixSVG" package
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
An n by 2 numeric matrix of two-dimensional spatial coordinates for n spots.

`X`: 
A n by p numeric matrix of p covariates for n spots.

`libsize_inc`: 
Whether to account for the library size. `libsize_inc=TRUE` by default.

`libsize`: 
A numeric vector of the library size for n spots. If `libsize_inc=TRUE`, then libsize will be the total expression counts of all genes on each spot by default. If `libsize_inc=FALSE`, then libsize=1 by default.

`vtest_zero_prop`: A numeric value between 0 and 1. The mixed effects (fixed effect and random effect) will be tested when the proportion of zero expression counts across all spots is less than `vtest_zero_prop`. Otherwise, only the fixed effect will be tested.

`ncore`:
Number of cores used for parallel computation

`n_perm`:	
Number of permutation replicates for testing the random effect.

`sig`:	
Significance level for detecting spatially variable genes based on the adjusted P-values of the Benjamini-Hochberg method.

## Value
`results`:	A list of results for all genes, including:

* `model0`: The estimation results under the null model, which contain the estimated coefficients of covariates (`beta`), the estimated variance of residual (`tau`), an indicator of the estimation convergence (`converge`), and the number of iterations (`iter`);

* `pval`: The combined P-value of mixSVG based on 13 kinds of transformations of spatial coordinates accounting for different spatial patterns;

* `pval_pat`: A 13 by 3 matrix of P-values of mixSVG, where 13 rows represent 13 transformations of spatial coordinates, and three columns represent the test for mixed effect, fixed effect and random effect, respectively.

`pval_all`:	
A matrix of P-values for all genes. The first column contains the original P-values of mixSVG, and the second column contains the P-values adjusted by the Benjamini-Hochberg method.

`pval_sig`:	
A matrix of P-values for the detected spatially variable genes (whose adjusted P-values > `sig`). The first column contains the original P-values of mixSVG, and the second column contains the P-values adjusted by the Benjamini-Hochberg method.

## Example
```
library(mixSVG)

# Example data of human breast cancer
data(example)
dim(rawcount) # the expression count matrix of 14,789 genes and 251 spots 
dim(rawcoord) # the matrix of two-dimensional spatial coordinates for 251 spots.

# Data preprocessing
## Select spots with at least one expression count
keep_spot <- which(colSums(rawcount) >= 1)
coord = rawcoord[keep_spot,]

## Select genes expressed in at least 2% of spots  
count = rawcount[rowMeans(rawcount>0) >= 0.02, keep_spot]

# Run mixSVG (It takes around 2.5 mins with 20 cores, and around 30 mins with a single core.)
mixSVG_output = mixSVG(count, coord, ncore = 20)

# P-values and adjusted P-values (by Benjamini-Hochberg method) for all genes
View(mixSVG_output$pval_all) 

# P-values and adjusted P-values (by Benjamini-Hochberg method) for the detected spatially variable genes (with adjusted P-values >0.05)
View(mixSVG_output$pval_sig)


# On the other hand, one may investigate the robustness of mixSVG under the null by:

## 1. Generate the null (i.e., all genes are non-spatially variable) by shuffling the locations of spots
n = ncol(count)
pos = colnames(count)
perm_index = sample(1:n, n, replace = F)
count0 = count[, perm_index]
colnames(count0) = pos

## 2. run mixSVG
mixSVG_output0 = mixSVG(count0, coord, ncore = 20)

## 3. Generate the Q-Q plot for P-values under the null
library(gaston)
qqplot.pvalues(mixSVG_output0$pval_all[,1], cex = 0.1)

## Note: If the Q-Q plot indicates that mixSVG inflates under the null, it is suggested to decrease vtest_zero_prop (e.g., 0.98) or increase the gene selection threshold, e.g., select genes expressed in at least 5% of spots.  


```









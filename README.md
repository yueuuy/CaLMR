# CaLMR

CaLMR (<u>Ca</u>usal analysis of <u>L</u>atent exposures using <u>M</u>endelian <u>R</u>andomization) is a Bayesian MR method to test the causal relationships between latent exposures and an outcome using GWAS summary-level association statistics of multiple traits co-regulated by the exposures. Both univariable and multivariable versions of CaLMR are available.
## Installation
``` R
# if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("yueuuy/CaLMR")
library(MASS)
library(LaplacesDemon)
```

## Examples
### CaLMR (Uni)
Assume there are six observable traits, named B1 to B6, and an outcome variable called Y. In the ``sumtable``, the columns containing the summary coefficients should be labeled as **"beta.B1"**, \..., **"beta.B6"**, and **"beta.Y"**, and the columns containing the summary standard error should be labeled as **"se.B1"**, \..., **"se.B6"**, and **"se.Y"**. The column and row names of the $1^{st}$ to $K^{th}$ variables in the input correlation matrix should match with the order in ``traitvec``.
``` R
# load the simulated summary data 
data(sample_data)
colnames(sample_data) # check how to label the column names
# load the simulated correlation matrix
data(samplecorr)
colnames(sample_data)
# check the column names and the rownames of the data match with the orders in traitvec
traitvec = paste0("B", 1:6)
identical(rownames(samplecorr), colnames(samplecorr), traitvec)
# assume the signs of thetak's are all positive
calmr_uni(sumtable=sample_data, Corr.mat=samplecorr, K = 6, traitvec, outcome="Y", sign=rep(1,K), T, burnin)
```

### CaLMR (Multi)
For the same six observable traits, B1 to B6, assume they are regulated by two latent exposures. Specifically, B1, B2, and B4 regulated by one latent exposure, while B2, B4, B5, and B6 regulated by the other. The sample codes for conducting CaLMR: 
``` R
# Assume all theta_kl are positive
calmr_multi(sumtable=sample_data, Corr.mat=samplecorr, grp=list(c("B1", "B2", "B4"), c("B2", "B4", "B5", "B6")),
            L=2, K=6, traitvec, outcome="Y", sign=list(c(1,1,0,1,0,0), c(0,1,0,1,1,1)), T, burnin)

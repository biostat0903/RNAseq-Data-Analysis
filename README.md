# RNAseq-Data-Analysis
The repository includes the algorithms for RNA-Seq data analysis, including differentially expressed (DE) and gene-pathway interactions (GPI). 

## [isoVCT](https://peerj.com/articles/3797/)

Variance Component Test with isoform (isoVCT) is an algorithm for identifying the DE genes for next generation seqencing (NGS) data. It can also automatically select a suitable distribution (Poisson or Negative Binomial) of each gene. 

### Getting Started 
First, you will need to install [R](https://cran.r-project.org/mirrors.html) as well as the packages listed under the requirements header below. All of the packages you can download from CRAN package. <br>
The packages are as following： `lme4, pipeR, MASS`.

### Support
Before contacting us, please try the following: <br>
* The wiki has tutorials on RNA-Seq data, generalized linear mixed model (GLMM), variance component test and score statistics.<br>
* The methods are described in the [papers](https://peerj.com/articles/3797/). <br>

If that doesn't work, you can get in touch with us via the email.

### Citation
If you use the software, please cite <br>
[Sheng Yang, Fang Shao, Weiwei Duan, Yang Zhao, Feng Chen. Variance component testing in identifying the differentially expressed genes from RNA-seq data. PeerJ. 5:e3797.](https://peerj.com/articles/3797/)

--------------------
## PEA
Permutaion based gEne-pAthway interactions identification in binary phenotype (PEA) is an algorithm for testing the gene-pathway interaction.
### Installation 
The code depends on `Rcpp` package.
```R
#install packages
install.packages("Rcpp")
#library Rcpp package
library(Rcpp)
#complie the cpp code
sourceCpp("yourpath/PEA.cpp")
```
### Example
We have upload the example data (example.Rdata).<br>
* `example[[1]]`: the binary phenotye.<br>
* `example[[2]]`: the covariate.<br>
* `example[[3]]`: the expresion level of genes in the core pathway.<br>
* `example[[4]]`: the expresionlevel of potential gene.<br>
### Testing of the GPI
```R
PEA_res <- PEA(example[[1]], example[[2]], example[[3]], example[[4]])
[[1]] 
[[2]] 
[[3]] 
```
### Support
Before contacting us, please try the following: <br>
* The methods are described in the [papers](https://peerj.com/articles/3797/). <br>

If that doesn't work, you can get in touch with us via the email.

### Citation
If you use the software, please cite <br>


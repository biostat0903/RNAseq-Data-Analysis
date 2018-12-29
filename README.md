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
The code depends on `Rcpp` package. If your computer do not install the `Rcpp` package, please install it first. 
For the linear algebra, we use the [Armadillo](http://arma.sourceforge.net/) library. 

```R
#install packages (without Rcpp package)
install.packages("Rcpp")
#library Rcpp package
library(Rcpp)
#complie the Cpp code in R environment (yourpath: the path of PEA in your computer )
sourceCpp("yourpath/PEA.cpp")
```
### Example
We have upload the example data (simData.Rdata).<br>
* `simData[[1]]`: vector (n by 1), the binary phenotye.<br>
* `simData[[2]]`: matrix (n by m), the covariate.<br>
* `simData[[3]]`: matrix (n by p), the expresion level of genes in the core pathway.<br>
* `simData[[4]]`: vector (n by 1), the expresionlevel of potential gene.<br>
### Testing of the GPI
```R
> PEA_res <- PEA(simData[[1]], simData[[2]], simData[[4]], simData[[3]], 10000)
> PEA_res
$U
[1] 41.57699
$P
[1] 0.0592
$betaHat
          [,1]
[1,] -2.400715

```
### Support
Before contacting us, please try the following: <br>
* The methods are described in the [papers](https://peerj.com/articles/3797/). <br>

If that doesn't work, you can get in touch with us via the email.

### Citation
If you use the software, please cite <br>


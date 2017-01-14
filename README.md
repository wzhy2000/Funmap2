# Funmap2

Functional Mapping Package for R (Version 2.4)

LSKAT is an R program that performs association testing between a set of SNPs and a longitudinal quantitative trait in population samples. The main reference for this program is

# Reference

[1] Ma, C. X., G. Casella, and R. L. Wu, 2003. Functional mapping of quantitative trait loci underlying the character process: A theoretical framework. Genetics 161(4): 1751-1762.

## Abstract:

The Funmap package is developed to identify quantitative trait loci (QTL) for a longitudinal, or vectorized, phenotypic trait as based on the Funmap model[1]. Version 2.2 has implemented 5 functional curves( Logistic Curve, Bi-Exponential Curve, Pharmacology Curve, Nonparametric method, Exponentiation Curve) and three cross types(backcross, F2 cross and RILs by selfing). the Function FM2.qtlscan is the easiest way to call the computational model by five arguments of the phenotype file, the genotype file, the marker file, the curve type and the cross type. The computational model automatically performs the hypothesis tests and permutation. The package can output a brief report for the raw data, results of the QTL scanning and permutation results. For some results, figures will be outputted to a PDF file and a Rdata file will be generated. For more details, please refer to the document. 

## Installation Instructions:

### Required software and packages
    
1. R (http://www.r-project.org/)
    
2. Package SKAT, mvtnorm, snpStats, snowfall.
    
Please install the required R package before you install LSKAT package. After the  installation of `SKAT`, `mvtnorm`, `snpStats` and `snowfall` package, please install the **LSKAT** as following steps.

 
### Install LSKAT on LUNIX or Mac OSX

```
git clone https://github.com/wzhy2000/Funmap2.git

cd Funmap2

R CMD INSTALL Funmap2

```

### Install LSKAT on Windows

1) Please download windows package from (https://github.com/wzhy2000/Funmap/blob/master/Funmap2-win.zip)

2) Install the package in R GUI by selecting the menu "Packages|Install package(s) from local zip files..."

##Usage instructions##

LSKAT is an R package which provides the main functions:

> 1) NULL model estimation to estimate the model parameters and the residuals.

> 2) Gene association test by LSKAT

> 3) LSKAT pipeline for PLINK data set with longitudial phenotype traits.


In general, two ways are common to do data analysis by LSKAT, one is test the genetic association between the longitudinal phenotype traits and SNP-set genotypic matrix, i.e, test LSKAT using the SNP-set matrix given the estimated NULL model which assumes no genetic effects contibute to phenotype traits.


```
## NULL model estimation
r.nodel <- longskat_est_model( phe.long.matrix, phe.cov.matrix, phe.time.matrix);

## Gene association test
r.lskat <- longskat_gene_test( r.model, snp.mat);
```

Other way is to run LSKAT on whole PLINK data set using the pipeline provided by the function `longskat_gene_plink`


```
r.lskat <- longskat_gene_plink( file.plink.bed,  file.plink.bim,  file.plink.fam,  
    file.phe.long,  file.phe.cov, NULL,  file.gene.set );
```

All functions and examples in the LSKAT are available in the manual (https://github.com/ZWang-Lab/LSKAT/blob/master/manual.pdf).

Sample Script:
The following source shows how to call the main function in R.

```
library(Funmap2);
FM2.qtlscan( "pheno.csv", "geno.csv", "marker.csv", CURVE_LC, CROSS_BC );
```

The following source shows how to call the main function in R.

```
library(Funmap2);
FM2.simu_test( par_NP, CURVE_NP, CROSS_F2 );
```

The following source shows how to call the main function in R.

```
library(Funmap2);
par <- FM2.param( par_LC, CURVE_LC, CROSS_BC );
dat <- FM2.simulate(par);
ret <- FM2.qtlmodel(dat);
FM2.report( "report.pdf", dat, ret );
```

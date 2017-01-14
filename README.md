# Funmap2 (Version 2.4)

Functional Mapping Package, is an R program that performs QTL mapping between a set of markers and a longitudinal quantitative trait in population samples. 

# Reference

[1] Ma, C. X., G. Casella, and R. L. Wu, 2003. Functional mapping of quantitative trait loci underlying the character process: A theoretical framework. Genetics 161(4): 1751-1762.

## Abstract:

The Funmap package is developed to identify quantitative trait loci (QTL) for a longitudinal, or vectorized, phenotypic trait as based on the Funmap model[1]. Version 2.4 has implemented 9 functional curves( Logistic Curve, Bi-Exponential Curve, Pharmacology Curve, Nonparametric method, Exponentiation Curve) , 13 covariance structures and 3 cross types (backcross, F2 cross and RILs by selfing). the Function FM2.pipe is the easiest way to call the computational model by five arguments of the phenotype file, the genotype file, the marker file, the cross type and the curve type. The computational model automatically performs the hypothesis tests and permutation. The package can output a brief report for the raw data, results of the QTL scanning and permutation results. For some results, figures will be outputted to a PDF file and a Rdata file will be generated. For more details, please refer to the document. 

## Document

> 1) Vinegette (https://github.com/wzhy2000/Funmap2/blob/master/manual.pdf)

> 2) Manual (https://github.com/wzhy2000/Funmap2/blob/master/manual.pdf)

## Installation Instructions:

### Required software and packages
    
1. R (http://www.r-project.org/)
    
2. Package mvtnorm, parallel (required in R <= 3.0).

Please install the required R package before you install Funmap package. After the  installation of `mvtnorm`, please install the **Funmap2** as following steps.

 
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

Funmap2 is an R package which provides the main functions:

> 1) QTL mapping for the experiment data .

> 2) Gene association test by LSKAT

> 3) LSKAT pipeline for PLINK data set with longitudial phenotype traits.


In general, two ways are common to do data analysis by LSKAT, one is test the genetic association between the longitudinal phenotype traits and SNP-set genotypic matrix, i.e, test LSKAT using the SNP-set matrix given the estimated NULL model which assumes no genetic effects contibute to phenotype traits.

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
All functions and examples in the LSKAT are available in the manual (https://github.com/ZWang-Lab/LSKAT/blob/master/manual.pdf).

Sample Report:

Link: http://statgen.psu.edu/software/funmap/report_demo.2.2.pdf


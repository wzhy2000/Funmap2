# Funmap2 (Version 2.4.2)

Functional Mapping Package, is an R program that performs QTL mapping between a set of markers and a longitudinal quantitative trait in population samples. 

The package has been moved from http://statgen.psu.edu/software/funmap.html to here.

# Reference

[1] Ma, C. X., G. Casella, and R. L. Wu, 2003. Functional mapping of quantitative trait loci underlying the character process: A theoretical framework. Genetics 161(4): 1751-1762.

## Abstract:

The Funmap2 package is developed to identify quantitative trait loci (QTL) for a longitudinal, or vectorized, phenotypic trait as based on the Funmap model[1]. Version 2.4 has implemented 9 functional curves( Logistic Curve, Bi-Exponential Curve, Pharmacology Curve, Nonparametric method, Exponentiation Curve) , 13 covariance structures and 3 cross types (backcross, F2 cross and RILs by selfing). the Function **FM2.pipe()** is the easiest way to call the computational model by five arguments of the phenotype file, the genotype file, the marker file, the cross type and the curve type. The computational model automatically performs the hypothesis tests and permutation. The package can output a brief report for the raw data, results of the QTL scanning and permutation results. For some results, figures will be outputted to a PDF file and a Rdata file will be generated. For more details, please refer to the document. 

## Workflow

![Hi](https://github.com/wzhy2000/Funmap2/raw/master/Funmap2-workflow.png?v=4&s=200 "Funmap2 workflow")

## Document

> 1) Vignette (https://github.com/wzhy2000/Funmap2/blob/master/Funmap2-vignette-v2.4.pdf)

> 2) Manual (https://github.com/wzhy2000/Funmap2/blob/master/Funmap2-manual.pdf)

## Installation Instructions:

### Required software and packages
    
> 1. R (http://www.r-project.org/)
    
> 2. Package mvtnorm, parallel (required in R <= 3.0).

Please install the required R package before you install Funmap2 package. After the  installation of `mvtnorm`, please install the **Funmap2** as following steps.

### Install Funmap2 on LINUX or Mac OSX

> 1. use install_github function in R console

```
library("devtools");
install_github("wzhy2000/Funmap2/Funmap2")
```
> 2. use R command in a command window

```
git clone https://github.com/wzhy2000/Funmap2.git
cd Funmap2
R CMD INSTALL Funmap2
```

### Install Funmap2 on Windows

> 1. use install_github function in R console

```
library("devtools");
install_github("wzhy2000/Funmap2/Funmap2")
```
> 2. use the menu "Packages|Install package"

1) Please download windows package from (https://github.com/wzhy2000/Funmap2/raw/master/windows/Funmap2.zip)

2) Rename the downloaded package to **Funmap2.zip** (It would be failed if long name package is installed directly)

3) Install the package in R GUI by selecting the menu "Packages|Install package(s) from local zip files..."

## Usage Instructions

Funmap2 is an R package which provides:

> 1) pipeline to do QTL mapping for the experiment data.

> 2) pipeline to do simulation test.

> 3) data loading, estimation, QTL scaning, permutation and report functions.

In general, pipeline is a common way to do data analysis by Funmap2, the following codes show how to call the pipeline in R.

```
library(Funmap2);
FM2.pipe( "pheno.csv", "time.csv", "geno.csv", "marker.csv", "BC", "Logistic" );
```

The following codes show how to call the simulation test pipeline in R.

```
library(Funmap2);
FM2.simu.pipe();
```

The following codes show how to call the main function in R.

```
library(Funmap2);
dat <- FM2.simulate();
ret <- FM2.qtlscan(dat);
FM2.report( "report.pdf", dat, ret );
```
All functions and examples in the Funmap2 are available in the manual (https://github.com/wzhy2000/Funmap2/blob/master/Funmap2-manual.pdf).

Sample Report:

Link: https://github.com/wzhy2000/Funmap2/blob/master/demo/report_demo.2.4.pdf


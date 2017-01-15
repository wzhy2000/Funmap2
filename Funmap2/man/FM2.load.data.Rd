\name{FM2.load.data}
\alias{FM2.load.data}

\title{ Loading data from the experiment files}

\description{ 
Load the experimental data from experiment files.
} 

\usage{ 
FM2.load.data( pheno.csv, time.csv, geno.csv, marker.csv, cross.type, 
    curve.type = NULL, 
    covar.type = NULL, 
    pdf.file = NULL, 
    log = FALSE ) 
}

\arguments{

\item{pheno.csv}{a CSV file of phenotypic traits . The format is described in the details section of \code{\link{FM2.load.data}}.}
\item{time.csv}{a CSV file of measured time. The format is described in the details section of \code{\link{FM2.load.data}}.}
\item{geno.csv}{a CSV file of genotype marker. The format is described in the details section of \code{\link{FM2.load.data}}.}
\item{marker.csv}{a CSV file of marker definition. The format is described in the details section of \code{\link{FM2.load.data}}}
\item{cross.type}{string indicating the cross type, optional values are "F2", "BC" and "RIL".}
\item{curve.type}{string indicating the curve type, optional values are "auto", "Logistic", "Exponential", "Power", "Legendre2", .., the full list of curve is described in the details section of \code{\link{FM2.get.curve}}. . 'auto' or NULL force the function to do curve fitting by calling \code{\link{FM2.estimate.data}}}
\item{covar.type}{string indicating the covariance type, including "auto", "AR1", "SAD1", "ARMA", "CS", .., the full list of covariance structure is described in the details section of \code{\link{FM2.get.covariance}}. 'auto' or NULL force the function to estimate the covariance structure by calling \code{\link{FM2.estimate.data}}}
\item{pdf.file}{string variable suggesting a PDF file name to illustrate the performance of curve fitting.}
\item{log}{boolean variable indicating whether logarithm is applied to the phenotype data.}
}

\value{
This function returns a S3 object with the class label of \code{FM2.dat} which structure is identical with the result from \code{\link{FM2.simulate}}.

\item{obj.cross}{the cross object.}
\item{obj.curve}{the curve object.}
\item{obj.covar}{the covariance structure(object).}
\item{obj.gen}{the genotype data, including \code{geno.csv}, \code{marker.csv}, \code{marker.obj}, \code{marker.table}, \code{genos.matrix} .}
\item{obj.phe}{the phenotype list, each item is explained in below section.}
\item{obj.phe$pheno.csv}{the phenotype file.} 
\item{obj.phe$time.csv}{the measured time file.}
\item{obj.phe$sample.obs}{sample size} 
\item{obj.phe$sample.times}{the measure times}
\item{obj.phe$log}{boolean value indicating whether log function is applied to the phenotype data}
\item{obj.phe$pheY}{matrix, longitudinal phenotype traits}
\item{obj.phe$pheX}{matrix, covariate data for all individuals.}
\item{obj.phe$pheT}{matrix, measured times for all individuals.}
\item{obj.phe$est.covar}{the estimation for covariance structure obtained from the call of \code{\link{FM2.estimate.data}}.}
\item{obj.phe$est.curve}{the estimation for curve object obtained from the call of \code{\link{FM2.estimate.data}}.}
\item{obj.phe$summary.curve}{the curve fitting results for all selected curves, if \code{curve.type} is NULL or 'auto', all available curves are estimated.}
\item{sobj.phe$ummary.covar}{the results of covariance estimation for all selected covariances, if \code{covar.type} is NULL or 'auto', all available covariances are estimated.}

}

\details{


The function returns a data object,which the structure is same as the simulation data, 
You can use \code{print} or \code{str} command to check the details in the data object.
The following conents are exported by \code{str} command


\preformatted{

List of 6
 $ obj.curve:Formal class 'fg.curve.log' [package "Funmap2"] with 2 slots
 $ obj.covar:Formal class 'fg.covariance.SAD1' [package "Funmap2"] with 3 slots
 $ obj.cross:List of 14
 $ obj.gen  :List of 5
  ..$ geno.csv    : chr "../populus.BC.geno.csv"
  ..$ marker.csv  : chr "../populus.BC.marker.csv"
  ..$ marker.obj  : NULL
  ..$ marker.table:'data.frame':        275 obs. of  4 variables:
  .. ..$ Marker : Factor w/ 275 levels "A/15-620D   ",..: 227 251 122 123 124 244 186 243 272 61 ...
  .. ..$ Dist   : num [1:275] 0 25.1 37 38.2 39.3 ...
  .. ..$ grp_idx: int [1:275] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ Group  : Factor w/ 22 levels "D1","D10","D11",..: 1 1 1 1 1 1 1 1 1 1 ...
  ..$ genos.matrix: int [1:78, 1:275] -1 -1 0 1 0 1 0 1 0 0 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:78] "1" "10" "11" "12" ...
  .. .. ..$ : chr [1:275] "marker1" "marker2" "marker3" "marker4" ...
 $ obj.phe  :List of 12
  ..$ pheno.csv    : chr "../populus.BC.pheno.csv"
  ..$ time.csv     : NULL
  ..$ log          : logi FALSE
  ..$ sample.obs   : int 78
  ..$ sample.times : int [1:11] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ pheY         : num [1:78, 1:11] 1.3 2.1 1 1.7 1.2 1.1 0.9 1 1.8 2.4 ...
  ..$ pheX         : NULL
  ..$ pheT         : int [1:78, 1:11] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ est.covar    :List of 2
  ..$ est.curve    :List of 7
  ..$ summary.curve:List of 5
  .. ..$ type   : chr "Logistic"
  .. ..$ par    : num [1:4] -6.401 32.937 5.131 0.455
  .. ..$ summary:'data.frame':  10 obs. of  9 variables:
  .. .. ..$ type: Factor w/ 10 levels "ABRK","Bi-Exponential",..: 8 3 1 9 4 2 10 5 6 7
  .. .. ..$ parm: num [1:10] 3 6 4 3 2 4 2 3 4 5
  .. .. ..$ AIC : num [1:10] 361 367 363 366 363 ...
  .. .. ..$ AICc: num [1:10] 5.66 5.75 5.69 5.72 5.68 ...
  .. .. ..$ BIC : num [1:10] 368 381 372 373 368 ...
  .. .. ..$ SSE : num [1:10] 7371 7371 7370 7879 7766 ...
  .. .. ..$ MSE : num [1:10] 8.59 8.59 8.59 9.18 9.05 ...
  .. .. ..$ RMSE: num [1:10] 2.93 2.93 2.93 3.03 3.01 ...
  .. .. ..$ R2  : num [1:10] 0.00491 0.00491 0.00486 0.07421 0.05876 ...
  ..$ summary.covar:List of 4
  .. ..$ type   : chr "SAD1"
  .. ..$ par    : num [1:2] 1.055 -0.819
  .. ..$ summary:'data.frame':  12 obs. of  4 variables:
  .. .. ..$ type: Factor w/ 12 levels "AR1","ARH1","ARMA(1,1)",..: 1 8 3 2 4 5 12 9 6 7 ...
  .. .. ..$ L   : num [1:12] -1249 -1132 -1227 -1190 -1850 ...
  .. .. ..$ AIC : num [1:12] 2502 2268 2460 2404 3703 ...
  .. .. ..$ BIC : num [1:12] 2506 2273 2468 2432 3708 ...
 - attr(*, "class")= chr "FM2.dat"
}


The phenotype file, measured time file,  genotype file and marker definition file must be a CSV file. The following sections illustrate the format of each data file.

1) \bold{The phenotype file.}\cr
The first column is individual ID and the rest columns are sample data 
for every measurement. It looks like the following file.

\preformatted{
ID.	1st	2nd	3rd	4th	5th	6th	7th
1	2.9033	4.118	6.1495	7.8161	9.8379	12.963	14.918
2	4.3306	5.3783	7.0647	9.3624	11.439	NA	15.701
3	2.3997	4.052	5.5431	7.6933	9.8471	NA	12.849
4	3.3044	4.154	5.8924	7.7133	9.2144	10.945	NA
...
}

Please note missing data is coded as space or NA in all four data files. 

2) \bold{The measurement time file.}
The first column is individual ID and the rest columns are sample data 
for each measurement. It looks like the following file.

\preformatted{
ID.	1st	2nd	3rd	4th	5th	6th	7th
1	1	2	3	4	5	6	7
2	1	2	3	4	5	NA	7
3	1	2	3	4	5	NA	7
4	1	2	3	4	5	6	NA
...
}


3) \bold{The genotype file.}
The first column is individual ID and the rest columns are 
genotype data for each marker. An example is shown in the following table.
Three genotypes (aa=0, Aa=1, AA=2) and missing data(coded as NA or -1)
are valid maker values.

\preformatted{
ID.	marker1	marker2	marker3	marker4	marker5	marker6
1	1	1	0	1	NA	0
2	1	1	1	1	0	0
3	1	1	1	0	1	1
...
}

4) \bold{The marker list file.}
The first column is marker's ID, the rest columns are the marker's name,
distance, group index, and group name for every marker.In the marker file,
the distance field is a distance (in cM) in one chromosome or linkage group. 
The header row should be included.

\preformatted{
id	Marker	Dist.	Grp_idx	Group
1	marker1	0	1	G1
2	marker2	20	1	G1
3	marker3	40	1	G1
...
}

}

\examples{

# Load the pre-installed data for this example
file.pheno.csv <- system.file("extdata","populus.BC.pheno.csv", package="Funmap2")
file.geno.csv <- system.file("extdata","populus.BC.geno.csv", package="Funmap2")
file.marker.csv <- system.file("extdata","populus.BC.marker.csv", package="Funmap2")

# Call the pipeline without permutation.
# Can't select QTL using pvalue due to the missing of permutation result.
dat <- FM2.load.data( file.pheno.csv, NULL, file.geno.csv, file.marker.csv, "BC",
    curve.type="auto", 
    covar.type="auto",
    pdf.file="FM2.test.load.pdf");

str( dat );
print( dat );
plot( dat, pdf.file="test.FM2.load.data.pdf")

}

\keyword{Data}

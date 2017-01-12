\name{FM2.pipe}
\alias{FM2.pipe}

\title{QTL mapping pipeline}

\description{ 
 Perform standard pipeline for the experiment data based on Functional Mapping framework.
} 

\usage{ 
FM2.pipe( pheno.csv, time.csv, geno.csv, marker.csv, cross.type, 
    curve.type = NULL, 
    covar.type = NULL, 
    model = "MLE", 
    grp.idx = NULL, 
    pdf.prefix = NULL, 
    threshold = 0.05, 
    threshold.type = "pvalue", 
    options = list() )
}

\arguments{
\item{pheno.csv}{a CSV file of phenotypic traits . The format is described in the details section of \code{\link{FM2.load.data}}.}
\item{time.csv}{a CSV file of measure time. The format is described in the details section of \code{\link{FM2.load.data}}.}
\item{geno.csv}{a CSV file of genotype marker. The format is described in the details section of \code{\link{FM2.load.data}}.}
\item{marker.csv}{a CSV file of marker definition. The format is described in the details section of \code{\link{FM2.load.data}}}
\item{cross.type}{string value indicating the cross type, including "F2", "BC" and "RIL".}
\item{curve.type}{string value indicating the curve type, including "auto", "logistic", "Exponential", "Power", "Legendre2", .., the full list is described in the details section of \code{\link{FM2.get.curve}}.}
\item{covar.type}{string value indicating the type of covariance structure, including "auto", "AR1", "SAD1", "ARMA", "CS", .., the full list is described in the details section of \code{\link{FM2.get.covariance}}.}
\item{model}{string value indicating the computation algorithm, currently only "MLE" is option.}
\item{grp.idx}{a numeric vector indicating which chromosomes or groups get involved in the QTL scaning.}
\item{pdf.prefix}{string value indicating the pefix name of pdf file exported by the pipeline}
\item{threshold}{a numeric value indicating the criteria of signicant QTLs.}
\item{threshold.type}{string value indicating the selection method of significant QTL, three optional values. 'pvalue', 'LR' and 'count' }
\item{options}{optional values for the pipeline, see the following details.} 
}

\value{
The data object and the result object of hypothesis tests are returned in a list object, 
\item{dat}{data object with the S3 class label of "FM2.dat", the structure is decribed in the details section of \code{\link{FM2.load.data}} }
\item{ret}{result object with the S3 class label of "FM2.qtl.mle", the structure is decribed in the details section of \code{\link{FM2.qtlscan}} }
}

\details{ 

The \code{options} paramaters can slightly adjust the results and greatly speed up the computational process. the below explains all items in the \code{options} list:

1) \strong{debug},  default=FALSE, indicating whether this function outputs the debug information.\cr
2) \strong{n.cores}, default=1, the cluster count or multiple cores for parallel permutation, used in \code{\link{FM2.qtlscan}} and \code{\link{FM2.permutation}} .\cr
3) \strong{scan.step}, default=2, an interval distance used to scan flanking marker, default is 1cm, used in \code{\link{FM2.qtlscan}}.\cr
4) \strong{peak.count}, default=5, a number determines how many significant QTLs will be selected, used in \code{\link{FM2.select.qtl}}.\cr
5) \strong{permu.loop}, default=100, the count of permutation loop, used in \code{\link{FM2.permutation}}.\cr
6) \strong{permu.filter.ratio},  default=1, indicating whether fast estimation algorithm on the basis of QTL filter is used or not in \code{\link{FM2.permutation}}. No any optimation for permutation under the default condition(=1). If 0.01 is specified, permutation is performed on top 0.01 QTLs which are highly associated with phenotypic traits. \cr

\strong{FM2.pipe} is a main pipeline of the Funmap2 package, it encapsulates a consecutive procedures, including:

1) Loading the phenotype, genotype and marker file, which is performed in the function \code{\link{FM2.load.data}}.\cr

2) Summarize the data object returned by the function \code{\link{FM2.load.data}} \cr

3) Performing the hypothesis tests on all chromosomes or specified chromosomes, which is performed in the function \code{\link{FM2.qtlscan}}.\cr

4) Summarize the result object returned by the function \code{\link{FM2.qtlscan}} \cr

5) Taking a long time to execute permutation parallelly or not, which is performed in the function \code{\link{FM2.permutation}}.\cr

6) Selecting the significant QTLs according to the selection method and threshold, which is performed in the function \code{\link{FM2.select.qtl}}.\cr

7) Outputing a PDF report which includes the summary of the data and QTL scaning results, which is performed in the function \code{\link{FM2.report}}\cr

}

\examples{

# Load the pre-installed data for the example
file.pheno.csv <- system.file("extdata","populus.pheno.csv", package="Funmap2")
file.geno.csv <- system.file("extdata","populus.geno.csv", package="Funmap2")
file.marker.csv <- system.file("extdata","populus.marker.csv", package="Funmap2")

# Call the pipeline without permutation.
# Can't select QTL using pvalue due to the missing of permutation result.
r <- FM2.pipe( file.pheno.csv, NULL, file.geno.csv, file.marker.csv, "BC",
    curve.type="logistic", 
    covar.type="auto", 
    grp.idx = c(1:5),
    threshold = 3,
    threshold.type = "count", 
    options=list(permu.loop=0) );

# Show the summary information of data object
show(r$dat);
# Show the summary information of result object
show(r$ret);

# Change the QTL criteria
r$ret <- FM2.select.qtl(r$ret, threshold = 40, threshold.type = "LR" );
show(r$ret);

# Make a report for the data analysis.
FM2.report("test.FM2.pipe.pdf", r$dat, r$ret );

}

\keyword{QTL Scaning}
\keyword{Pipeline}

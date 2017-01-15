\name{FM2.qtlscan}
\alias{FM2.qtlscan}

\title{ QTL scanning}

\description{ 
Perform QTL scaning for all QTLs to detect the significant ones based on the hypothesis test.
} 

\usage{ FM2.qtlscan( dat, model="MLE", grp.idx=NULL, options=list() ) }

\arguments{
\item{dat}{a data object returned by \code{\link{FM2.simulate}} or \code{\link{FM2.load.data}}}
\item{model}{string value indicateing which method will be used to test hypothesis, one optional value currently.}
\item{grp.idx}{a numeric vector indicating which chromosomes or groups get involved in the QTL scaning.}
\item{options}{optional list for QTL scanning, including:\cr
\strong{debug}: default=FALSE, indicating whether this function outputs the debug information.\cr
\strong{n.cores}: default=1, a number of cpu cores for parallel computation.\cr
\strong{scan.step}: default=1, an interval distance used to scan flanking marker, default is 1cm.\cr
\strong{peak.count}: default=5, a number indicating how many top(or significant) QTLs will be selected.\cr
}
}

\value{
This function returns a S3 object with the class label of \code{FM2.qtl.mle}, including the following items:

\item{param}{ list recording the parameters for this function call.}
\item{obj.phe}{the phenotype data, copied from the data object (\code{dat}).}
\item{obj.gen}{the genotype data, copied from the data object (\code{dat}).}
\item{obj.curve}{the curve object, copied from the data object (\code{dat}).}
\item{obj.covar}{the covariance structure, copied from the data object (\code{dat}).}
\item{obj.cross}{the cross object, copied from the data object (\code{dat}).}
\item{cross.type}{the cross type, copied from the data object (\code{dat}).}
\item{covar.type}{the type of covariance structure, copied from the data object (\code{dat}).}
\item{curve.type}{the curve type, copied from the data object (\code{dat}).}
\item{full.res}{a matrix for all QTLs with the postion, likelihood ratio, curve parameters of different genes and covariance parameters.}
\item{threshold.type}{the selection method of significant QTLs, used in \code{\link{FM2.select.qtl}} }
\item{threshold}{the criteria of significant QTLs, used in \code{\link{FM2.select.qtl}} }
\item{qtl.peaks}{a numeric vector indicating the row index of significant QTLs in above matrix. This item is available after the calling of \code{\link{FM2.select.qtl}}}
\item{obj.permu}{the permutation result, obtained from the calling of \code{\link{FM2.permutation}}}

}

\details{

This function returns a result object which can be inspected by the following method: \cr

1) \code{str} command. \cr

2) \code{print} command. \cr

3) \code{\link[=summary.FM2.qtl.mle]{summary}} command. \cr

4) \code{\link[=plot.FM2.qtl.mle]{plot}} command. \cr

\bold{Hypothesis Test:} \strong{For different genotypes, all parameters are identical.} \cr

The hypothesis testing scans every marker by the specified step (1cm). It maybe take a long time, so the Funmap2 package displays its progress after each chromosome (linkage group) has been calculated.\cr

After the QTL scanning, the package identifies the 5 top QTLs. At most one significant QTL is selected within each choromosome (group). The top QTLs are strongly displayed at the head of report. \cr

}

\examples{

# Simulate the data and do QTL scaning
dat <- FM2.simulate();
ret <- FM2.qtlscan(dat);
show(ret);

# Load the example data and do QTL scaning
file.pheno.csv <- system.file("extdata","populus.BC.pheno.csv", package="Funmap2")
file.geno.csv <- system.file("extdata","populus.BC.geno.csv", package="Funmap2")
file.marker.csv <- system.file("extdata","populus.BC.marker.csv", package="Funmap2")
dat <- FM2.load.data( file.pheno.csv, NULL, file.geno.csv, file.marker.csv, "BC",
    curve.type="logistic", 
    covar.type="AR1")

# Call the QTL scaning process
ret <- FM2.qtlscan(dat, grp.idx = c(1:2) );
show(ret);
plot(ret, pdf.file="test.FM2.qtlscan.pdf");


}

\keyword{QTL Scaning}

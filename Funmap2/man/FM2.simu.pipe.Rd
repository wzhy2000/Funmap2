\name{FM2.simu.pipe}
\alias{FM2.simu.pipe}

\title{Pipeline for simulation test}

\description{ 
Demostrate the simulation test using the pipeline function.
} 

\usage{ FM2.simu.pipe( cross.type = "BC", curve.type="Logistic", covar.type="AR1", 
    simu.mrkdist = rep(20,10), 
    simu.qtlpos = 95, 
    simu.obs = 800, 
    simu.times = 8, 
    par.X = NULL, 
    par0 = NULL, 
    par1 = NULL, 
    par2 = NULL, 
    par.covar = NULL, 
    phe.missing = 0.01, 
    marker.missing = 0.01,  
    threshold = 0.05, 
    threshold.type = "pvalue", 
    model = "MLE", 
    pdf.prefix = NULL, 
    options = list() ) }

\arguments{

\item{cross.type}{string value indicating the cross type, including "F2", "BC" and "RIL".}
\item{curve.type}{string value indicating the curve type, including "logistic", "Exponential", "Power", "Legendre2", .., the full list is described in the details section of \code{\link{FM2.get.curve}}.}
\item{covar.type}{string value indicating the covariance type, including "AR1", "SAD1", "ARMA", "CS", .., the full list is described in the details section of \code{\link{FM2.get.covariance}}.}
\item{simu.mrkdist}{numeric vector indicating the distance between the genomic marker.  }
\item{simu.qtlpos}{numeric value indicating the significant QTL position.}
\item{simu.obs}{numeric value indicating the sample size.}
\item{simu.times}{numeric value indicating the measured times. }
\item{par.X}{numeric vector indicating covariate parameters. }
\item{par0}{numeric vector indicating curve parameters for gene QQ, default value is retrived from the curve object. }
\item{par1}{numeric vector indicating curve parameters for gene Qq, default value is retrived from the curve object. }
\item{par2}{numeric vector indicating curve parameters for gene qq, default value is retrived from the curve object. }
\item{par.covar}{numeric vector indicating covariance parameters, default value is retrived from the covariance structure. }
\item{phe.missing}{numeric value indicating the missing rate of phenotypic traits. }
\item{marker.missing}{numeric value indicating the missing rate of genomic markers.}
\item{threshold}{numeric value indicating the criteria of signicant QTLs.}
\item{threshold.type}{string value indicating the selection method of significant QTL, three optional values. 'pvalue', 'LR' and 'count' }
\item{model}{string value indicating the computation algorithm, currently only "MLE" is option.}
\item{pdf.prefix}{string value indicating the pefix name of pdf file exported by the pipeline}
\item{options}{optional list for the pipeline, see the details in \code{\link{FM2.pipe}}.} 
}

\value{
A list including the data object(\strong{dat}) and the result object of QTL scaning (\strong{ret}).
\item{dat}{ data object describled in \code{\link{FM2.simulate}} }
\item{ret}{ result object of QTL scaning with permutation cutoff table describled in \code{\link{FM2.qtlscan}} }
}

\details{

The \code{options} paramater is described in \code{\link{FM2.pipe}}. \cr  

\strong{FM2.simu.pipe} demostrates how to use the Funmap2 to do a simulation test, which includes the following steps:

1) Simulate a raw data object on the basis of the parameters.\cr

2) Perform QTL scaning on all QTLs based on the hypothesis test.\cr

3) Execute permutation to get a cutoff for significant QTLs.\cr

4) Summarize all objects and plot all figures.\cr

5) Export a PDF report including all summary information and figures.\cr

}

\examples{

# Call the pipeline for the simulation test. This test doesn't 
# call permutation to determine the p-value for QTL peaks.
r <- FM2.simu.pipe("RIL", "Logistic", "SAD1", simu.obs=1000, simu.times = 7,  
    threshold = 1,  threshold.type = "count", options = list(permu.loop=0) );

# Summarize the data object
show(r$dat);

# Summarize the result object
show(r$ret);
}

\keyword{QTL Scaning}
\keyword{Simulation}
\keyword{Pipeline}

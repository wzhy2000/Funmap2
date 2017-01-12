\name{FM2.simulate}
\alias{FM2.simulate}

\title{Data simulation}

\description{ 
  Create a simulation data object for pipeline demonstration.
} 

\usage{ FM2.simulate( cross.type = "BC", 
    curve.type = "Logistic", 
    covar.type = "AR1", 
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
    pdf.file = NULL ) }

\arguments{

\item{cross.type}{string value indicating the cross type, including "F2", "BC" and "RIL".}
\item{curve.type}{string value indicating the curve type, including "logistic", "Exponential", "Power", "Legendre2", .., the curve list is described in the details section of \code{\link{FM2.get.curve}}.}
\item{covar.type}{string value indicating the covariance type, including "AR1", "SAD1", "ARMA", "CS", .., the covariance list is described in the details section of \code{\link{FM2.get.covariance}}.}
\item{simu.mrkdist}{numeric vector indicating the distance between the genomic marker.  }
\item{simu.qtlpos}{numeric value indicating the significant QTL position.}
\item{simu.obs}{numeric value indicating the sample size.}
\item{simu.times}{numeric value indicating the measured times. }
\item{par.X}{numeric vector indicating covariate parameters. }
\item{par0}{numeric vector indicating curve parameters for gene QQ, default value is retrived from the curve object.}
\item{par1}{numeric vector indicating curve parameters for gene Qq, default value is retrived from the curve object.}
\item{par2}{numeric vector indicating curve parameters for gene qq, default value is retrived from the curve object.}
\item{par.covar}{numeric vector indicating covariance parameters, default value is retrived from the covariance object.}
\item{phe.missing}{numeric value indicating the missing rate of phenotypic traits. }
\item{marker.missing}{numeric value indicating the missing rate of genomic markers.}
\item{pdf.file}{string variable suggesting a PDF file name to illustrate the performance of curve fitting.}
}

\value{
This function returns a S3 object with the class label of \code{FM2.dat}, which structure is same as the experiment data object obtained from the function \code{\link{FM2.load.data}}. \cr
}

\details{
The structure of simulation data is identical to experiment data object. The different points are listed below:

1) The items of pheno_file, geno_file and marker_file are 
made up by the Funmap2 package and will be used to assign 
the output filename as a filename prefix. \cr

2) In the genotype is coded by 1=Qq 2=QQ for backcross 
, 0=qq, 1=Qq, 2=QQ for F2 intercross and 0=qq, 2=QQ for RILs intercross. \cr

}

\examples{
dat <- FM2.simulate("RIL", "Logistic", "SAD1", simu.obs=1000, simu.times = 7 );
#summarize the data information.
summary( dat );
plot(dat, pdf.file="test.FM2.simulate.pdf");
str(dat);
}

\keyword{Data}
\keyword{Simulation}


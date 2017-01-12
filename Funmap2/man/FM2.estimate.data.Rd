\name{FM2.estimate.data}
\alias{FM2.estimate.data}
\title{Phenotype estimation}
\description{ 
Identifying or estimating the parameters of curve and covariance structure.
} 

\usage{ FM2.estimate.data(dat, curve.type=NULL, covar.type=NULL, pdf.file=NULL ) }

\arguments{
\item{dat}{ a data object returned by \code{\link{FM2.simulate}} or \code{\link{FM2.load.data}} }
\item{curve.type}{string value indicating the specific curve type. optional values are listed in \code{\link{FM2.get.curve}}. The curve fitting is performed using the least-squares if 'auto' or NULL is assigned.}
\item{covar.type}{string value indicating the specific type of covariance structure, optional values are listed in \code{\link{FM2.get.covariance}}. the MLE process identifies covariance matrix if 'auto' or NULL is assigned.  }
\item{pdf.file}{string value suggesting a PDF file name to illustrate the performance of curve fitting.}
}

\value{
A data object with the results of curve fitting and covariance identifying. 
This function updates or adds the estimation of curve fitting and covariance structure.
\code{\link{FM2.load.data}} illustrates the structure of data object. 

}

\examples{
# data simulation using the default parameters
dat <- FM2.simulate();
dat;

# estimate the parameter of curve object and covariance structure
dat <- FM2.estimate.data(dat);
dat;

}

\keyword{Data}

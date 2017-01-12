\name{FM2.report}
\alias{FM2.report}

\title{PDF report for data and result object}

\description{ 
Output a PDF report including the summary information and figures for the data and result object.
} 

\usage{ FM2.report( file.report.pdf, dat, res=NULL, options=list( debug=F ) ) }

\arguments{
\item{ file.report.pdf}{PDF file name.}
\item{dat}{a data object returned by \code{\link{FM2.simulate}} or \code{\link{FM2.load.data}} }
\item{res}{a result object returned by \code{\link{FM2.qtlscan}} }
\item{options}{option list including whether debug information is outputted.}
}

\value{
No return value.
}

\details{

This function don't use the HaruPDF package anymore!!!
It outputs the summary information and figures into a PDF file for the data object and the result object.

The following link is an example for \code{FM2.report}.

\url{http://statgen.psu.edu/software/funmap/report_demo.pdf}

}

\examples{
dat <- FM2.simulate();
res <- FM2.qtlscan(dat);
FM2.report("test.FM.report.pdf", dat, ret);
}

\keyword{Report}

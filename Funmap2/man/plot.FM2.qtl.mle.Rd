\name{plot.FM2.qtl.mle}
\alias{plot.FM2.qtl.mle}

\title{ Plotting figures of QTL scanning }

\description{ 
Plot the figures based on the results of hypothesis test.
} 

\usage{\method{plot}{FM2.qtl.mle}(x, plot.type=NULL, pdf.file=NULL, ... ) }

\arguments{
\item{x}{a result object of hypothesis tests returned by \code{\link{FM2.qtlscan}}.}
\item{plot.type}{a number indicating which figure is plotted.}
\item{pdf.file}{a pdf file name for the figure output, if no pdf file is specified, the plot command can output this figure in the R console. }
\item{...}{additional arguments affecting the plot produced.}
}

\details{

The result object of QTL scaning can output three kinds of figure according to the parameter 'plot.type', including:

1) The LR profile for all chromosomes.\cr
2) The LR profile for QTL postion. \cr
3) The curve for QTL postion. \cr

The examples can be viewed in the following url.

\url{http://statgen.psu.edu/software/funmap/plot.t10-1.jpg}.\cr

\url{http://statgen.psu.edu/software/funmap/plot.t10-2.jpg}.\cr

\url{http://statgen.psu.edu/software/funmap/plot.t10-3.jpg}.\cr

}

\examples{
#check the codes in FM2.qtlscan()
}

\keyword{Plot}


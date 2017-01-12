\name{plot.FM2.dat}
\alias{plot.FM2.dat}
\title{Plotting figures of data object}

\description{ 
  Draw figures for a data object.
} 

\usage{ \method{plot}{FM2.dat}( x, plot.type=NULL, pdf.file=NULL, ... )  }

\arguments{
\item{x}{a data object returned by \code{\link{FM2.simulate}} or\code{\link{FM2.load.data}}   .}
\item{plot.type}{number, the plot type, 1 is for tiled curves and 2 is for overlapping curves.}
\item{pdf.file}{a pdf file name for the figure output, if no pdf file is specified, the plot command can output this figure in the R console. }
\item{...}{additional arguments affecting the plot produced.}
}

\details{


Two figures can be outputted to R console.\cr

1) tiled curves for every individuals. \cr

2) overlapping curves for every individuals \cr

An example of this command is available in the following URL. \cr

\url{http://statgen.psu.edu/software/funmap/plot.data1.jpg}. \cr

\url{http://statgen.psu.edu/software/funmap/plot.data2.jpg}. \cr

}

\examples{
#check the codes in FM2.simulate() or FM2.load.data()
}

\keyword{Data}
\keyword{Plot}

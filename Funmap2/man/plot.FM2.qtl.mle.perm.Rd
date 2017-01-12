\name{plot.FM2.qtl.mle.perm}
\alias{plot.FM2.qtl.mle.perm}

\title{Plotting figure of permutation result}

\description{ 
Draw a cutoff profile on the basis of the permutation result.
} 

\usage{ \method{plot}{FM2.qtl.mle.perm}(x, pdf.file=NULL, ... ) }
 
\arguments{                                                  
\item{x}{an object of permutation result returned by \code{\link{FM2.permutation}}.}
\item{pdf.file}{a pdf file name for the figure output, if no pdf file is specified, the plot command can output this figure in the R console. }
\item{...}{additional arguments affecting the plot produced.}
}

\details{

This summary exports a figure based on the cutoff table in the permutation result. \cr
An example of this command is available in the following URL.\cr

\url{https://raw.githubusercontent.com/wzhy2000/Funmap2/master/img/plot.perm.jpg}.

}

\examples{
#check the example in the FM2.permutation()
}

\keyword{Permutation}
\keyword{Plot}

\name{FM2.select.qtl}
\alias{FM2.select.qtl}

\title{Selecting significant QTLs}

\description{ 
 Select significant QTLs according to the threshold and method.
} 

\usage{ FM2.select.qtl( res,  threshold=0.05, threshold.type="pvalue" ) }

\arguments{
\item{res}{a result object returned by \code{\link{FM2.qtlscan}}}
\item{threshold}{numeric value indicating the criteria of signicant QTLs.}
\item{threshold.type}{string value indicating the selection method of significant QTL, three optional values. 'pvalue', 'LR' and 'count'. }
}

\value{
A results with updated significant QTLs is returned. 
}

\details{

Three methods can be used to show the significant QTLs: \cr

1) use \code{show} or \code{print} command to check the significant QTLs;  \cr

2) use the \code{\link[=plot.FM2.qtl.mle]{plot}} command to show the significant QTLs and genetic curve at significant QTLs; \cr

3) access the items of result object \cr

e.g.
\preformatted{
> res <- FM2.select.qtl(res, threshold=5, threshold.type="count")
> cat("The significant QTL list:\n");
> show( res$full.res[ res$qtl.peaks, 1:3] );

}

}

\examples{

dat <- FM2.simulate();
ret <- FM2.qtlscan(dat);
ret <- FM2.select.qtl( ret,  threshold=40, threshold.type="LR" )
plot(ret, pdf.file="test.FM2.select.qtl.pdf");

}

\keyword{QTL Scaning}

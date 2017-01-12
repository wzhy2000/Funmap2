\name{summary.FM2.qtl.mle}
\alias{summary.FM2.qtl.mle }

\title{Summary of the result object.}

\description{ 
Summarize information for the QTL scanning results based on the hypothesis test.
} 

\usage{ \method{summary}{FM2.qtl.mle}( object, ... ) }

\arguments{
\item{object}{a result object returned by \code{\link{FM2.qtlscan}} which stores the results of hypothesis tests.}
\item{...}{additional arguments affecting the summary produced.}
}

\details{

The following sections demostrate the context of summary report. including: \cr

1) Estimated parameters. \cr
2) The significant QTL postions. \cr

\preformatted{
Hypothesis test 10: 
    a2=a1 and b2=b1 and r2=r1
------------------------------------
          Model: Logistic Curve
          Cross: Backcross
       QTL pos.: 50.1  (Group:8)
         QTL LR: 66.516   
    QTL p-value: 0.000    
Grwoth para(Qq): a2=  30.615, b2=  10.776, r2=   0.538
Grwoth para(qq): a1=  23.707, b1=   9.449, r1=   0.615
            rho: 0.953   
         sigma2: 8.637   
------------------------------------
 No.  Grp      Pos.       LR        a1        b1        r1        a0        b0        r0
   1    8   50.100    66.516    30.615    10.776     0.538    23.707     9.449     0.615
   2   12  113.100    55.190    29.865     9.736     0.528    25.207     8.725     0.586
   3   13   12.000    50.963    29.518     9.723     0.526    24.926     8.906     0.602
   4   18   10.000    25.684    29.236     9.550     0.536    25.523     8.963     0.584
   5    1  151.300    24.162    25.998     8.520     0.575    28.801     9.672     0.536

}

}

\examples{
#check the codes in FM2.qtlscan()
}

\keyword{Summary}


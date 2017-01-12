\name{summary.FM2.qtl.mle.perm}
\alias{summary.FM2.qtl.mle.perm}

\title{Summary of permutation result}

\description{ 
Summarize the permutation result.
} 

\usage{ \method{summary}{FM2.qtl.mle.perm}( object, ...) }

\arguments{
\item{object}{an object of permutation result returned by \code{\link{FM2.permutation}}.}
\item{...}{additional arguments affecting the summary produced.}
}

\details{

The summary command gives a table of cutoff values 
which starts at 90% and stops at 1% or even smaller.
If the p-value of x% want to be estimated, the 
permutation count should be greater than 100/x.

The following gives an example of this summary command.

\preformatted{
Permutation result:   

------------------------------------
          Curve: Logistic Curve
          Cross: BC
           Loop: 100
------------------------------------

p-value    Cutoff
0.90000    0.76583
0.80000    1.45845
0.70000    2.21657
0.60000    3.09488
0.50000    4.16930
0.40000    5.52038
0.30000    7.33172
0.20000    9.97031
0.10000   14.86003
0.09000   15.63818
0.08000   16.50743
0.07000   17.51136
0.06000   18.68964
0.05000   20.09660
0.04000   21.77278
0.03000   23.96467
0.02000   26.98845
0.01000   32.80592

}

}
\examples{
#check the example in the FM2.permutation()
}

\keyword{Summary}
\keyword{Permutation}


\name{FM2.permutation}
\alias{FM2.permutation}

\title{Permutation}

\description{ 
  Execute permutation to get the cutoff value for significance levels p=0.05 and 0.01.
} 

\usage{ FM2.permutation( dat, res, grp.idx=NULL, options=list() ) }

\arguments{
\item{dat}{a data object returned by \code{\link{FM2.simulate}} or \code{\link{FM2.load.data}} }
\item{res}{a result object returned by \code{\link{FM2.qtlscan}}}
\item{grp.idx}{a numeric vector indicating which chromosomes or groups get involved in the permutation test.}
\item{options}{optional value for permutation control, including.\cr
\bold{debug}: default=FALSE, indicating whether this function outputs the debug information.\cr
\bold{n.cores} : default=1, the number of cpu cores for parallel computation.\cr
\bold{scan.step}: default=1, an interval distance used to scan flanking marker, default is 1cm.\cr
\bold{permu.loop}: default=100, the count of permutation loop.\cr
\bold{permu.filter.ratio}: default=1, indicating whether fast estimation algorithm on the basis of QTL filter is used or not. No any optimation for permutation with the default value(=1). If 0.01 is specified, permutation is performed on top 0.01 QTLs which are highly associated with phenotypic traits. \cr
}
}

\value{

This function returns a new result object with the update of permutation results. The result objects is described in \code{\link{FM2.qtlscan}}. \cr 

Here we explain the updated results only in the item of \code{obj.perm}, a S3 object mainly including a p-value matrix(\code{pv.table}). \cr

\item{cross.type}{the cross type}
\item{curve.type}{the curve type}
\item{covar.type}{the type of covariance structure}
\item{permu.loop}{the permutation count}
\item{param}{a list recording the paramters in this functrion calling, including \code{permu.loop}, \code{permu.filter.ratio}, \code{scan.step}, \code{n.cores} }
\item{full.res}{matrix recording all the permutation results} 
\item{pv.table}{matrix which has two columns, the first column is significance leve and nd the second column is cutoff value.} 
}


\details{
If permutation count is greater than 100, more precise cutoff will be obtained.
For example, 10000 times permutation can give the significance table which
looks like the following table.

\tabular{rrrrrrrrr}{
0.9\tab0.8\tab0.7\tab0.6\tab0.5\tab0.4\tab0.3\tab0.2\tab0.1\cr
0.09\tab0.08\tab0.07\tab0.06\tab0.05\tab0.04\tab0.03\tab0.02\tab0.01\cr
0.009\tab0.008\tab0.007\tab0.006\tab0.005\tab0.004\tab0.003\tab0.002\tab0.001\cr
0.0009\tab0.0008\tab0.0007\tab0.0006\tab0.0005\tab0.0004\tab0.0003\tab0.0002\tab0.0001\cr
}

If clusters or multiple CPU cores are available, the permutaion can use \strong{snowfall} package to do parallel computation.
In order to do that, the following is necessary.

1) \strong{snowfall} is installed.\cr

2) The cluster count should be specified in the \code{options} parameter\cr

e.g. \code{options=list(n.cores=10)}

}

\examples{
dat <- FM2.simulate();
ret <- FM2.qtlscan(dat);
ret <- FM2.permutation(dat, ret, options=list(n.cores=10, permu.filter.ratio = 0.02, scan.step=2));

# only print the permutation part in the result object
show(ret$obj.permu);

# draw cutoff curve based on permutation results.
plot(ret$obj.permu, pdf.file="test.FM2.permu.pdf");

}

\keyword{Permutation}

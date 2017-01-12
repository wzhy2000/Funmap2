\name{summary.FM2.dat}
\alias{summary.FM2.dat}

\title{Summary of the data object}

\description{ 
  Summarize information for the data object.
} 

\usage{ \method{summary}{FM2.dat}( object, ... ) }

\arguments{
\item{object}{a data object return by \code{\link{FM2.simulate}} or \code{\link{FM2.load.data}} }
\item{...}{additional arguments affecting the summary produced.}
}

\value{
No return values, only output the summary information on the R console.
}

\details{

The data object is described in \code{\link{FM2.load.data}}.

The following example demonstrates summary command for a data object.

\preformatted{
Data set for FunMap model:
------------------------------------
          Date: 2010-03-19 03:49:36
          Model: Logistic Curve
          Cross: F2
    Pheno. file: simu.pheno.LC.F2
     Geno. file: simu.geno.LC.F2
     Maker file: simu.marker.LC.F2
    Sample size: 100       
   Sample times: 7         
   Marker count: 6         
          LC  a: 19.83678  
              b: 8.96370   
              r: 0.47202   
            rho: 0.75430   
         sigma2: 0.58849   
------------------------------------

}

}


\examples{
#check the codes in FM2.simulate() or FM2.load.data()
}

\keyword{Data}
\keyword{Summary}

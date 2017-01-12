\name{FM2.get.covariance}
\alias{FM2.get.covariance}

\title{Retriving covariance structure}

\description{ 
Retrive a covariance structure which characterizes the correlation between the measured phenotype.
} 

\usage{ FM2.get.covariance(covar.type) }

\arguments{

\item{covar.type}{string value indicating the type of covariance structure, such as "AR1", "SAD1", full list is described in the "details" section }

}

\value{
This functions returns a S4 object of covariance structure. 
You can use \code{show} or \code{print} command to check the summary information.
}

\details{

13 covariance structures are implemented in current version, including:

 \tabular{rll}{
[1] \tab "AR1"        \tab First-order Autoregressive \cr
[2] \tab "SAD1"       \tab First-order Structured Antedependence \cr
[3] \tab "ARMA(1,1)"  \tab First-order Autoregressive Moving Average \cr
[4] \tab "ARH1"       \tab Heterogeneous Autoregressive \cr
[5] \tab "CS"         \tab Compound Symmetry \cr
[6] \tab "CSH"        \tab Heterogeneous Compound Symmetry \cr
[7] \tab "VS"         \tab Variance Components \cr
[8] \tab "SI"         \tab Scaled Identity \cr
[9] \tab "FA1"        \tab Factor Analytic - First-order \cr
[10] \tab "FAH1"       \tab Heterogeneous Factor Analytic - First-order \cr
[11] \tab "TOEP"       \tab Toeplitz \cr
[12] \tab "TOEPH"      \tab Heterogeneous Toplitz \cr
[13] \tab "HF"         \tab Huynh-Feldt \cr
}	

The following summarize the parameters of covariance structure 'AR1'. 

\preformatted{
> x<-FM2.get.covariance("AR1");
> show(x);
     Class : fg.covariance.AR1 
Covar.Type : AR1 
Parameters : rho sigma2 

}

 
}
\examples{
x<-FM2.get.covariance("SAD1");
x;
show(x);
}

\keyword{Covariance}

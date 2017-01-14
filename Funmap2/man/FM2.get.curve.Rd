\name{FM2.get.curve}
\alias{FM2.get.curve}

\title{Retriving a curve object.}

\description{ 
  Retrive a curve object.
} 

\usage{ FM2.get.curve(curve.type) }

\arguments{
\item{curve.type}{string value indicating the curve type, full list in the details.}
}

\value{
This function reurns a S4 object of curve. The structure is described in the details section.
You can use \code{show} or \code{print} command to check the summary information.
}

\details{

9 curves have been implemented in current version, including:

\bold{1) "Logistic" }
        \deqn{g(t) = \frac{a}{1+b*e^{-r*t}} }{g(t)=a/(1+b*exp(-rt))}


\bold{2) "Bi-Logistic"}
        \deqn{g(t) = \frac{a1}{1+b1*e^{-r1*t}} + \frac{a2}{1+b2*e^{-r2*t}} }{ g(t)=a1/(1+b1*exp(-r1*t)) + a2/(1+b2*exp(-r2*t)) }

\bold{3) "Pharmacology"}
        \deqn{g(t) = \frac{E_{max}*t}{Ec_{50}+t} + E_{0} }{ g(t)=Emax*t/(Ec50+t) + E0 }
        
\bold{4) "Exponential"}
        \deqn{g(t)= a*e^{-r*t} }{ g(t)=a1*exp(-r*t) }

\bold{5) "Bi-Exponential"}
        \deqn{g(t) = a_{1}*e^{-r_{1}*t} + a_{2}*e^{-r_{2}*t} }{ g(t)=a1*exp(-r1*t) + a2*exp(-r2*t) }

\bold{6) "Power"}
        \deqn{g(t)= a*t^b }{ g(t)=a*t^b }

\bold{7) "Legendre2"}, Legendre Polynomial(2nd-order)
        \deqn{ g(t) = u_{0}+u_{1}*t + u_{2}*(3*t^2-1)/2 }{ y = u0 + u1*t + u2*1/2*(3*t^2-1) }

\bold{8) "Legendre3"}, Legendre Polynomial(3rd-order)
        \deqn{ g(t) = u_{0}+u_{1}*t + u_{2}*(2*t^2-1)/2 + u_{3}*(5*t^3-3t)/2 }{ y = u0 + u1*t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) }

\bold{9) "Legendre4"}, Legendre Polynomial(4th-order)
        \deqn{ g(t) = u_{0}+u_{1}*t + u_{2}*(2*t^2+1)/2 + u_{3}*(5*t^3-3t)/2 + ... }{ y = u0 + u1*t + u2*1/2*(3*t^2-1) + u3*1/2*(5t^3-3t) + u4*1/8*(35*ti^4-30*ti^2+3) }

The following introduces the summary infomation of a curve object. 

\preformatted{
> x<-FM2.get.curve("Logistic");
> show(x);
     Class : fg.curve.log 
Curve Type : Logistic 
Parameters : a b r 
   Formula : y = a/(1+b*exp(-r*t)) 

}
 
}

\examples{

curve <- FM2.get.curve("Logistic");
show(curve);

}

\keyword{Curve}
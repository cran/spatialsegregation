\name{shannonF}
\alias{shannonF}
\alias{shannon.index}
\docType{package}

\title{Spatial Shannon Index}

\description{
Compute the spatial and aspatial Shannon index for a given multitype point pattern.
}


\details{
The form of Shannon index is  \var{H = 1 - E(o)/E(N)}, where \var{E(N)} is the global entropy and \var{E(o)} is the local entropy 
calculated as \var{E(o)= - sum pi_tau log(pi_tau)}, where the sum is over the different types present in the pattern, and \var{pi_tau} is 
the expected frequency of type \var{tau} points in a neighbourhood of a typical point of the pattern. 

The function \code{shannonF} is the calculation function. Uses function \code{\link{segregationFun}}. 

The function \code{shannon.index} is a shortcut to get the non-spatial Shannon index. 
}

\usage{
  
  shannonF(X, r=NULL, v2=FALSE, ...)
  shannon.index(X, spatial=FALSE, ...)
}

\arguments{
\item{X}{Multitype point pattern of class \code{ppp} (see package 'spatstat')}
\item{r}{Vector of sizes for neighbourhoods, e.g. \code{geometric} graph with different ranges.}
\item{spatial}{If FALSE, return the classical aspatial index value.}
\item{v2}{If TRUE, use the real number of types in neighbourhoods as the log-base instead of total population type count.}
\item{...}{Further parametes for the function \code{\link{segregationFun}}.}
}

\value{
 Returns an \code{fv}-object, see \code{spatstat} for more information.
}
\author{
Tuomas Rajala
University of Jyvaskyla, Finland
tuomas.rajala@iki.fi
}
\references{
Rajala, Illian: A family of spatial biodiversity measures based on graphs, Env. Ecol. Stat. 2012

Reardon, O'sullivan: Measures of spatial segregation. Sociological methodology, 34:121-162, 2004.
}

\keyword{ package }

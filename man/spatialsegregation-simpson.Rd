\name{simpsonF}
\alias{simpsonF}
\alias{simpson.index}
\docType{package}

\title{Spatial Simpson Index}
\description{
Compute the spatial and aspatial Simpson index for a given multitype point pattern.
}


\details{
The form of Simpson index is  \var{S = 1 - sum pi_tau}, where the sum is over the types of the pattern, and \var{pi_tau} is like
in Shimatani\& Kubota 2004.

The function \code{simpsonF} is the main calculation function. Uses function \code{\link{segregationFun}}.

The function \code{simpson.index} is a shortcut to get a single value for the pattern. Uses 4-nn graph by default. 

}

\usage{
simpsonF(X, r=NULL, ...)
simpson.index(X, spatial=FALSE, ...)
}

\arguments{
\item{X}{Multitype point pattern of class \code{ppp} (see package 'spatstat')}
\item{r}{Vector of sizes for neighbourhoods, e.g. \code{geometric} graph with different ranges.}
\item{spatial}{If FALSE, return the classical aspatial index value.}
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

Shimatani, Kubota: Quantitative assesment of multispecies spatial pattern with high species diversity. Ecological Research, 19, 2004.
}

\keyword{ package }

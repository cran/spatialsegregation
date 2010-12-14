\name{biomassF}
\alias{biomassF}

\title{Individual Species Area Relationship}

\description{
Compute the biomass around individuals points.
}


\details{
Computes the neighbourhood for each point and then sums up the biomass in that neighbourhood. 
}

\usage{
biomassF(X, r=NULL, target=NULL, v2=FALSE, ...)
}
\arguments{
\item{X}{Multitype point pattern of class \code{ppp} (see package 'spatstat')}
\item{r}{Vector of sizes for neighbourhoods, e.g. \code{geometric} graph with different ranges.}
\item{target}{Default NULL. Calculate only for target type. If NULL compute mean over all types.}
\item{v2}{Logical. Return the average biomass instead of just sum.}
\item{...}{Further parameters for the function \code{segregationFun}.}
}
\value{
 Returns an \code{fv}-object, see \code{spatstat} for more information.
}

\author{
Tuomas Rajala
University of Jyvaskyla, Finland
tuomas.a.rajala@jyu.fi
}
\references{
Rajala, Illian: Graph-based description of mingling and segregation in multitype spatial point patterns. To appear 2010 
}

\keyword{ package }

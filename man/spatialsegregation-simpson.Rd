\name{spatialsegregation-simpson}
\alias{simpsonF}
\alias{simpson_index}
\docType{package}

\title{Spatial Simpson Index}
\description{
Compute the spatial and aspatial Simpson index for a given multitype point pattern.
}


\details{
  \tabular{ll}{
    Date: \tab 2009-03-09\cr
    License: \tab GPL v2 or later\cr
  }


The form of Simpson index is  \var{S = 1 - sum pi_tau}, where the sum is over the types of the pattern, and \var{pi_tau} is like
in Shimatani\& Kubota 2004.

The function \code{simpsonF} is the main calculation function. Uses function \code{\link{segregationFun}} with parameter \code{ftype=3}. 

The function \code{simpson_index} is a shortcut to get a single value for the pattern. Uses 4-nn graph by default. 

}

\usage{
simpsonF(pp, parvec = 1:20, graph_type = "knn", ...)
simpson_index(pp, graph_type = "knn", graph_par = 4, spatial = TRUE, ...)
}

\arguments{
\item{pp}{Multitype point pattern (see package \code{spatstat})}
\item{parvec}{Default 1:20. Vector for the neighbourhood defining graph, e.g. \code{geometric} graph with different \code{r}.}
\item{graph_type}{Default \code{knn}. Type of the neighbourhood graph. Accepts: \code{knn}, \code{geometric}, \code{delauney}, \code{gabriel}.}
\item{graph_par}{The single parameter to calculate the value in.}
\item{spatial}{If FALSE, return the classical aspatial index value.}
\item{...}{Further tuning parametes for the function segregationFun.}
}

\author{
Tuomas Rajala
University of Jyvaskyla, Finland
tarajala@maths.jyu.fi
}
\references{
Rajala, Illian: Graph-based description of mingling and segregation in multitype spatial point patterns. To appear 2009.

Shimatani, Kubota: Quantitative assesment of multispecies spatial pattern with high species diversity. Ecological Research, 19, 2004.
}

\keyword{ package }

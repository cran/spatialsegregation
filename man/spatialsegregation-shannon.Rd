\name{spatialsegregation-shannon}
\alias{shannonF}
\alias{shannon_index}
\docType{package}

\title{Spatial Shannon Index}

\description{
Compute the spatial and aspatial Shannon index for a given multitype point pattern.
}


\details{
  \tabular{ll}{
    Date: \tab 2009-03-09\cr
    License: \tab GPL v2 or later\cr
  }

The form of Shannon index is  \code{$H = 1 - E(o)/E(N)$, where $E(N)$} is the global entropy and \code{$E(o)$} is the local entropy 
calculated as \code{$E(o)= - \sum \pi\_\tau log(\pi\_\tau)$}, where the sum is over the types of the pattern, and \code{$\pi\_\tau$} is 
the expected relative frequency of type \code{$\tau$} points in a neighbourhood of a typical point of the pattern. 

The function \code{shannonF} is the main calculation function. Uses function \code{\link{segregationFun}} with parameter \code{ftype=2}. 

The function \code{shannon_index} is a shortcut to get a single value for the pattern. Uses 4-nn graph by default. 
}

\usage{
  shannonF(pp, parvec = 1:20, graph_type = "knn", ...)
  shannon_index(pp, graph_type = "knn", graph_par = 4, spatial = TRUE, ...)
}

\arguments{
\item{pp}{Multitype point pattern (see package 'spatstat')}
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

Reardon, O'sullivan: Measures of spatial segregation. Sociological methodology, 34:121-162, 2004.
}

\keyword{ package }

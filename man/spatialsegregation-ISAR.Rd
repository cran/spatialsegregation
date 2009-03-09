\name{spatialsegregation-ISAR}
\alias{isarF}
\alias{isar_index}

\title{Individual Species Area Relationship}

\description{
Compute the Individual Species Area Relationship ( ISAR ) for a given multitype point pattern.
}


\details{

  \tabular{ll}{
    Date: \tab 2009-03-04\cr
    License: \tab GPL v2 or later\cr
  }

Extension of ISAR-function introduced in WGGH 2007. In effect calculates the expected amount of different types present
in the neighbourhood of a typical point.

The function \code{isarF} is the main calculation function. Uses function 'segregationFun' with parameter \code{ftype=4}. 

The function \code{isar_index} is a shortcut to get a single value for the pattern. Uses 4-nn graph by default.

}

\usage{
isarF(pp, parvec = 1:20, graph_type = "knn", type = NULL, ...)
isar_index(pp, graph_type = "knn", graph_par = 4, type = NULL, ...)
}
\arguments{
\item{pp}{Multitype point pattern (see package 'spatstat')}
\item{parvec}{Default 1:20. Vector for the neighbourhood defining graph, e.g. \code{geometric} graph with different \code{r}.}
\item{graph_type}{Default \code{knn}. Type of the neighbourhood graph. Accepts: \code{knn}, \code{geometric}, \code{delauney}, \code{gabriel}.}
\item{type}{Default NULL. Calculate only for type 'type'. If NULL compute mean over all types.}
\item{graph_par}{The single parameter to calculate the value in.}
\item{...}{Further tuning parametes for the function segregationFun.}
}

\author{
Tuomas Rajala
University of Jyvaskyla, Finland
tarajala@maths.jyu.fi
}
\references{
Rajala, Illian: Graph-based description of mingling and segregation in multitype spatial point patterns. To appear 2009.

Wiegand, Gunatilleke, Gunatilleke, Huth: How individual species structure diversity in tropical forests. PNAS, nov 16, 2007. 
}

\keyword{ package }

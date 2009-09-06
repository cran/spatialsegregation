\name{spatialsegregation-ISAR}
\alias{isarF}
\alias{isar.index}

\title{Individual Species Area Relationship}

\description{
Compute the Individual Species Area Relationship ( ISAR ), or Local Species Richness, for a given multitype point pattern.
}


\details{

  \tabular{ll}{
    Date: \tab 2009-09-06\cr
    License: \tab GPL v2 or later\cr
  }

Extension of ISAR-function introduced in WGGH07. In effect calculates the expected amount of different types present
in the neighbourhood of a point in the pattern.

The function \code{isarF} is the calculation function for different neighbourhoods. Uses function \code{\link{segregationFun}}.

The function \code{isar.index} is a shortcut to get a single value for the pattern. Uses 4-nn graph by default.

}

\usage{
isarF(X, r=NULL, target=NULL, v2=FALSE, ... )
isar.index(X, r=4, ntype="knn", ...)
}
\arguments{
\item{X}{Multitype point pattern of class \code{ppp} (see package 'spatstat')}
\item{r}{Vector of sizes for neighbourhoods, e.g. \code{geometric} graph with different ranges.}
\item{target}{Default NULL. Calculate only for target type. If NULL compute mean over all types.}
\item{v2}{Logical. Estimate species/neighbours -ratios instead of just species/1.}
\item{ntype}{Set the default in isar.index to \code{knn} neighbourhood type.}
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
Rajala, Illian: Graph-based description of mingling and segregation in multitype spatial point patterns. To appear 2009.

Wiegand, Gunatilleke, Gunatilleke, Huth: How individual species structure diversity in tropical forests. PNAS, nov 16, 2007. 
}

\keyword{ package }

\name{spatialsegregation-helpers}
\alias{clean_up_data}
\alias{freqs}
\alias{minusID}
\alias{shake}


\title{Functions for the aid of segregation measures}

\description{
Small functions included in package \code{spatialsegregation}, used for manipulation of forest datasets which
have dbh-values (pp with a component dbh).
}


\details{

  \tabular{ll}{
    Date: \tab 2009-04-03\cr
    License: \tab GPL v2 or later\cr
  }

Small functions to manipulate multitype point patterns.

\code{clean_up_data}: Returns a subsample fullfilling the given constrains.

\code{freqs}: Returns the abundance vector.

\code{minusID}: Returns a 0-1-vector indicating inclusion in a simple minus-correction.

\code{shake}: Shakes the pattern, i.e. adds a random displacement shift to each point.
}

\usage{
 clean_up_data(pp, dbh = 10, atleast = 10)
 freqs(pp)
 minusID(pp, minusR, dbh = 0, atleast = 0)
 shake(pp, a = 0.001)
}
\arguments{
\item{pp}{Multitype point pattern (see package 'spatstat')}
\item{atleast}{Include specii with abundance atleast \code{atleast}.}
\item{dbh}{Include only those points with dbh atleast \code{dbh}.}
\item{minusR}{Range from the border withing which to exluce points (used for correction of estimates).}
\item{a}{Size of displacement: x+Unif(-a,a), y+Unif(-a,a). }
}

\author{
Tuomas Rajala
University of Jyvaskyla, Finland
tarajala@maths.jyu.fi
}

\keyword{ spatial }
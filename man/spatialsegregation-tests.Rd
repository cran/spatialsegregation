\name{Dixon's 2-type contingency table tests}
\alias{dixon}
\alias{print.segtest}


\title{Dixon's 2-type contingency table tests}

\description{
Computes the tests of segregation using nearest neighbour contingency tables introduced 
by Philip Dixon in his paper "Testing spatial segregation using a nearest-neighbor contingency table",
Ecology, 75, p.1940-1948 (1994). The tests are an improvement on the Pielou's test of segregation.

The test is defined only for two-type spatial pattern.
}

\usage{
   dixon(X, prepR=0)
}
\arguments{
\item{X}{Bivariate i.e. 2-type point pattern (see package 'spatstat')}
\item{prepR}{Computes first the geometric neighbours with this distance, and then finds the nearest neighbours.}
}

\details{
See the paper by Dixon for details.
}

\author{
Tuomas Rajala
University of Jyvaskyla, Finland
tuomas.rajala@iki.fi
}

\references{
Philip Dixon: "Testing spatial segregation using a nearest-neighbor contingency table",
Ecology, 75, p.1940-1948 (1994).
}
\keyword{ spatial }

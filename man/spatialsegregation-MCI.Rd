\name{mciF}
\alias{mciF}
\title{Mean Composite Information}
\description{
Compute the Mean Composite Information for a given multitype point pattern. See Podani\&Czaran 1997.
}


\details{
The function \code{mciF} is the main calculation function. Uses function \code{\link{segregationFun}}. 
}

\usage{
mciF(X, r=NULL, target=NULL, ...)
}

\arguments{
\item{X}{Multitype point pattern of class \code{ppp} (see package 'spatstat')}
\item{r}{Vector of sizes for neighbourhoods, e.g. \code{geometric} graph with different ranges.}
\item{target}{If given, look at the surroundings of this type only.}
\item{...}{Further parameters for the function \code{segregationFun}.}
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
Podani, Czaran: Individual-centered analysis of mapped point patterns representing multi-species assemblages. J. Veg. Sci. 8: 259-270, 1997.
}

\keyword{ package }

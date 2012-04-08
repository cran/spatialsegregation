\name{exposurepps}
\alias{exposurepps}

\title{Example datasets for package spatialsegregation}
\description{Example datasets with 9 different scattering+exposure combinations.}
\usage{
   data(exposurepps)
}
\format{
A list with nine elements of class \code{ppp}.
}
\details{
A list of 9 point patterns with different degrees of intra-species clustering and 
inter-species mingling.

The patterns are synthetically produced using a combination of the functions of this
package as an energy function of a Gibbs model. 
}


\author{
Tuomas Rajala
University of Jyvaskyla, Finland
tuomas.rajala@iki.fi
}

\examples{
 data(exposurepps)
 par(mfrow=c(3,3), mar=c(2,2,2,2))
 for(i in 1:9)plot(exposurepps[[i]])
}
\keyword{ spatial }

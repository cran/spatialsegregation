\name{spatialsegregation-package}
\alias{spatialsegregation-package}
\alias{spatialsegregation}
\docType{package}

\title{Spatial Segregation Measures}

\description{
Collection of measures of spatial multitype exposure: segregation vs. mingling.
}


\details{
  
  \tabular{ll}{
    Date: \tab 2009-03-09\cr
    License: \tab GPL v2 or later\cr
  }

This is a collection of summaries for multitype spatial point patterns (check package 'spatstat' for more info about them). 
The package is developed for the article Rajala\&Illian 2009, and consists of summaries for detecting inter-type forces in the 
pattern.

See the helps of one of the function for further information how to use the package.  


\preformatted{
Functions
------------------
segregationFun   -   General calculation function, please use one of the following:
minglingF        -   Mingling index
shannonF         -   spatial Shannon index
simpsonF         -   spatial Simpson index
isarF            -   ISAR function
mingling\_index   -   shortcut for a single value
shannon\_index    -   -''-
simpson\_index    -   -''-
isar\_index       -   -''-
}

}

\author{
Tuomas Rajala
University of Jyvaskyla, Finland
tarajala@maths.jyu.fi
}

\references{
Graz: The behaviour of the species mingling index \code{$m_{sp}$} in relation to species dominance and dispersion. Eur. J. forest research. 123:87-92, 2004.

Lewandowski, Pommerening: Zur Beschreibung der Waldstruktur - Erwartete und beobachtete Arten-Durchmischung. Forstwiss Centralbl, 116:129-139, 1997.

Rajala, Illian: Graph-based description of mingling and segregation in multitype spatial point patterns. To appear 2009.

Reardon, O'sullivan: Measures of spatial segregation. Sociological methodology, 34:121-162, 2004.

Shimatani, Kubota: Quantitative assesment of multispecies spatial pattern with high species diversity. Ecological Research, 19, 2004.

Wiegand, Gunatilleke, Gunatilleke, Huth: How individual species structure diversity in tropical forests. PNAS, nov 16, 2007. 
}

\keyword{ package }

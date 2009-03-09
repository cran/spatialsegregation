\name{spatialsegregation-mingling}
\alias{minglingF}
\alias{mingling_index}

\title{Spatial Mingling Index}
\description{
Compute the Mingling index for a given multitype point pattern.
}


\details{
  \tabular{ll}{
    Date: \tab 2009-03-09\cr
    License: \tab GPL v2 or later\cr
  }


Extension of Mingling index introduced by Lewandowski\&Pommerening 1997. Measures the proportion of alien points in
the neighbourhood of a specific type typical point of the pattern.

If no specific type is given, the function means over all types. A typewise value is more useful.

The function \code{minglingF} is the main calculation function. Uses function \code{\link{segregationFun}} with parameter \code{ftype=1}. 

The function \code{mingling_index} is a shortcut to get a single value for the pattern. Uses 4-nn graph by default, which is the original 
Mingling index used by Lewandowski\&Pommerening 1997 and Graz 2004. 
}

\usage{
 minglingF(pp, parvec = 1:20, graph_type = "knn", type = NULL, ratio = FALSE, ...)
 mingling_index(pp, graph_type = "knn", graph_par = 4, type = NULL, ...)
}

\arguments{
\item{pp}{Multitype point pattern (see package 'spatstat')}
\item{parvec}{Default 1:20. Vector for the neighbourhood defining graph, e.g. \code{geometric} graph with different \code{r}.}
\item{graph_type}{Default \code{knn}. Type of the neighbourhood graph. Accepts: \code{knn}, \code{geometric}, \code{delauney}, \code{gabriel}.}
\item{type}{Default NULL. Calculate only for type 'type'. If NULL compute mean over all types.}
\item{ratio}{Default FALSE. If TRUE, scale the typewise values \code{$M\_t$} using formula \code{$(1-M\_\tau)/\lambda\_\tau$} which equals 1 for Poisson.}
\item{graph_par}{The single parameter to calculate the value in.}
\item{...}{Further tuning parametes for the function segregationFun.}
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
}

\keyword{ package }

\name{spatialsegregation-segregationFun}
\alias{segregationFun}
\alias{minusID_gf}
\alias{plot.segfcl} 
\alias{print.segfcl} 
\alias{segfcl}
\alias{spatialsegregation_function_explanations}
\docType{package}

\title{Spatial Segregation Function}

\description{
Compute the spatial exposure (segregation vs. mingling) features from a given multitype point pattern. Usage of 
shortcuts \code{minglingF}, \code{isarF}, \code{shannonF} and \code{simpsonF} highly recommended.
}

\usage{
     segregationFun(pp, fpar = NULL, graph_parvec = 1:20, graph_type = "knn",
    toroidal = FALSE, dbg = FALSE, funtype = 1, doDists = FALSE,
    prepR = 0, prepGraph = NULL, included = NULL, minusR = NULL,
    relative = FALSE)
}

\arguments{
\item{pp}{Multitype point pattern (see package 'spatstat')}
\item{fpar}{Default NULL. Parameter(s) for the measure. Mingling: c(i,j), where i= only for type i (0 for all), j=1 -> ratio version. ISAR: i, i=type (integer). Others: NULL.}
\item{graph_parvec}{Default 1:20 (\code{parvec} in shortcuts). Vector for the neighbourhood defining graph, e.g. "geometric" graph with different "r".}
\item{graph_type}{Default "knn". Type of the neighbourhood graph. Accepts: "knn", "geometric", "delauney", "gabriel".}
\item{toroidal}{Default FALSE. If TRUE, use a toroidal correction in distance calculation. Works only for rectangular windows and "geometric" or "knn" graph.}
\item{dbg}{Default FALSE. Print additional runtime texts.}
\item{funtype}{Default 1. 1=Mingling index, 2=Shannon, 3=Simpson, 4=ISAR.}
\item{doDists}{Default TRUE. Precalculate distances for speed. Be aware of memory requirements, O(n*n)!}
\item{prepR}{Default 0. If >0, shrink the search space for connections in the graph by searching only points within distance R (precalculates geometric graph).}
\item{prepGraph}{Precalculated graph for the point pattern. If not NULL, The graph\_par, graph\_parvec, toroidal and prepR are ignored and calculations are carried 
      using only prepGraph. Useful for huge datasets.}
\item{included}{boolean-vector of length |pp|. included[i]==TRUE => pp[i] included in calculations. See the next parameter for quick minus-sampling.}
\item{minusR}{If given included-vector is created with points with distance atleast minusR from the border.}
\item{relative}{Default FALSE. If TRUE, scale the parvec to a unit distance using intensity. Only for "geometric" graph. Useful for comparison of different scale pp's.}
}

\details{
General function for computing the spatial exposure (segregation/mingling) features. Used by minglingF, shannonF, simpsonF and isarF,
which should be preferred for better (and nicer) outcome.

Possible neighbourhood relations for the spatial version include geometric, k-nearest neighbours, Delauney, and Gabriel.
Delauney and Gabriel are parameter free, so given \code{parvec} has no meaning. In geometric graph, \code{parvec} is a vector
 of distances (sizes of the surrounding 'disc') and for k-nn \code{parvec} is the vector of neighbourhood abundances for each point 
to consider in the calculation of the spatial exposure measures.

The basic type of spatial summary uses range 'r', or 'geometric' graph with varying neighbourhood parameter 'r'. We default to 'k-nn'
graph as the exposure effect should be free of spatial scattering.

For \code{geometric} and \code{knn}, the calculations are done by shrinking the graph given by the largest value of \code{parvec}. If dealing with large datasets, 
it is advisable to give preprocessing range, \code{prepR}. The algorithm first calculates a geometric graph with parameter
\code{prepR}, and uses this as basis for finding the needed neighbourhoods. Speeds up calculations. 
\code{prepGraph}, if given, functions like the preprocess geometric graph. But make sure \code{prepR} is large enough (e.g. in \code{geometric}, \code{prepR}>max(\code{parvec})).

The \code{doDists} option speeds up calculations by precomputing the pairwise distances but takes O(n\^2) memory! 

For border correction, use \code{minusR} for reduced border correction (for rectangular windows only). If using \code{geometric} or \code{knn} neighbourhoods,
the option \code{toroidal} for toroidal correction is also available. The vector \code{included} can be given for more specific \code{minus}-correction,
only those points with TRUE (1) value are used in calculation. However, the neighbourhoods are calculated with all points.

  \tabular{ll}{
    Date: \tab 2009-03-09\cr
    License: \tab GPL v2 or later\cr
  }
}

\author{
Tuomas Rajala
University of Jyvaskyla, Finland
tarajala@maths.jyu.fi
}

\keyword{ spatial }

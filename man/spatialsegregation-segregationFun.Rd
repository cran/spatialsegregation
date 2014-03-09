\name{spatialsegregation-segregationFun}
\alias{segregationFun}
\alias{minusID.gf}
\alias{kGraphs} 
\alias{sg.modify.pp}
\docType{package}

\title{Spatial Segregation Function}

\description{
Compute the spatial exposure (segregation vs. mingling) features from a given multitype point pattern. Usage of 
shortcuts \code{minglingF}, \code{isarF}, \code{shannonF} and \code{simpsonF} highly recommended.
}

\usage{
segregationFun(X, fun="isar", r=NULL, ntype="geometric", funpars=NULL, 
toroidal=FALSE, minusRange=TRUE,  included=NULL, dbg=FALSE, 
doDists=FALSE, prepRange=0.0, prepGraph=NULL, prepGraphIsTarget=FALSE,
weightMatrix=NULL, translate=FALSE)
}

\arguments{
\item{X}{Multitype point pattern of class \code{ppp} (see package 'spatstat')}
\item{fun}{Default "isar". Takes "isar","mingling","shannon" and "simpson", see below.}
\item{r}{Vector for the neighbourhood defining graph, e.g. "geometric" graph with different ranges. See below.}
\item{ntype}{Default "geometric". Type of the neighbourhood graph. Accepts: "knn", "geometric", "delauney", "gabriel".}
\item{funpars}{Default NULL. Parameter(s) for the measure. Mingling: c(i,j), where i= only for type i (0 for all), j=1 -> ratio version. ISAR: i, i=type (integer). Shannon: 0 or 1, see \code{v2} in \code{shannonF}. Simpson: none.}
\item{toroidal}{Default FALSE. If TRUE, use a toroidal correction in distance calculation. Works at the moment only for rectangular windows and "geometric" or "knn" graph.}
\item{minusRange}{If TRUE, adaptive minus-sampling is employed. Overrides \code{included}-vector. If given as a positive number, \code{included}-vector is created with points with distance atleast minusRange from the border.}
\item{included}{boolean-vector of length |pp|. included[i]==TRUE => pp[i] included in calculations. Used for minus-sampling border correction.}
\item{dbg}{Default FALSE. Print additional runtime texts.}
\item{doDists}{Default TRUE. Precalculate distances for speed. Be aware of memory requirements, n*(n-1)!}
\item{prepRange}{Default 0. If >0, shrink the search space for neighbourhoods by searching only points within distance R i.e. precalculates a geometric graph.}
\item{prepGraph}{Precalculated graph for the point pattern. If given, The \code{prepRange}, \code{dodists} and \code{toroidal} are ignored and calculations are carried 
      using the prepGraph as a starting point. Useful for huge datasets.}
\item{prepGraphIsTarget}{If TRUE, precalculated graph \code{prepGraph} is used to calculate a single function value directly, all other neighbourhood parameters are ignored.}
\item{weightMatrix}{See \code{isarF} for this.}
\item{translate}{Use translation correction (see e.g. documentation of \code{spatstat::Kest} for details). Used only in mingling index.}
}
\details{
This is the general function for computing the spatial exposure (segregation/mingling) features. Used by minglingF, shannonF, simpsonF and isarF,
which should be preferred for better (and nicer) outcome.

Possible neighbourhood relations for the spatial version include geometric, k-nearest neighbours, Delauney, and Gabriel.
Delauney and Gabriel are parameter free, so given \code{r} has no meaning. In geometric graph, \code{r} is a vector
 of distances (sizes of the surrounding 'disc') and for k-nn \code{r} is the vector of neighbourhood abundances for each point 
to consider in the calculation of the spatial exposure measures. The basic type of spatial summary uses range, or 'geometric' 
graph connections with varying neighbourhood parameter. 

For \code{geometric} and \code{knn}, the calculations are done by shrinking the graph given by the largest value of \code{r}. If dealing with large datasets, 
it is advisable to give preprocessing range, \code{prepRange}. The algorithm first calculates a geometric graph with parameter
\code{prepRange}, and uses this as basis for finding the needed neighbourhoods. Speeds up calculations. \code{prepGraph}, if given, works as the 
preprocessed geometric graph. But make sure \code{prepRange} is large enough (e.g. in \code{geometric}, \code{prepRange}>max(\code{r})).

The \code{doDists} option speeds up calculations by precomputing the pairwise distances but takes n*(n-1) memory! 

For border correction, use \code{minusRange} for reduced border correction (for rectangular windows only). If using \code{geometric} or \code{knn} neighbourhoods,
the option \code{toroidal} for toroidal correction is also available. The vector \code{included} can be given for more specific \code{minus}-correction,
only those points with TRUE (1) value are used in calculation. However, the neighbourhoods are calculated with all points.

}

\value{
	Returns an object of class \code{fv}, see spatstat for more details. Basically a list with the computed values and parameter values.
}

\author{
Tuomas Rajala
University of Jyvaskyla, Finland
tuomas.rajala@iki.fi
}

\keyword{ spatial }

#the generic function for segregation measures
# isarF, minglingF, shannonF, simpsonF and mciF
#  are the wrappers that should be used.
#
#Tuomas Rajala <tuomas.rajala@jyu.fi>
# last change : 280410 
#################################################
#constant: the supported neighbourhoods
kGraphs<-c("geometric","knn","gabriel","delauney")

segregationFun<-function(X, fun="isar", r=NULL, ntype="geometric", funpars=NULL, 
		                 toroidal=FALSE, minusRange=0,  included=NULL, dbg=FALSE, 
						 doDists=FALSE, prepRange=0.0, prepGraph=NULL, prepGraphIsTarget=FALSE)
# function types:
#	1 mingling
#   2 shannon
#   3 simpson
#   4 ISAR
#	5 MCI

{
	# a note about the neighbourhoods
	note<-NULL
	# turn the neighbourhood type into an integer
	ntypei<-charmatch(ntype, kGraphs) - 1 #minus1 for c-side...
	# check 
	if(is.na(ntypei))stop("Error segregationFun: Wrong neighbourhood type.")
	
	# init neighbourhood parameters if not given: geometric range taken from Kest in spatstat
	# r not given
	if(is.null(r))
	{
		# snip
		if(ntypei==0)
		{
			
			W <- X$window
			npoints <-X$n
			lambda <- npoints/area.owin(W)
			rmaxdefault <- rmax.rule("G", W, lambda)
			breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
			rvals <- breaks$r
			rmax  <- breaks$max
			parvec<-rvals
		}
		# end snip
		# defaults for knn, gabriel and delauney. TODO: k-nn not yet clear what is a good range
		else
			parvec <- switch(ntype, "knn"=1:20, 0)
	}
	else parvec <- r # r given
	
	# TODO: for the c++-module coherence with 'spatgraphs', skip mass geometric
	if(ntypei>1) ntypei <- ntypei+1  
	
	# if a precalculated Graph is given, check its from spatgraphs
	if(class(prepGraph)!="sg" & !is.null(prepGraph) )stop("Error segregationFun: Prepared graph is not of class 'sg'.")
	
	#add a note about prepGraph
	if(!is.null(prepGraph)) note<-c(note,paste("PrepGraph given, type ",prepGraph$gtype,", par ",prepGraph$par,sep=""))
	
	# turn the wanted function type into an integer
	funi <- charmatch(fun, (kFuns<-c("mingling","shannon","simpson","isar","mci")))
	if(is.na(funi)) stop("Error segregationFun: wrong function type.")
	
	# TODO: modify the pp, see below in the function	
	X<-sg.modify.pp(X)
	
	# if a minus  (border) correction is to be used, compute the ones to exclude
	if(minusRange>0)
	{
		included<-minusID.gf(X, minusRange)
		note<-c(note,paste("Minus correction, radius=", minusRange, ";", sep=""))
	}
	# if we accept that all points are good for computation, included vector is all 1's
	if(is.null(included) | length(included)!=X$n) included<-rep(1,X$n)
	
	# check if the prepGraph is given when used as the target neighbourhood configuration
	if(prepGraphIsTarget && is.null(prepGraph)) stop("Error segregationFun: prepGraph not given but needed for calculation.")
		
	# this is for the c-side check to know if we are giving a prepGraph
	prepGraph$'isnull'<- as.integer(is.null(prepGraph))
	
	# set parameters according to prepGraph
	if(prepGraphIsTarget)
	{
		parvec<-prepGraph$parameters
		ntype<-prepGraph$type
	}
		
	
	# include lambdas
	X$area<-as.numeric(area.owin(X$window))
	
	# the main call 
	res<-.External("fun_c", as.integer(dbg), X, as.numeric(funpars), 
					as.integer(ntypei), as.numeric(parvec), 
					as.integer(funi), as.integer(toroidal), 
					as.numeric(prepRange), as.integer(doDists), 
					as.integer(included), prepGraph, 
					as.integer(prepGraphIsTarget),
					PACKAGE="spatialsegregation")

	# turn the result list into a matrix: 1 col per valuetype, 1 row per parameter
	a<-t(matrix(unlist(res),ncol=length(parvec)))

	# name the rows par1, par2, ...
	#rownames(a)<-paste("par",1:length(parvec),sep="")#no need for this anymore
 
	# the unit name: different from Kest etc, reports the meaning of the unit depending on -
	#   the neighbourhood type
	unitname<-switch(ntypei,"1"="range","2"="neighbour",NULL)
	
	list(v=a, included=as.logical(included), 
		 parvec=parvec, unitname=unitname, 
		 ntype=ntype, note=note
        )
}


#####################

minusID.gf<- function (pp0, minusRange)
{
	id <- (pp0$x < (pp0$window$x[2] - minusRange)) & (pp0$x > (pp0$window$x[1] +
					minusRange)) & (pp0$y < (pp0$window$y[2] - minusRange)) & (pp0$y >
				(pp0$window$y[1] + minusRange))
	id
}

#####################
sg.modify.pp<-function(pp)
{
	n<-length(pp[["x"]])
	
	if(length(pp[["mass"]]) != n ) # set the masses
	{
		if(length(pp[["marks"]])< n | !is.numeric(pp[["marks"]])) pp$mass<-rep(1.0,n)
		else pp$mass<-pp$marks
	}
	pp$mass<-as.numeric(pp$mass)
	
	if(length(pp[["types"]]) < n) # set the types
	{
		if( (is.factor(pp$marks) | (is.integer(pp$marks)) & length(pp[["marks"]])==n) ) x<-as.factor(pp$marks) 
		else x<-as.factor(rep(1,n))
	}
	else x<-as.factor(pp$types)
	y<-rep(1,n)
	m<-levels(x)
	for(i in 1:length(m))
	{
		y[ which(x==m[i]) ]<-i
		
	}
	pp$types<-as.integer(y)
	pp$marks<-NULL
	
	if(is.null(pp[["z"]]) || length(pp[["z"]])!=length(pp[["x"]])) pp$z<-rep(0.0,n) # if 2D only
	if(is.null(pp[["window"]][["z"]])) pp$window$z<-as.numeric(c(0.0,1.0)) # if 2D only
	
	pp$window$x<-as.numeric(pp$window$x)
	pp$window$y<-as.numeric(pp$window$y)
	pp$window$z<-as.numeric(pp$window$z)
	pp$x<-as.numeric(pp$x)
	pp$y<-as.numeric(pp$y)
	pp$z<-as.numeric(pp$z)	
	pp
}

#####################


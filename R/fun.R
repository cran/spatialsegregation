#the generic function for segregation measures
#

segregationFun<-function(pp, fpar=NULL, graph_parvec=1:20, graph_type="knn", 
		                 toroidal=FALSE, dbg=FALSE, funtype=1, 
						 doDists=FALSE, prepR=0.0, prepGraph=NULL, prepGraphIsTarget=FALSE, included=NULL, minusR=NULL, relative=FALSE)
# funtypes:
#	1 mingling
#   2 shannon
#   3 simpson
#   4 ISAR

{
	nGRAPHS<-c("geometric","knn","gabriel","delauney")
	if(!graph_type%in%nGRAPHS)stop("Error segregationFun: Wrong graphtype.")
	typei<-which(graph_type==nGRAPHS)-1
	
	if(typei>2)graph_parvec=0
	if(typei>1) typei=typei+1 # for the c++-module coherence with 'spatgraphs' 
	if(class(prepGraph)!="sg" & !is.null(prepGraph) )stop("Error segregationFun: Prepared graph is not of type 'sg'.")
	if(is.null(funtype))stop("Error segregationFun: wrong function type.")
	pp<-sg_modify_pp(pp)
	
	note<-""
	if(relative & typei==0){sum0<-summary(pp);note<-"r scaled with sqrt(pi*lambda)."; graph_parvec<- graph_parvec/sqrt(pi*sum0$int)}
	
	if(!is.null(minusR)) included<-minusID_gf(pp, minusR)
	if(is.null(included) | length(included)!=pp$n) included<-rep(1,pp$n)
	
	if(prepGraphIsTarget && is.null(prepGraph)) stop("Error segregationFun: prepGraph not given but needed for calculation.")
	prepGraph$'isnull'<- as.integer(is.null(prepGraph))
	
		
	res<-.External("fun_c", as.integer(dbg), pp, as.numeric(fpar), 
					as.integer(typei), as.numeric(graph_parvec), 
					as.integer(funtype), as.integer(toroidal), 
					as.numeric(prepR), as.integer(doDists), 
					as.integer(included), prepGraph, 
					as.integer(prepGraphIsTarget),
					PACKAGE="spatialsegregation")
	
			
	a<-matrix(unlist(res),ncol=length(graph_parvec))
	colnames(a)<-paste("par",1:length(graph_parvec),sep="")
	rownames(a)<-paste("type",1:dim(a)[1],sep="")
	list(v=a, included=as.logical(included), note=note, parvec=graph_parvec)
}
#####################

minusID_gf<- function (pp0, minusR)
{
	id <- (pp0$x < (pp0$window$x[2] - minusR)) & (pp0$x > (pp0$window$x[1] +
					minusR)) & (pp0$y < (pp0$window$y[2] - minusR)) & (pp0$y >
				(pp0$window$y[1] + minusR))
	id
}

#####################
sg_modify_pp<-function(pp)
{
	n<-length(pp[["x"]])
	if(length(pp[["mass"]]) < n ) # set the masses
	{
		if(length(pp[["marks"]])< n | !is.numeric(pp[["marks"]])) pp$mass<-rep(1.0,n)
		else pp$mass<-pp$marks
	}
	if(length(pp[["types"]]) < n) # set the types
	{
		if( (is.factor(pp$marks) | is.integer(pp$marks)) & length(pp[["marks"]])==n ) pp$types<-pp$marks 
		else pp$types<-rep(1,n)
	}
	pp$mass<-as.numeric(pp$mass)
	pp$types<-as.integer(pp$types)
	
	if(is.null(pp[["z"]]) || length(pp[["z"]])!=length(pp[["x"]])) pp$z<-rep(0.0,n) # if 2D only
	if(is.null(pp[["window"]][["z"]])) pp$window$z<-as.numeric(c(0.0,1.0)) # if 2D only
	pp$marks<-NULL
	pp$window$x<-as.numeric(pp$window$x)
	pp$window$y<-as.numeric(pp$window$y)
	pp$window$z<-as.numeric(pp$window$z)
	pp$x<-as.numeric(pp$x)
	pp$y<-as.numeric(pp$y)
	pp$z<-as.numeric(pp$z)
	pp
}

#####################


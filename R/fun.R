#the generic function for segregation measures
#

segregationFun<-function(pp, fpar=NULL, graph_parvec=1:20, graph_type="knn", 
		                 toroidal=FALSE, dbg=FALSE, funtype=1, 
						 doDists=FALSE, prepR=0.0, prepGraph=NULL, included=NULL, minusR=NULL, relative=FALSE)
# funtypes:
#	1 mingling
#   2 shannon
#   3 simpson
#   4 ISAR

{
	nGRAPHS<-c("geometric","knn","gabriel","delauney")
	if(!graph_type%in%nGRAPHS)return("Error segregationFun: Wrong graphtype.")
	typei<-which(graph_type==SG_SUPPORTED_GRAPHS)-1
	
	if(typei>1)graph_parvec=0
	
	if(class(prepGraph)!="sg" & !is.null(prepGraph) )return("Error segregationFun: Prepared graph is not of type 'sg'.")
	if(is.null(funtype))return("Error segregationFun: wrong function type.")
	pp<-sg_modify_pp(pp)
	
	note<-""
	if(relative & typei==0){sum0<-summary(pp);note<-"r scaled with sqrt(pi*lambda)."; graph_parvec<- graph_parvec/sqrt(pi*sum0$int)}
	
	if(!is.null(minusR)) included<-minusID_gf(pp, minusR)
	if(is.null(included) | length(included)!=pp$n) included<-rep(1,pp$n)
	prepGraph$'isnull'<- as.integer(is.null(prepGraph))
	
	
	res<-.External("fun_c", as.integer(dbg), pp, as.numeric(fpar), 
					as.integer(typei), as.numeric(graph_parvec), 
					as.integer(funtype), as.integer(toroidal), 
					as.numeric(prepR), as.integer(doDists), 
					as.integer(included), prepGraph, 
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
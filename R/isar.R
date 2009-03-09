# mingling.R
# 
# ISAR function for multitype spatial point pattern and graphs
# see: T. Wiegand et al.: How individual species structure diversity in tropical forests. PNAS, November 16 2007

# Author: Tuomas Rajala <tarajala@maths.jyu.fi>
###############################################################################
# 


isarF<-function(pp, parvec=1:20, graph_type="knn", type=NULL, ...)
#ISAR for various graphs
#If target type is not given (type=NULL) take over all types mean
#type as integer 1,...,S=number of types

{
	
	note2<-paste("ISAR for type",type)
	if(is.null(type))
	{
		type<-0
		note2<-"ISAR over all types."
	}

	res<-segregationFun(pp, fpar=type, graph_type=graph_type, graph_parvec=parvec, funtype=4, ...)
	
#   the poisson values

	if(graph_type=="geometric")mdeg<-function(l,k)pi*l*k^2
	if(graph_type=="knn")mdeg<-function(l,k)k
	if(graph_type=="delauney")mdeg<-function(l,k) 6
	if(graph_type=="gabriel")mdeg<-function(l,k) 4
	
	poisson<-NULL
	sum0<-summary(pp)
	l<-sum0$marks[,3]
	for(i in 1:length(parvec)) poisson<-c(poisson, sum(1-(1-l/sum(l))^mdeg(sum(l),parvec[i]) ) )
	
	if(pp$markformat=="none")pp$markformat<-"vector" # TODO: fix this bug
	f<-freqs(pp[res$included])
	w<-f/sum(f)
	
	segfcl(list(I=apply(res$v,2,mean),sd=apply(res$v,2,sd),typewise=res$v,
					Iw=apply(res$v,2,weighted.mean,w=w),
					gtype=graph_type,par=res$parvec,note=res$note,note2=note2,poisson=poisson))
}

###############################################################################

###############################################################################
isar_index<-function(pp, graph_type="knn", graph_par=4, type=NULL, ...)
#if type=NULL compute the arithmetic mean over all types
#type as integer 1,...,S=number of types
{
	if(is.null(type)){type<-0;nvec<-paste("Mean ISAR over all types.")}
	else nvec<-paste("ISAR for type",type)
	I0<-isarF(pp=pp, graph_type=graph_type, parvec=graph_par[1], type=type, ...)
	I<-I0$I
	names(I)<-nvec
	I
}
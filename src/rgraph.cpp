#include <R.h>
#include "Graph.h"
#include "dists.h"
#include "Pp.h"


extern "C" {

SEXP graph_c(SEXP Args)
{
//	e<-.Extern("graph_c", pp, typei, par , preprocessR, as.integer(toroidal), as.integer(doDists), dbg)

	Pp pp;
	double *prepR, *par;
	int *gtype, *doDists, *toroidal, *dbg;
	Graph graph;

//start parsing the args
	Args = CDR(Args);
	pp.Init(CAR(Args)); // init pp

	Args = CDR(Args);
	gtype = INTEGER(CAR(Args)); //what type of graph

	Args = CDR(Args);
	par = REAL(CAR(Args)); // graph par

	Args = CDR(Args);
	prepR = REAL(CAR(Args)); // if preprocessing 
	
	Args = CDR(Args);
	toroidal = INTEGER(CAR(Args)); // if toroidal correction
	
	Args = CDR(Args);
	doDists = INTEGER(CAR(Args)); // if the distances are precalculated and stored
	
	Args = CDR(Args);
	dbg = INTEGER(CAR(Args)); // if debug messages
				
//	void Init(Pp *pp0, double *par, double *prepR, int *doDists, int *toroidal, int *dbg );
	graph.Init(&pp, gtype, par, prepR, doDists, toroidal, dbg);
	graph.sg_calc();
	
	if(*dbg)printf("\n");
	return graph.toSEXP();
}

	

}

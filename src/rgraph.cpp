#include <R.h>
#include "Graph.h"
#include "Pp.h"


extern "C" {

//SEXP graph_c(SEXP Args)
//{
//	Pp pp;
//	double *prepR, *par, *d0;
//	int *gtype, *doDists, *toroidal, *dbg;
//	Graph graph;
//	d0 = new double;
//	d0[0] =0.0;
////start parsing the args
//	Args = CDR(Args);
//	pp.Init(CAR(Args)); // init pp
//
//	Args = CDR(Args);
//	gtype = INTEGER(CAR(Args)); //what type of graph
//
//	Args = CDR(Args);
//	par = REAL(CAR(Args)); // graph par
//
//	Args = CDR(Args);
//	prepR = REAL(CAR(Args)); // if preprocessing
//
//	Args = CDR(Args);
//	toroidal = INTEGER(CAR(Args)); // if toroidal correction
//
//	Args = CDR(Args);
//	doDists = INTEGER(CAR(Args)); // if the distances are precalculated and stored
//
//	Args = CDR(Args);
//	dbg = INTEGER(CAR(Args)); // if debug messages
//
//	graph.Init(&pp, gtype, par, prepR , doDists, d0, toroidal,  *gtype, dbg);
//	graph.sg_calc();
//
//	if(*dbg)printf("\n");
//	return graph.toSEXP();
//
//}


}

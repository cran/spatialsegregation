#include <R.h>
#include "Fun.h"



extern "C" {

//	res<-.External("fun_c", as.integer(dbg), pp, as.numeric(fpar),
//						as.integer(typei), as.numeric(graph_parvec),
//						as.integer(funtype), as.integer(toroidal),
//						as.numeric(prepR), as.integer(doDists),
//						as.integer(included), prepGraph,
//						PACKAGE="spatialsegregation")

SEXP fun_c(SEXP Args)
{

	Pp pp;
	Graph graph;
	Fun fun;
	double *prepR, *fpar, *par, *parvec;
	int *gtype, *doDists, *toroidal, *ftype, *dbg, parn, *incl, prepG=0, *prepGraphIsTarget;
	SEXP prepGraph;

//start parsing the args
	Args = CDR(Args);
	dbg = INTEGER(CAR(Args)); // if debug messages
	if(*dbg)printf("Parsing parameter:");

	Args = CDR(Args);
	pp.Init(CAR(Args)); // init pp

	if(*dbg)printf(".");
	Args = CDR(Args);
	fpar = REAL(CAR(Args)); // additional function parameters.

	if(*dbg)printf(".");
	Args = CDR(Args);
	gtype = INTEGER(CAR(Args)); //what type of graph

	if(*dbg)printf(".");
	Args = CDR(Args);
	parvec = REAL(CAR(Args)); // graph parvec
	parn = length(CAR(Args));

	if(*dbg)printf(".");
	Args = CDR(Args);
	ftype = INTEGER(CAR(Args)); //what type of function

	if(*dbg)printf(".");
	Args = CDR(Args);
	toroidal = INTEGER(CAR(Args)); // if toroidal correction

	if(*dbg)printf(".");
	Args = CDR(Args);
	prepR = REAL(CAR(Args)); // if preprocessing R given

	if(*dbg)printf(".");
	Args = CDR(Args);
	doDists = INTEGER(CAR(Args)); // if the distances are precalculated and stored

	if(*dbg)printf(".");
	Args = CDR(Args);
	incl = INTEGER(CAR(Args)); // the inclusion vector

	if(*dbg)printf(".");
	Args = CDR(Args);
	prepGraph = CAR(Args); // possibly precalculated graph
	prepG = 1- INTEGER(getListElement(prepGraph,"isnull"))[0];

	if(*dbg)printf(".");
	Args = CDR(Args);
	prepGraphIsTarget = INTEGER(CAR(Args)); // use only the precalculated graph to get a value


	if(*dbg)printf("done.\n");
	par = &parvec[parn-1];



	//	void Init(Pp *pp0, double *par, double *prepR, int *doDists, int *toroidal, int *dbg );
	if(*dbg)printf("Init graph...");
	graph.Init(&pp, gtype, par, prepR , doDists, toroidal, dbg);
	if(prepG)// if precalculated graph, set the edges, overrule the prepR-parameter
	{
		if(*dbg)printf("loading precalculated edges...");
		std::vector<std::vector<int> > prepNodelist;
		VectsxpToVector(getListElement(prepGraph,"edges"), prepNodelist);
		graph.setNodelist(&prepNodelist);
		graph.prepR = REAL(getListElement(prepGraph,"parameters"));
		*graph.oldpar = *graph.prepR;
		*graph.prepDone = 1;
	}



	//	void Init(Graph *g0, double *par0, int *parn, int *gt, int *ft, double *fpar, int *dbg0);
	if(*dbg)printf("Init fun...");
	fun.Init(&graph, parvec, &parn, gtype, ftype, fpar, incl, dbg);
	if(*dbg)printf("done.\n");
	if(*dbg)printf("Calculating:\n");

	if(*prepGraphIsTarget){ fun.re_calculate(); }
	else{  fun.calculate();  }

	if(*dbg)printf("done.\n");
	return fun.toSEXP();
}



} //extern

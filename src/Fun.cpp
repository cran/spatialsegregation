#include "Fun.h"

Fun::Fun()
{
}

Fun::~Fun()
{
}


void Fun::Init(Graph *g0, double *par0, int *parn, int *gt, int *ft, double *fpar0, int *trans0, int *included0, int *dbg0)
{
	int i;
	graph = g0;
	std::vector<double> * pvec;
	for(i=0;i<*parn;i++)
	{
		parvec.push_back(par0[i]);
		pvec = new std::vector<double>;
		value.push_back(*pvec);
		value.at(i).push_back(0.0);
		delete pvec;
	}
	included = included0;
	gtype = gt;
	ftype = ft;
	fpar = fpar0;
	dbg = dbg0;

	trans = trans0;
	if(*trans>0){
		graph->pp->calcTransWeights();
	}
}


SEXP Fun::toSEXP(SEXP pp)
//transform a std::vector<double > to SEXP
{
	SEXP res, *onevalue;
	double *p, *m;

	PROTECT(res = allocVector(VECSXP, this->value.size()));
	int i,j;
	for(i=0;i < (int)this->value.size();i++)
	{
		onevalue = new SEXP;
		PROTECT(*onevalue = allocVector(REALSXP, this->value.at(i).size()));
		p = REAL(*onevalue);

		for(j=0;j<(int)this->value.at(i).size();j++)
			p[j] = this->value.at(i).at(j);
		SET_VECTOR_ELT(res, i, *onevalue);
		UNPROTECT(1);
		delete onevalue;
	}
	UNPROTECT(1);
//	set the mass element of pp
	m = REAL(getListElement(pp,"mass"));
	for(int i=0; i< graph->pp->size();i++)
		m[i] = graph->pp->getMass2(&i);
	return res;
}


void Fun::calculate()
{
	int i;
	std::vector<double> resvec;
	for(i=parvec.size()-1 ; i >= 0 ; i--)
	{
		if(*this->dbg)printf("Fun %i/%i: graph[",(int)parvec.size()-i,(int)parvec.size());

		// update graph
		graph->par = &parvec[i];
		graph->sg_calc();
		*graph->oldpar = *graph->par;

		if(*this->dbg)printf("] Value[ ");
		// calc index
		if(*ftype == 1)
			resvec = mingling(graph, fpar, dbg, included);
		if(*ftype == 2)
			resvec = shannon(graph, fpar, dbg, included);
		if(*ftype == 3)
			resvec = simpson(graph, fpar, dbg, included);
		if(*ftype == 4)
			resvec = isar(graph, fpar, dbg, included);
		if(*ftype == 5)
			resvec = mci(graph, fpar, dbg, included);
		if(*ftype == 6)
			resvec = biomass(graph, fpar, dbg, included);
				value.at(i) = resvec;
		if(*this->dbg)printf(" ]\n");
	}
}

void Fun::re_calculate()
{
		std::vector<double> resvec;
		if(*this->dbg)printf(" Value[ ");
		// calc index
		if(*ftype == 1)
			resvec = mingling(graph, fpar, dbg, included);
		if(*ftype == 2)
			resvec = shannon(graph, fpar, dbg ,included);
		if(*ftype == 3)
			resvec = simpson(graph, fpar, dbg, included);
		if(*ftype == 4)
			resvec = isar(graph, fpar, dbg, included);
		if(*ftype == 5)
			resvec = mci(graph, fpar, dbg, included);
		value.at(0) = resvec;
		if(*this->dbg)printf(" ]\n");
}


#include "Fun.h"

Fun::Fun()
{
}

Fun::~Fun()
{
}


void Fun::Init(Graph *g0, double *par0, int *parn, int *gt, int *ft, double *fpar0, int *included0, int *dbg0)
{
	int i,j,empty;
	graph = g0;
	std::vector<double> * pvec;
	for(i=0;i<*parn;i++)
	{
		pvec = new std::vector<double>;
		parvec.push_back(par0[i]);
		value.push_back(*pvec);
		value.at(i).push_back(0.0);
	}
	included = included0;
	gtype = gt;
	ftype = ft;
	fpar = fpar0;
	dbg = dbg0;

	for(i=0; i< *graph->pp->S;i++)//need to put 0 to some lambdas based on the inclusion vector
	{
		empty=1;
		for(j=0;j<*graph->pp->n;j++)
			if(graph->pp->type[j]==i+1 && included[j])
			{
				empty=0;
				break;
			}
		if(empty)graph->pp->lambdas[i]=0.0;

	}
}


SEXP Fun::toSEXP()
//transform a std::vector<double > to SEXP
{
	SEXP res, *onevalue;
	double *p;
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
	}
	UNPROTECT(1);
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
		*graph->oldpar = *graph->par;
		graph->par = &parvec[i];
		graph->sg_calc();

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
			value.at(0) = resvec;
		if(*this->dbg)printf(" ]\n");
}


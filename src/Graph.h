#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <math.h>
#include <vector>
#include "Pp.h"

#ifndef GRAPH_H_
#define GRAPH_H_

class Graph
{

public:
	Pp *pp; //the point pattern
	double *par;
	double opar;
	double *oldpar;
	int	   *doDists;
	int    *dbg;
	double *prepR;
	int    *gtype;
	int    *inc;
	double  mdeg;
	std::vector<std::vector<int> > nodelist;
	std::vector<int> typeIncluded;
	Graph();
	virtual ~Graph();

	void Init(Pp *pp0, int *gtype0, double *par, double *prepR, int *doDists, double *preDists, int *toroidal, int *inc, int *dbg );
	void setNodelist(std::vector<std::vector<int> > *nodelist_new);
	void setNodelist(SEXP);
	void addNew(int , int);
	void sg_calc();

	void sg_geometric();
	void sg_geometric(double *);
	void sg_shrink_geometric(double *);
	void sg_mass_geometric();
	void sg_knn();
	void sg_shrink_knn();
	void sg_gabriel();
	void sg_delaunay();
	void sg_MST();
	void sg_markcross();
	void sg_SIG();
	void sg_RST();
	void sg_RNG();
	void sg_CCC();
	void sg_STIR();

	void sg_cut(double *R);
	void sg_prune(double *lev);


	void remove_duplicates();
	SEXP toSEXP();
};

int compare_doubles(const void *a, const void *b);
double Attenuate(double r, double alpha); // used by STIR graph

#endif /*GRAPH_H_*/

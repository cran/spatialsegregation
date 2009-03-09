#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <math.h>
#include <vector>
#include "dists.h"
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
	int    *toroidal;
	int    *dbg;
	double *prepR;
	int    *gtype;
	std::vector<std::vector<int> > nodelist;
	std::vector<double> distTriangle;
	std::vector<double> * pdists;
	double (Graph::*Dist)(int*, int*);
	double Dist1(int*, int*);
	double Dist2(int*, int*);
	Graph();
	virtual ~Graph();

	void Init(Pp *pp0, int *gtype0, double *par, double *prepR, int *doDists, int *toroidal, int *dbg );
	void setNodelist(std::vector<std::vector<int> > *nodelist_new);
	void addNew(int , int);
	void sg_calc();

	void sg_geometric();
	void sg_geometric(double *);
	void sg_shrink_geometric(double *);
	void sg_mass_geometric();
	void sg_knn();
	void sg_shrink_knn();
	void sg_gabriel();
	void sg_delauney();
	void sg_MST();
	void sg_markcross();
	void sg_SIG();
	void sg_RST();
	void sg_RNG();
	void sg_CCC();

	void sg_cut(double *R);
	void sg_prune(double *lev);


	void remove_duplicates();
	SEXP toSEXP();
};

#endif /*GRAPH_H_*/

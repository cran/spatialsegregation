/* Functional class for the segregation measures
 * ftype=
 *        1   mingling
 *        2   shannon
 *        3   simpson
 *        4   ISAR
 *
 * Supports Geometric and k-nn graphs, and toroidal correction
 * TODO: border correction
 * by: Tuomas Rajala
 *
 * 	300508
 *
 */
#include <R.h>
#include <vector>
#include "Graph.h"
#include "mingling.h"
#include "shannon.h"
#include "simpson.h"
#include "isar.h"
#include "mean_sd.h"
#ifndef FUN_H_
#define FUN_H_


class Fun
{
	Graph *graph;
	std::vector<std::vector <double> > value;
	std::vector<double> parvec;
	int *gtype; // 0 = geometric, 1 = knn
	int *ftype;
	int *included;
	double *fpar;
	int *dbg;

public:
	Fun();
	virtual ~Fun();
	void Init(Graph *g0, double *par0, int *parn, int *gt, int *ft, double *fpar, int *included0, int *dbg0);
	void calculate();
	void re_calculate();
	SEXP toSEXP();
};

#endif /*FUN_H_*/

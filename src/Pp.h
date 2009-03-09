#include <Rdefines.h>
#include <Rinternals.h>
#include <R.h>
#include <vector>
#include "Rextras.h"
#ifndef PP_H_
#define PP_H_

class Pp
{
public:
	double *x;
	double *y;
	double *z;
	int m;
	int *n;
	int *type;
	double *mass;
	int s;
	int *S;
	std::vector<double > lambdas;
	double *xlim;
	double *ylim;
	double *zlim;
	Pp();
	virtual ~Pp();
	void Init(double *x, double *y, double *z, int *type, double *mass0, int *n, double *xlim, double *ylim, double *zlim);
	void Init(SEXP);
};

#endif /*PP_H_*/

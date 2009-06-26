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
	int toroidal;
	int *S;
	std::vector<double > lambdas;
	double *xlim;
	double *ylim;
	double *zlim;
	std::vector<double> distTriangle;
	std::vector<double> * pdists;

	Pp();
	virtual ~Pp();



	void Init(double *x, double *y, double *z, int *type, double *mass0, int *n, double *xlim, double *ylim, double *zlim);
	void Init(SEXP);

	void toggleToroidal();
	double getDist(int*,int*);
	double (Pp::*getDistp)(int*, int*);
	double Dist1(int*, int*);
	double Dist2(int*, int*);
	void calcDists();
};

#endif /*PP_H_*/

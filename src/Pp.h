#include <Rdefines.h>
#include <Rinternals.h>
#include <R.h>
#include <vector>
#include "Point.h"
#include "Rextras.h"
#ifndef PP_H_
#define PP_H_

class Pp
{
	std::vector<Point> points;
	int m;
	int ntypes;
	int tor;
	int *toroidal;


	double (Pp::*dist)(int*, int*);
	double distEuclidian(int*, int*);
	double distPrecalculated(int*, int*);
	std::vector<double> distTriangle;
	std::vector<double> * pdists;
	std::vector<int> typevec;

public:
	std::vector<double > lambdas;
	double lambda;
	double *xlim;
	double *ylim;
	double *zlim;

	Pp();
	virtual ~Pp();

	void Init(double *x, double *y, double *z, int *type, double *mass0, int *n, double *xlim, double *ylim, double *zlim);
	void Init(SEXP);

	double getX(int *);
	double getY(int *);
	double getZ(int *);
	int    getT(int *);

	int	   getTypevec(int *);
	void   setToroidal(int *);
	int    size();
	int    getNtypes();


	double getDist(int *, int *);
	void   calcDists();
	void   setDists(double *);

	int    getCluster(int *);
	int    nsize(int*); // neighbours
	void   setMass(int *, double );
	double getMass(int *);
	void   setMass2(int *, double );
	double getMass2(int *);

	int Empty(int *, int *, int *); //voronoi

	void addNeighbour(int *i, int *j);
	void removeNeighbour(int *i, int *j);
	void clearNeighbourhood(int *i);
	int  getNeighbour(int *i, int *j);
	void setCluster(int *i, int *j);
	void movePoint(int *, double *, double *);
	SEXP toSEXP();
};

#endif /*PP_H_*/

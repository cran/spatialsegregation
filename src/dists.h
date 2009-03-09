#import <math.h>
#import <R.h>
#import <Rinternals.h>
#import <Rdefines.h>
#import <vector>
#import "Pp.h"

//double getDist0(double *x, double *y, double *z, int *n, int *i, int *j, std::vector<double> *dist);
double getDist(Pp *pp, int *i, int *j, int *toroidal);
double getDist(int *i, int *j, int *n, std::vector<double> *dist);
void calcDists(Pp *pp, std::vector<double> *dist, int *toroidal);
int compare_doubles(const void *a, const void *b);
int Empty(double *x, double *y, int *n, int i, int j, int k, double *x00, double *y00, double *R20);

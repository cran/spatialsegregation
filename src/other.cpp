/********************************************************
 *  functions for spatgraphs
 *
 * Tuomas Rajala  <tuomas@sokkelo.net>
 * 260609
 *
 * ****************/

#import "other.h"

/**********************************************************************************/
int compare_doubles(const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

/**********************************************************************************/

int Empty(double *x, double *y, int *n, int i, int j, int k, double *x00, double *y00, double *R20)
// check if the circumcircle of three point triangle is empty of other points.
// puts to (x00,y00,R20) the circumcircle
// See: http://mathworld.wolfram.com/Circumcircle.html
{
	int m;
	double x0,y0,R2,bx,by,a,c,d2,xxyy1,xxyy2,xxyy3,x13,x23,x21,y13,y21,y23;
	xxyy1 = x[i]*x[i]+y[i]*y[i];
	xxyy2 = x[j]*x[j]+y[j]*y[j];
	xxyy3 =	x[k]*x[k]+y[k]*y[k];
	y23 = y[j]-y[k];
	y13 = y[i]-y[k];
	y21 = y[j]-y[i];
	x23 = x[j]-x[k];
	x13 = x[i]-x[k];
	x21 = x[j]-x[i];
	bx = -( xxyy1*y23-xxyy2*y13-xxyy3*y21 );
	by =  ( xxyy1*x23-xxyy2*x13-xxyy3*x21 );
	a  = x[i]*y23-x[j]*y13-x[k]*y21;
	c  = -(xxyy1*(x[j]*y[k]-y[j]*x[k])-xxyy2*(x[i]*y[k]-y[i]*x[k])-xxyy3*(x[j]*y[i]-x[i]*y[j]));
	R2 = (bx*bx+by*by-4.0*a*c)/(4.0*a*a);
	x0 = -bx/(2.0*a);
	y0 = -by/(2.0*a);
	x00[0] = x0;
	y00[0] = y0;
	R20[0] = R2;
	for(m=0;m<*n;m++)
	{
		if(m!=i & m!=j & m!=k)
		{
			d2 = pow(x0-x[m],2)+pow(y0-y[m],2);
			if(d2<R2) return 0;
		}
	}
	return 1;
}
/**********************************************************************************/

#include "Pp.h"
/********************************************************************************************/
Pp::Pp()
{
}
/********************************************************************************************/
Pp::~Pp()
{
}
/********************************************************************************************/
void Pp::Init(SEXP Argspp)
{
	int i,j,old;
	double Area;
	m = length(getListElement(Argspp, "x"));
	n = &m;
	x = REAL(getListElement(Argspp, "x"));
	y = REAL(getListElement(Argspp, "y"));
	z = REAL(getListElement(Argspp, "z"));
	type = INTEGER(getListElement(Argspp, "types"));
	mass = REAL(getListElement(Argspp,"mass"));
	typevec.clear();
	for(i=0;i<m;i++)
	{
		old = 0;
		for(j=0;j<(int)typevec.size();j++)
			if(typevec.at(j)==type[i]){ old = 1;break;}
		if(!old)
			typevec.push_back(type[i]);
	}
	s = typevec.size();
	S = &s;
	toroidal = 0;
	xlim = REAL(getListElement(getListElement(Argspp, "window") ,"x"));
	ylim = REAL(getListElement(getListElement(Argspp, "window") ,"y"));
	zlim = REAL(getListElement(getListElement(Argspp, "window") ,"z"));
	Area = (xlim[1]-xlim[0])*(ylim[1]-ylim[0])*(zlim[1]-zlim[0]);
	for(i=0;i<s;i++)
	{
		lambdas.push_back(0.0);
		for(j=0;j<m;j++)
			if(type[j]==typevec.at(i))
				lambdas[i]=lambdas[i]+1.0;
		lambdas[i]=lambdas[i]/Area;
	}
	getDistp = &Pp::Dist1;
}
/********************************************************************************************/
void Pp::Init(double *x0, double *y0, double *z0, int *type0, double *mass0, int *n0, double *xlim0, double *ylim0, double *zlim0)
{
	int i,j,old;
	toroidal = 0;
	n = n0; m = *n;
	x = x0;
	y = y0;
	z = z0;
	type =type0;
	mass =mass0;
	typevec.clear();
	for(i=0;i<m;i++)
	{
		old = 0;
		for(j=0;j<(int)typevec.size();j++)
			if(typevec.at(j)==type[i]){ old = 1;break;}
		if(!old)
			typevec.push_back(type[i]);
	}
	s = typevec.size();
	S = &s;

	xlim = xlim0;
	ylim = ylim0;
	zlim = zlim0;
	getDistp = &Pp::Dist1;
}
/********************************************************************************************/
void Pp::toggleToroidal(){this->toroidal = 1-this->toroidal;}
double Pp::getDist(int *i, int *j) {return (this->*getDistp)(i,j);}
/********************************************************************************************/

void Pp::calcDists()
{
	int i,j,k, *n;
	double d;
	n = this->n;
	for(i=0;i<*n-1;i++)
	{
		k = (int) i*(*n)-i*(i+1)/2;
		for(j=i+1;j<*n;j++)
		{
			d = this->Dist1(&i, &j);
			pdists->push_back(d);
		}
	}
	getDistp = &Pp::Dist2;
}

/********************************************************************************************/
double Pp::Dist1(int *i, int *j)
{
	if(*i==*j) return 0.0;
	if(*i>*j) return Dist1(j, i);
	if(this->toroidal)
		return	sqrt(
					pow( fmin( this->xlim[1]-this->xlim[0]-fabs(this->x[*i]-this->x[*j]) , fabs(this->x[*i]-this->x[*j]) ) ,2.0) +
					pow( fmin( this->ylim[1]-this->ylim[0]-fabs(this->y[*i]-this->y[*j]) , fabs(this->y[*i]-this->y[*j]) ) ,2.0) +
					pow( fmin( this->zlim[1]-this->zlim[0]-fabs(this->z[*i]-this->z[*j]) , fabs(this->z[*i]-this->z[*j]) ) ,2.0)   );
	else
		return 	sqrt(
				pow( this->x[*i]- this->x[*j]  ,2.0) +
				pow( this->y[*i]- this->y[*j]  ,2.0) +
				pow( this->z[*i]- this->z[*j]  ,2.0)   );
}
/**********************************************************************************/


double Pp::Dist2(int *i, int *j)
{
	if(*i==*j) return 0.0;
	if(*i>*j) return Dist2(j, i);
	return this->pdists->at( *j-*i -1 + (int)((*i)*(*this->n)-(*i)*(*i+1)/2) ) ;
}

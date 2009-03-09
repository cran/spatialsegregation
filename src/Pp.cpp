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
	std::vector<int> temp;
	for(i=0;i<m;i++)
	{
		old = 0;
		for(j=0;j<(int)temp.size();j++)
			if(temp.at(j)==type[i]){ old = 1;break;}
		if(!old)
			temp.push_back(type[i]);
	}
	s = temp.size();
	S = &s;

	xlim = REAL(getListElement(getListElement(Argspp, "window") ,"x"));
	ylim = REAL(getListElement(getListElement(Argspp, "window") ,"y"));
	zlim = REAL(getListElement(getListElement(Argspp, "window") ,"z"));
	Area = (xlim[1]-xlim[0])*(ylim[1]-ylim[0])*(zlim[1]-zlim[0]);
	for(i=0;i<s;i++)
	{
		lambdas.push_back(0.0);
		for(j=0;j<m;j++)
			if(type[j]==i+1)
				lambdas[i]=lambdas[i]+1.0;
		lambdas[i]=lambdas[i]/Area;
	}
}
/********************************************************************************************/
void Pp::Init(double *x0, double *y0, double *z0, int *type0, double *mass0, int *n0, double *xlim0, double *ylim0, double *zlim0)
{
	int i,j,old;
	n = n0; m = *n;
	x = x0;
	y = y0;
	z = z0;
	std::vector<int> temp;
	type =type0;
	mass =mass0;
	for(i=0;i<m;i++)
	{
		old = 0;
		for(j=0;j<(int)temp.size();j++)
			if(temp.at(j)==type[i]){ old = 1;break;}
		if(!old)
			temp.push_back(type[i]);
	}
	s = temp.size();
	S = &s;

	xlim = xlim0;
	ylim = ylim0;
	zlim = zlim0;
}
/********************************************************************************************/

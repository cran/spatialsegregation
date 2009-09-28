#include "mingling.h"
/**********************************************************************************/

std::vector<double> mingling(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)printf("mingling[");
	int target_type;
	int i,j,n=0,m,neq, dbg0, incS=0;
	double a1[2];
	std::vector<double> value;
	value.clear();
	if(*dbg)printf("(type=%i,ratio=%i)",(int)fpar[0],(int)fpar[1]);
	if((int)fpar[0]==0)
	{
		dbg0 = *dbg;
		*dbg = 0;
		a1[1] = fpar[1];
		for(i=0;i< *graph->pp->S;i++)
		{
			if(graph->pp->lambdas[i]>0)
			{
				incS++;
				a1[0] = (double) graph->pp->typevec.at(i);
				value.push_back(mingling(graph, a1, dbg ,included).at(0));
			}

		}
		//M = M/((double) incS);
		*dbg = dbg0;
	}
	else // target type given
	{
		target_type = (int) fpar[0];
		value.push_back(0.0);
		for(i=0;i< (int)graph->nodelist.size() ;i++)
			if(included[i] &&  graph->pp->type[i] == target_type)
			{
				n=n+1;
				neq = 0;
				m=graph->nodelist[i].size();
				if(m>0)
				{
					for(j=0;j<m;j++)
						if(target_type != graph->pp->type[graph->nodelist[i][j]-1]) neq=neq+1;
					value.at(0) = value.at(0) + (double) neq/(double) m;
				}
			}
		if(n>0) value.at(0) = value.at(0)/(double)n;
		if(fpar[1]>0) // ratio version (1-M)/ (lambda_t/lambda)
		{
			if(*dbg)printf("M=%1.3f -> ",value.at(0));
			double lambda=0.0,ala;
			for(i=0; i<*graph->pp->S;i++) lambda = lambda + graph->pp->lambdas[i];
			ala = graph->pp->lambdas[target_type-1] / lambda;
			value.at(0) = (1.0-value.at(0))/(double)ala;
			if(*dbg)printf("%1.3f",value.at(0));
		}
	}
	if(*dbg)printf("]");
	return value;
}

/**********************************************************************************/

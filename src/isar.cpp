#include "isar.h"
/**********************************************************************************/

std::vector<double> isar(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)printf("isar[ %f, ",*fpar);
	int i, j, n, k, m, dbg0;
	double a1[1];
	std::vector<double> value;
	value.clear();

	if((int)fpar[0]<=0)
	{
		dbg0 = *dbg;
		*dbg = 0;
		for(i=0;i< *graph->pp->S;i++)
		{
			if(graph->pp->lambdas[i]>0)
			{
				a1[0] = (double) i+1;
				value.push_back(isar(graph,a1,dbg,included).at(0));
			}
		}
//		I = I/((double) incS);
		*dbg = dbg0;
	}
	else
	{
		int target_type = (int) fpar[0];
		value.push_back(0.0);
		for(j=0; j< *graph->pp->S;j++) // sum over all species
		{

			n=0;
			m=0;
			for(i=0;i< (int)graph->nodelist.size();i++)
				if(included[i] && graph->pp->type[i]==target_type) //sum over target species
				{
					n++;
					for(k=0;k < (int)graph->nodelist[i].size();k++)  // check if type j present
					{
						if(graph->pp->type[graph->nodelist[i][k]-1]==j+1)
						{
							m++;
							break;
						}
					}
				}
			value.at(0)=value.at(0)+(double)m/(double)n;
		}
	}
	if(*dbg)printf("]");
	return value;
}

/**********************************************************************************/

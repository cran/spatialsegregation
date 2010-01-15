#include "simpson.h"

//std::vector<double> simpson0(Graph *graph, double *fpar, int *dbg, int *included)
//{
//	if(*dbg)printf("Simpson[");
//
//	int i;
//	double value0=0.0;
//	std::vector<double> value;
//
//	value = simpson_typewise(graph, fpar, dbg, included);
//
//	//and finally the estimate
//	for(i=0;i< *graph->pp->S;i++)
//		value0 = value0 + value.at(i);
//	value.clear();
//	value.push_back( 1.0 - value0 );
//
//	if(*dbg)printf(" ]");
//	return value;
//
//}

std::vector<double> simpson(Graph *graph, double *fpar, int *dbg, int *included)
{
	return simpson_typewise(graph, fpar, dbg, included);
}



std::vector<double> simpson_typewise(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)printf("typewise[");

	int n=graph->nodelist.size(),i,j,k,m,*S, N=0;
	S = graph->pp->S;
	double *piitautau = new double[*S], piitautaui, lambda0=0.0;
	double ptau, bardeg=0.0, *bardegtautau = new double[*S], bardegtautaui;
	std::vector<double> value;

	for(i=0;i<*S;i++)// compute overall intensity&format barpii_tau
	{
		lambda0= lambda0 + graph->pp->lambdas[i];
		piitautau[i]=0.0;
		bardegtautau[i]=0.0;
		m = 0;
		for(j=0;j<n;j++) // first the piitau for each tau node...
		{
		  if(included[j] && graph->pp->type[j]==i+1)// in our minus-sample and i of tau type
		  {
			  bardeg = bardeg + (double) graph->nodelist.at(j).size();
			  N = N+1;
			  piitautaui = 0.0;
			  bardegtautaui = 0.0;
			  m=m+1;//number of tau types
			  if(graph->nodelist.at(j).size()>0)
			  {
				  for(k=0;k<(int)graph->nodelist.at(j).size(); k++)
				  {
					  if(graph->pp->type[graph->nodelist.at(j).at(k)-1]==i+1)
					  {
						  piitautaui = piitautaui + 1;               // TODO: add the mass
						  bardegtautaui = bardegtautaui + 1;
					  }
				  }
				  piitautaui = piitautaui / (double)graph->nodelist.at(j).size();
			  }
			  bardegtautau[i] = bardegtautau[i] + bardegtautaui;
			  piitautau[i] = piitautau[i] + piitautaui; // sum for the arith. mean...
		  }
		}
		if(m>0)
		{
			bardegtautau[i] = bardegtautau[i]/(double)m;
			piitautau[i] = piitautau[i]/(double) m;//then the arith. mean over all sampled nodes
		}

	}
	bardeg = bardeg/(double)N;
//	value = 0.0;
	value.resize(0);
	//and finally the estimate
	for(i=0;i<*S;i++)
	{
		if(graph->pp->lambdas[i]>0)
		{
			ptau = graph->pp->lambdas[i]/lambda0;
			value.push_back( ptau*bardegtautau[i]/bardeg );
		}

		//value = value + barpiitau*piitautau[i]; //~pi_tau. Arith. mean, maybe use geometric?
	}
//	value.at(0) = 1 - value.at(0);

	if(*dbg)printf(" ]");
	return value;
}
//EOF

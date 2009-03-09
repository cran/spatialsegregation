#include "shannon.h"


std::vector<double> shannon0(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)printf("shannon[");
	int i,S=*graph->pp->S;
	double Eglobal=0.0, Elocal=0.0, ratio, lambda0=0, X;
	std::vector<double> value;


	for(i=0;i<S;i++)// compute overall intensity&format piitau
		lambda0= lambda0 + graph->pp->lambdas[i];

	for(i=0;i<S;i++)// compute the global entropy
	{
		ratio = graph->pp->lambdas[i]/lambda0;
		X = 0.0;
		if(ratio>0) X =  graph->pp->lambdas[i] * (log(ratio) / log(S) );
		Eglobal= Eglobal + X;
	}
	Eglobal = - Eglobal/lambda0;//=~E(N)

	value = piitauf(graph, fpar, dbg, included); // compute the pi_tau's, most demanding part

	for(i=0;i<S;i++)//~pi_tau
	{
		ratio = value.at(i);
		if(ratio>0)
			Elocal = Elocal + ratio* (log(ratio)/log(S));
	}
	//and finally the estimate
	value.clear();
	value.push_back(0.0);

	value.at(0) =  -Elocal /(double) Eglobal;

	if(*dbg)printf("]");
	return value;
}

std::vector<double> shannon(Graph *graph, double *fpar, int *dbg, int *included)
{
	return piitauf(graph, fpar, dbg, included);
}


std::vector<double> piitauf(Graph *graph, double *fpar, int *dbg, int *included)
{
	if(*dbg)printf("piitau[");
	int i,j,k,m,S=*graph->pp->S;
	double piitauj, piitau[S];
	std::vector<double> value;
	value.clear();

	for(i=0;i<S;i++)//compute the pii_tau's
	{
//		if(graph->pp->lambdas[i]>0) // not right, should not exclude any targets
		{
			m = 0;
			piitau[i]=0.0;

			for(j=0;j< *graph->pp->n;j++) // first the piitau for each node...
			 if(included[j]) // exclude the minus-sampling sources
			 {
				m=m+1;
				piitauj = 0.0;
				if(graph->nodelist[j].size()>0)
				{
					for(k=0;k<(int)graph->nodelist[j].size(); k++)
					{
						if(graph->pp->type[graph->nodelist[j][k]-1]==i+1)
							piitauj = piitauj + 1;
					}
					piitauj = piitauj / (double)graph->nodelist[j].size();
				}
				piitau[i] = piitau[i] + piitauj; // sum for the arith. mean...
			 }
			if(m>0) piitau[i] = piitau[i]/(double) m;//then the arith. mean over all sampled nodes
			value.push_back(piitau[i]);
		}
	}

	if(*dbg)printf("]");
	return value;
}
//EOF

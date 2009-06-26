#include "Graph.h"

#ifndef MAX_DOUBLE
const double MAX_DOUBLE = 9999999;
#endif

Graph::Graph()
{
}
/********************************************************************************************/
Graph::~Graph()
{
}
/********************************************************************************************/
double Graph::Dist(int *i, int *j)
{
    return this->pp->getDist(i, j);
}
/********************************************************************************************/
/********************************************************************************************/
void Graph::Init(Pp *pp0, int *gtype0, double *par0, double *prepR0, int *doDists0, int *toroidal0, int *dbg0)
{
	if(*dbg0)printf("intializing graph-object... ");

	pp = pp0;
	par=par0;
	prepR=prepR0;
	opar = *par;
	oldpar = &opar;
	doDists=doDists0;dbg=dbg0;
	toroidal= toroidal0;
	nodelist.resize(*pp->n);

	gtype = gtype0;

	if(toroidal)pp->toggleToroidal();
	if(*doDists)  // distance triangle
	{
		if(*dbg)printf("Precalculating distances...");
		pp->calcDists();
		if(*dbg)printf("ok. ");
	}
	pdone = 0;
	prepDone = &pdone;

	if(*dbg)printf(" done.\n");

}

/********************************************************************************************/
void Graph::setNodelist(std::vector<std::vector<int> > *nodelist_new)
{
	int i;
	nodelist.clear();nodelist.resize(0);
	for(i=0;i<(int)nodelist_new->size();i++)
		nodelist.push_back(nodelist_new->at(i));
}
/********************************************************************************************/
SEXP Graph::toSEXP()
//transform a std::vector<std::vector<int> > to SEXP, desctructive
{
	SEXP graph, *node;

	PROTECT(graph = allocVector(VECSXP, this->nodelist.size()));
	int i,j, *p, n;
	for(i=0;i< (int) this->nodelist.size();i++)
	{
		node = new SEXP;
		PROTECT(*node = allocVector(INTSXP, this->nodelist[i].size()) );
		p = INTEGER(*node);
		n = (int) this->nodelist[i].size();
		if(n<1) p[0]=NULL;
		else
			for(j=0;j<n;j++)
			{
				p[j] = (int) this->nodelist[i][j];
			};
		this->nodelist[i].clear();
		SET_VECTOR_ELT(graph, i, *node);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return graph;
}
/********************************************************************************************/
void Graph::remove_duplicates()
// remove duplicates in edge lists, especially the simple delauney method creates them
{
	int i,j,k,isnew;
	std::vector<int> *node;
	for(i=0; i < (int) this->nodelist.size() ;i++)
	{
		node = new std::vector<int>;
		node->resize(0);
		for(j=0; j < (int) this->nodelist.at(i).size() ; j++)
		{
			isnew=1;
			for(k=0;k<(int)node->size();k++)
			{
				if(node->at(k) == this->nodelist.at(i).at(j))
				{
					isnew=0;
					break;
				}
			}
			if(isnew)
				node->push_back(nodelist.at(i).at(j));
		}
		nodelist.at(i).swap(*node);
		delete node;
	}
}
/********************************************************************************************/
void Graph::addNew(int i, int j)
// add j to i
{
	int k,isnew=1;
	for(k=0; k< (int) nodelist.at(i).size();k++)
	{
		if(nodelist.at(i).at(k) == j)
		{
			isnew=0;
			break;
		}
	}
	if(isnew)
		nodelist.at(i).push_back(j);
}
/********************************************************************************************/
//The graph methods
/********************************************************************************************/
void Graph::sg_calc()
{
	// preprocess if requested
	if( *prepR > 0 & *prepDone == 0 )
	{
		if(*dbg)printf("Preprocessing[");
		this->sg_geometric(prepR);
		if(*dbg)printf("] ok.\n ");
		*prepDone = 1;
	}
	//start the calculation
	if(*gtype==0) //geometric
	{
		if(*oldpar > *par)
			this->sg_shrink_geometric(par);
		else
			this->sg_geometric();
	}
	else if(*gtype==1) //knn
	{
		if(*oldpar > *par)
			this->sg_shrink_knn();
		else
			this->sg_knn();
	}
	else if(*gtype==2) this->sg_mass_geometric();
	else if(*gtype==3) this->sg_gabriel();
	else if(*gtype==4) this->sg_delauney();
	else if(*gtype==5) this->sg_MST();
	else if(*gtype==6) this->sg_markcross();
	else if(*gtype==7) this->sg_SIG();
	else if(*gtype==8) this->sg_RST();
	else if(*gtype==9) this->sg_RNG();
	else if(*gtype==10) this->sg_CCC();
}

/********************************************************************************************/
void Graph::sg_geometric()
{
 Graph::sg_geometric(par);
}

void Graph::sg_geometric(double *R)
{
	if(*dbg)printf("Geometric (R=%f):",*R);
	int i,j;
	double dist;
	for(i=0;i<(*pp->n-1);i++)
		for(j=i+1;j<*pp->n;j++)
		{
			dist = Dist(&i,&j);
			if(dist<*R){
				nodelist[i].push_back(j+1);
				nodelist[j].push_back(i+1);
			}
		}
	if(*dbg)printf(" Ok.");
}

void Graph::sg_shrink_geometric(double *R)
{
	if(*dbg)printf("Geometric (R=%f) (shrinking):",*R);
	int i,j,j0;
	double dist;
	std::vector<int> *node;
	for(i=0; i < *pp->n ; i++)
	{
		node = new std::vector<int>;
		for(j=0;j < (int) this->nodelist[i].size() ; j++)
		{
			j0 = nodelist[i][j]-1;
			dist = Dist(&i,&j0);
			if(dist<*R)
				node->push_back(j0+1);
		}
		nodelist[i].clear();nodelist[i].resize(0);
		for (j = 0; j < (int)node->size(); ++j) this->nodelist[i].push_back(node->at(j));
		delete node;
	}
	if(*dbg)printf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_mass_geometric()
{
	if(*dbg)printf("Mass-geometric:");
	int i,j;
	double dist;
	for(i=0;i<(*pp->n-1);i++)
		for(j=i+1;j<*pp->n;j++)
		{
			dist = Dist(&i,&j);
			if(dist< pp->mass[i]){
				nodelist[i].push_back(j+1);
				nodelist[j].push_back(i+1);
			}
		}
	if(*dbg)printf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_knn()
{

	int i,j,l,*k,kk;
	kk = (int) par[0];
	k = &kk;
	std::vector<int> *node;
	if(*prepR==0)// if not preprocessed
	{
		if(*dbg)printf("%i-nn:",*k);
		double dists2_i[*pp->n], dists2_i2[*pp->n];
		for(i=0;i<*pp->n;i++) //for each point
		{
			for(j=0;j<*pp->n;j++) dists2_i2[j]=dists2_i[j]= Dist(&i,&j); //gather the distances to others
			qsort( dists2_i, *pp->n, sizeof(double),compare_doubles); // sort distances, rising
			for(j=1;j<=*k;j++) // find the k nearest
				for(l=0;l<*pp->n;l++)
					if( dists2_i[j] == dists2_i2[l] ) //with distance comparison
					{
						nodelist[i].push_back(l+1);
						break;
					}

		}
	}
	else{ //preprocessed
		if(*dbg)printf("%i-nn (shrinking):",*k);
		double *dists2_i, *dists2_i2;
		for(i=0;i<*pp->n;i++) //for each point
		{
			node = new std::vector<int>;
			dists2_i = new double [nodelist[i].size()];
			dists2_i2 = new double [nodelist[i].size()];
			if( (int) nodelist[i].size()<*k){ printf("\n preprocessing R too small, not enough neighbours (point #%i)!!\n",i+1);}
			else// in case there is more than k precalc neighbours
			{
				for(l =	0;l < (int) nodelist[i].size();l++)
				{
					j = nodelist[i][l]-1;
					dists2_i2[l]=Dist(&i,&j); //gather the distances to others, given preprocessing
					dists2_i[l]=dists2_i2[l];
				}
				qsort( dists2_i, nodelist[i].size() , sizeof(double),compare_doubles); // sort distances, rising
				for(j=0;j<*k;j++) // find the k nearest
					for(l=0;l< (int) nodelist[i].size();l++)
						if( dists2_i[j] == dists2_i2[l] ) //with distance comparison
						{
							node->push_back(nodelist[i][l]);
							break;
						}
				nodelist[i].clear();nodelist[i].resize(0);
				for(j=0;j < (int) node->size();j++) nodelist[i].push_back( (*node)[j] );
				delete node;
				delete[] dists2_i;
				delete[] dists2_i2;
			}
		}
	}
	 if(*dbg)printf(" Ok.");
}

void Graph::sg_shrink_knn()
{
	double *R0=prepR, R=1;
	prepR = &R;
	this->sg_knn();
	prepR = R0;
}

/********************************************************************************************/
void Graph::sg_gabriel()
{
	int kk = (int) par[0];
	if(*dbg & kk>0)printf("%i-",kk);
	if(*dbg)printf("Gabriel:");
	int i,j,k, empty,m,l,h;
	double x0,y0,R2, d;
	std::vector<int> *node;

	if(*this->prepR == 0) // no preprocessing done,a heavy looping
	  for(i=0;i<(*pp->n-1);i++)
	  {
		  for(j=i+1;j<*pp->n;j++)
		  {
			  x0 = fabs(pp->x[i]-pp->x[j])/2.0+fmin(pp->x[i],pp->x[j]);
			  y0 = fabs(pp->y[i]-pp->y[j])/2.0+fmin(pp->y[i],pp->y[j]);
			  R2 = ( pow(pp->x[i]-pp->x[j],2) + pow(pp->y[i]-pp->y[j],2) )/4.0;
			  //		brute force
			  empty = 1+kk;
			  for(k=0;k<*pp->n;k++)
			  {
				  if(k != i)
					  if( k != j)
					  {
						  d = pow(x0-pp->x[k],2) + pow(y0-pp->y[k],2);
						  if( d<R2 )
						  {
							  empty = empty - 1;
							  if(empty == 0) break;
						  }
					  }
			  }
			  if(empty)
			  {
				  this->nodelist[i].push_back(j+1);this->nodelist[j].push_back(i+1);
			  }
		  }
	  }
	else{ // preprocessed: nodelist has the restricted neighbourhoods to look trough
		if(*dbg)printf("(prepd): ");
		for(i = 0 ; i< *pp->n ;  i++)
		{
			node = new std::vector<int>;
			for( l=0 ; l < (int) this->nodelist[i].size(); l++ )
			{
				j = this->nodelist[i][l]-1;
				x0 = fabs(this->pp->x[i]-this->pp->x[j])/2.0+fmin(this->pp->x[i],this->pp->x[j]);
				y0 = fabs(this->pp->y[i]-this->pp->y[j])/2.0+fmin(this->pp->y[i],this->pp->y[j]);
				R2 = (pow(this->pp->x[i]-this->pp->x[j],2) + pow(this->pp->y[i]-this->pp->y[j],2) )/4.0;
				empty = 1+kk;
				for(m=0; m < (int)this->nodelist[i].size();m++) // the small ball is included in the preprocessing ball
				{
					k =  (int) this->nodelist[i][m]-1;
					if(k != i)
						if( k != j)
						{
							d = pow(x0-pp->x[k],2) + pow(y0-pp->y[k],2);
							if( d<R2 )
							{
								empty = empty - 1;
								if(empty == 0) break;
							}
						}
				}
				if(empty)
				{
					node->push_back(j+1);
				}
			}
			nodelist[i].clear();nodelist[i].resize(0);
			for(h=0;h<(int)node->size();h++) nodelist[i].push_back(node->at(h));
			delete node;
		}
	}
	  if(*dbg)printf(" Ok.");
}

/********************************************************************************************/
void Graph::sg_delauney()
{
//Naive algorithm, checks the interiors of triangle circumcircles.
//For 2D patterns

	if(*dbg)printf("Delauney: ");
	int i,j,k,l,h;
	double dummy[2];
	std::vector<int> *node;
	if(*this->prepR==0) // no preprocessing done, heavy looping
	{
		if(*dbg)printf("(raw):");
		for(i = 0 ; i< *pp->n-2 ; i++ )
			for(j = i+1 ; j < *pp->n-1 ; j++ )
				for(k = j+1 ; k < *pp->n ; k++ )
					if( Empty(pp->x,pp->y,pp->n,i,j,k, dummy, dummy, dummy) )
					{
						addNew(i,j+1); addNew(i,k+1);
						addNew(j,i+1); addNew(j,k+1);
						addNew(k,i+1); addNew(k,j+1);
					}
	}
	else{ // preprocessed: nodelist has the restricted neighbourhoods to look trough for triangles
		if(*dbg)printf("(prepd): ");
		for(i = 0 ; i< *pp->n ;  i++)
		{
			node = new std::vector<int>;
			for( l=0 ; l < (int)nodelist[i].size()-1; l++ )
			{
				j = nodelist[i][l]-1;
				for(h=l+1; h < (int)nodelist[i].size();h++ )
				{
					k = nodelist[i][h]-1;
					if( Empty(pp->x,pp->y,pp->n,i,j,k,dummy, dummy, dummy) )
					{
						node->push_back(j+1);
						node->push_back(k+1);
					}
				}
			}
			nodelist[i].clear();nodelist[i].resize(0);
			for(h=0;h<(int)node->size();h++) nodelist[i].push_back((*node).at(h));
			delete node;
		}
		this->remove_duplicates();
	}
	if(*dbg)printf(" Ok.");

}
/********************************************************************************************/
void Graph::sg_MST()
{
  if(*this->dbg) printf("MST:");
  int i,j,k=0,l=0,zz,k0=0,l0=0;
  int done[*this->pp->n],dn;
  double apu0,apu1,apu2;
  done[0] = 0;
  dn = 1;
  int left=*this->pp->n-dn;
//   double sum;
  while( left > 0 )
  {
    //if(*print_sg) if((left+1)%100==0) printf("MST: %i/%i         \r",left,*n);
    apu2 = MAX_DOUBLE;
    for(i=1; i<*this->pp->n;i++){
      zz = 1;
      apu0=apu2;
      for(j=0; j<dn;j++){
        if(i == done[j] ) {
          zz=0;
          break;
        }
        apu1 = Dist(&i,&done[j]); //dists[i*(*n)+done[j]];
        if( apu1<apu0 ){
          apu0=apu1;
          k0=i;
          l0=done[j];
        }
      }
      if(zz){
        if(apu0<apu2){
          apu2=apu0;
          k=k0;
          l=l0;
        }
      }
    }
//     sum+=apu2;
    done[dn] = k;
    dn++;
    left--;
    this->nodelist[l].push_back(k+1);
    //e[l*(*n)+k] = 1;
//     printf("\r");
  }
  if(*this->dbg)printf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_markcross()
{
	if(*dbg)printf("Markcross: ");
	int i,j;
	double dist;
	for(i=0;i<(*pp->n-1);i++)
		for(j=i+1;j<*pp->n;j++)
		{
			dist = Dist(&i,&j);
			if(dist< pp->mass[i]+pp->mass[j]){
				nodelist[i].push_back(j+1);
				nodelist[j].push_back(i+1);
			}
		}
	if(*dbg)printf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_SIG()
{
	if(*dbg)printf("Spheres-of-Influence:");
	int i,j,dbg0=*dbg;
	double dist;
	for(i=0;i<*pp->n;i++)
	{
		dist = MAX_DOUBLE;
		for(j=0;j<*pp->n;j++)
			if(i!=j) dist = fminf(dist, Dist(&i,&j));
		pp->mass[i]=dist;
	}
	*dbg=0;
	sg_markcross();
	*dbg=dbg0;
	if(*dbg)printf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_RST()
{
  if(*dbg) printf("Radial Spanning Tree (o=(%f,%f,%f)): ",par[0],par[1],par[2]);
  nodelist.resize(*pp->n-1);
  int i,j,k,foc_i=*pp->n-1;
  double apu0,apu1,apu2,apu3;

  for(i=0;i<*pp->n-1;i++)
  {
    apu0 = Dist(&i,&foc_i);//dists[i*(*n+1)+*n];
    apu3=MAX_DOUBLE;
    k=-1;
    for(j=0;j<*pp->n-1;j++)
    {
      if(j!=i)
      {
        apu1 = Dist(&j,&foc_i);//dists[j*(*n+1)+*n];
        if(apu1 < apu0 )
        {
          apu2 = Dist(&i,&j);//dists[i*(*n+1)+j];
          if( apu2 < apu3 )
          {
            apu3 = apu2;
            k = j;
          }
        }
      }
    }
    if(k>-1) addNew(k,i+1);//e[k*(*n)+i] = 1;
  }
  if(*dbg) printf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_RNG()
{
	if(*dbg) printf("Relative neighbourhood: ");
	int i,j,k,isempty;
    for(i=0;i<(*pp->n-1);i++)
    {
        for(j=i+1;j<*pp->n;j++)
        {
        	isempty = 1;
        	for(k=0;k<*pp->n;k++)
        		if(k!=i&k!=j)
        			if(Dist(&i,&k) < Dist(&i,&j))
        				if(Dist(&j,&k) < Dist(&j,&i))
        				{isempty=0;break;}
        	if(isempty)
        	{
        		addNew(i,j+1);
        		addNew(j,i+1);
        	}
        }
    }
    if(*dbg) printf(" Ok.");
}
/********************************************************************************************/
void Graph::sg_CCC()
{
	if(*dbg) printf("Class Cover Catch for type=%i: ",(int)par[0]);
	int i,j, type0=(int)par[0];
	for(i=0; i<*pp->n;i++)
	{
		pp->mass[i]=-MAX_DOUBLE;
		if(pp->type[i]==type0)
		{
			pp->mass[i]=MAX_DOUBLE;
			for(j=0;j<*pp->n;j++)
				if(j!=i & pp->type[j]!=type0) pp->mass[i]=fminf(pp->mass[i],Dist(&i,&j));
		}
	}
	for(i=0;i<*pp->n;i++) //TODO: optimize this
		if(pp->type[i]==type0)
			for(j=0;j<*pp->n;j++)
				if(i!=j)
					if(pp->type[j]==type0)
						if(Dist(&i,&j)< pp->mass[i])
							addNew(i,j+1);
	if(*dbg) printf(" Ok.");
}
/********************************************************************************************/
// cut (remove) the graph edges longer than R
void Graph::sg_cut(double *R)
{
	int i,j,k, count=0;
	if(*dbg)printf("Cutting the graph (R=%f):",*R);
	std::vector<int > *pnode;
	for(i=0;i < *pp->n;i++)
	{
		pnode = new std::vector<int>;
		pnode->resize(0);
		for(j=0; j < (int)nodelist.at(i).size();j++)
		{
			k = nodelist.at(i).at(j)-1;
			if( Dist(&i, &k) < *R )
				pnode->push_back(k+1);
			else
				count++;
		}
		nodelist.at(i).swap(*pnode);
		delete pnode;
	}
	if(*dbg)printf(" ok (%i edges cut). ",count);
}
/********************************************************************************************/
// prune branches less than lev hops long

void Graph::sg_prune(double *lev)
{
	int level = (int) *lev , i, leaf, count=0, prev, next;
	std::vector<int> left;
	std::vector<int> branch, *pnode;
	left.resize(0);
	branch.resize(0);
	if(*dbg)printf("Pruning the graph (level=%i):",level);

	for(i=0; i < (int)nodelist.size(); i++ ) // get the leaves
	{
		if( (int)nodelist.at(i).size() == 1 )
			left.push_back(i+1);
	}
	if(*dbg)printf("found %i leaves, pruning...",(int)left.size());

	while(!left.empty()) // go each branch trough starting from the leaf
	{
		leaf = left.back();
		branch.push_back(leaf);
		prev = leaf;
		next = nodelist.at(leaf-1).at(0);
		while((int) nodelist.at(next-1).size()==2)
		{
			branch.push_back(next);
			if(nodelist.at(next-1).at(0) != prev)
			{
				prev = next;
				next = nodelist.at(next-1).at(0);
			}
			else
			{
				prev = next;
				next = nodelist.at(next-1).at(1);
			}
//			for(i=0;i < branch.size();i++)
//				if(j == branch.at(i)){notnew=1; break;}
//			if(notnew)break;
		}
//		printf("leaf:%i, branch length: %i\n",left.back(),branch.size());
		if((int)branch.size() <= level) // if short enough branch, cut it.
		{
			pnode = new std::vector<int>;
			pnode->resize(0);
			for(i=0; i < (int)branch.size();i++)
				nodelist.at(branch.at(i)-1).clear();

			for(i=0; i < (int)nodelist.at(next-1).size();i++)
				if(nodelist.at(next-1).at(i) != prev)
					pnode->push_back(nodelist.at(next-1).at(i));
			nodelist.at(next-1).swap(*pnode);
			delete pnode;
			count++;
		}
		left.pop_back();
		branch.clear();
		branch.resize(0);
	}
	if(*dbg)printf(" Ok (%i branches pruned).",count);
}
/********************************************************************************************/
// EOF

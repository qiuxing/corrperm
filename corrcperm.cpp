//#include "Python.h"
#include "corrcperm.h"
#include "math.h"
#include <pthread.h>
#include <iostream>
//#include <iomanip>
#include <algorithm>
#include <vector>

using namespace std;
/* matrix functions*/

data2d::data2d(int ni, int nj) {
  int i,j;
  row = ni;
  column = nj;
  data  = new double *[row];
  for (i = 0; i < row; i++)
    data[i] = new double[column];

  for (i = 0; i < row; i++)
    for (j = 0; j < column; j++)
      data[i][j] = 0.0;
}

data2d::~data2d() {
  for (int i = 0; i < row; i++){
    delete data[i];
  }
  delete data;
}


inline double fishertransform(double x){
  return 0.5 * log((1.0+x)/(1.0-x));
}

double corrcperm::transform(double x){
  if (fisher) {
    return fishertransform(x);
  }else {
    return x;
  }
}

//ni is the number of genes
//nj is the number of columns in each sample (double that total)
//
//note:  client code is responsible for ensuring that permutvec is valid before all calls
corrcperm::corrcperm(int ni, int nj, int groupsize, bool fisherenabled){
  this->fisher = fisherenabled;

  genes = ni;
  cols = nj;
  permutvec = new int[nj*2];
  datamat = new data2d(ni, nj*2);

  nstats = new double[genes];

  ignoredgenes = new bool[genes];
  for (int i = 0; i < genes; i++){
    ignoredgenes[i] = false;
  }

  groupsumsallocated = false;

  numgroups = (int)ceil((double)cols/(double)groupsize);
  gsize = groupsize;
}

corrcperm::~corrcperm(){

  delete permutvec;
  delete datamat;

  if(groupsumsallocated){
    for (int i = 0; i < numgroups*2; i++){
      delete sums[i];
      delete sumsqs[i];
    }
    delete sums;
    delete sumsqs;  
  }

  delete ignoredgenes;

  delete nstats;
}

void rearrange(double * source, double * dest, int * permut, int length){
  /*used for permuting a single row of the data matrix */
  for (int i = 0; i < length; i++){
    dest[i] = source[permut[i]];
  }
}

void corrcperm::rearrangeall(){
  /*rearrange the data matrix according to the permutation vector*/
  double * newrow = new double[cols*2];
  for (int i = 0; i < genes; i++){
    rearrange(datamat->data[i], newrow, permutvec, cols*2);
    double * temp = datamat->data[i];
    datamat->data[i] = newrow;
    newrow = temp;
  }
  delete newrow;
}

void corrcperm::data_set(int i, int j, double val){
  datamat->data[i][j] = val;
}


void corrcperm::permutevecset(int n, int val){
  /*set the values of the permutation vector*/
  permutvec[n] = val;
}

void corrcperm::makepvec(int seed){
  srand(seed);
  std::vector<int> p;
  std::vector <int>::iterator it;

  int j = 0;

  for (int j = 0; j < cols*2; j++){
    p.push_back(j);
  }
  random_shuffle(p.begin(), p.end());

  int i = 0;
  for (it = p.begin(); it < p.end(); it++){
    permutvec[i] = *it;
    //cout << setw(3) << *it << " ";
    i++;
  }			
  //cout << endl;
}

void corrcperm::makerandomdata(int seed, double mu, double std){
  srand(seed);
  double x1, x2, w, rand1, rand2;
  for (int i = 0; i < genes; i++){
    for (int j = 0; j < cols*2; j++){
      do {
	rand1 = ((double)rand())/( (double)(RAND_MAX)+(double)(1) );
	rand2 = ((double)rand())/( (double)(RAND_MAX)+(double)(1) );
	x1 = 2.0 * rand1 - 1.0;
	x2 = 2.0 * rand1 - 1.0;
	w = x1 * x1 + x2 * x2;
      } while ( w >= 1.0);
      w = sqrt ((-2.0 * log(w))/w);
      datamat->data[i][j] = mu + std*x1*w;
      //      cout << datamat->data[i][j] << " ";
    }
    //cout << endl;
  }
}  

//splits the slides into multiple groups, and computes the a
//correlation vectors for gene i for each group
//the vectors are stored in vecmat
//vecmat is structured as vecmat[groupnum][gene]
void corrcperm::corvec(int gi, data2d * vecmat)
{
  double sx1,sx2,sxs1,sxs2,sxy1,sxy2,sy1,sy2,sys1,sys2,corTemp1,corTemp2;
/*sx1,sx2: sum of xi in group 1 or 2
  sxs1,sxs2: sum of xi^2 first, and then sqrt(n*(sum of xi^2)-sx^2)
  sxy1,sxy2: sum of xiyi
  sy1,sy2: sum of yi
  sys1,sys2: sum of yi^2 */
  int gj,gjtemp,itemp,islides,isub;
  //  isub = (int)ceil((double)cols/(double)numgroups);
  isub = gsize;
  double ** d = datamat->data;
  for (islides=0; islides < numgroups; ++islides)
  {
    corTemp1 = 0.0;
    corTemp2 = 0.0;

    sx1 = 0.0;
    sx2 = 0.0;
    sxs1 = 0.0;
    sxs2 = 0.0;
    gjtemp = 0;

    itemp = 0;
    while(itemp < isub && itemp+islides*isub < cols)
    { 
      //double val1 = d[gi][permutvec[itemp+islides*isub]];
      double val1 = d[gi][itemp+islides*isub];//***
      sx1 += val1;
      sxs1 += val1*val1;

      //double val2 = d[gi][permutvec[itemp+islides*isub+cols]];
      double val2 = d[gi][itemp+islides*isub+cols];//***
      sx2 += val2;
      sxs2 += val2*val2;

      itemp++;                    
    } 

    double xterm1 = sqrt(itemp*sxs1-sx1*sx1);
    double xterm2 = sqrt(itemp*sxs2-sx2*sx2);

    for (int gj = 0; gj < genes; gj++)
    {
      if (gj != gi)
      {
        sxy1 = 0.0;
        sy1 = 0.0;
	sys1 = 0.0;

        sy2 = 0.0;
        sys2 = 0.0;
	sxy2 = 0.0;

        int itemp2 = 0;
	while(itemp2 < isub && itemp2+islides*isub < cols)
	  {
	    //int c1 = permutvec[itemp2+islides*isub];
	    int c1 = itemp2+islides*isub;//***
	    double datagj1 = d[gj][c1];
	    sxy1 += d[gi][c1] * datagj1;
	    sy1 += datagj1;
	    sys1 += datagj1*datagj1;

	    //int c2 = permutvec[itemp2+islides*isub+cols];
	    int c2 = itemp2+islides*isub+cols;//***
	    double datagj2 = d[gj][c2];
	    sxy2 += d[gi][c2]*datagj2;
	    sy2 += datagj2;
	    sys2 += datagj2*datagj2;

	    itemp2++;                    
	  } 

        corTemp1 = (itemp2*sxy1-sx1*sy1)/(xterm1*sqrt(itemp2*sys1-sy1*sy1));
        corTemp2 = (itemp2*sxy2-sx2*sy2)/(xterm2*sqrt(itemp2*sys2-sy2*sy2));
       
	//apply the transformation
	vecmat->data[islides][gjtemp] = transform(corTemp1);
	vecmat->data[islides+numgroups][gjtemp] = transform(corTemp2);
       
        gjtemp++; 
      }
    }
  }
}

double corrcperm::Nstat(int igen, int kern)
{

  data2d * vecmat = new data2d(2*numgroups, genes);

  double Ltemp1,Ltemp2;
  double N1 = 0.0;
  double N2 = 0.0;
  double N;
  int i,j,k;

  double ** vd = vecmat->data;

  corvec(igen, vecmat);

  if (kern == 2)
  {
    for (i=0;i<numgroups;++i)
    {
      for (j=0;j<numgroups;++j)
      {
        Ltemp1 = 0.0;
        for (k=0;k<genes-1;++k)
        {
	  double val = vd[i][k] - vd[j+numgroups][k];
	  Ltemp1 += val*val;
        }
        Ltemp1 = sqrt(Ltemp1);
        N1 = N1 + Ltemp1;      
      }
      for (j=0;j<i;++j)
      {
        Ltemp1 = 0.0;
        Ltemp2 = 0.0;
        for (k=0;k<genes-1;++k)
        {
	  double val = vd[i][k] - vd[j][k];
	  Ltemp1 += val*val;
 
	  double val2 = vd[i+numgroups][k] - vd[j+numgroups][k];
	  Ltemp2 += val2*val2;
        }
        Ltemp1 = sqrt(Ltemp1); 
        Ltemp2 = sqrt(Ltemp2);    
        N2 = N2 + Ltemp1+ Ltemp2;      
      }
    }
  }
  else
  {
    for (i=0;i<numgroups;++i)
    {
      for (j=0;j<numgroups;++j)
      {
        Ltemp1 = 0.0;
        for (k=0;k<genes-1;++k)
        {
	  Ltemp1 += fabs(vd[i][k] - vd[j+numgroups][k]);
        }
        N1 = N1 + Ltemp1;      
      }
      for (j=0;j<i;++j)
      {
        Ltemp1 = 0.0;
        Ltemp2 = 0.0;
        for (k=0;k<genes-1;++k)
        {
	  Ltemp1 += fabs(vd[i][k] - vd[j][k]);       
	  Ltemp2 += fabs(vd[i+numgroups][k] - vd[j+numgroups][k]);
        }
        N2 = N2 + Ltemp1+ Ltemp2;      
      }
    }
  }
  N = N1-N2;
  N = 2 * N/(numgroups*numgroups);

  delete vecmat;
  return N;
}

/*calculategroupsums calculates various intermediate values used by other functions and stores them in sums and sumsqs*/
/* if these values are not calculated, the functions will not give the correct results*/
//gsize and numgroups should be correctly set before calling this function
void corrcperm::calculategroupsums(){

  //init sums and sumsqs
  //possible problem here if the number of groups changes.  Don't change the number of groups!
  if (!groupsumsallocated){
    sums = new double*[numgroups*2];
    sumsqs = new double*[numgroups*2];
    for (int i = 0; i < numgroups*2; i++){
      sums[i] = new double[genes];
      sumsqs[i] = new double[genes];
    }
    groupsumsallocated = true;
  }


  double ** data = datamat->data;

  for (int g = 0; g < numgroups; g++){

    for (int i = 0; i < genes; i++){
      int j = 0;
    double sum1 = 0.0;
    double sumsq1 = 0.0;
    double sum2 = 0.0;
    double sumsq2 = 0.0;
      while(j < gsize && j+g*gsize < cols){
	//double val1 = data[i][permutvec[j+g*gsize]];
	double val1 = data[i][j+g*gsize];
	sum1 += val1;
	sumsq1 += val1*val1;

	//double val2 = data[i][permutvec[j+g*gsize+cols]];
	double val2 = data[i][j+g*gsize+cols];
	sum2 += val2;
	sumsq2 += val2*val2;

	j++;
      }
      sums[g][i] = sum1;
      sumsqs[g][i] = sumsq1;
      sums[g+numgroups][i] = sum2;
      sumsqs[g+numgroups][i] = sumsq2;
    }
  }
}

/*calculates the correlation between genei and genej in the specified group*/
/* calculategroupsums should be called before this function is called*/
double corrcperm::correlation(int genei, int genej, int group){
  double ** data = datamat->data;
  //  int groupsize = (int)ceil((double)cols/(double)numgroups);
  int groupstart = (group<numgroups)?(group*gsize):((group-numgroups)*gsize+cols);
  int groupend;
  if(group == numgroups*2-1){
    groupend = cols*2;
  }
  else if (group == numgroups-1){
    groupend = cols;
  }else if (group < numgroups){
    groupend = (group+1)*gsize;
  }else{
    groupend = (group+1-numgroups)*gsize+cols;
    }
  double sum = 0.0;

  for (int i = groupstart; i < groupend; i++){
    //sum += data[genei][permutvec[i]]*data[genej][permutvec[i]];
    sum += data[genei][i]*data[genej][i];
  }
  double n = groupend-groupstart;
  double top = ((n*sum)-(sums[group][genei]*sums[group][genej]));
  double bottomleft  = sqrt(n*sumsqs[group][genei]-sums[group][genei]*sums[group][genei]); 
  double bottomright = sqrt(n*sumsqs[group][genej]-sums[group][genej]*sums[group][genej]);


  double result = top/(bottomleft*bottomright);
  return result;
}

/*equivalent to corvec*/
void corrcperm::corvec2(int gene, data2d * vecmat)
{

  for (int group = 0; group < numgroups*2; group++){
    for ( int i = 0; i < genes; i++){
       
      if ( i > gene){
        vecmat->data[group][i-1] = transform(correlation(gene, i, group));
      }else if (i < gene){
	vecmat->data[group][i] = transform(correlation(gene, i, group));
      }
    }
  }
}


/*equivalent to corvec, but slightly faster*/
/*note that calculategroupsums should be called before calling corvec3*/
//gsize and numgroups should be set correctly as well
void corrcperm::corvec3(int gene, data2d * vecmat){
  
  double ** data = datamat->data;
  for (int group = 0; group < numgroups; group++){
    int groupstart = group*gsize;
    int groupend;

    if (group == numgroups-1){
      groupend = cols;
    }else if (group < numgroups){
      groupend = (group+1)*gsize;
    }    

    double * sg1 = sums[group];
    double * sg2 = sums[group+numgroups];
    double * sqg1 = sumsqs[group];
    double * sqg2 = sumsqs[group+numgroups];
    
    double n = groupend-groupstart;
    double sx1 = sg1[gene];
    double sx2 = sg2[gene];
    double bottomleft1 = sqrt(n*sqg1[gene] - sx1*sx1);
    double bottomleft2 = sqrt(n*sqg2[gene] - sx2*sx2); 
   
    for (int genej = 0; genej < genes; genej++){
      double sum1 = 0.0;
      double sum2 = 0.0;
      

      for (int i = groupstart; i < groupend; i++){
	//int index1 = permutvec[i];
	int index1 = i;
	sum1 += data[gene][index1]*data[genej][index1];
	//int index2 = permutvec[i+cols];
	int index2 = i+cols;
	sum2 += data[gene][index2]*data[genej][index2];
      }

      double sy1 = sg1[genej];
      double sy2 = sg2[genej];
      double top1 = ((n*sum1)-(sx1*sy1));
      double top2 = ((n*sum2)-(sx2*sy2));
      double bottomright1 = sqrt(n*sqg1[genej]-sy1*sy1);
      double bottomright2 = sqrt(n*sqg2[genej]-sy2*sy2);
      
      double result1 = transform(top1/(bottomleft1*bottomright1));
      double result2 = transform(top2/(bottomleft2*bottomright2));

      if (genej > gene){
        vecmat->data[group][genej-1] = result1;
	vecmat->data[group+numgroups][genej-1] = result2;
      }else if (genej < gene){
	vecmat->data[group][genej] = result1;
	vecmat->data[group+numgroups][genej] = result2;
      }
    }
  }
}

/*calculates the nstatistic using corvec3*/
/*calculategroupsums should be called before calling this function*/
double corrcperm::Nstat3(int gene, int kern)
{
  data2d * vecmat = new data2d(2*numgroups, genes);

  double Ltemp1,Ltemp2;
  double N1 = 0.0;
  double N2 = 0.0;
  double N;
  int i,j,k;
  double ** vd = vecmat->data;

  corvec3(gene, vecmat);
//after calling corvec, vecmat now stores the correlation vectors between "gene" and all other genes in each group

  if (kern == 2)//the kernel function is f(n) = n^2
  {
    for (i=0;i<numgroups;++i)
    {
      for (j=0;j<numgroups;++j)
      {
        Ltemp1 = 0.0;
        for (k=0;k<genes-1;++k)
        {
	  double val = vd[i][k] - vd[j+numgroups][k];
	  Ltemp1 += val*val;
        }
        Ltemp1 = sqrt(Ltemp1);
        N1 = N1 + Ltemp1;      
      }
      for (j=0;j<i;++j)
      {
        Ltemp1 = 0.0;
        Ltemp2 = 0.0;
        for (k=0;k<genes-1;++k)
        {
	  double val = vd[i][k] - vd[j][k];
	  Ltemp1 += val*val;
 
	  double val2 = vd[i+numgroups][k] - vd[j+numgroups][k];
	  Ltemp2 += val2*val2;
        }
        Ltemp1 = sqrt(Ltemp1); 
        Ltemp2 = sqrt(Ltemp2);    
        N2 = N2 + Ltemp1+ Ltemp2;      
      }
    }
  }
  else//use the kernel function f(n) = abs(n)
  {
    for (i=0;i<numgroups;++i)
    {
      for (j=0;j<numgroups;++j)
      {
        Ltemp1 = 0.0;
        for (k=0;k<genes-1;++k)
        {
	  Ltemp1 += fabs(vd[i][k] - vd[j+numgroups][k]);
        }
        N1 = N1 + Ltemp1;      
      }
      for (j=0;j<i;++j)
      {
        Ltemp1 = 0.0;
        Ltemp2 = 0.0;
        for (k=0;k<genes-1;++k)
        {
	  Ltemp1 += fabs(vd[i][k] - vd[j][k]);       
	  Ltemp2 += fabs(vd[i+numgroups][k] - vd[j+numgroups][k]);
        }
        N2 = N2 + Ltemp1+ Ltemp2;      
      }
    }
  }
  N = N1-N2;
  N = 2 * N/(numgroups*numgroups);

  delete vecmat;
  return N;
}


void corrcperm::allNstats(int kern){
  /*calculate all the nstatistics and store them in nstats*/

  calculategroupsums();

  for (int i = 0; i < genes; i++){
    nstats[i] = Nstat3(i, kern);
  }
}

double corrcperm::getNstat(int gene){
  /*return a single value from the stored nstatistics*/
  /*needed for the swig interface*/
  return nstats[gene];
}

/*a structure passed to each thread containing its arguments*/
struct thread_data{
  int thread_id;
  int start;
  int end;
  int kern;
  corrcperm * c;
};

/*this is the function each thread runs to compute the nstatistics*/
/*it computes nstatistics from the specified start to end (not including end)*/
void *threadfun(void * arg){
  struct thread_data * thread_args = (struct thread_data *)arg;
  
  int thread_id = thread_args->thread_id;
  int start = thread_args->start;
  int end = thread_args->end;
  int kern = thread_args->kern;
  corrcperm * c = thread_args->c; 

  //  cout << thread_id << " " << start << " " << end << endl;
  for (int i = start; i < end; i++){
    if (!c->ignoredgenes[i]){
      c->nstats[i] = c->Nstat3(i, kern);
    }else{
      c->nstats[i] = 0.0;//not quite sure if this is mathematically valid...
    }
  }
}

/*the same as allNstats() only multithreaded*/
void corrcperm::threadedallNstats(int kern, int numthreads){
  //numgroups = groups;
  //  gsize = (int)ceil((double)cols/(double)numgroups);

  calculategroupsums();

  pthread_t threads[numthreads];
  struct thread_data threadargs[numthreads];//contains the arguments to each thread
  
  int rc;
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);//make sure to create the threads as joinable
  
  //create the threads
  for (int t = 0; t < numthreads; t++){
    //set up each threads arguments
    threadargs[t].thread_id = t;
    threadargs[t].start = (genes * t)/numthreads;
    threadargs[t].end = (genes * (t + 1))/numthreads;
    threadargs[t].kern = kern;
    threadargs[t].c = this;
    
    //create the threads
    rc = pthread_create(&threads[t], &attr, threadfun, (void *)&threadargs[t]);
    if(rc){
      cout << "In main: creating thread " << t << endl;
      exit(-1);
    }
  }
  
  //wait for the threads to finish
  void * status;
  pthread_attr_destroy(&attr);
  for(int t=0; t< numthreads; t++){
    rc = pthread_join(threads[t], &status);
    if (rc)
      {
	cout << "ERROR; return code from pthread_join() is " << rc << endl;
	exit(-1);
      }
    //cout << "Completed join with thread " << t << " status = " << (long)status << endl;
  }
}

void corrcperm::ignoregene(int gene){
  ignoredgenes[gene] = true;
}



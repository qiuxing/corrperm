
#ifndef DATA2D
#define DATA2D
struct data2d{
data2d(int ni, int nj);
~data2d();
double **data;
  int row;
  int column;
};
#endif

#ifndef CORRCPERM
#define CORRCPERM
class corrcperm
{
private:
  void corvec(int gi, data2d * vecmat);
  void corvec2(int gene, data2d * vecmat);
  //  void corvec3(int gene, int numgroups, data2d * vecmat);
  void corvec3(int gene, data2d * vecmat);

  double transform(double x);

  double ** sums;
  double ** sumsqs;
  int numgroups;
  int gsize;
  bool groupsumsallocated;



public:
    int genes, cols;
    bool fisher;  //fisher transform enabled
    //    corrcperm(int ni, int nj, bool fisherenabled = true);
    corrcperm(int ni, int nj, int groupsize, bool fisherenabled = true);


    ~corrcperm();
    data2d *datamat;

    void data_set(int i, int j, double val);
    

    int * permutvec; 


    void makerandomdata(int seed, double mu, double std);
    
    void permutevecset(int n, int val);
    void makepvec(int seed);
 
    double Nstat(int igen, int kern);

    void calculategroupsums();
    double correlation(int genei, int genej, int group);
    double Nstat3(int gene, int kern);


    void allNstats(int kern);
    double getNstat(int gene);


    double * nstats;

    void rearrangeall();

    void threadedallNstats(int kern, int numthreads);

    void ignoregene(int gene);
    bool * ignoredgenes;
};
#endif
 

%module corrcperm
%{
#include "corrcperm.h"
%}

//struct data2d{
//data2d(int ni, int nj);
//~data2d();
//double **data;
//  int row;
//  int column;
//};

//%include corrcperm.h
class corrcperm
{
public:
    %immutable;
    int genes, cols;

    %nothread;
    corrcperm(int ni, int nj, int groupsize, bool fisher);
    ~corrcperm();
    void data_set(int i, int j, double val);

    void makerandomdata(int seed, double mu, double std);

    void permutevecset(int n, int val);

    void makepvec(int seed);

    double Nstat(int igen, int kern);	

    void calculategroupsums();
//    double correlation(int genei, int genej, int group, bool fisher=true);
    double Nstat3(int gene, int kern);

    void allNstats(int kern);
    double getNstat(int gene);

  

    void rearrangeall();

    %thread;	
    void threadedallNstats(int kern, int numthreads);

    void ignoregene(int gene);
};


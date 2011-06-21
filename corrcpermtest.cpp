#include "corrcperm.h"
#include "math.h"
#include <iostream>
#include <iomanip>
#include <time.h>

using namespace std;

void printmat(data2d *datamat){
  int rows = datamat->row;
  int cols = datamat->column;
  cout << rows << " rows by " << cols << " columns\n";
  //cout << datamat->data << "\n";
  cout << setiosflags(ios::left);
  for (int i = 0; i < rows; i++){
    for (int j = 0; j < cols; j++){
      cout << setw(5) << datamat->data[i][j] << " ";
    }
    //    cout << datamat->data[i];
    cout << endl;
  }
  cout << resetiosflags(ios::left);
}

void printmatpermed(corrcperm &c){
  data2d * datamat = c.datamat;
  int rows = datamat->row;
  int cols = datamat->column;
  cout << rows << " rows by " << cols << " columns\n";
  //cout << datamat->data << "\n";
  cout << setiosflags(ios::left);
  for (int i = 0; i < rows; i++){
    for (int j = 0; j < cols; j++){
      cout << setw(5) << datamat->data[i][c.permutvec[j]] << " ";
    }
    //    cout << datamat->data[i];
    cout << endl;
  }
  cout << resetiosflags(ios::left);
}

void printtrianglemat(int size, float ** mat){
  for (int i = 0; i < size; i++){
    for (int j = 0; j <= i; j++){
      cout << mat[i][j] << " ";
    }
    cout << "\n";
  }
}



  



#define TIME(f)					\
  start = time(NULL);				\
  f;						\
  end = time(NULL);				\
  cout << end-start << " seconds" << endl



/*void checkcorr(int seed, int rows, int cols, bool fisher){
  corrcperm c = corrcperm(rows, cols);
  c.makerandomdata(seed, 0.0, 2.0);
  c.makepvec(seed);
 
  time_t start;
  time_t end;

  TIME(c.corrcoefthread(1, fisher));

  double diftotal = 0.0;
  double maxdif = 0.0;
  TIME( 
       c.calculategroupsums(1);
 
       for (int i = 0; i < rows; i++){
	 for (int j = 0; j < i; j++){
	   double result1 = c.correlation(i, j, 0, fisher);
	   double result2 = c.correlation(i, j, 1, fisher);
	   double difference1 = fabs(result1 - c.correlations1[i][j]);
	   double difference2 = fabs(result2 - c.correlations2[i][j]);

	   diftotal += difference1 + difference2;
	   maxdif = (maxdif > difference1)?maxdif:difference1;
	   maxdif = (maxdif > difference2)?maxdif:difference2;
	   if (difference1 > 0.0000001){
	     cout << "result1 different " << i << " " << j << " " << result1 << " " << c.correlations1[i][j] << endl;
	   }
	   if (difference2 > 0.0000001){
	     cout << "result2 different " << i << " " << j << " " << result1 << " " << c.correlations2[i][j] << endl;
	   }
	 }
       });
  cout << "total difference " << diftotal << endl;
  cout << "maximum difference " << maxdif << endl;
  }*/

int main(){
  //  int rows = 22200;
  int rows = 80;
  int cols = 800;
  int groupsize = 100;
  int kern = 2;
  corrcperm c = corrcperm(rows, cols, groupsize);
  //  bool fisher = true;


  int seed = 1234; 
  time_t start;
  time_t end;


  //checkcorr(seed, rows,cols, true);
  //checkcorr(seed, rows, cols, false);




    c.makerandomdata(seed, 0.0, 2.0);
    c.makepvec(seed);
    //printmat(c.datamat);
    //printmatpermed(c);
    int numthreads = 4;
    int numperms = 1;

    
    /*TIME(
	 for (int i =0; i < numperms; i++){
	   c.makegroupsortedpvec(i, groups);
	   c.threadedallNstats(groups, kern, numthreads);
	 }
	 );*/
    c.threadedallNstats(kern, numthreads);
    for (int i = 0; i < rows; i++){
      cout << c.getNstat(i) << " ";
    }
    cout << endl;

    cout << endl;
    c.rearrangeall();
    for (int i = 0; i < cols*2; i++){
      c.permutevecset(i, i);
    }
    c.threadedallNstats(kern, numthreads);
    for (int i = 0; i < rows; i++){
      cout << c.getNstat(i) << " ";
    }
    cout << endl;
    //c.makepvec(seed);
    //c.calculategroupsums(groups);
    /*   for (int i = 0; i < rows; i++){
      double difference = c.Nstat3(i, 1) - c.Nstat(i, groups, 1);
      if (difference > 0.0000000001){
	cout << c.Nstat3(i, 1) << " " << c.Nstat(i, groups, 1) << endl;
	}
    }
    
    for (int i = 0; i < rows; i++){
      double difference = c.Nstat3(i, 2) - c.Nstat(i, groups, 2);
      if (difference > 0.0000000001){
	cout << c.Nstat3(i, 2) - c.Nstat(i, groups, 2) << endl;
      }
      }
    */

  
    /*TIME(
	 for (int i = 0; i < rows; i++){
	   c.Nstat(i, groups, kern);
	 }
	 );

    TIME(c.calculategroupsums(groups);
	 for (int i = 0; i < rows; i++){
	   c.Nstat3(i, kern);
	 }
	 );
    */

    

    /*
   c.makegroupsortedpvec(seed, groups);
     TIME(
	 for (int i = 0; i < rows; i++){
	   c.Nstat(i, groups, kern);
	 }
	 );
    */
}



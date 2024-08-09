/* File: main.cc
 * --------------
 * Calculate Fidelity of d-wave pair field for the presence of V
 *
 */ 


#include <iostream>
#include <fstream>
#include <bitset>
#include <vector>
#include <algorithm>
#include <cmath>
#include "para.h"
#include "cluster.h"
#include "hamiltonian.h"
#include "arpack.h"

#define NNZPCOL 60


using namespace std;

int choose(int n, int m) {
   if(n<=m) return 1;
   int temp=1;
   int number=n;
   if(m>n/2){
      while(number > m) temp *= number--;	
      number = n-m;
      while(number) temp = temp/number--; 
   } else {
      while(number > n-m) temp *= number--;	
      while(m) temp = temp/m--; 
   }
   return temp;
}

int main(int argc, char *argv[]){

//read in parameters
   Para para;
   ifstream infile("input");
   if(!infile){
    //  cerr << "error: enable to open input file:"
    //       << infile << endl;
      return -1;
   } else {
      infile >> para.nup >> para.ndn  
          >> para.U >> para.t;  
      cout << "   Parameters:" << endl;
      cout << "\tCluster size:\t" << N << endl;
      cout << "\t# of spin up(down) electrons:\t" << para.nup << " , " << para.ndn << endl;     
      cout << "\ton-site Coulomb repulsion U:\t" <<  para.U  << endl;
      cout << "\tnearest-neighbor hopping t:\t" << para.t  << endl;     
   }

   //build up Hilbert space states
   int Hsize = choose(N, para.nup)*choose(N, para.ndn);
   vector< bitset<2*N> > states;
   states.reserve(Hsize);
   genState(para, states);
   cout << "   Generating Hilbert space finished!" << endl;
   cout << "\tHsize = "<< states.size() << endl;

   vector<HEle> hami;
   hami.reserve(NNZPCOL*Hsize);
   hami.clear();
   genHami(Hsize, states, para, &hami);
   cout << "   Generating Hamiltonian matrix done!" << endl;
   cout << "\tSize of non-zero elements = " << hami.size() << endl;

   //int nev = NEV;
   //double eval[nev], *evec1 = new double[Hsize*nev];
   //EDarpack(Hsize, &hami, eval, evec1, nev);
   //cout << "print the two leading eigenvalues: " << eval[0] << " " << eval[1] << endl;

   double *evec1 = new double[Hsize*Hsize];
   double *eval = new double[Hsize];
   convertMatrix(Hsize, &hami, evec1);
   for(int i=0; i<Hsize; i++){
   //   for(int j=0; j<Hsize; j++){
   	 cout << i << " " << 0 << " " << evec1[i*Hsize+0] << endl;
   //   }
   }
   EDlapack(Hsize, evec1, eval);
   cout << "print the two leading elements: " << eval[Hsize-1] << " " << eval[Hsize-2] << endl;

   return 0;
}

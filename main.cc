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
      cerr << "error: enable to open input file:"
           << infile << endl;
      return -1;
   } else {
      infile >> para.nup >> para.ndn  
          >> para.U >> para.V  
          >> para.t >> para.tp >> para.tpp >> para.alpha
          >> para.UsStart >> para.UsEnd  
          >> para.UsStep >> para.lambda; 
      cout << "   Parameters:" ;
      cout << endl;
      cout << "\tCluster size:\t" << N << endl;
      cout << "\t# of spin up(down) electrons:\t" << para.nup << " , " << para.ndn << endl;     
      cout << "\ton-site Coulomb repulsion U:\t" <<  para.U  << endl;
      cout << "\tnearest-neighbor Coulomb repulsion V:\t" << para.V << endl;
      cout << "\tnearest-neighbor hopping t, t', t'':\t" << para.t  << " , " << para.tp << " , " << para.tpp << endl;     
      cout << "\talpha:\t" << para.alpha << endl;     
      cout << "\tScan of s-wave pair field strength Us:" << endl;
      cout << "\t   start:\t" << para.UsStart << "\tend:\t"<< para.UsEnd << endl;     
      cout << "\t   step length:\t" << para.UsStep << "\tD-wave delta Us:\t"<< para.lambda << endl;     
   }

   //build up Hilbert space states
   int Hsize = choose(N, para.nup)*choose(N, para.ndn);
   vector< bitset<2*N> > states;
   states.reserve(Hsize);
   genState(para, states);
   cout << "   Generating Hilbert space finished!" << endl;
   cout << "\tHsize = "<< states.size() << endl;

   //generate nloc, startH and endH
   ofstream outfile("Extended_s-wave_ABD.out"); 
   ofstream outfile2("eigenvalue.out"); 
   ofstream outfile3("Nq_Sq.out"); 

   vector<HEle> hami;
   hami.reserve(NNZPCOL*Hsize);
   vector<HEle> Oextswave;
   Oextswave.reserve(NNZPCOL*3*Hsize);
   genExtendedswave(Hsize, states, para, &Oextswave);
   cout << "   Generating extended s-wave matrix done!" << endl;
   cout << "\tSize of non-zero elements = " << Oextswave.size() << endl;
   vector<HEle> BdaggerB;
   BdaggerB.reserve(NNZPCOL/2*Hsize);
   genBdaggerB(Hsize, states, para, &BdaggerB);
   cout << "   Generating BdaggerB matrix done!" << endl;
   cout << "\tSize of non-zero elements = " << BdaggerB.size() << endl;
   vector<HEle> DdaggerD;
   DdaggerD.reserve(NNZPCOL*3*Hsize);
   genDdaggerD(Hsize, states, para, &DdaggerD);
   cout << "   Generating DdaggerD matrix done!" << endl;
   cout << "\tSize of non-zero elements = " << DdaggerD.size() << endl;

   vector<HEle> MatrixSq;
   MatrixSq.reserve(Hsize);
   genSq(Hsize, states, para, &MatrixSq);
   cout << "   Generating Sq matrix done!" << endl;

   vector<HEle> MatrixNq;
   MatrixNq.reserve(Hsize);
   genNq(Hsize, states, para, &MatrixNq);
   cout << "   Generating Nq matrix done!" << endl;

   for(double Us=para.UsStart; Us<=para.UsEnd; Us+=para.UsStep){
      cout << endl << "************************" << endl;

      // calculate groundstate at Us
      para.Us=Us;
      if(para.Us < 0.0) para.Us = para.Us/N;
      cout << " Calculating Us: " << para.Us << endl;

      hami.clear();
      genHami(Hsize, states, para, &hami);
      cout << "   Generating Hamiltonian matrix done!" << endl;
      cout << "\tSize of non-zero elements = " << hami.size() << endl;

      int nev = NEV;
      double eval[nev], *evec1 = new double[Hsize*nev];
      EDarpack(Hsize, &hami, eval, evec1, nev);

      outfile2 << para.Us << "\t" << eval[0] << endl;

      // Calculate v' * O_extended_s-wave * v
      double tsumAA = 0.0, *evect = new double[Hsize], tsumDD = 0.0, tsumBB = 0.0;
      double tsumSq = 0.0, tsumNq = 0.0;

      // Calculate A^\dagger A
      for(int i=0; i<Hsize; i++){
          evect[i] = 0.0;
      } 
      for(vector<HEle>::iterator Hrow_itr=Oextswave.begin();
          Hrow_itr!=Oextswave.end(); Hrow_itr++) {
          int mycol=Hrow_itr->getx();
          int myrow=Hrow_itr->gety();
          double myVal=Hrow_itr->getVal();
          evect[mycol]+=myVal*evec1[myrow];
      }
      for(int i=0; i<Hsize; i++){
          tsumAA += evect[i]*evec1[i];
      } 

      // Calculate B^\dagger B
      for(int i=0; i<Hsize; i++){
          evect[i] = 0.0;
      } 
      for(vector<HEle>::iterator Hrow_itr=BdaggerB.begin();
          Hrow_itr!=BdaggerB.end(); Hrow_itr++) {
          int mycol=Hrow_itr->getx();
          int myrow=Hrow_itr->gety();
          double myVal=Hrow_itr->getVal();
          evect[mycol]+=myVal*evec1[myrow];
      }
      for(int i=0; i<Hsize; i++){
          tsumBB += evect[i]*evec1[i];
      }
 
      // Calculate D^\dagger D
      for(int i=0; i<Hsize; i++){
          evect[i] = 0.0;
      } 
      for(vector<HEle>::iterator Hrow_itr=DdaggerD.begin();
          Hrow_itr!=DdaggerD.end(); Hrow_itr++) {
          int mycol=Hrow_itr->getx();
          int myrow=Hrow_itr->gety();
          double myVal=Hrow_itr->getVal();
          evect[mycol]+=myVal*evec1[myrow];
      }
      for(int i=0; i<Hsize; i++){
          tsumDD += evect[i]*evec1[i];
      }

      // Calculate Nq
      for(int i=0; i<Hsize; i++){
          evect[i] = 0.0;
      } 
      for(vector<HEle>::iterator Hrow_itr=MatrixNq.begin();
          Hrow_itr!=MatrixNq.end(); Hrow_itr++) {
          int mycol=Hrow_itr->getx();
          int myrow=Hrow_itr->gety();
          double myVal=Hrow_itr->getVal();
          evect[mycol]+=myVal*evec1[myrow];
      }
      for(int i=0; i<Hsize; i++){
          tsumNq += evect[i]*evec1[i];
      }

      // Calculate Sq
      for(int i=0; i<Hsize; i++){
          evect[i] = 0.0;
      } 
      for(vector<HEle>::iterator Hrow_itr=MatrixSq.begin();
          Hrow_itr!=MatrixSq.end(); Hrow_itr++) {
          int mycol=Hrow_itr->getx();
          int myrow=Hrow_itr->gety();
          double myVal=Hrow_itr->getVal();
          evect[mycol]+=myVal*evec1[myrow];
      }
      for(int i=0; i<Hsize; i++){
          tsumSq += evect[i]*evec1[i];
      }

 
      outfile << para.Us << "\t" << tsumAA/N << "\t" << tsumBB/N << "\t" << tsumDD/N << endl;

 
      // calculate d-wave pair density matrix 
      int Nsquare = N*N;
      double *Dmat = new double[Nsquare], *Bmat = new double[Nsquare], *Amat = new double[Nsquare]; 
      double *wD = new double[N], *wB = new double[N], *wA = new double[N]; 

      genDmatrix(Hsize, states, para, evec1, Dmat);
      // genDdaggerD(Hsize, states, para, &DdaggerD); 
      EDlapack(N, Dmat, wD);
      cout << "print the two leading elements: " << wD[N-1] << " " << wD[N-2] << endl;

      genBmatrix(Hsize, states, para, evec1, Bmat);
      EDlapack(N, Bmat, wB);
      cout << "print the two leading elements: " << wB[N-1] << " " << wB[N-2] << endl;

      outfile3 << para.Us << "\t" << tsumNq << "\t" << tsumSq << "\t" << wD[N-1]/wD[N-2] << endl;
   
   }
   outfile.close();
   outfile2.close();
   outfile3.close();
   return 0;
}



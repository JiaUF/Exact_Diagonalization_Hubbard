/* File: hamiltonian.h
 * ----------------
 * Define the class HEle, which saves one nnz with its col, row and value. 
 * The function genHami prototype is defined here 
 */

#ifndef _H_hmlt
#define _H_hmlt

#include <iostream>
#include "para.h"
#include <vector>
#include <bitset>

using namespace std;

class HEle {
   int indx, indy;
   double val;

public:
   HEle(int i, int j, double v){
      indx=i; indy=j; val=v;
   }
   HEle(): indx(0), indy(0), val(0.0) {}   

   int getx(){ return indx;}
   int gety(){ return indy;}
   double getVal(){ return val;}
};

void genState(Para, vector<bitset<2*N> > &);
void genHami(int, vector<bitset<2*N> >, Para, vector<HEle> *);
void genExtendedswave(int, vector<bitset<2*N> >, Para, vector<HEle> *);
void genBdaggerB(int, vector<bitset<2*N> >, Para, vector<HEle> *);
void genDdaggerD(int, vector<bitset<2*N> >, Para, vector<HEle> *);
void genSq(int, vector<bitset<2*N> >, Para, vector<HEle> *);
void genNq(int, vector<bitset<2*N> >, Para, vector<HEle> *);
void genDmatrix(int, vector<bitset<2*N> >, Para, double *, double *);
void genBmatrix(int, vector<bitset<2*N> >, Para, double *, double *);
vector<bitset<2*N> >::const_iterator dichotomy_search(vector<bitset<2*N> >::const_iterator, vector<bitset<2*N> >::const_iterator, bitset<2*N>);
bool operator<(bitset<2*N>, bitset<2*N> &);

#endif

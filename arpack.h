/* File: parpack.h
 * ----------------
 * Define a set of important numbers and the EDparpack function prototype 
 */

#ifndef _H_arpack
#define _H_arpack

#include "hamiltonian.h"

#define NEV 20

using namespace std;

void EDlapack(const int, double *, double *);

void EDarpack(int, vector<HEle> *, 
               double *, double *, int);


#endif

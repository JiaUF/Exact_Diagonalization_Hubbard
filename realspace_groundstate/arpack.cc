/* File: parpack.cc
 * -----------------
 * Implementation of EDparpack function
 */


#include <iostream>
#include <fstream>
#include <bitset>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include "hamiltonian.h"
#include "arpack.h"

#define MAXITERATION 1000
#define INFINITISIMAL 0.0001  

using namespace std;

/* DSYEV prototype */
extern "C" {
    void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );

    void dsaupd_( int *ido, const char *bmat, int *n, const char *which,
                  int *nev, double *tol, double *resid,
                  int *ncv, double *V, int *ldv,
                  int *iparam, int *ipntr, double *workd,
                  double *workl, int *lworkl, int *info);

    void dseupd_( int *rvec, char *HowMny, int *select,
                  double *d, double *Z, int *ldz,
                  double *sigma,  const char *bmat, int *n,
                  const char *which, int *nev, double *tol,
                  double *resid, int *ncv, double *V,
                  int *ldv, int *iparam, int *ipntr,
                  double *workd, double *workl,
                  int *lworkl, int *info);
}



inline void AlgVectorInit(double *v, int size) {for(int i=0; i<size; i++) v[i]=0;}
inline void printVector(double *v, int size) {for(int i=0; i<size; i++) cout << v[i] << endl;}       
void convertMatrix(int thsize, vector<HEle> * phami, double * evec1){

   for(vector<HEle>::iterator Hrow_itr=phami->begin();
       Hrow_itr!=phami->end(); Hrow_itr++) {
       int mycol=Hrow_itr->getx();
       int myrow=Hrow_itr->gety();
       double myVal=Hrow_itr->getVal();
       //if(count==1) cout << "\t  " << mycol << " " << myrow << " " << myVal << endl;
       evec1[mycol*thsize + myrow]+=myVal;
   }//for temp_pt (j)

}


void EDlapack(int hhsize, double *a, double *w){

   int lda2 = hhsize;
   int lda = hhsize, info, lwork;
   double wkopt;
   /* Executable statements */
   cout << endl << " DSYEV Example Program Results" << endl;
   /* Query and allocate the optimal workspace */
   lwork = -1;
   dsyev_( "Vectors", "Upper", &lda2, a, &lda, w, &wkopt, &lwork, &info );
   lwork = (int)wkopt;
   double work[lwork];
   /* Solve eigenproblem */
   dsyev_( "Vectors", "Upper", &lda2, a, &lda, w, work, &lwork, &info );
   if(info > 0){
      cout << "The algorithm failed to compute eigenvalues." << endl;
      exit(0);
   }
   cout << "eigenvalues are: " << w[0] << "  "
                               << w[1] << "  "
                               << w[2] << endl;
   return;
}


/*
void EDarpack( int Hsize, vector<HEle> *phami, 
               double *eval, double *evec, int nev){

//                 PARPACK Initialization                 
!     %---------------------------------------------------%^M
!     | This program uses exact shifts with respect to    |^M
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |^M
!     | IPARAM(3) specifies the maximum number of Arnoldi |^M
!     | iterations allowed.  Mode 1 of PSSAUPD is used    |^M
!     | (IPARAM(7) = 1).  All these options may be        |^M
!     | changed by the user. For details, see the         |^M
!     | documentation in PSSAUPD.                         |^M
!     %---------------------------------------------------%^M
    if(phami->empty()){
       cout << "The Hamiltonian matrix is empty!" << endl;
       return;
    }
    if(nev<=0 || nev>Hsize){
       cout << "Please set a different nev value" << endl;
       return;
    }
//    int tag=0;                    //MPI tag mark
    int rvec=0;                   //The tag of post process, 1 for process
    int iparam[11]={0};
    int ipntr[11]={0};
    char bmat[2]="I";
    int ido=0;
    char which[3]="SA";             //Keep the smallest part of the eignvalue spectrum.
    double tol=INFINITISIMAL;          //The error tolorance
    int ldv=Hsize;
    int info = 0;
    int ncv=min(Hsize, max(2*nev, nev+3));
    int lworkl=ncv*(ncv+8);         //The length of work space
    double *resid_pt=(double *)malloc(ldv*sizeof(double));
    double *V_pt=(double *)malloc(ldv*ncv*sizeof(double));
    double *workd_pt=(double *)malloc(3*ldv*sizeof(double));
    double *workl_pt=(double *)malloc(lworkl*sizeof(double));
    for(int i=0;i<lworkl;i++) workl_pt[i]=0;
    for(int i=0;i<3*ldv;i++) workd_pt[i]=0;
//    int ndigit = -3;
//    int logfil = 6;
//    int msaupd = 1;
    int ishfts=1;
    int mode=1;
    iparam[0] = ishfts;
    iparam[2] = MAXITERATION;
    iparam[6] = mode;
    //The x vector in y=OP*x multiplication. It is used to receive data from processors
    double *vector_in=(double *)malloc(Hsize*sizeof(double));           
    int count=0;
!        %---------------------------------------------%^M
!        | Repeatedly call the routine PSSAUPD and take| ^M
!        | actions indicated by parameter IDO until    |^M
!        | either convergence is indicated or maxitr   |^M
!        | has been exceeded.                          |^M
!        %---------------------------------------------%^M
 
    dsaupd_( &ido, bmat, &Hsize, which, &nev, &tol,
             resid_pt, &ncv, V_pt, &ldv, iparam, ipntr, workd_pt, 
             workl_pt, &lworkl, &info); //run the first iteration
    count++;
!         %-------------------------------------------%^M
!         | M A I N   L O O P (Reverse communication) |^M
!         %-------------------------------------------%^M
^M
!           %--------------------------------------%^M
!           | Perform matrix vector multiplication |^M
!           |              y <--- OP*x             |^M
!           | The user should supply his/her own   |^M
!           | matrix vector multiplication routine |^M
!           | here that takes workd(ipntr(1)) as   |^M
!           | the input, and return the result to  |^M
!           | workd(ipntr(2)).                     |^M
!           %--------------------------------------%                    ^M

    for(;ido==1||ido==-1;count++) {

       AlgVectorInit(workd_pt+ipntr[1]-1, Hsize);     //prepare the output space
        
       //for(i=0;i<nloc;i++){
          for(vector<HEle>::iterator Hrow_itr=phami->begin(); 
              Hrow_itr!=phami->end(); Hrow_itr++) {
              int mycol=Hrow_itr->getx();
              int myrow=Hrow_itr->gety();
              double myVal=Hrow_itr->getVal();
              //if(count==1) cout << "\t  " << mycol << " " << myrow << " " << myVal << endl;
              workd_pt[ipntr[1]+mycol-1]+=myVal*workd_pt[ipntr[0]+myrow-1];
          }//for temp_pt (j)
       //}//for i

       dsaupd_( &ido, bmat, &Hsize, which, &nev, 
                &tol, resid_pt, &ncv, V_pt, &ldv, iparam, ipntr,
                workd_pt, workl_pt, &lworkl, &info);
        
       cout << "\t\tThe " << count << "-th iteration of dsaupd_" << endl;
        
    }//for ido iteration
!     %----------------------------------------%^M
!     | Either we have convergence or there is |^M
!     | an error.                              |^M
!     %----------------------------------------%^M
    if(info<0){      //Error message. Check the documentation in PSSAUPD
!        %--------------------------%^M
!        | Error message. Check the |^M
!        | documentation in PSSAUPD.|^M
!        %--------------------------%^M
        cout << " Print Log: Error appears in _saupd, info = " << info << endl 
             << " Check the documentation of _saupd. " << endl
             << " Print Log: nev=" << nev << "    ncv=" << ncv;
        //     << " nloc=" << nloc << endl << endl;
        return;
    } else {
!        %-------------------------------------------%^M
!        | No fatal errors occurred.                 |^M
!        | Post-Process using PSSEUPD.               |^M
!        |                                           |^M
!        | Computed eigenvalues may be extracted.    |  ^M
!        |                                           |^M
!        | Eigenvectors may also be computed now if  |^M
!        | desired.  (indicated by rvec = .true.)    | ^M
!        %-------------------------------------------%^M
        rvec=1;
        int ldz=Hsize;
        double sigma;
        char HowMny[4]="All";
        int *select_pt=(int *)malloc(ncv*sizeof(double));
        double *d_pt=(double *)malloc(nev*sizeof(double));
        //        (*eigenvectors_pt)=(double *)malloc(nev*ldz*sizeof(double));  
        //record the eigenvectors space pointor by its address
        
        cout << "\tPrint Log: Ready for post process pdseupd, ncv=" 
             << ncv << endl << endl;
        
        dseupd_( &rvec, HowMny, select_pt, d_pt, evec, 
                 &ldz, &sigma, bmat, &Hsize, which, &nev, &tol, resid_pt, &ncv, 
                 V_pt, &ldv, iparam, ipntr, workd_pt, workl_pt, &lworkl, &info);
!        %----------------------------------------------%^M
!        | Eigenvalues are returned in the first column |^M
!        | of the two dimensional array D and the       |^M
!        | corresponding eigenvectors are returned in   |^M
!        | the first NEV columns of the two dimensional |^M
!        | array V if requested.  Otherwise, an         |^M
!        | orthogonal basis for the invariant subspace  |^M
!        | Corresponding to the eigenvalues in D is     |^M
!        | returned in V.                               |^M
!        %----------------------------------------------%^M
        if(info){     //Error condition
            cout << "Print Log: Error appears in _seupd, infor = " << info << endl 
                 << "Check the documentation of _seupd. " << endl 
                 << "Print Log: nev=" << nev << "    ncv=" << ncv;
            return;
        } else {
            cout << "\tPrint Log: ************************" << endl << endl;
            cout << "\tPrint Log: Eigenvalues and eigenvectors are calculated. Start to copy the data." << endl << endl;
            cout << "\tPrint Log: " << endl << endl << "Eigenvalues:" <<endl << endl;
            for(int i=0;i<nev;i++){
                eval[i]=d_pt[i]; 
                //cout << evec[0+Hsize*i] << endl;                            
                //cout << "\tPrint Log: eigens No. " << i << ":   " << d_pt[i] << endl; 
                //cout << "\tPrint Log: eigens No. " << i << ":   " << *(d_pt+i) << endl; 
                cout << "\tPrint Log: eigens No. " << i << ":   " << eval[i] << endl; 
            }
            free(select_pt);            
        }
!        %------------------------------------------%^M
!        | Print additional convergence information |^M
!        %------------------------------------------%^M
        if(info==1){
           cout << endl << "Maximum number of iterations reached."<< endl;
        } else if(info==3){
           cout << endl << "No shifts could be applied during implicit Arnoldi update"
                        << " try increasing NCV." << endl;
        }    
        cout << endl << "_SDRV1" << endl << " ======" << endl;
        cout << "\t Size of the matrix is " << Hsize << endl;
        //cout << "\t The number of processors is " << nprocs << endl;
        cout << "\t The number of Ritz values requested is " << nev << endl;
        cout << "\t The number of Arnoldi vectors generated (NCV) is " << ncv << endl;
        cout << "\t What portion of the spectrum: " << which << endl;
        // cout << " The number of converged Ritz values is " << nconv << endl;
        cout << "\t The number of Implicit Arnoldi update iterations taken is " << iparam[2] << endl;
        cout << "\t The number of OP*x is " << iparam[8] << endl;
        cout << "\t The convergence criterion is " << tol << endl;
        cout << endl;
    } //end else if(info<0)
    
    free(resid_pt);
    free(V_pt);
    free(workd_pt);
    free(workl_pt);
    free(vector_in);

    return;
}//EDarpack
*/







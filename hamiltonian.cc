/* File: hamiltonian.cc
 * ----------------
 * Implementation of genHami with many auxiliaryauxiliary  functions
 */

#include <iostream>
#include <bitset>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "para.h"
#include "hamiltonian.h"
#include "cluster.h"

using namespace std;


Cluster cl(N); //N is defined elsewhere
typedef pair<bitset<2*N>, double> statType;

//function prototype
vector<statType> delta(char, int, int, statType, double);
vector<statType> deltaDagger(char, int, int, statType, double);
vector<statType> dwave(int, statType);
vector<statType> dwaveDagger(int, vector<statType>);
vector<statType> swaveExtend(statType, double);
vector<statType> swaveExtendDagger(vector<statType>, double);
vector<statType> compress(vector<statType>);
statType Boper(int, statType);
statType BoperDagger(int, statType);
void signCount(int, statType &);
void signCount(int, statType &);
bool crt(char, int, statType &);
bool anhl(char, int, statType &);
bool tst(char, int, statType);
bool count(unsigned long, int, int);
bool countSingle(unsigned long, int);
double rho(char, int, statType &);


void genState(Para para, vector<bitset<2*N> > &states){
   int nup = para.nup, ndn = para.ndn;
   for(unsigned long i=0; i<(1UL<<(N)); i++){
      //if(i % 10000 == 0) cout << "genState at i = " << i << endl;
      if(countSingle(i, nup)){
         for(unsigned long j=0; j<(1UL<<(N)); j++){
            if(countSingle(j, ndn)){
               bitset<2*N> bitvec(i*(1UL<<N)+j);
               states.push_back(bitvec);
            }
         }
      }
   }
   return;
}

void genNq(int Hsize, vector<bitset<2*N> > states, Para pa, 
           vector<HEle> *phami){

   for(int i=0; i<Hsize; i++){
      bitset<2*N> stat = states[i];
      statType sp(stat,1); 
      double tsum_pi= 0.0, tsum_0 = 0.0;
      double rho0 = rho('u', 0, sp) + rho('d', 0, sp);
      for(int nj=0; nj<N; nj++){
         double rhoj = rho('u', nj, sp) + rho('d', nj, sp);
         tsum_pi += cos(3.1415926 * cl.xcoord(nj) + 3.1415926 * cl.ycoord(nj)) * rhoj * rho0;
         tsum_0  += rhoj * rho0;
      }
      HEle he(i, i, tsum_pi);
      phami->push_back(he);
   }
   return;
}

void genSq(int Hsize, vector<bitset<2*N> > states, Para pa, 
           vector<HEle> *phami){

   for(int i=0; i<Hsize; i++){
      bitset<2*N> stat = states[i];
      statType sp(stat,1); 
      double tsum_pi= 0.0, tsum_0 = 0.0;
      double rho0 = rho('u', 0, sp) - rho('d', 0, sp);
      for(int nj=0; nj<N; nj++){
         double rhoj = rho('u', nj, sp) - rho('d', nj, sp);
         tsum_pi += cos(3.1415926 * cl.xcoord(nj) + 3.1415926 * cl.ycoord(nj)) * rhoj * rho0;
         tsum_0  += rhoj * rho0;
      }
      HEle he(i, i, tsum_pi);
      phami->push_back(he);
   }
   return;
}

void genBmatrix(int Hsize, vector<bitset<2*N> > states, 
             Para pa, double *evect1, double *mat) {

   // initialize mat
   for(int i=0; i<N; i++)
      for(int j=0; j<N; j++)
          mat[i*N + j] = 0.0;
   cout << "initialize genDmatrix done " << endl;

   //iterate over all states
   for(int i=0; i<Hsize; i++){
      //if(i % 1000 == 0) 
      //cout << "genHami at i = " << i << endl;
      bitset<2*N> stat = states[i];

      //s-wave pair field term
      for(int nj=0; nj<N; nj++){
         statType sp(stat,1); 
         statType sp1(Boper(nj,sp));
         if(sp1.second == 0.0) continue;
         for(int ni=0; ni<N; ni++){
             statType sp2(BoperDagger(ni,sp1));
             if(sp2.second != 0.0){

               vector<bitset<2*N> >::const_iterator itr2
                      = dichotomy_search(states.begin(), states.end(), sp2.first);
               int j=itr2-states.begin();
               mat[ni*N + nj] += (sp2.second) * evect1[j] * evect1[i]; 
             }
         }
      }

   }//iterate over all states
}//genBmatrix

void genBdaggerB(int Hsize, vector<bitset<2*N> > states, 
             Para pa, vector<HEle> *phami) {

   //iterate over all states
   for(int i=0; i<Hsize; i++){
      //if(i % 1000 == 0) 
      //cout << "genHami at i = " << i << endl;
      bitset<2*N> stat = states[i];
      vector<statType> iset;
      iset.reserve(100);
      iset.clear();

      //s-wave pair field term
      for(int nj=0; nj<N; nj++){
         statType sp(stat,1); 
         statType sp1(Boper(nj,sp));
         if(sp1.second == 0.0) continue;
         for(int ni=0; ni<N; ni++){
             statType sp2(BoperDagger(ni,sp1));
             if(sp2.second != 0.0)
               iset.push_back(make_pair(sp2.first,sp2.second));
         }
      }

      vector<statType> inset;
      inset=compress(iset);
      vector<statType>::const_iterator itr = inset.begin();
      while(itr != inset.end()){
         vector<bitset<2*N> >::const_iterator itr2 
               // = find(states.begin(), states.end(), itr->first);
                = dichotomy_search(states.begin(), states.end(), itr->first);
         int j=itr2-states.begin();
         HEle he(i, j, itr->second);
         phami->push_back(he);
         itr++;
      }

   }//iterate over all states
}//genBdaggerB

void genDmatrix(int Hsize, vector<bitset<2*N> > states, 
             Para pa, double *evect1, double *mat) {

   // initialize mat
   for(int i=0; i<N; i++)
      for(int j=0; j<N; j++)
         mat[i*N + j] = 0.0;
   cout << "initialize genDmatrix done " << endl;

   //iterate over all states
   for(int i=0; i<Hsize; i++){
      //if(i % 1000 == 0) 
      //cout << "genDmatrix at i = " << i << endl;
      bitset<2*N> stat = states[i];

      //d-wave pair field term
      for(int nj=0; nj<N; nj++){
         statType sp(stat,1); 
         vector<statType> tset(dwave(nj,sp));
         for(int ni=0; ni<N; ni++){
             vector<statType> myset(dwaveDagger(ni,tset));
             vector<statType>::const_iterator itrin = myset.begin();
             while(itrin != myset.end()){
               //iset.push_back(make_pair(itrin->first,1*(itrin->second)));
               //
               vector<bitset<2*N> >::const_iterator itr2
                      = dichotomy_search(states.begin(), states.end(), itrin->first);
               int j=itr2-states.begin();
               mat[ni*N + nj] += (itrin->second) * evect1[j] * evect1[i];
               itrin++;
             }
         }
      }

   }//iterate over all states
}//genDmatrix


void genDdaggerD(int Hsize, vector<bitset<2*N> > states, 
             Para pa, vector<HEle> *phami) {

   //iterate over all states
   for(int i=0; i<Hsize; i++){
      //if(i % 1000 == 0) 
      //cout << "genHami at i = " << i << endl;
      bitset<2*N> stat = states[i];
      vector<statType> iset;
      iset.reserve(100);
      iset.clear();

      //s-wave pair field term
      for(int nj=0; nj<N; nj++){
         statType sp(stat,1); 
         vector<statType> tset(dwave(nj,sp));
         for(int ni=0; ni<N; ni++){
             vector<statType> myset(dwaveDagger(ni,tset));
             vector<statType>::const_iterator itrin = myset.begin();
             while(itrin != myset.end()){
               iset.push_back(make_pair(itrin->first,1*(itrin->second)));
               itrin++;
             }
         }
      }

      vector<statType> inset;
      inset=compress(iset);
      vector<statType>::const_iterator itr = inset.begin();
      while(itr != inset.end()){
         vector<bitset<2*N> >::const_iterator itr2 
               // = find(states.begin(), states.end(), itr->first);
                = dichotomy_search(states.begin(), states.end(), itr->first);
         int j=itr2-states.begin();
         HEle he(i, j, itr->second);
         phami->push_back(he);
         itr++;
      }

   }//iterate over all states
}//genDdaggerD

void genAmatrix(int Hsize, vector<bitset<2*N> > states, 
             Para pa, double *evect1, double *mat) {

   // initialize mat
   for(int i=0; i<N; i++)
      for(int j=0; j<N; j++)
         mat[i*N + j] = 0.0;
   cout << "initialize genDmatrix done " << endl;

   //iterate over all states
   for(int i=0; i<Hsize; i++){
      //if(i % 1000 == 0) 
      //cout << "genDmatrix at i = " << i << endl;
      bitset<2*N> stat = states[i];

      //d-wave pair field term
      for(int nj=0; nj<N; nj++){
         statType sp(stat,1); 
         vector<statType> tset(dwave(nj,sp));
         for(int ni=0; ni<N; ni++){
             vector<statType> myset(dwaveDagger(ni,tset));
             vector<statType>::const_iterator itrin = myset.begin();
             while(itrin != myset.end()){
               //iset.push_back(make_pair(itrin->first,1*(itrin->second)));
               //
               vector<bitset<2*N> >::const_iterator itr2
                      = dichotomy_search(states.begin(), states.end(), itrin->first);
               int j=itr2-states.begin();
               mat[ni*N + nj] += (itrin->second) * evect1[j] * evect1[i];
               itrin++;
             }
         }
      }

   }//iterate over all states
}//genAmatrix




void genExtendedswave(int Hsize, vector<bitset<2*N> > states, 
             Para pa, vector<HEle> *phami) {

   //iterate over all states
   for(int i=0; i<Hsize; i++){
      //if(i % 1000 == 0) 
      //cout << "genHami at i = " << i << endl;
      bitset<2*N> stat = states[i];
      //vector<statType> iset;
      //iset.reserve(100);
      //iset.clear();

      //Extended s-wave term
      statType sp(stat,1); 
      vector<statType> tset(swaveExtend(sp, pa.alpha));
      vector<statType> myset(swaveExtendDagger(tset, pa.alpha));
      //vector<statType>::const_iterator itrin = myset.begin();
      //while(itrin != myset.end()){
      //  iset.push_back(make_pair(itrin->first,1*(itrin->second)));
      //  itrin++;
      //}

      /*
      for(int nj=0; nj<N; nj++){
         statType sp(stat,1); 
         vector<statType> tset(swaveExtend(nj,sp));
         for(int ni=0; ni<N; ni++){
             vector<statType> myset(swaveExtendDagger(ni,tset));
             vector<statType>::const_iterator itrin = myset.begin();
             while(itrin != myset.end()){
               iset.push_back(make_pair(itrin->first,1*(itrin->second)));
               itrin++;
             }
         }
      }
      */
      

      vector<statType> inset;
      inset=compress(myset);
      vector<statType>::const_iterator itr = inset.begin();
      while(itr != inset.end()){
         vector<bitset<2*N> >::const_iterator itr2 
               // = find(states.begin(), states.end(), itr->first);
                = dichotomy_search(states.begin(), states.end(), itr->first);
         int j=itr2-states.begin();
         HEle he(i, j, itr->second);
         phami->push_back(he);
         itr++;
      }

   }//iterate over all states
}//genExtendedswave


void genHami(int Hsize, vector<bitset<2*N> > states, 
             Para pa, vector<HEle> *phami) {

   //iterate over all states
   for(int i=0; i<Hsize; i++){
      //if(i % 1000 == 0) 
      //cout << "genHami at i = " << i << endl;
      bitset<2*N> stat = states[i];
      vector<statType> iset;
      iset.reserve(100);
      iset.clear();

      //Hubbard U
      statType sp(stat,1);
      double sum=0;
      for(int ni=0; ni<N; ni++){
         if(tst('u',ni,sp) && tst('d',ni,sp)){
            sum += 1;
         }
      }
      if(sum){
         iset.push_back(make_pair(stat,sum*pa.U));
      }
     
      //hoppint t
      //for(int ni=0; ni<N; ni++){
      //   vector<int> mynei(cl.neighbor(ni)); 
      //   vector<int>::const_iterator itr = mynei.begin();
      //   while(itr != mynei.end()){
      //      string str("ud");
      //      for(int inds=0; inds<2; inds++){
      //         statType sp(stat,1);
      //         if(anhl(str[inds],ni,sp) && crt(str[inds], *itr, sp)){
      //   	  iset.push_back(make_pair(sp.first, -1*(sp.second)*pa.t));
      //         }
      //      }           
      //      itr++;
      //   }
      //}
      for(int ni=0; ni<N; ni++){
         string str("ud");
         
         int nx = cl.xneighbor(ni);
         double phx = cl.xphase(ni);
         for(int inds=0; inds<2; inds++){
            statType sp(stat,1);
            if(anhl(str[inds],ni,sp) && crt(str[inds], nx, sp)){
               iset.push_back(make_pair(sp.first, -1*(sp.second)*pa.t*phx));
            }
            statType spp(stat,1);
            if(anhl(str[inds],nx,spp) && crt(str[inds], ni, spp)){
               iset.push_back(make_pair(spp.first, -1*(spp.second)*pa.t*phx));
            }
         }           

         int ny = cl.yneighbor(ni);
         double phy = cl.yphase(ni);
         for(int inds=0; inds<2; inds++){
            statType sp(stat,1);
            if(anhl(str[inds],ni,sp) && crt(str[inds], ny, sp)){
               iset.push_back(make_pair(sp.first, -1*(sp.second)*pa.t*pa.alpha*phy));
            }
            statType spp(stat,1);
            if(anhl(str[inds],ny,spp) && crt(str[inds], ni, spp)){
               iset.push_back(make_pair(spp.first, -1*(spp.second)*pa.t*pa.alpha*phy));
            }
         }           

      }

      //nearest neighbor hopping t'
      for(int ni=0; ni<N; ni++){
         vector<int> mynei(cl.nnneighbor(ni)); 
         vector<int>::iterator itr = mynei.begin();
         while(itr != mynei.end()){
            string str("ud");
            for(int inds=0; inds<2; inds++){
               statType sp(stat,1);
               if(anhl(str[inds],ni,sp) && crt(str[inds], *itr, sp)){
                  iset.push_back(make_pair(sp.first, -1*(sp.second)*pa.tp));
               }
            }
            itr++;
         }        
      }

      //nearest neighbor hopping t''
      for(int ni=0; ni<N; ni++){
         vector<int> mynei(cl.nnnneighbor(ni)); 
         vector<int>::iterator itr = mynei.begin();
         while(itr != mynei.end()){
            string str("ud");
            for(int inds=0; inds<2; inds++){
               statType sp(stat,1);
               if(anhl(str[inds],ni,sp) && crt(str[inds], *itr, sp)){
                  iset.push_back(make_pair(sp.first, -1*(sp.second)*pa.tpp));
               }
            }
            itr++;
         }        
      }

      //nearest neighbor Coulomb interaction V
      sum=0;
      for(int ni=0; ni<N; ni++){
         vector<int> mynei(cl.neighbor(ni)); 
         vector<int>::iterator itr = mynei.begin();
         while(itr != mynei.end()){
            string str("ud");
            for(int inds1=0; inds1<2; inds1++){
               for(int inds2=0; inds2<2; inds2++){
                  statType sp(stat,1);
                  if(tst(str[inds1],ni,sp) && tst(str[inds2],*itr,sp)){
                     sum++;
                  }
               }
            }
            itr++;
         }        
      }
      if(sum){
         iset.push_back(make_pair(stat,sum*pa.V*0.5));
      }

      //s-wave pair field term B_dagger B
      for(int nj=0; nj<N; nj++){
         statType sp(stat,1); 
         statType sp1(Boper(nj,sp));
         if(sp1.second == 0.0) continue;
         for(int ni=0; ni<N; ni++){
             statType sp2(BoperDagger(ni,sp1));
             if(sp2.second != 0.0)
               iset.push_back(make_pair(sp2.first,sp2.second*pa.lambda));
         }
      }

      //d-wave pair field term D_dager D
      for(int nj=0; nj<N; nj++){
         statType sp(stat,1); 
         vector<statType> tset(dwave(nj,sp));
         for(int ni=0; ni<N; ni++){
             vector<statType> myset(dwaveDagger(ni,tset));
             vector<statType>::const_iterator itrin = myset.begin();
             while(itrin != myset.end()){
               iset.push_back(make_pair(itrin->first,1*(itrin->second)*pa.lambda));
               itrin++;
             }
         }
      }

      //extended s-wave term A_dagger A
      statType sp1(stat,1);
      vector<statType> tset(swaveExtend(sp1, pa.alpha));
      vector<statType> myset(swaveExtendDagger(tset, pa.alpha));
      vector<statType>::const_iterator itrin = myset.begin();
      while(itrin != myset.end()){
         iset.push_back(make_pair(itrin->first,1*(itrin->second)*pa.Us));
         itrin++;
      }

      vector<statType> inset;
      inset=compress(iset);
      vector<statType>::const_iterator itr = inset.begin();
      while(itr != inset.end()){
         vector<bitset<2*N> >::const_iterator itr2 
               // = find(states.begin(), states.end(), itr->first);
                = dichotomy_search(states.begin(), states.end(), itr->first);
         int j=itr2-states.begin();
         HEle he(i, j, itr->second);
         phami->push_back(he);
         itr++;
      }

   }//iterate over all states
}//genHami


//to check if iterator works when delete
vector<statType> compress(vector<statType> mset){
   map<unsigned long, double> mymap;
   vector<statType>::const_iterator itr1 = mset.begin();
   while(itr1 != mset.end()){
      unsigned long stat = (itr1->first).to_ulong(); 
      double val = itr1->second;
      mymap[stat] = mymap[stat]+val; 
      itr1++;
   }
 
   vector<statType> myvec;
   map<unsigned long, double>::const_iterator itr2 = mymap.begin();
   while(itr2 != mymap.end()){
      if(itr2->second != 0.0){
         bitset<2*N> bs(itr2->first);
         pair< bitset<2*N>, double> myp(bs, itr2->second);
         myvec.push_back(myp);
      }
      itr2++;
   }
   return myvec;
}

statType Boper(int ni, statType st){
   statType myst1 = st;
   if(anhl('d',ni,myst1) && anhl('u',ni,myst1)) {}
   else { myst1.second = 0.0;}
   return myst1;
}

statType BoperDagger(int ni, statType st){
   statType myst1 = st;
   if(crt('u',ni,myst1) && crt('d',ni,myst1)) {}
   else { myst1.second = 0.0;}
   return myst1;
}

vector<statType> delta(char a, int ni, int nj, statType st, double alpha = 1.0){
   vector<statType> myset;
   statType myst = st;
   switch(a){
      case 'p': break;
      case 'n': myst.second *= -1; break;
      default: cout << "wrong input for delta!" << endl;
   }
   if(alpha != 1.0) myst.second *= alpha;
   statType myst1 = myst;
   if(anhl('u',nj,myst1) && anhl('d',ni,myst1))
      myset.push_back(myst1);
   statType myst2 = myst;
   if(anhl('u',ni,myst2) && anhl('d',nj,myst2))
      myset.push_back(myst2);
   return myset;
}

vector<statType> deltaDagger(char a, int ni, int nj, statType st, double alpha = 1.0){
   vector<statType> myset;
   statType myst = st;
   switch(a){
      case 'p': break;
      case 'n': myst.second *= -1; break;
      default: cout << "wrong input for delta!" << endl;
   }
   if(alpha != 1.0) myst.second *= alpha;
   statType myst1 = myst;
   if(crt('d',nj,myst1) && crt('u',ni,myst1))
      myset.push_back(myst1);
   statType myst2 = myst;
   if(crt('d',ni,myst2) && crt('u',nj,myst2))
      myset.push_back(myst2);
   return myset;
}

vector<statType> swaveExtend(statType st, double alpha){
   vector<statType> myset;
   for(int ni=0; ni<N; ni++){
      int nx = cl.xneighbor(ni), ny = cl.yneighbor(ni);
      double phx = cl.xphase(ni), phy = cl.yphase(ni);
      vector<statType> tset(delta('p',ni, nx, st, phx));
      myset.insert(myset.begin(), tset.begin(), tset.end());
      // setting alpha in y direction
      vector<statType> tset2(delta('p',ni, ny, st, alpha*phy));
      myset.insert(myset.begin(), tset2.begin(), tset2.end());
   }
   return myset;
} 

vector<statType> swaveExtendDagger(vector<statType> st, double alpha){
   vector<statType> myset;
   for(int ni=0; ni<N; ni++){
      int nx = cl.xneighbor(ni), ny = cl.yneighbor(ni);
      double phx = cl.xphase(ni), phy = cl.yphase(ni);
      vector<statType>::const_iterator itr = st.begin();
      while(itr != st.end()){
         vector<statType> tset(deltaDagger('p',ni, nx, *itr, phx));
         myset.insert(myset.begin(), tset.begin(),tset.end());
         // setting alpha in y direction
         vector<statType> tset2(deltaDagger('p',ni, ny, *itr, alpha*phy));
         myset.insert(myset.begin(), tset2.begin(),tset2.end());
         itr++;
      }
   }
   return myset;
}

vector<statType> dwave(int n, statType st){
   int nx = cl.xneighbor(n), ny = cl.yneighbor(n);
   vector<statType> myset(delta('p',n, nx, st));
   vector<statType> tset(delta('n',n, ny, st));
   myset.insert(myset.begin(), tset.begin(), tset.end());
   return myset;
} 

vector<statType> dwaveDagger(int n, vector<statType> st){
   vector<statType> myset;
   int nx = cl.xneighbor(n), ny = cl.yneighbor(n);
   vector<statType>::const_iterator itr = st.begin();
   while(itr != st.end()){
      vector<statType> tset(deltaDagger('p',n, nx, *itr));
      myset.insert(myset.begin(), tset.begin(),tset.end());
      vector<statType> tset2(deltaDagger('n',n, ny, *itr));
      myset.insert(myset.begin(), tset2.begin(),tset2.end());
      itr++;
   }
   return myset;
}

// c^{\dagger}_{n,\sigma} c_{n,\sigma}
double rho(char a, int n, statType &st){
   switch(a){   
      case 'u':
         if(st.first.test(n+N)) return 1.0;
         else return 0.0;
         break;
      case 'd':
         if(st.first.test(n)) return 1.0; 
         else return 0.0;
         break;
      default:
         cout << "state::rho has wrong input!" << endl;
         return 0.0;
   }
   return 0.0;
}

//high bits for spin up, low bits for spin down
bool crt(char a, int n, statType &st){
   switch(a){   
      case 'u':
         if(st.first.test(n+N)) return false;
         st.first.set(n+N);
         signCount(n+N, st);
         break;
      case 'd':
         if(st.first.test(n)) return false; 
         st.first.set(n);
         signCount(n, st);
         break;
      default:
         cout << "state::crt has wrong input!" << endl;
         return false;
   }
   return true;
}

bool anhl(char a, int n, statType &st){
   switch(a){   
      case 'u':
         if(!st.first.test(n+N)) return false;
         st.first.reset(n+N);
         signCount(n+N, st);
         break;
      case 'd':
         if(!st.first.test(n)) return false; 
         st.first.reset(n);
         signCount(n, st);
         break;
      default:
         cout << "state::anhl has wrong input!" << endl;
         return false;
   }
   return true;
}

bool tst(char a, int n, statType st){
   switch(a){
      case 'u':
         return st.first.test(n+N);
      case 'd':
         return st.first.test(n);
      default:
         return false;
   }
}

void signCount(int n, statType &st){
   unsigned long mask = (1UL<<n) - 1;
   unsigned long init = mask & st.first.to_ulong();
   int sum=0;
   while(init){
      sum += (init & 1UL);
      init >>= 1;
   }
   if(sum & 1UL == 1UL) st.second *= -1;
   //cout << "signCount" << s << " " << n << " " << sum << endl;
   return;
}

//high bits for spin up, low bits for spin down
bool count(unsigned long s, int up, int dn){
   int sumUp=0, sumDn=0;
   for(int i=0; i<N; i++){
      if(s & 1UL<<i) sumDn++;   
      if(s & 1UL<<(i+N)) sumUp++;
   }
   if(sumUp == up && sumDn == dn) return true;
   else return false;
}

bool countSingle(unsigned long s, int n){
   int sum=0;
   for(int i=0; i<N; i++){
      if(s & 1UL<<i) sum++;
   }
   return sum==n;
}

vector<bitset<2*N> >::const_iterator dichotomy_search(vector<bitset<2*N> >::const_iterator start_itr, vector<bitset<2*N> >::const_iterator end_itr, bitset<2*N> content)
{
    vector<bitset<2*N> >::const_iterator initialenditr=end_itr;
 //       vector<bitset<2*N> >::const_iterator initialstitr=start_itr;
    vector<bitset<2*N> >::const_iterator pitr=start_itr;
  
    if((*start_itr)==content) return start_itr;
    
    for(; start_itr<end_itr-1;)
    {
        if((*pitr)==content)
        {
//            cout << "find the state No." << pitr-initialstitr << endl;
            
            return pitr;
        }
        
        if((*pitr)<content) start_itr=pitr;
        else end_itr=pitr;
        pitr=start_itr+(end_itr-start_itr+1)/2;
    }
    
    return initialenditr;
   
}

//my new implemented < for bitset<2*N>
bool  operator<(bitset<2*N> input1, bitset<2*N> &input2)
{
   return input1.to_ulong() < input2.to_ulong();
}


bool  operator>(bitset<2*N> input1, bitset<2*N> &input2)
{
   return input1.to_ulong() > input2.to_ulong();
}

bool  operator==(bitset<2*N> input1, bitset<2*N> &input2)
{
   return input1.to_ulong() == input2.to_ulong();
}

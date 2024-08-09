#ifndef _H_para
#define _H_para

const int N=10;
const int dim=1;

class Para {
public:
   int nup, ndn;
   double U, V, t, tp, tpp, alpha, lambda, Us,ltol, htol;
   double lStart, lEnd, lStep, dl;
   double UsStart, UsEnd, UsStep, dUs;
};

#endif

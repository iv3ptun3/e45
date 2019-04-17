#ifndef CHI2FCN_H
#define CHI2FCN_H

//#include "globals.hh"
//#include <math.h>
#include "minuit2.hh"
#include "Minuit2/FCNBase.h"

class chi2fcn : public FCNBase{
  explicit chi2fcn();
  ~chi2fcn();
public:
  double operator()(int np, double *g, double*u, int iflag); 
};

#endif

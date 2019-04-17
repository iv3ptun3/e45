///let's make a RKtracking code by S. Hwang
#ifndef RKTRACKING_HH
#define RKTRACKING_HH

#include "minuit2.hh"

class RKTracking{
private:
  double RKPar[5];
  double StepPar[5];

public:
  RKTracking();
  ~RKTracking();

  double RKTracking(double *xin,double *yin,double *zin, int npoints,double *rkpar, int iflag);

};

#endif

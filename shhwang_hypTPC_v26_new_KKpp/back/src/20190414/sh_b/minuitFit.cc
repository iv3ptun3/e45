#include "track.hh"
#include "minuitFit.hh"

//#include "cfortran/minuit.h"
//#include "minuit.h"

//#include "chi2FCN.hh"
#include "minuit2.hh"

void chi2(int *npar, double *grad, double *fval,
	  double *xval, int *iflag, void (*futil)() );
MinuitFCN fcn =  chi2 ;

extern double RKChi2[MAX_ITERATION];
extern double RKPara[MAX_ITERATION][NUM_PARA_RK];
extern double RKPadHitTmp[MAX_HIT_IN_TRACK+1][3];
extern double RKPadErrTmp[MAX_HIT_IN_TRACK+1][3];

extern int RKNumHitTmp;

int minuitFit(double* initPara, float* finalChi2, double* para,
	      Track* aTrack){
  int i,j;
  const int nPara = NUM_PARA_RK;
  char* paraName[NUM_PARA_RK]={"s","t","p","theta","phi"};
  
  double error[NUM_PARA_RK];
  double stepPara[NUM_PARA_RK];
  double lowLimit[NUM_PARA_RK];
  double upLimit[NUM_PARA_RK];


  //  minuitInit(-1);

  stepPara[0]=0.05;
  stepPara[1]=0.05;
  stepPara[2]=0.005;
  stepPara[3]=0.1;
  stepPara[4]=1.;
  ///create chi2 function based on helix function
  //  chi2FCN chi2fcn(aTrack);
  return 0;
}

//double chi2Fcn::operator(int *npar, double *grad, double *fval, 
//			 double *xval, int * iflag){
//}



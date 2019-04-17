#include <math.h>
#include <float.h>

#include "minuitFit.h"
#include "minuit.h"
//#include "traceInSolenoid.h"
//#include "arraySize.h"
#include "track.h"
//#include "rungeKuttaTrack.h"
//#include <switch.h>

void chi2(int *npar, double *grad, double *fval,
        double *xval, int *iflag, void (*futil)() );
MinuitFCN fcn =  chi2 ;
extern double RKChi2[MAX_ITERATION];
extern double RKPara[MAX_ITERATION][NUM_PARA_RK];
/*extern struct RKVirtualPlane;*/
/*
extern double RKPadHit[MAX_HIT_IN_TRACK+1][3];
*/
extern double RKPadHitTmp[MAX_HIT_IN_TRACK+1][3];
extern double RKPadHitlocalTmp[MAX_HIT_IN_TRACK+1][3];
extern double RKPadErrTmp[MAX_HIT_IN_TRACK+1][3];

/*
extern int    RKNumHit;
*/
extern int    RKNumHitTmp;
extern int    RKPadSector[MAX_HIT_IN_TRACK+1];
extern int    RKPadLay[MAX_HIT_IN_TRACK+1];

/*************************************************
 *************************************************/
void minuitInit(double printLevel)
/* Control PrintOut LEVEL -1:no, 0:min., 1:default, 3:max. */
{
  //  double argument[10];
  //  int nArg;
  //  int errFlag;

  /* read, write, SAVE */
  //  MNINIT( 5, 6, 7 );

  //  argument[0] = printLevel; nArg=1;
  //  MNEXCM( fcn, "set print",      argument, nArg, & errFlag, 0 );

  //  argument[0]= 1; nArg=1;
  //  MNEXCM( fcn, "set nowarnings", argument, nArg, & errFlag, 0 );

}
/*************************************************
 *************************************************/



/*************************************************
 *************************************************/
int minuitFit(double* initPara, float* finalChi2, double* para,
              Track* aTrack, int sector, int lay)
{
  G4cout<<sector<<":"lay<<G4endl;

  /*
  int i,j;
  const int nPara = NUM_PARA_RK;
  double argument[10];
  int nArg;
  int errFlag;
  int flag;
  int ivarbl;
  int numStep;
  double fedm, errdef;
  double low, up;
  int npari, nparx, istat;
  char name[256];
  double tmpChi2;
  double rKPadHit[MAX_HIT_IN_TRACK+1][3] = {{0.}};
  double rKPadHitlocal[MAX_HIT_IN_TRACK+1][3] = {{0.}};
  double rKPadErr[MAX_HIT_IN_TRACK+1][3] = {{0.}};
  int rKLay[MAX_HIT_IN_TRACK];
  int rKSec[MAX_HIT_IN_TRACK];
  /// test com-out //
  extern Track* ATrack;
  char* paraName[NUM_PARA_RK] = {
    "s", "t", "P", "theta", "phi"
  };
  // s, t defines a starting point on a plane //

  double error[NUM_PARA_RK];      // Errors for minuit parmaeters   //
  double stepPara[NUM_PARA_RK];   // Step size in Minuit            //
  double lowLimit[NUM_PARA_RK];   // Lower limit of para. in Minuit //
  double upLimit[NUM_PARA_RK];    // Upper limit of para. in Minuit //

  // Init Minuit function //
  */

  G4cout<<"minuit initialization"<<G4endl;
  ////  //  minuitInit(-1);
  
  /* Set step size */


  /*
  stepPara[0] = 0.05;     // [mm]   //
  stepPara[1] = 1.;  // [mm]  //
  stepPara[2] = 0.005;  // [GeV/c]//
  stepPara[3] = 0.1; // [rad] //
  stepPara[4] = 1.;     // [rad]  //

  for(i=0; i<NUM_PARA_RK; i++){
    upLimit[i]  = lowLimit[i] = 0.;
    para[i]     = initPara[i];
  }

  for( i=0; i< NUM_PARA_RK; i++ ) {  // Set parameters //
    MNPARM( i+1, paraName[i], initPara[i], stepPara[i],
            lowLimit[i], upLimit[i], & flag );
    if(flag>0) {
      fprintf(stderr,"Minimize: set parameter error: %d exit.\n",i);
      fflush(stdout);
      exit(1);
    }
  }

  //////////////////////////////  traceInSolenoid(nPara, para, &tmpChi2, ATrack,2, tpcPar4minuit);

  //  if(sw->resolution==1){
    // resolution check  -- from TPCana0.47.res//
    j = 0; 
    for( i = 0; i < aTrack->numHits; i++){
      // if( RKPadLay[i] != lay ){//
      //also rejecting lay0,8//
      if( RKPadLay[i] != lay && RKPadLay[i] != 0 && RKPadLay[i] != 8 ){
	rKPadHit[j][0] = aTrack->x[i][0];
	rKPadHit[j][1] = aTrack->x[i][1];
	rKPadHit[j][2] = aTrack->x[i][2];

	rKPadHitlocal[j][0] = aTrack->xlocal[i][0];
	rKPadHitlocal[j][1] = aTrack->xlocal[i][1];
	rKPadHitlocal[j][2] = aTrack->xlocal[i][2];
      
	rKPadErr[j][0] = aTrack->err[i][0];
	rKPadErr[j][1] = aTrack->err[i][1];
	rKPadErr[j][2] = aTrack->err[i][2];

	rKLay[j] = aTrack->lay[i];
	rKSec[j] = aTrack->sector[i];

	j++;
	
      }
     
    }
    aTrack->numHits = j;

    if(aTrack->numHits<4){
      return 1;
    }
    
    for( i = 0; i < aTrack->numHits; i++){
      aTrack->x[i][0] = rKPadHit[i][0];
      aTrack->x[i][1] = rKPadHit[i][1];
      aTrack->x[i][2] = rKPadHit[i][2];

      aTrack->xlocal[i][0] = rKPadHitlocal[i][0];
      aTrack->xlocal[i][1] = rKPadHitlocal[i][1];
      aTrack->xlocal[i][2] = rKPadHitlocal[i][2];

      aTrack->err[i][0] = rKPadErr[i][0];
      aTrack->err[i][1] = rKPadErr[i][1];
      aTrack->err[i][2] = rKPadErr[i][2];

      aTrack->lay[i] = rKLay[i];
      aTrack->sector[i] = rKSec[i];
    }
    ////////////
    //  }//sw

  MNSETI("final track fit");

  argument[0]=0;  nArg=1;
  MNEXCM( fcn, "set str", argument, nArg, & errFlag, 0 );

  argument[0]=1;  nArg=1;
  MNEXCM( fcn, "call fcn", argument, nArg, & errFlag, 0 );
  
  // call chi2() many times.//
  argument[0]=MAX_ITERATION; argument[1]=MINUIT_TOLERANCE;  nArg=2;
  MNEXCM( fcn, "migrad",   argument, nArg, & errFlag, 0 );

  // Get status //
  MNSTAT( &tmpChi2, & fedm, & errdef, & npari, & nparx, & istat );
  *finalChi2 = tmpChi2;
  //
  if(istat != 3){
    printf("minuit didn't completely converge istat:%d\n",istat);
  }else{
    printf("minuit converged istat:%d\n",istat);
  }
  //

  // get fit results //
  for( i=0; i<nPara; i++ ) {
    MNPOUT( i+1, name, &( para[i] ), &( error[i] ), & low, & up, & ivarbl );
  }

  //  exec: init. of fcn //
  MNEXCM( fcn, "stop", argument, nArg, & errFlag, 0 );

  //  if(sw->resolution==1){
    //resolution check mode//
    //return to original value//
    for( j=0; j < MAX_HIT_IN_TRACK+1; j++){
      aTrack->x[j][0] = RKPadHitTmp[j][0];
      aTrack->x[j][1] = RKPadHitTmp[j][1];
      aTrack->x[j][2] = RKPadHitTmp[j][2];

      aTrack->xlocal[j][0] = RKPadHitlocalTmp[j][0];
      aTrack->xlocal[j][1] = RKPadHitlocalTmp[j][1];
      aTrack->xlocal[j][2] = RKPadHitlocalTmp[j][2];
    
      aTrack->err[j][0] = RKPadErrTmp[j][0];
      aTrack->err[j][1] = RKPadErrTmp[j][1];
      aTrack->err[j][2] = RKPadErrTmp[j][2];

      aTrack->lay[j] = RKPadLay[j];
      aTrack->sector[j] = RKPadSector[j];
    }
    aTrack->numHits = RKNumHitTmp;

    //  }//sw

  //////
  
  // calculate residual for all hits ...
//     Track parameters were determined
//     without "lay"th hit for residual mode.//
    ////////////////////  aTrack->numStep = traceInSolenoid(nPara, para, &tmpChi2, ATrack, 1, tpcPar4minuit);

  *finalChi2 = tmpChi2;

  return 0;
*/
}
/////////////////////////////////



/*************************************************
 *************************************************/
void chi2(int *npar, double *grad, double *fval,
        double *xval, int *iflag, void (*futil)() ) {
  int errFlag;
  extern int RKNumIter;
  extern int eveNum;
  extern Track* ATrack;

  //////////////////////////////  errFlag = traceInSolenoid(*npar,xval,fval, ATrack, 0, tpcPar4minuit);

  if( isinf(*fval) != 0) {

    /*    *fval = DBL_MAX-1; */
    *fval = 1000000000000000000000.;
  }

  if( *iflag==4 ) {   /* set *fval   */
    if(errFlag == -1){
      /*      *iflag=3;*/
      *fval = 1000000000000000000000.;
    }
    if(RKNumIter < MAX_ITERATION ){
      int i;

      for(i=0; i < NUM_PARA_RK; i++){
        RKPara[RKNumIter][i] = xval[i];
      }

      RKChi2[RKNumIter] = *fval;

    }

    RKNumIter++;
  }
  
  else if(*iflag==1) { RKNumIter=0; } /* for initialize : read input data */
  else if(*iflag==2) { ; } /* calculate grad 1st derivatives of
                              fval(: option) */
  else if(*iflag==3) { ; } /* for finish.(after fit is finished) */

  /*
  printf("%d %lf, %lf %lf %lf %lf\n",RKNumIter
	 ,xval[0],xval[1],xval[2],xval[3],xval[4]);
  */
}
/*************************************************
 *************************************************/


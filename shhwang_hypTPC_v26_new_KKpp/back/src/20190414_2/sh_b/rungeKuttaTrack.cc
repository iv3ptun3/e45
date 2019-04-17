#include "rungeKuttaTrack.hh"
#include "minuitFit.hh"
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "TMath.h"
#include "common.hh"
#define SQ(a) ((a)*(a))

int RKNumHitTmp = 0;               /* Number of hits in the track       */
double RKPadHitTmp[MAX_HIT_IN_TRACK+1][3] = {{0.}}; /* x,y,z of hits in the track*/
double RKPadHitlocalTmp[MAX_HIT_IN_TRACK+1][3] = {{0.}}; /* x,y,z of hits in the track*/
double RKPadErrTmp[MAX_HIT_IN_TRACK+1][3] = {{0.}}; /* error of hit pos in the track*/
int    RKPadSector[MAX_HIT_IN_TRACK+1] = {0.}; /* Pad sector of hits in the track */
int    RKPadLay[MAX_HIT_IN_TRACK+1] = {0.};    /* Pad layer of hits in the track  */
double RKResXYZ[MAX_HIT_IN_TRACK+1][4] = {{0.}};/* residual in (x,y,z,sqrt(x^2+y^2+z^2)) */
double RKInitRes[MAX_HIT_IN_TRACK+1][4] = {{0.}};


int    RKNumIter=0;
double RKChi2[MAX_ITERATION] = {0.};
double RKPara[MAX_ITERATION][NUM_PARA_RK] = {{0.}};
Track* ATrack;



int  rungeKuttaTrack(Track* aTrack)
{
  int i,j;
  static const std::string funcname = "[RungeKuttaTrack]";
  std::cout<<"rungeKuttaTrack"<<std::endl;
  std::cout<<"aTrack->numHits:  "<<aTrack->numHits<<std::endl;
  std::cout<<"aTrack->trkQual:  "<<aTrack->trkQual<<std::endl;

  int ndf=2*(aTrack->numHits-aTrack->nout)-5;
  setInitPara(aTrack,aTrack->rKInitPara);

  ///initiallization of final RK parameter
  aTrack->RKPFinal[0]=  aTrack->RKPFinal[1]=  aTrack->RKPFinal[2]= -100;

  if(ndf>0){
    minuitFit(aTrack->rKInitPara, &(aTrack->chi2),aTrack->rKFinalPara, aTrack);
    ndf=2*(aTrack->ngood-aTrack->nout)-5;
    if(ndf<1){
      aTrack->chi2Prob=-1.;
    }else{
      aTrack->chi2Prob = TMath::Prob(aTrack->chi2,ndf);
      for(j=0;j<RKNumIter && j < MAX_ITERATION; j++){
	aTrack->rKChi2[j]=RKChi2[j];
	aTrack->rKNumIter=RKNumIter;
      }
    }
  } 



  return 0;
  
}
/***************************************
 ***************************************/
int setInitPara(Track* aTrack, double* initPara){

  int i,j;
  //  int setVirtualPlane(Track* aTrack);
  //  setVirtualPlane(Track* aTrack);
  //0->x0
  //1->z0
  //2->mom
  //3->theta
  //4->phi
  initPara[0]=initPara[1]=0;
  initPara[2]=aTrack->mom[3];
  initPara[3]=atan2(aTrack->mom[1],sqrt(aTrack->mom[0]*aTrack->mom[0]+aTrack->mom[2]*aTrack->mom[2]));//atan2(py,Pt)
  initPara[4]=atan2(aTrack->mom[0],aTrack->mom[2]);
  if(initPara[4] < 0.0){
    initPara[4] += 2.0*PI;
  }


  ///input rkparameter
  for( i = 0 ; i < NUM_PARA_RK ; i++){
    aTrack->rKInitPara[i]=initPara[i];
  }

  ///initiallization of resolution
  for(i = 0 ; i < MAX_HIT_IN_TRACK ; i++){
    for(j = 0 ; j < 4 ; j++){
      RKResXYZ[i][j] = -1000.;
      RKInitRes[i][j] = -1000.;
    }
  }
  
  return 0;
}


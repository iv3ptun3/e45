#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include "filter.h"
#include "kmf.h"
int i,j,k;
float Rerr;

int kmfinit1_()
{
  
  for(i = 0; i < PARASIZE; i++){
    for(j = 0; j < PARASIZE; j++){
      for(k= 0; k < PLANESIZE; k++){
		Cmats[i][j][k] = 0.;
      }
    }
  }

  Cmats[0][0][0] = 1.0;
  Cmats[1][1][0] = 1.0;
  Cmats[2][2][0] = 0.02;
  Cmats[3][3][0] = 0.02;
  Cmats[4][4][0] = 0.02; 

  return 0;
  
}

int kmfinit2_()
{
double Cms,reso,meas,measz,zcoor;


  KMFpar0[0][0] = kmf_.KMFpar0[0][0]; 
  KMFpar0[1][0] = kmf_.KMFpar0[1][0];
  KMFpar0[2][0] = kmf_.KMFpar0[2][0];
  KMFpar0[3][0] = kmf_.KMFpar0[3][0];
  KMFpar0[4][0] = kmf_.KMFpar0[4][0];
  
  reso  = kmf_.reso[0];
  Cms   = pow(kmf_.mulsth[0],2);
  meas  = kmf_.meas[0];
  measz = kmf_.measz[0];
  zcoor = kmf_.zcoor[0]; 
  
  for(i = 0; i < PARASIZE; i++){
    for(j = 0; j < PLANESIZE; j++){
      KMFpar0[i][j] = 0.;
      KMFpar1[i][j] = 0.;
      SMTpar[i][j]  = 0.;
    }
  }
  
  for(i = 0; i < PARASIZE; i++){
    for(j = 0; j < PARASIZE; j++){
      for(k= 0; k < PLANESIZE; k++){
		Fmat[i][j][k]  = 0.;
		Cmat0[i][j][k] = 0.;
		Cinv0[i][j][k] = 0.;
		Cmat1[i][j][k] = 0.;
		Qmat[i][j][k]  = 0.;
      }
    }
  }
  
  for(i = 0; i < PLANESIZE; i++){
    KMFchi2[i]    = 0.;
    KMFchi2new[i] = 0.;
    SMTchi2[i]    = 0.;
    SMTchi2new[i] = 0.;
  }

  for(i = 0; i < PLANESIZE; i++){
	KMFresi[i] = 0.;
	SMTresi[i] = 0.;
  }

  KMFpar0[0][0] = kmf_.KMFpar0[0][0]; 
  KMFpar0[1][0] = kmf_.KMFpar0[1][0];
  KMFpar0[2][0] = kmf_.KMFpar0[2][0];
  KMFpar0[3][0] = kmf_.KMFpar0[3][0];
  KMFpar0[4][0] = kmf_.KMFpar0[4][0];

  KMFpar1[0][0] = kmf_.KMFpar1[0][0]; 
  KMFpar1[1][0] = kmf_.KMFpar1[1][0];
  KMFpar1[2][0] = kmf_.KMFpar1[2][0];
  KMFpar1[3][0] = kmf_.KMFpar1[3][0];
  KMFpar1[4][0] = kmf_.KMFpar1[4][0];

  for(i=0;i<PARASIZE;i++){
	for(j=0;j<PARASIZE;j++){
	  Fmat[i][j][0] = kmf_.Fplane[j][i];
	  Cmat1[i][j][0] = Cmats[i][j][0];
	}
  }

  if(kmf_.idwi == 1){
	KMFresi[0] =  meas - KMFpar1[0][0];
  }else if(kmf_.idwi == 2){
	KMFresi[0] =  meas - KMFpar1[1][0];
  }

  kmf_.kmfresi[0] =  KMFresi[0];
  Rerr = pow(0.035,2)- Cmat1[0][0][0];
  if(Rerr < 0){
	kmf_.kmferr[0] = 1.;
  }else{
	kmf_.kmferr[0] = sqrt(Rerr);
  }
  
  /*******
  Rerr = 0.2;
  
  if(Rerr == 0){
	KMFchi2new[0] = 0;
	KMFchi2[0] = 0;
  }else{
	KMFchi2new[0] = pow(KMFresi[0],2)/Rerr;
	KMFchi2[0]    = KMFchi2new[0];
  }
  *******/

	KMFchi2new[0] = 0;
	KMFchi2[0] = 0;

	kmf_.kmfchi2       = KMFchi2[0];
	kmf_.kmfchi2new[0] = KMFchi2new[0];
	kmf_.smtchi2 = 0.;

	return 0;
  
}


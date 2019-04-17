//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h> 
//#include <math.h>
#include "kmf.h"
#include "filter.h"

//int filter_Q();
//int filter_C0();
//int filter_C0_1();
//int filter_C1();
//int filter_par1();
//int filter_resi();
//int filter_chi2();

int i,j,k,idwi;
int plane;

double Cms,reso,meas,measz,zcoor,path;
double c30,s30;
double c45,s45;

int kmfilter_()
{
  iplane = kmf_.iplane - 1;
  
  c30 = cos(30.0*3.14159/180.0);
  s30 = sin(30.0*3.14159/180.0);	
  c45 = cos(45.0*3.14159/180.0);
  s45 = sin(45.0*3.14159/180.0);


  /* track parameter = 5 vector KMFpar0 or KMFpar1  */
  /* x,y,tx,ty,lam */
  /* KMFpar0 = par{k|k-1} */
  KMFpar0[0][iplane] = kmf_.KMFpar0[0][iplane]; 
  KMFpar0[1][iplane] = kmf_.KMFpar0[1][iplane];
  KMFpar0[2][iplane] = kmf_.KMFpar0[2][iplane];
  KMFpar0[3][iplane] = kmf_.KMFpar0[3][iplane];
  KMFpar0[4][iplane] = kmf_.KMFpar0[4][iplane];
  
  Cms   = pow(kmf_.mulsth[iplane],2);
  meas  = kmf_.meas[iplane];
  measz = kmf_.measz[iplane];
  zcoor = kmf_.zcoor[iplane]; 
  path  = kmf_.path[iplane]; 
  idwi  = kmf_.idwi; 
  reso  = kmf_.reso[iplane];
  
  /* track model F{k} */
  for(i=0;i<PARASIZE;i++){
	for(j=0;j<PARASIZE;j++){
	  Fmat[i][j][iplane] = kmf_.Fplane[j][i];
	}
  }
  
  filter_Q(iplane - 1);    /* MS error matrix */
  
  filter_C0(iplane - 1);   /* cov matrix C{k|k-1} */
  
  filter_C0_1(iplane - 1);   /* inverse of cov matrix C{k|k-1} */
  
  filter_C1(iplane - 1);  /* cov matrix C{k|k} */

  filter_par1(iplane - 1);  /*  KMFpar1 = par{k|k} */
  
  filter_resi(iplane - 1); /* residual */
  
  filter_chi2(iplane - 1);  /* chi2 */


  return 0;
}

int filter_Q(int iplane0)
{
  int iplane1;
  double tx,ty,lam;
  double acont,bcont,ccont;
  double Cxx,Cyy,Cxy;
  
  iplane1 = iplane0 + 1;
  
  for(i = 0; i < PARASIZE; i++){
    for(j = 0; j < PARASIZE; j++){
      Qmat[i][j][iplane1] = 0.0;
    }
  }

  tx  = KMFpar1[2][iplane0];
  ty  = KMFpar1[3][iplane0];
  lam = KMFpar1[4][iplane0];

  acont = 1+pow(tx,2)+pow(ty,2);
  bcont = 1+pow(tx,2);
  ccont = 1+pow(ty,2);
  
  Cxx = bcont*acont*Cms;
  Cyy = ccont*acont*Cms;
  Cxy = (tx*ty)*acont*Cms;

  /* thick */
  Qmat[0][0][iplane1] = Cxx*pow(path,2)/3;
  Qmat[0][1][iplane1] = Cxy*pow(path,2)/3;
  Qmat[0][2][iplane1] = Cxx*pow(path,1)/2;
  Qmat[0][3][iplane1] = Cxy*pow(path,1)/2;
  
  Qmat[1][0][iplane1] = Qmat[0][1][iplane1];
  Qmat[1][1][iplane1] = Cyy*pow(path,2)/3;
  Qmat[1][2][iplane1] = Cxy*pow(path,1)/2;
  Qmat[1][3][iplane1] = Cyy*pow(path,1)/2;

  Qmat[2][0][iplane1] = Qmat[0][2][iplane1];
  Qmat[2][1][iplane1] = Qmat[1][2][iplane1];
  Qmat[2][2][iplane1] = Cxx;
  Qmat[2][3][iplane1] = Cxy;

  Qmat[3][0][iplane1] = Qmat[0][3][iplane1];
  Qmat[3][1][iplane1] = Qmat[1][3][iplane1];
  Qmat[3][2][iplane1] = Qmat[2][3][iplane1];
  Qmat[3][3][iplane1] = Cyy;

  /* thin
  Qmat[2][2][iplane1] =  bcont*acont*Cms;
  Qmat[3][3][iplane1] =  ccont*acont*Cms;
  Qmat[2][3][iplane1] =(tx*ty)*acont*Cms;
  Qmat[3][2][iplane1] = Qmat[2][3][iplane1];
  Qmat[4][4][iplane1] = pow(ty,4)*pow(tx,2)/(acont*bcont)*Cms/pow(lam,2);
  */
  /*
  Qmat[4][4][iplane1] = kmf_.eloss[iplane1];
  */
  return 0;
  
}

int filter_C0(int iplane0)
{
  int iplane1;
  double Emat[PARASIZE][PARASIZE][PLANESIZE];
  int aba_t();

  iplane1 = iplane0 + 1;

  aba_t(Fmat,Cmat1,Emat,iplane1,iplane0);

  for(i = 0; i < PARASIZE; i++){
    for(j = 0; j < PARASIZE; j++){
      Cmat0[i][j][iplane1] = Emat[i][j][iplane1] + Qmat[i][j][iplane1];
	}
  }

  return 0;
  
}

int filter_C0_1(int iplane0)
{
  int iplane1;
  int gaussj();
  float vect[PARASIZE];
  double Matrix[6][6];
  double B[PARASIZE][PARASIZE];
  
  iplane1 = iplane0 + 1;

  for(i = 0; i < PARASIZE;i ++){
	for(j = 0; j < PARASIZE; j++){
	  Matrix[i+1][j+1] = Cmat0[i][j][iplane1];
	}
	vect[i] = 0.0;
  }
	  
  gaussj(Matrix,PARASIZE,vect,PARASIZE);
	  
    
  for(i = 0; i < PARASIZE; i++){
	for(j = 0; j < PARASIZE; j++){
	  Cinv0[i][j][iplane1] = Matrix[i+1][j+1];
	}
  }

/*** for test ***/
  /***
  for(i = 0; i < PARASIZE; i++){
	for(j = 0; j < PARASIZE; j++){      
	  B[i][j] = 0.0;
	  for(k = 0; k < PARASIZE; k++){      
		B[i][j] += Cinv0[i][k][iplane1]*Cmat0[k][j][iplane1];
	  }
	}
  }
  ***/
  return 0;
}

int filter_C1(int iplane0)
{
  int iplane1;
  int gaussj();
  float vect[PARASIZE];
  double Matrix[6][6];
  double B[PARASIZE][PARASIZE];
  double C[PARASIZE][PARASIZE];

  iplane1 = iplane0 + 1;

  for(i = 0; i < PARASIZE; i++){
    for(j = 0; j < PARASIZE; j++){
      Matrix[i+1][j+1] = Cinv0[i][j][iplane1];
    }
  }
  
  if(idwi == 1){     /* x */
	Matrix[1][1] += 1/pow(reso,2);
  }else if(idwi == 2){     /* y */
	Matrix[2][2] += 1/pow(reso,2);
  }else if(idwi == 3){     /* u */
	Matrix[1][1] += pow(c30,2)/pow(reso,2);
	Matrix[1][2] +=   -c30*s30/pow(reso,2);
	Matrix[2][1] +=   -c30*s30/pow(reso,2);
	Matrix[2][2] += pow(s30,2)/pow(reso,2);
  }else if(idwi == 4){     /* v */
	Matrix[1][1] += pow(c30,2)/pow(reso,2);
	Matrix[1][2] +=    c30*s30/pow(reso,2);
	Matrix[2][1] +=    c30*s30/pow(reso,2);
	Matrix[2][2] += pow(s30,2)/pow(reso,2);
  }else if(idwi == 5){     /* u */
	Matrix[1][1] += pow(c45,2)/pow(reso,2);
	Matrix[1][2] +=   -c45*s45/pow(reso,2);
	Matrix[2][1] +=   -c45*s45/pow(reso,2);
	Matrix[2][2] += pow(s45,2)/pow(reso,2);
  }else if(idwi == 6){     /* v */
	Matrix[1][1] += pow(c45,2)/pow(reso,2);
	Matrix[1][2] +=    c45*s45/pow(reso,2);
	Matrix[2][1] +=    c45*s45/pow(reso,2);
	Matrix[2][2] += pow(s45,2)/pow(reso,2);
  }else{
	printf("1 error in odd or even plane\n");
  }
  
  for(i = 0; i < PARASIZE; i++){
	for(j = 0; j < PARASIZE; j++){
	  C[i][j] = Matrix[i+1][j+1];
	}
  }

  gaussj(Matrix,PARASIZE,vect,PARASIZE);
	 	  
  for(i = 0; i < PARASIZE; i++){
	for(j = 0; j < PARASIZE; j++){
	  Cmat1[i][j][iplane1] = Matrix[i+1][j+1];
	}
  }
  
  for(i = 0; i < PARASIZE; i++){
	for(j = 0; j < PARASIZE; j++){      
	  B[i][j] = 0.0;
	  for(k = 0; k < PARASIZE; k++){      
		B[i][j] += C[i][k]*Cmat1[k][j][iplane1];
	  }
	}
  }
 
  return 0;
}

int filter_par1(int iplane0)
{
  int iplane1;
  double par[PARASIZE][PLANESIZE];
  double par1[PARASIZE][PLANESIZE];
  double temp;
  int MaxV();
  
  iplane1 = iplane0 + 1;
  for(i = 0; i < PARASIZE; i++){
    for(j = 0; j < PLANESIZE; j++){
	  par[i][j] = 0.0;
	  par1[i][j] = 0.0;
    }
  }
    
  MaxV(Cinv0,KMFpar0,par,iplane1);

  if(idwi == 1){     /* x */
	par[0][iplane1] += meas/pow(reso,2);
  }else if(idwi == 2){     /* y */
	par[1][iplane1] += meas/pow(reso,2);
  }else if(idwi == 3){     /* u */
	par[0][iplane1] +=  c30*meas/pow(reso,2);
	par[1][iplane1] += -s30*meas/pow(reso,2);
  }else if(idwi == 4){     /* v */
	par[0][iplane1] +=  c30*meas/pow(reso,2);
	par[1][iplane1] +=  s30*meas/pow(reso,2);
  }else if(idwi == 5){     /* u */
	par[0][iplane1] +=  c45*meas/pow(reso,2);
	par[1][iplane1] += -s45*meas/pow(reso,2);
  }else if(idwi == 6){     /* v */
	par[0][iplane1] +=  c45*meas/pow(reso,2);
	par[1][iplane1] +=  s45*meas/pow(reso,2);
  }else{
	printf("2 error in odd or even plane\n");
  }

  MaxV(Cmat1,par,par1,iplane1);

  for(i = 0; i < PARASIZE; i++){
    KMFpar1[i][iplane1] = par1[i][iplane1];  
	kmf_.KMFpar1[i][iplane1] = par1[i][iplane1];  
  }

  if(idwi == 1){
	temp = (meas - KMFpar0[0][iplane1])*Cmat0[0][0][iplane1]
	  /(Cmat0[0][0][iplane1]+reso)+KMFpar0[0][iplane1];
	  }

  return 0;
  
}

int filter_resi(int iplane0)
{
   int iplane1;
  
  iplane1 = iplane0 + 1;
  
  if(idwi == 1){     /* x */
	KMFresi[iplane1] = meas - KMFpar1[0][iplane1];
  }else if(idwi == 2){     /* y */
	KMFresi[iplane1] = meas - KMFpar1[1][iplane1];
  }else if(idwi == 3){     /* u */
	KMFresi[iplane1] = meas - 
	  (c30*KMFpar1[0][iplane1] - s30*KMFpar1[1][iplane1]);
  }else if(idwi == 4){     /* v */
	KMFresi[iplane1] = meas - 
	  (c30*KMFpar1[0][iplane1] + s30*KMFpar1[1][iplane1]);
  }else if(idwi == 5){     /* u */
	KMFresi[iplane1] = meas - 
	  (c45*KMFpar1[0][iplane1] - s45*KMFpar1[1][iplane1]);
  }else if(idwi == 6){     /* v */
	KMFresi[iplane1] = meas - 
	  (c45*KMFpar1[0][iplane1] + s45*KMFpar1[1][iplane1]);
  }else{
	printf("4 error in odd or even plane\n");
  }

  kmf_.kmfresi[iplane1] = KMFresi[iplane1];
  
return 0;
}

int filter_chi2(int iplane0)
{
  int iplane1;
  double Rerr;
  
  iplane1 = iplane0 + 1;

  if(idwi == 1){     /* x */
	Rerr = pow(reso,2) - Cmat1[0][0][iplane1];
  }else if(idwi == 2){     /* y */
	Rerr = pow(reso,2) - Cmat1[1][1][iplane1];
  }else if(idwi == 3){     /* u */
    Rerr = pow(reso,2) - 
	      (pow(c30,2)*Cmat1[0][0][iplane1] + 
	       pow(s30,2)*Cmat1[1][1][iplane1] -
	       c30*s30*(Cmat1[0][1][iplane1] + Cmat1[1][0][iplane1]));
  }else if(idwi == 4){     /* v */
    Rerr = pow(reso,2) - 
	      (pow(c30,2)*Cmat1[0][0][iplane1] + 
	       pow(s30,2)*Cmat1[1][1][iplane1] +
	       c30*s30*(Cmat1[0][1][iplane1] + Cmat1[1][0][iplane1]));
  }else if(idwi == 5){     /* u */
    Rerr = pow(reso,2) - 
	      (pow(c45,2)*Cmat1[0][0][iplane1] + 
	       pow(s45,2)*Cmat1[1][1][iplane1] -
	       c45*s45*(Cmat1[0][1][iplane1] + Cmat1[1][0][iplane1]));
  }else if(idwi == 6){     /* v */
    Rerr = pow(reso,2) - 
	      (pow(c45,2)*Cmat1[0][0][iplane1] + 
	       pow(s45,2)*Cmat1[1][1][iplane1] +
	       c45*s45*(Cmat1[0][1][iplane1] + Cmat1[1][0][iplane1]));
  }else{
	printf("error in odd or even plane\n");
  }

  KMFchi2new[iplane1] = pow(KMFresi[iplane1],2)/Rerr;
  KMFchi2[iplane1]    = KMFchi2[iplane0] + KMFchi2new[iplane1];

  kmf_.kmfchi2             = KMFchi2[iplane1];
  kmf_.kmfchi2new[iplane1] = KMFchi2new[iplane1];
  kmf_.kmferr[iplane1]     = sqrt(Rerr);
  
  return 0;
}


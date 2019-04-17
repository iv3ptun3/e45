//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h> 
//#include <math.h>
#include "filter.h"
#include "kmf.h"
int i,j,k;

int ab_tc(INmat1,INmat2,INmat3,OUTmat,plane0,plane1)
double INmat1[PARASIZE][PARASIZE][PLANESIZE];
double INmat2[PARASIZE][PARASIZE][PLANESIZE];
double INmat3[PARASIZE][PARASIZE][PLANESIZE];
double OUTmat[PARASIZE][PARASIZE][PLANESIZE];
int plane0,plane1;
{
  double Xmat[PARASIZE][PARASIZE];

  for(i=0;i<PARASIZE;i++){
    for(j=0;j<PARASIZE;j++){      
	  for(k=0;k<PLANESIZE;k++){      
		OUTmat[i][j][k]=0.;
	  }
	  Xmat[i][j]=0.;
	}
  }
  
  for(i=0;i<PARASIZE;i++){
    for(j=0;j<PARASIZE;j++){
      for(k=0;k<PARASIZE;k++){
	Xmat[i][j]=Xmat[i][j]+INmat1[i][k][plane0]*INmat2[j][k][plane1];
      }
    }
  }
  
  for(i=0;i<PARASIZE;i++){
    for(j=0;j<PARASIZE;j++){
      for(k=0;k<PARASIZE;k++){
	OUTmat[i][j][plane0]=OUTmat[i][j][plane0]+
	  Xmat[i][k]*INmat3[k][j][plane1];
      }
    }
  }
  return 0;
  
  
}

int aba_t(INmat1,INmat2,OUTmat,plane1,plane0)
double INmat1[PARASIZE][PARASIZE][PLANESIZE];
double INmat2[PARASIZE][PARASIZE][PLANESIZE];
double OUTmat[PARASIZE][PARASIZE][PLANESIZE];
int plane0,plane1;
{
  double Xmat[PARASIZE][PARASIZE];

  for(i=0;i<PARASIZE;i++){
    for(j=0;j<PARASIZE;j++){      
	  for(k=0;k<PLANESIZE;k++){      
		OUTmat[i][j][k]=0.;
	  }
	  Xmat[i][j]=0.;
	}
  }


  for(i=0;i<PARASIZE;i++){
    for(j=0;j<PARASIZE;j++){
      for(k=0;k<PARASIZE;k++){
		Xmat[i][j]=Xmat[i][j]+INmat1[i][k][plane1]*INmat2[k][j][plane0];
      }
    }
  }
  
  for(i=0;i<PARASIZE;i++){
    for(j=0;j<PARASIZE;j++){
      for(k=0;k<PARASIZE;k++){
		OUTmat[i][j][plane1]=OUTmat[i][j][plane1]
		  +Xmat[i][k]*INmat1[j][k][plane1];
      }
    }
  }

  return 0;
  
  
}

int MaxV(INmat,INvec,OUTvec,plane)
double INmat[PARASIZE][PARASIZE][PLANESIZE];
double INvec[PARASIZE][PLANESIZE];
double OUTvec[PARASIZE][PLANESIZE];
int plane;
{

  for(i=0;i<PARASIZE;i++){
	for(k=0;k<PLANESIZE;k++){
	  OUTvec[i][k]=0.;
	}
  }

  for(i=0;i<PARASIZE;i++){
    for(k=0;k<PARASIZE;k++){
      OUTvec[i][plane]=OUTvec[i][plane]
		+INmat[i][k][plane]*INvec[k][plane];
    }
  }
  
  return 0;
  
}

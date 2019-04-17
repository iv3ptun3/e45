/* MINUIT functions
 *						by tsugu
 */

#include <stdio.h>
#include <string.h>
#include "minuit.hh"


#define min( a, b ) ( ( a < b ) ? a : b )

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* CERN Library functions written in FORTRAN */
void mintio_( int*, int*, int* );
void mncomd_( MinuitFCN, char*, int*, MinuitFUtil, int );
void mncont_( MinuitFCN, int*, int*, int*,
	double*, double*, int*, MinuitFUtil );
void mnemat_( double*, int* );
void mnerrs_( int*, double*, double*, double*, double* );
void mnexcm_( MinuitFCN, char*, double*, int*, int*, MinuitFUtil, int );
void mninpu_( int*, int* );
void mnintr_( MinuitFCN, MinuitFUtil );
void mnpars_( char*, int*, int );
void mnpout_( int*, char*, double*, double*, double*, double*, int*, int );
void mnstat_( double*, double*, double*, int*, int*, int* );
void mninit_( int*, int*, int* );
void mnparm_( int*, char*, double*, double*, double*, double*, int*, int );
void mnseti_( char*, int );

void MinuitChopString( char* str, int len ){
  if ( ( len < 1 ) || ( NULL == str ) ) return;
  while ( len -- ) if ( FORTRAN_DEFAULT_CHAR != str[len] ) break;
  str[ ++ len ] = ( char )( NULL );
}

/* Functions call original CERN Library functions */
void MINTIO( int iread, int iwrite, int isave ){
  mintio_( & iread, & iwrite, & isave );
}

void MNCOMD( MinuitFCN fcn, char* ichstr, int* icondn, MinuitFUtil futil ){
  int lchstr = min( MAX_CH_LEN, strlen( ichstr ) );
  char chstr[MAX_CH_LEN];
  strncpy( chstr, ichstr, lchstr );
  mncomd_( fcn, chstr, icondn, futil, lchstr );
}

void MNCONT( MinuitFCN fcn, int num1, int num2, int npt,
	double* xpt, double* ypt, int* nfound, MinuitFUtil futil ){
  mncont_( fcn, & num1, & num2, & npt, xpt, ypt, nfound, futil );
}

void MNEMAT( double* emat, int ndim ){ mnemat_( emat, & ndim ); }

void MNERRS( int num, double* eplus,
	double* eminus, double* eparab, double* globcc ){
  mnerrs_( & num, eplus, eminus, eparab, globcc );
}

void MNEXCM( MinuitFCN fcn, char* ichcom,
	double* arglis, int narg, int* ierflg, MinuitFUtil futil ){
  int lchcom = min( MAX_CH_LEN, strlen( ichcom ) );
  char chcom[MAX_CH_LEN];
  strncpy( chcom, ichcom, lchcom );
  mnexcm_( fcn, chcom, arglis, & narg, ierflg, futil, lchcom );
}

void MNINIT( int ird, int iwr, int isav ){
  mninit_( & ird, & iwr, & isav );
}

void MNINPU( int nunit, int* ierr ){ mninpu_( & nunit, ierr ); }
void MNINTR( MinuitFCN fcn, MinuitFUtil futil ){ mnintr_( fcn, futil ); }

void MNPARM( int num, char* ichnam,
	double stval, double step, double bnd1, double bnd2, int* ierflg ){
  int lchnam = min( MAX_CH_LEN, strlen( ichnam ) );
  char chnam[MAX_CH_LEN];
  strncpy( chnam, ichnam, lchnam );
  mnparm_( & num, chnam, & stval, & step, & bnd1, & bnd2, ierflg, lchnam );
}

void MNPARS( char* ichstr, int* icondn ){
  int lchstr = min( MAX_CH_LEN, strlen( ichstr ) );
  char chstr[MAX_CH_LEN];
  strncpy( chstr, ichstr, lchstr );
  mnpars_( chstr, icondn, lchstr );
}

void MNPOUT( int num, char* chnam, double* val,
	double* error, double* bnd1, double* bnd2, int* ivarbl ){
  mnpout_( & num, chnam, val, error, bnd1, bnd2, ivarbl, MINUIT_VAR_LEN );
  MinuitChopString( chnam, MINUIT_VAR_LEN );
  chnam[MINUIT_VAR_LEN] = ( char )( NULL );
}

void MNSETI( char* ictitle ){
  int lctitle = min( MAX_CH_LEN, strlen( ictitle ) );
  char ctitle[MAX_CH_LEN];
  strncpy( ctitle, ictitle, lctitle );
  mnseti_( ctitle, lctitle );
}

void MNSTAT( double* fmin, double* fedm,
	double* errdef, int* npari, int* nparx, int* istat ){
  mnstat_( fmin, fedm, errdef, npari, nparx, istat );
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

/* End of file */


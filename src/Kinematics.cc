// -*- C++ -*-

#include "Kinematics.hh"

#include <Randomize.hh>
#include <cmath>
#include <cstdio>

#include <TF2.h>
#include "FuncName.hh"
#include <TMath.h>

namespace Kinematics
{
const G4double b = 1.62 * CLHEP::GeV * CLHEP::fermi / CLHEP::hbarc;
static int gNumOfTracks;
static std::vector<double> gX0;
static std::vector<double> gY0;
static std::vector<double> gU0;
static std::vector<double> gV0;
static void fcn_vertex(int &npar, double *gin, double &f, double *par, int iflag){
  double chisqr=0.;
  for(int i=0; i<gNumOfTracks; ++i){
    chisqr += TMath::Power(par[0]-gX0[i]-gU0[i]*par[2], 2);
    chisqr += TMath::Power(par[1]-gY0[i]-gV0[i]*par[2], 2);
  };
  f = chisqr;
};


//______________________________________________________________________________
G4ThreeVector
HarmonicFermiMomentum( G4int type )
{
  /*      THIS ROUTINE GENERATES FERMI MOMENTUM KF(3) BY HARMONIC */
  /*      OSCILATOR MODEL FOR 1S AND 1P */
  /*      INPUT; L=0 OR 1 */
  /*      OUTPUT; KF(3) */
  ///what is unit?? mev? gev?
  ////in the g3LEPS, b=1.64 (fm)
  ////in the g3LEPS, b=1.73 hicks
  ////in the g3LEPS, b=1.93 fitting
  ////hbar = 197.327053 MeV*fm --. 0.197 --> GeV*fm

  G4double x = 0.;
  switch( type ){
  case 0: {
    static const G4double ymax = std::exp(-1);
    while( true ){
      x = G4RandFlat::shoot();
      G4double y = x * x * std::exp(-b * b * x * x);
      G4double yy = G4RandFlat::shoot() * ymax;
      // G4cout << "iterate " << x << " " << b << " "
      // 	     << y << " " << yy << G4endl;
      if ( yy < y ) break;
    }
    break;
  }
  case 1: {
    static const G4double ymax = std::exp(-2) * 4 / (b * b);
    while( true ){
      x = G4RandFlat::shoot();
      G4double y = b * b * x * x * x * x * std::exp(-b * b * x * x);
      G4double yy = G4RandFlat::shoot() * ymax;
      // G4cout << "iterate " << x << " " << b << " "
      // 	     << y << " " << yy << G4endl;
      if ( yy < y ) break;
    }
    break;
  }
  default:
    G4Exception( FUNC_NAME, " ", RunMustBeAborted, "" );
    break;
  }

  G4double theta = std::acos( G4RandFlat::shoot( -1., 1. ) );
  G4double phi   = G4RandFlat::shoot( -CLHEP::pi, CLHEP::pi );

  return G4ThreeVector( x * std::sin(theta) * std::cos(phi),
			x * std::sin(theta) * std::sin(phi),
			x * std::cos(theta) );
}

//______________________________________________________________________________
G4int
HarmonicFermiMomentumDeuteron( G4double* Kf )
{
  //deut_pro : 0: deuteron
  //deut_pro : 1: proton
  //// for deuteron E27
  /////// Copy of H. Takahashi-san's code
  /////// revised to geant4 by S.Hwang
  G4double ymax;
  G4double x, y, yy, theta, phi;
  /* THIS ROUTINE GENERATES FERMI MOMENTUM KF(3) BY HARMONIC */
  /* OSCILATOR MODEL FOR 1S AND 1P */
  /* INPUT; L=0 OR 1 */
  /* OUTPUT; KF(3) */
  ///what is unit?? mev? gev?
  ////in the g3LEPS, b=1.64 (fm)
  ////in the g3LEPS, b=1.73 hicks
  ////in the g3LEPS, b=1.93 fitting
  // G4double b = 1.62;
  // G4double b = 1.93;
  // G4double b = 2.1;

  ymax = 20.;
  do {
    x = G4RandFlat::shoot();
    y = (x*x) * (13.* std::exp(-8.*x)+30000.*std::exp(-37.*x));
    yy = G4RandFlat::shoot() * ymax;
  } while (yy > y);

  theta = std::acos( G4RandFlat::shoot( -1., 1. ) );
  phi   = ( G4RandFlat::shoot( -1., 1. ) )*CLHEP::pi;
  // IsotropicAngle(&theta, &phi);
  Kf[0] = x * std::sin(theta) * std::cos(phi);
  Kf[1] = x * std::sin(theta) * std::sin(phi);
  Kf[2] = x * std::cos(theta);
  // G4cout<<x<<":"<<Kf[0]<<":"<<Kf[1]<<":"<<Kf[2]<<":"<<G4endl;
  return 0;
}

//______________________________________________________________________________
G4double
Legendre( G4int order, G4double x )
{
  switch( order ){
  case 0:
    return 1.;
  case 1:
    return x;
  case 2:
    return 1./2.*( 3.*x*x - 1. );
  case 3:
    return 1./2.*( 5.*x*x*x - 3.*x );
  case 4:
    return 1./8.*( 35.*x*x*x*x - 30.*x*x + 3. );
  case 5:
    return 1./8.*( 63.*x*x*x*x*x - 70.*x*x*x + 15.*x );
  case 6:
    return 1./16.*( 231.*x*x*x*x*x*x - 315.*x*x*x*x + 105.*x*x - 5. );
  case 7:
    return 1./16.*( 429.*x*x*x*x*x*x*x -
		    693.*x*x*x*x*x +
		    315.*x*x*x -
		     35.*x );
  case 8:
    return 1./128.*(  6435.*x*x*x*x*x*x*x*x -
		     12012.*x*x*x*x*x*x +
		      6930.*x*x*x*x -
		      1260.*x*x +
		        35. );
  case 9:
    return 1./128.*( 12155.*x*x*x*x*x*x*x*x*x -
		     25740.*x*x*x*x*x*x*x +
		     18018.*x*x*x*x*x -
		      4620.*x*x*x +
		       315.*x );
  case 10:
    return 1./256.*(  46189.*x*x*x*x*x*x*x*x*x*x -
		     109395.*x*x*x*x*x*x*x*x +
		      90090.*x*x*x*x*x*x -
		      30030.*x*x*x*x +
		       3465.*x*x -
		         63. );
  case 11:
    return 1./256.*(  88179.*x*x*x*x*x*x*x*x*x*x*x -
		     230945.*x*x*x*x*x*x*x*x*x +
		     218790.*x*x*x*x*x*x*x -
		      90090.*x*x*x*x*x +
		      15015.*x*x*x -
		        693.*x );
  case 12:
    return 1./1024.*(  676039.*x*x*x*x*x*x*x*x*x*x*x*x -
		      1939938.*x*x*x*x*x*x*x*x*x*x +
		      2078505.*x*x*x*x*x*x*x*x -
		      1021020.*x*x*x*x*x*x +
		       225225.*x*x*x*x -
		        18018.*x*x +
		          231. );
  default:
    G4cout << "#E " << FUNC_NAME << " invalid order : " << order << G4endl;
    return 0.;
  }
}

G4ThreeVector
MultitrackVertex(int ntrack, double *x0, double *y0, double *u0, double *v0){

  gNumOfTracks = ntrack;
  gX0.clear();
  gY0.clear();
  gU0.clear();
  gV0.clear();
  for(int i=0; i<gNumOfTracks; ++i){
    gX0.push_back(x0[i]);
    gY0.push_back(y0[i]);
    gU0.push_back(u0[i]);
    gV0.push_back(v0[i]);
  }

  double par[3] = {0, 0, 0};
  double err[3] = {999., 999., 999.};

  TMinuit *minuit = new TMinuit(3);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn_vertex);

  double arglist[10];
  int ierflg = 0;

  arglist[0] = 1; //error level for ch2 minimization
  minuit->mnexcm("SET ERR", arglist, 1, ierflg);
  minuit->mnexcm("SET NOW", arglist, 1, ierflg);

  TString name[3] = {"x", "y" ,"z"};
  const double FitStep[3] = {0.001, 0.001, 0.001};
  const double LowLimit[3] = {-50., -50., -100};
  const double UpLimit[3] = {50., 50., 100};
  for(int i=0; i<3; i++){
    minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }
  minuit->Command("SET STRategy 0");

  //arglist[0] = 500.;
  //arglist[1] = 1;
  arglist[0] = 1000.;
  arglist[1] = 0.1;

  int Err;
  double bnd1, bnd2;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  //minuit->mnimpr();
  //minuit->mnexcm("MINOS", arglist, 0, ierflg);
  //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

  double amin, edm, errdef;
  int nvpar, nparx, icstat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  for(int i=0; i<3; i++){
    minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
  }
  delete minuit;

  return G4ThreeVector(par[0], par[1], par[2] - 143);
}


G4ThreeVector SphericalRandom(){
  G4double theta = std::acos( G4RandFlat::shoot( -1., 1. ) );
  G4double phi   = G4RandFlat::shoot( -CLHEP::pi, CLHEP::pi );
  return G4ThreeVector( std::sin(theta) * std::cos(phi),
      std::sin(theta) * std::sin(phi),
      std::cos(theta) );
}
G4ThreeVector FermiGasMomentum(double p_f){
  double u = G4RandFlat::shoot(0.,1.);
  u = TMath::Power(u, 1./3.);
  double p = p_f * u;
  G4ThreeVector dir = SphericalRandom();
  return p * dir;
}

G4ThreeVector RotateAlongBeam(G4ThreeVector Beam, G4ThreeVector Vect, G4double phi){
	G4double BeamMag = Beam.mag();
	auto z_Beam = Beam* (1./BeamMag);
	auto y_Beam = z_Beam.cross(Vect);
	y_Beam = y_Beam * (1./y_Beam.mag());
	auto x_Beam = y_Beam.cross(z_Beam);
	
	auto Pz_Beam = Vect * z_Beam;
	auto Pt_Beam = Vect * x_Beam;//X axis is defined as the transverse component direction. 

	auto x_rot = x_Beam * cos(phi) + y_Beam * sin(phi);
	auto y_rot = -x_Beam * sin(phi) + y_Beam * cos(phi);

	auto Pt_rot = Pt_Beam * x_rot; 
	auto Pz_rot = Pz_Beam * z_Beam;

	auto P_rot = Pt_rot + Pz_rot;

	return P_rot;
}

}

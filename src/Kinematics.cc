// -*- C++ -*-

#include "Kinematics.hh"

#include <Randomize.hh>

#include "FuncName.hh"

namespace Kinematics
{

const G4double b = 1.62 * CLHEP::GeV * CLHEP::fermi / CLHEP::hbarc;


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

}

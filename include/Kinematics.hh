// -*- C++ -*-

#ifndef KINEMATICS_HH
#define KINEMATICS_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include "TMinuit.h"

namespace Kinematics
{
  G4ThreeVector HarmonicFermiMomentum( G4int type );
  G4int         HarmonicFermiMomentumDeuteron( G4double* Kf );
  G4double      Legendre( G4int order, G4double x );
  G4ThreeVector SphericalRandom( void );
  G4ThreeVector FermiGasMomentum(double p_f = 250);//MeV
	G4ThreeVector MultitrackVertex(int ntrack, double *x0, double *y0,
			    double *u0, double *v0);

  inline G4String ClassName( void )
  {
    static G4String s_name("Kinematics");
    return s_name;
  }
}

#endif

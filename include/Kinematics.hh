// -*- C++ -*-

#ifndef KINEMATICS_HH
#define KINEMATICS_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>

namespace Kinematics
{
  G4ThreeVector HarmonicFermiMomentum( G4int type );
  G4int         HarmonicFermiMomentumDeuteron( G4double* Kf );
  G4double      Legendre( G4int order, G4double x );

  inline G4String ClassName( void )
  {
    static G4String s_name("Kinematics");
    return s_name;
  }
}

#endif

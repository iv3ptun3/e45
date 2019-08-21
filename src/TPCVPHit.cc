// -*- C++ -*-

#include "TPCVPHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCVPHit> TPCVPHitAllocator;

//_____________________________________________________________________________
TPCVPHit::TPCVPHit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCVPHit::~TPCVPHit( void )
{
}

//_____________________________________________________________________________
void
TPCVPHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCVPHit::Print( void )
{
  VHitInfo::Print();
}

// -*- C++ -*-

#include "TPCTargetVPHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCTargetVPHit> TPCTargetVPHitAllocator;

//_____________________________________________________________________________
TPCTargetVPHit::TPCTargetVPHit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCTargetVPHit::~TPCTargetVPHit( void )
{
}

//_____________________________________________________________________________
void
TPCTargetVPHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCTargetVPHit::Print( void )
{
  VHitInfo::Print();
}

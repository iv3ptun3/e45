// -*- C++ -*-

#include "TPCVC1Hit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCVC1Hit> TPCVC1HitAllocator;

//_____________________________________________________________________________
TPCVC1Hit::TPCVC1Hit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCVC1Hit::~TPCVC1Hit( void )
{
}

//_____________________________________________________________________________
void
TPCVC1Hit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCVC1Hit::Print( void )
{
  VHitInfo::Print();
}

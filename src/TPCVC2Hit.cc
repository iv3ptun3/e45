// -*- C++ -*-

#include "TPCVC2Hit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCVC2Hit> TPCVC2HitAllocator;

//_____________________________________________________________________________
TPCVC2Hit::TPCVC2Hit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCVC2Hit::~TPCVC2Hit( void )
{
}

//_____________________________________________________________________________
void
TPCVC2Hit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCVC2Hit::Print( void )
{
  VHitInfo::Print();
}

// -*- C++ -*-

#include "TPCBH2Hit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCBH2Hit> TPCBH2HitAllocator;

//_____________________________________________________________________________
TPCBH2Hit::TPCBH2Hit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCBH2Hit::~TPCBH2Hit( void )
{
}

//_____________________________________________________________________________
void
TPCBH2Hit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCBH2Hit::Print( void )
{
  VHitInfo::Print();
}

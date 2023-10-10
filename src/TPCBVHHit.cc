// -*- C++ -*-

#include "TPCBVHHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCBVHHit> TPCBVHHitAllocator;

//_____________________________________________________________________________
TPCBVHHit::TPCBVHHit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCBVHHit::~TPCBVHHit( void )
{
}

//_____________________________________________________________________________
void
TPCBVHHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCBVHHit::Print( void )
{
  VHitInfo::Print();
}

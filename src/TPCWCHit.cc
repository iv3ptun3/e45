// -*- C++ -*-

#include "TPCWCHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCWCHit> TPCWCHitAllocator;

//_____________________________________________________________________________
TPCWCHit::TPCWCHit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCWCHit::~TPCWCHit( void )
{
}

//_____________________________________________________________________________
void
TPCWCHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCWCHit::Print( void )
{
  VHitInfo::Print();
}

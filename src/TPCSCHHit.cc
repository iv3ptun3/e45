// -*- C++ -*-

#include "TPCSCHHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCSCHHit> TPCSCHHitAllocator;

//_____________________________________________________________________________
TPCSCHHit::TPCSCHHit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCSCHHit::~TPCSCHHit( void )
{
}

//_____________________________________________________________________________
void
TPCSCHHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCSCHHit::Print( void )
{
  VHitInfo::Print();
}

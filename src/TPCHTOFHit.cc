// -*- C++ -*-

#include "TPCHTOFHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCHTOFHit> TPCHTOFHitAllocator;

//_____________________________________________________________________________
TPCHTOFHit::TPCHTOFHit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCHTOFHit::~TPCHTOFHit( void )
{
}

//_____________________________________________________________________________
void
TPCHTOFHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCHTOFHit::Print( void )
{
  VHitInfo::Print();
}

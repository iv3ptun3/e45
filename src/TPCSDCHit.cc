// -*- C++ -*-

#include "TPCSDCHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCSDCHit> TPCSDCHitAllocator;

//_____________________________________________________________________________
TPCSDCHit::TPCSDCHit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCSDCHit::~TPCSDCHit( void )
{
}

//_____________________________________________________________________________
void
TPCSDCHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCSDCHit::Print( void )
{
  VHitInfo::Print();
}

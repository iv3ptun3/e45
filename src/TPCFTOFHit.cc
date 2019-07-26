// -*- C++ -*-

#include "TPCFTOFHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCFTOFHit> TPCFTOFHitAllocator;

//_____________________________________________________________________________
TPCFTOFHit::TPCFTOFHit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCFTOFHit::~TPCFTOFHit( void )
{
}

//_____________________________________________________________________________
void
TPCFTOFHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCFTOFHit::Print( void )
{
  VHitInfo::Print();
}

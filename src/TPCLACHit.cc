// -*- C++ -*-

#include "TPCLACHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TPCLACHit> TPCLACHitAllocator;

//_____________________________________________________________________________
TPCLACHit::TPCLACHit( const G4String& name, G4Step* step )
  : G4VHit(),
    VHitInfo( name, step )
{
}

//_____________________________________________________________________________
TPCLACHit::~TPCLACHit( void )
{
}

//_____________________________________________________________________________
void
TPCLACHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCLACHit::Print( void )
{
  VHitInfo::Print();
}

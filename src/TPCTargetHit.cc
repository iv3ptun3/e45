// -*- C++ -*-

#include "TPCTargetHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>

G4Allocator<TPCTargetHit> TPCTargetHitAllocator;

//_____________________________________________________________________________
TPCTargetHit::TPCTargetHit( G4ThreeVector& axyz, G4ThreeVector& apxyz,
			    G4int tid, G4int pid, G4int Parentid,
			    G4double Mass, G4int Charge, G4ThreeVector& Vtxpos,
			    G4ThreeVector& Vtxmom, G4double avtxene,
			    G4double Kinene )
  : xyz(axyz), pxyz(apxyz), trackID(tid), particleID(pid),
    parentID(Parentid),
    mass(Mass), charge(Charge), vtxposi(Vtxpos),
    vtxmome(Vtxmom), vtxene(avtxene), kinene(Kinene)
{
}

//_____________________________________________________________________________
TPCTargetHit::~TPCTargetHit( void )
{
}

//_____________________________________________________________________________
void
TPCTargetHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCTargetHit::Print( void )
{
  G4cout << "Hit in Counter:" << xyz*(1./CLHEP::cm) << " cm, "
	 << tof/CLHEP::ns << " ns" << G4endl;
}

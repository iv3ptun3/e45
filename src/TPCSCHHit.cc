// -*- C++ -*-

#include "TPCSCHHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>

G4Allocator<TPCSCHHit> TPCSCHHitAllocator;

//_____________________________________________________________________________
TPCSCHHit::TPCSCHHit( void )
  : xyz(0.,0.,0.), pxyz(0.,0.,0.), tof(0.)
{
}

//_____________________________________________________________________________
TPCSCHHit::TPCSCHHit( G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t )
  : xyz(axyz), pxyz(apxyz), tof(t)
{
}

//_____________________________________________________________________________
TPCSCHHit::TPCSCHHit( G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
		    G4int tid, G4int pid, G4int did,G4double mass, G4int qq,
		    G4int parentid,
		    G4ThreeVector& Vtxpos, G4ThreeVector& Vtxmom,
		    G4double avtxene, G4double tlength )
  : xyz(axyz), pxyz(apxyz), tof(t), trackID(tid), particleID(pid),
    detectorID(did), massSH(mass), qqSH(qq), parentID(parentid),
    vtxmome(Vtxmom), vtxposi(Vtxpos), vtxene(avtxene), flength(tlength)
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
  G4cout << "Hit in Counter:" << xyz*(1./CLHEP::cm) << " cm, "
	 << tof/CLHEP::ns << " ns" << G4endl;
}

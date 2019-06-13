// ====================================================================
//   TPCACHit.cc
//
// ====================================================================
#include "TPCACHit.hh"

// allocator
G4Allocator<TPCACHit> TPCACHitAllocator;

//////////////////////////////
TPCACHit::TPCACHit()
  : xyz(0.,0.,0.), pxyz(0.,0.,0.), tof(0.)
//////////////////////////////
{
}

/////////////////////////////////////////////////////////////
TPCACHit::TPCACHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t)
  : xyz(axyz), pxyz(apxyz), tof(t)
/////////////////////////////////////////////////////////////
{
}

/////////////////////////////////////////////////////////////
TPCACHit::TPCACHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
			 G4int tid, G4int pid, G4int did,G4double mass, G4int qq, G4int parentid, G4ThreeVector& Vtxpos, G4ThreeVector& Vtxmom, G4double avtxene, G4double tlength)
  : xyz(axyz), pxyz(apxyz), tof(t), trackID(tid), particleID(pid),detectorID(did),massSH(mass),qqSH(qq), parentID(parentid),vtxmome(Vtxmom), vtxposi(Vtxpos), vtxene(avtxene), flength(tlength)
/////////////////////////////////////////////////////////////
{
}

///////////////////////////////
TPCACHit::~TPCACHit()
///////////////////////////////
{
}

//////////////////////////
void TPCACHit::Draw()
//////////////////////////
{
}

///////////////////////////
void TPCACHit::Print()
///////////////////////////
{
  // (global cordinate!!)
  G4cout << "Hit in Counter:" << xyz*(1./cm) << " cm, " 
	 << tof/ns << " ns" << G4endl;
}



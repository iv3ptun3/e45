// ====================================================================
//   TPCTargetHit.cc
//
// ====================================================================
#include "TPCTargetHit.hh"

// allocator
G4Allocator<TPCTargetHit> TPCTargetHitAllocator;
//TPCTargetHit(pos, mom, tid, pid, parentID, mass, charge, VertexPosition, VertexMomentum, VertexEnergy, kinEnergy);
TPCTargetHit::TPCTargetHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4int tid, 
			   G4int pid, G4int Parentid,G4double Mass, G4int Charge,
			   G4ThreeVector& Vtxpos, G4ThreeVector& Vtxmom, G4double avtxene, G4double Kinene)
  : xyz(axyz), pxyz(apxyz), trackID(tid), particleID(pid),
    parentID(Parentid),
    mass(Mass), charge(Charge), vtxposi(Vtxpos), vtxmome(Vtxmom), vtxene(avtxene), kinene(Kinene)
{
}


TPCTargetHit::~TPCTargetHit()
{
}


void TPCTargetHit::Draw()
{
}

void TPCTargetHit::Print()
{
  // (global cordinate!!)
  G4cout << "Hit in Counter:" << xyz*(1./cm) << " cm, " 
	 << tof/ns << " ns" << G4endl;
}



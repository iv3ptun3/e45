// ====================================================================
//   TPCPadHit.cc
//
// ====================================================================
#include "TPCPadHit.hh"

// allocator
G4Allocator<TPCPadHit> TPCPadHitAllocator;

//TPCPadHit::TPCPadHit()
//  : xyz(0.,0.,0.), pxyz(0.,0.,0.), tof(0.),vtxpos(0.,0.,0.),vtxmom(0.,0.,0.), edep(0.), parentID(0.),tLength(0.),beta(0.)
//{
//}
/*
TPCPadHit::TPCPadHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t)
  : xyz(axyz), pxyz(apxyz), tof(t)
{
}

TPCPadHit::TPCPadHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
		     G4int tid, G4int pid)
  : xyz(axyz), pxyz(apxyz), tof(t), trackID(tid), particleID(pid)
{
}


TPCPadHit::TPCPadHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
		     G4int tid, G4int pid,
		     G4int ilay, G4int irow)
  : xyz(axyz), pxyz(apxyz), tof(t), trackID(tid), particleID(pid),
    iLay(ilay),iRow(irow)
{
}
*/
/*
TPCPadHit::TPCPadHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
		     G4int tid, G4int pid,
		     G4int ilay, G4int irow, 
		     G4double b, G4double ed, G4int parentid,
		     G4double tlength, G4double Mass, G4int Charge)
  : xyz(axyz), pxyz(apxyz), tof(t), trackID(tid), particleID(pid),
    iLay(ilay),iRow(irow), beta(b),edep(ed), parentID(parentid),
    tLength(tlength), mass(Mass), charge(Charge)
{
}
*/

TPCPadHit::TPCPadHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
		     G4int tid, G4int pid,
		     G4int ilay, G4int irow, G4double Beta, G4double Ed, G4int Parentid,
		     G4double tlength, G4double Mass, G4int Charge,
		     G4ThreeVector& Vtxpos, G4ThreeVector& Vtxmom, G4double avtxene, G4double slength,
		     G4int Parentid_pid)
  : xyz(axyz), pxyz(apxyz), tof(t), trackID(tid), particleID(pid),
    iLay(ilay),iRow(irow), beta(Beta),edep(Ed), parentID(Parentid),
    Length(tlength), mass(Mass), charge(Charge), vtxmome(Vtxmom), vtxposi(Vtxpos), vtxene(avtxene), SLength(slength), parentID_pid(Parentid_pid)
{
}


TPCPadHit::~TPCPadHit()
{
}


void TPCPadHit::Draw()
{
}

void TPCPadHit::Print()
{
  // (global cordinate!!)
  G4cout << "Hit in Counter:" << xyz*(1./cm) << " cm, " 
	 << tof/ns << " ns" << G4endl;
}



// ====================================================================
//   TPCTargetHit.hh
//
// ====================================================================
#ifndef TPC_TARGET_HIT_H
#define TPC_TARGET_HIT_H
 
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class TPCTargetHit : public G4VHit {
private:
  G4ThreeVector xyz;
  G4ThreeVector pxyz;
  G4double tof;
  G4int trackID;
  G4int particleID;
  G4int iLay;
  G4int iRow;
  G4double beta;
  G4double edep; // Energy deposit
  G4int parentID;
  G4double Length;
  G4double mass;
  G4int charge;
  G4ThreeVector vtxmome;
  G4ThreeVector vtxposi;
  G4double vtxene;
  G4double kinene;
public:
  //  TPCTargetHit();
  TPCTargetHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4int tid, 
	       G4int pid, G4int Parentid,G4double Mass, G4int Charge,
	       G4ThreeVector& Vtxpos, G4ThreeVector& Vtxmom, 
	       G4double avtxene, G4double Kinene);
  TPCTargetHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
	    G4int tid, G4int pid, 
	    G4int ilay, G4int irow, 
	    G4double b, G4double ed, G4int parentid, G4double Length, 
	    G4double Mass, G4int Charge,
	       G4ThreeVector& vtxmom, G4ThreeVector& vtxpos, G4double vtxene, G4double Kinene);
  
  virtual ~TPCTargetHit();
  
  // copy constructor & assignment operator
  TPCTargetHit(const TPCTargetHit& right);
  const TPCTargetHit& operator=(const TPCTargetHit& right);
  
  // new/delete operators
  void* operator new(size_t);
  void operator delete(void* aHit);
  
  // set/get functions
  const G4ThreeVector& GetPosition() const { return xyz; }
  const G4ThreeVector& GetMomentum() const { return pxyz; }

  const G4ThreeVector& GetVtxPosition() const { return vtxposi; }
  const G4ThreeVector& GetVtxMomentum() const { return vtxmome; }
  G4double GetVtxEnergy() const { return vtxene; }
  G4double GetKinEnergy() const { return kinene; }

  G4double GetTOF() const { return tof; }
  G4int GetTrackID() const { return trackID; }
  G4int GetParticleID() const { return particleID; }
  //  G4int GetTargetLay() const { return iLay; }
  //  G4int GetTargetRow() const { return iRow; }
  void SetEdep(G4double aedep) { edep = aedep; }
  G4double GetEdep() const { return edep; }
  G4double GetBeta() const { return beta; }
  G4double GetMass() const { return mass; }
  G4int GetCharge() const { return charge; }
  G4int GetParentID() const { return parentID; }
  G4double GettLength() const { return Length; }


  // methods
  virtual void Draw();
  virtual void Print();
};

// ====================================================================
// inline functions
// ====================================================================
inline TPCTargetHit::TPCTargetHit(const TPCTargetHit& right)
  : G4VHit()
{
  xyz= right.xyz;
  pxyz= right.pxyz;
  tof= right.tof;
  edep = right.edep;
}

inline const TPCTargetHit& TPCTargetHit::operator=
(const TPCTargetHit& right)
{
  xyz= right.xyz;
  pxyz= right.pxyz;
  tof= right.tof;
  edep = right.edep;
  return *this;
}

// externally instanciated.
extern G4Allocator<TPCTargetHit> TPCTargetHitAllocator; 

inline void* TPCTargetHit::operator new(size_t)
{
  void* aHit= (void*)TPCTargetHitAllocator.MallocSingle();
  return aHit;
}

inline void TPCTargetHit::operator delete(void* aHit)
{
  TPCTargetHitAllocator.FreeSingle((TPCTargetHit*) aHit);
}

#endif

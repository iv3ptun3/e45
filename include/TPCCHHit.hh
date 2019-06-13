// ====================================================================
//   TPCCHHit.hh
//
// ====================================================================
#ifndef TPC_CH_HIT_H
#define TPC_CH_HIT_H
 
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class TPCCHHit : public G4VHit {

private:
  G4ThreeVector xyz;
  G4ThreeVector pxyz;
  G4double tof;
  G4double edep; 
  G4int trackID;
  G4int particleID;
  G4int detectorID;
  G4double massSH;
  G4double qqSH;
  G4int parentID;
  G4ThreeVector vtxmome;
  G4ThreeVector vtxposi;
  G4double vtxene;
  G4double flength;
public:
  TPCCHHit();
  //  TPCCHHit(G4ThreeVector& axyz, G4double t);
  TPCCHHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t);
  TPCCHHit(G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
	      G4int tid, G4int pid, G4int did,G4double mass, G4int qq, G4int parentid,
	      G4ThreeVector& vtxmom, G4ThreeVector& vtxpos, G4double vtxene, G4double tlength);

  virtual ~TPCCHHit();

  // copy constructor & assignment operator
  TPCCHHit(const TPCCHHit& right);
  const TPCCHHit& operator=(const TPCCHHit& right);
  
  // new/delete operators
  void* operator new(size_t);
  void operator delete(void* aHit);
  
  const G4ThreeVector& GetVtxPosition() const { return vtxposi; }
  const G4ThreeVector& GetVtxMomentum() const { return vtxmome; }
  G4double GetVtxEnergy() const { return vtxene; }

  // set/get functions
  const G4ThreeVector& GetPosition() const { return xyz; }
  const G4ThreeVector& GetMomentum() const { return pxyz; }
  G4double GetTOF() const { return tof; }

  G4int GetParentID() const { return parentID; }
  G4int GetTrackID() const { return trackID; }
  G4int GetParticleID() const { return particleID; }
  G4int GetDetectorID() const { return detectorID; }
  G4double GetParticleMassID() const { return massSH; }
  G4double GetParticleQqID() const { return qqSH; }
  void SetEdep(G4double aedep) { edep = aedep; }
  G4double GetEdep() const { return edep; }
  G4double GetLength() const { return flength; }

  // methods
  virtual void Draw();
  virtual void Print();
};

// ====================================================================
// inline functions
// ====================================================================
inline TPCCHHit::TPCCHHit(const TPCCHHit& right)
  : G4VHit()
{
  xyz= right.xyz;
  pxyz= right.pxyz;
  tof= right.tof;
  edep = right.edep;
}

inline const TPCCHHit& TPCCHHit::operator=
(const TPCCHHit& right)
{
  xyz= right.xyz;
  pxyz= right.pxyz;
  tof= right.tof;
  edep = right.edep;
  return *this;
}

// externally instanciated.
extern G4Allocator<TPCCHHit> TPCCHHitAllocator; 

inline void* TPCCHHit::operator new(size_t)
{
  void* aHit= (void*)TPCCHHitAllocator.MallocSingle();
  return aHit;
}

inline void TPCCHHit::operator delete(void* aHit)
{
  TPCCHHitAllocator.FreeSingle((TPCCHHit*) aHit);
}

#endif

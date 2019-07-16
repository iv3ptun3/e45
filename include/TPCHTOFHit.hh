// -*- C++ -*-

#ifndef TPC_HTOF_HIT_HH
#define TPC_HTOF_HIT_HH

#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//_____________________________________________________________________________
class TPCHTOFHit : public G4VHit
{
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
  TPCHTOFHit( void );
  //  TPCHTOFHit(G4ThreeVector& axyz, G4double t);
  TPCHTOFHit( G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t );
  TPCHTOFHit( G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
	      G4int tid, G4int pid, G4int did,G4double mass, G4int qq,
	      G4int parentid, G4ThreeVector& vtxmom, G4ThreeVector& vtxpos,
	      G4double vtxene, G4double tlength );
  virtual ~TPCHTOFHit( void );

  // copy constructor & assignment operator
  TPCHTOFHit(const TPCHTOFHit& right);
  const TPCHTOFHit& operator=(const TPCHTOFHit& right);

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

//_____________________________________________________________________________
inline TPCHTOFHit::TPCHTOFHit( const TPCHTOFHit& right )
  : G4VHit()
{
  xyz= right.xyz;
  pxyz= right.pxyz;
  tof= right.tof;
  edep = right.edep;
}

//_____________________________________________________________________________
inline const TPCHTOFHit&
TPCHTOFHit::operator =( const TPCHTOFHit& right )
{
  xyz= right.xyz;
  pxyz= right.pxyz;
  tof= right.tof;
  edep = right.edep;
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCHTOFHit> TPCHTOFHitAllocator;

//_____________________________________________________________________________
inline void*
TPCHTOFHit::operator new( size_t )
{
  void* aHit = (void*)TPCHTOFHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
TPCHTOFHit::operator delete( void* aHit )
{
  TPCHTOFHitAllocator.FreeSingle((TPCHTOFHit*) aHit);
}

#endif

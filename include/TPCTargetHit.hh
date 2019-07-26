// -*- C++ -*-

#ifndef TPC_TARGET_HIT_HH
#define TPC_TARGET_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCTargetHit : public G4VHit, public VHitInfo
{
public:
  TPCTargetHit( const G4String& name, G4Step* step );
  virtual ~TPCTargetHit( void );

  TPCTargetHit( const TPCTargetHit& right );
  const TPCTargetHit& operator=( const TPCTargetHit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCTargetHit::TPCTargetHit( const TPCTargetHit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCTargetHit&
TPCTargetHit::operator =( const TPCTargetHit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCTargetHit> TPCTargetHitAllocator;

//_____________________________________________________________________________
inline void*
TPCTargetHit::operator new( size_t )
{
  void* aHit = (void*)TPCTargetHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
TPCTargetHit::operator delete( void* aHit )
{
  TPCTargetHitAllocator.FreeSingle( (TPCTargetHit*) aHit );
}

#endif

// -*- C++ -*-

#ifndef TPC_TGT_VP_HIT_HH
#define TPC_TGT_VP_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCTargetVPHit : public G4VHit, public VHitInfo
{
public:
  TPCTargetVPHit( const G4String& name, G4Step* step );
  virtual ~TPCTargetVPHit( void );

  TPCTargetVPHit( const TPCTargetVPHit& right );
  const TPCTargetVPHit& operator=( const TPCTargetVPHit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCTargetVPHit::TPCTargetVPHit( const TPCTargetVPHit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCTargetVPHit&
TPCTargetVPHit::operator =( const TPCTargetVPHit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCTargetVPHit> TPCTargetVPHitAllocator;

//_____________________________________________________________________________
inline void*
TPCTargetVPHit::operator new( size_t )
{
  return TPCTargetVPHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
TPCTargetVPHit::operator delete( void* aHit )
{
  TPCTargetVPHitAllocator.FreeSingle( static_cast<TPCTargetVPHit*>( aHit ) );
}

#endif

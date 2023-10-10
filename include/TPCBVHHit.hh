// -*- C++ -*-

#ifndef TPC_BVH_HIT_HH
#define TPC_BVH_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCBVHHit : public G4VHit, public VHitInfo
{
public:
  TPCBVHHit( const G4String& name, G4Step* step );
  virtual ~TPCBVHHit( void );

  TPCBVHHit( const TPCBVHHit& right );
  const TPCBVHHit& operator=( const TPCBVHHit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCBVHHit::TPCBVHHit( const TPCBVHHit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCBVHHit&
TPCBVHHit::operator =( const TPCBVHHit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCBVHHit> TPCBVHHitAllocator;

//_____________________________________________________________________________
inline void*
TPCBVHHit::operator new( size_t )
{
  return TPCBVHHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
TPCBVHHit::operator delete( void* aHit )
{
  TPCBVHHitAllocator.FreeSingle( static_cast<TPCBVHHit*>( aHit ) );
}

#endif

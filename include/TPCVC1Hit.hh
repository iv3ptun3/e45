// -*- C++ -*-

#ifndef TPC_VC1_HIT_HH
#define TPC_VC1_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCVC1Hit : public G4VHit, public VHitInfo
{
public:
  TPCVC1Hit( const G4String& name, G4Step* step );
  virtual ~TPCVC1Hit( void );

  TPCVC1Hit( const TPCVC1Hit& right );
  const TPCVC1Hit& operator=( const TPCVC1Hit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCVC1Hit::TPCVC1Hit( const TPCVC1Hit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCVC1Hit&
TPCVC1Hit::operator =( const TPCVC1Hit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCVC1Hit> TPCVC1HitAllocator;

//_____________________________________________________________________________
inline void*
TPCVC1Hit::operator new( size_t )
{
  return TPCVC1HitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
TPCVC1Hit::operator delete( void* aHit )
{
  TPCVC1HitAllocator.FreeSingle( static_cast<TPCVC1Hit*>( aHit ) );
}

#endif

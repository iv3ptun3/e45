// -*- C++ -*-

#ifndef TPC_VC2_HIT_HH
#define TPC_VC2_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCVC2Hit : public G4VHit, public VHitInfo
{
public:
  TPCVC2Hit( const G4String& name, G4Step* step );
  virtual ~TPCVC2Hit( void );

  TPCVC2Hit( const TPCVC2Hit& right );
  const TPCVC2Hit& operator=( const TPCVC2Hit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCVC2Hit::TPCVC2Hit( const TPCVC2Hit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCVC2Hit&
TPCVC2Hit::operator =( const TPCVC2Hit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCVC2Hit> TPCVC2HitAllocator;

//_____________________________________________________________________________
inline void*
TPCVC2Hit::operator new( size_t )
{
  return TPCVC2HitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
TPCVC2Hit::operator delete( void* aHit )
{
  TPCVC2HitAllocator.FreeSingle( static_cast<TPCVC2Hit*>( aHit ) );
}

#endif

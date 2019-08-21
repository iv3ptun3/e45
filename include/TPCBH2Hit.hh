// -*- C++ -*-

#ifndef TPC_BH2_HIT_HH
#define TPC_BH2_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCBH2Hit : public G4VHit, public VHitInfo
{
public:
  TPCBH2Hit( const G4String& name, G4Step* step );
  virtual ~TPCBH2Hit( void );

  TPCBH2Hit( const TPCBH2Hit& right );
  const TPCBH2Hit& operator=( const TPCBH2Hit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCBH2Hit::TPCBH2Hit( const TPCBH2Hit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCBH2Hit&
TPCBH2Hit::operator =( const TPCBH2Hit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCBH2Hit> TPCBH2HitAllocator;

//_____________________________________________________________________________
inline void*
TPCBH2Hit::operator new( size_t )
{
  return TPCBH2HitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
TPCBH2Hit::operator delete( void* aHit )
{
  TPCBH2HitAllocator.FreeSingle( static_cast<TPCBH2Hit*>( aHit ) );
}

#endif

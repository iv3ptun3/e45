// -*- C++ -*-

#ifndef TPC_SCH_HIT_HH
#define TPC_SCH_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCSCHHit : public G4VHit, public VHitInfo
{
public:
  TPCSCHHit( const G4String& name, G4Step* step );
  virtual ~TPCSCHHit( void );

  TPCSCHHit( const TPCSCHHit& right );
  const TPCSCHHit& operator=( const TPCSCHHit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCSCHHit::TPCSCHHit( const TPCSCHHit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCSCHHit&
TPCSCHHit::operator =( const TPCSCHHit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCSCHHit> TPCSCHHitAllocator;

//_____________________________________________________________________________
inline void*
TPCSCHHit::operator new( size_t )
{
  return TPCSCHHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
TPCSCHHit::operator delete( void* aHit )
{
  TPCSCHHitAllocator.FreeSingle( static_cast<TPCSCHHit*>( aHit ) );
}

#endif

// -*- C++ -*-

#ifndef TPC_HTOF_HIT_HH
#define TPC_HTOF_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCHTOFHit : public G4VHit, public VHitInfo
{
public:
  TPCHTOFHit( const G4String& name, G4Step* step );
  virtual ~TPCHTOFHit( void );

  TPCHTOFHit( const TPCHTOFHit& right );
  const TPCHTOFHit& operator=( const TPCHTOFHit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCHTOFHit::TPCHTOFHit( const TPCHTOFHit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCHTOFHit&
TPCHTOFHit::operator =( const TPCHTOFHit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCHTOFHit> TPCHTOFHitAllocator;

//_____________________________________________________________________________
inline void*
TPCHTOFHit::operator new( size_t )
{
  return TPCHTOFHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
TPCHTOFHit::operator delete( void* aHit )
{
  TPCHTOFHitAllocator.FreeSingle( static_cast<TPCHTOFHit*>( aHit ) );
}

#endif

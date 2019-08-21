// -*- C++ -*-

#ifndef TPC_VP_HIT_HH
#define TPC_VP_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCVPHit : public G4VHit, public VHitInfo
{
public:
  TPCVPHit( const G4String& name, G4Step* step );
  virtual ~TPCVPHit( void );

  TPCVPHit( const TPCVPHit& right );
  const TPCVPHit& operator=( const TPCVPHit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCVPHit::TPCVPHit( const TPCVPHit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCVPHit&
TPCVPHit::operator =( const TPCVPHit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCVPHit> TPCVPHitAllocator;

//_____________________________________________________________________________
inline void*
TPCVPHit::operator new( size_t )
{
  return TPCVPHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
TPCVPHit::operator delete( void* aHit )
{
  TPCVPHitAllocator.FreeSingle( static_cast<TPCVPHit*>( aHit ) );
}

#endif

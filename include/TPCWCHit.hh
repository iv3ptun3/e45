// -*- C++ -*-

#ifndef TPC_WC_HIT_HH
#define TPC_WC_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCWCHit : public G4VHit, public VHitInfo
{
public:
  TPCWCHit( const G4String& name, G4Step* step );
  virtual ~TPCWCHit( void );

  TPCWCHit( const TPCWCHit& right );
  const TPCWCHit& operator=( const TPCWCHit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCWCHit::TPCWCHit( const TPCWCHit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCWCHit&
TPCWCHit::operator =( const TPCWCHit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCWCHit> TPCWCHitAllocator;

//_____________________________________________________________________________
inline void*
TPCWCHit::operator new( size_t )
{
  void* aHit = (void*)TPCWCHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
TPCWCHit::operator delete( void* aHit )
{
  TPCWCHitAllocator.FreeSingle( (TPCWCHit*) aHit );
}

#endif

// -*- C++ -*-

#ifndef TPC_LAC_HIT_HH
#define TPC_LAC_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCLACHit : public G4VHit, public VHitInfo
{
public:
  TPCLACHit( const G4String& name, G4Step* step );
  virtual ~TPCLACHit( void );

  TPCLACHit( const TPCLACHit& right );
  const TPCLACHit& operator=( const TPCLACHit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCLACHit::TPCLACHit( const TPCLACHit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCLACHit&
TPCLACHit::operator =( const TPCLACHit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCLACHit> TPCLACHitAllocator;

//_____________________________________________________________________________
inline void*
TPCLACHit::operator new( size_t )
{
  void* aHit = (void*)TPCLACHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
TPCLACHit::operator delete( void* aHit )
{
  TPCLACHitAllocator.FreeSingle( (TPCLACHit*) aHit );
}

#endif

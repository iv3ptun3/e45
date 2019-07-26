// -*- C++ -*-

#ifndef TPC_SDC_HIT_HH
#define TPC_SDC_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCSDCHit : public G4VHit, public VHitInfo
{
public:
  TPCSDCHit( const G4String& name, G4Step* step );
  virtual ~TPCSDCHit( void );

  TPCSDCHit( const TPCSDCHit& right );
  const TPCSDCHit& operator=( const TPCSDCHit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCSDCHit::TPCSDCHit( const TPCSDCHit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCSDCHit&
TPCSDCHit::operator =( const TPCSDCHit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCSDCHit> TPCSDCHitAllocator;

//_____________________________________________________________________________
inline void*
TPCSDCHit::operator new( size_t )
{
  void* aHit = (void*)TPCSDCHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
TPCSDCHit::operator delete( void* aHit )
{
  TPCSDCHitAllocator.FreeSingle( (TPCSDCHit*) aHit );
}

#endif

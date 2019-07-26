// -*- C++ -*-

#ifndef TPC_FTOF_HIT_HH
#define TPC_FTOF_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCFTOFHit : public G4VHit, public VHitInfo
{
public:
  TPCFTOFHit( const G4String& name, G4Step* step );
  virtual ~TPCFTOFHit( void );

  TPCFTOFHit( const TPCFTOFHit& right );
  const TPCFTOFHit& operator=( const TPCFTOFHit& right );

  void* operator new( size_t );
  void operator delete( void* aHit );

public:
  virtual void Draw( void );
  virtual void Print( void );
};

//_____________________________________________________________________________
inline
TPCFTOFHit::TPCFTOFHit( const TPCFTOFHit& right )
  : G4VHit( right ),
    VHitInfo( right )
{
}

//_____________________________________________________________________________
inline const TPCFTOFHit&
TPCFTOFHit::operator =( const TPCFTOFHit& right )
{
  VHitInfo::operator =( right );
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCFTOFHit> TPCFTOFHitAllocator;

//_____________________________________________________________________________
inline void*
TPCFTOFHit::operator new( size_t )
{
  void* aHit = (void*)TPCFTOFHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
TPCFTOFHit::operator delete( void* aHit )
{
  TPCFTOFHitAllocator.FreeSingle( (TPCFTOFHit*) aHit );
}

#endif

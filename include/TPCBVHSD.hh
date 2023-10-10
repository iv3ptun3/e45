// -*- C++ -*-

#ifndef TPC_BVH_SD_HH
#define TPC_BVH_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCBVHHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCBVHSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName( void );
  TPCBVHSD( const G4String& name );
  virtual ~TPCBVHSD( void );

private:
  G4THitsCollection<TPCBVHHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

//_____________________________________________________________________________
inline G4String
TPCBVHSD::ClassName( void )
{
  static G4String s_name("TPCBVHSD");
  return s_name;
}

#endif

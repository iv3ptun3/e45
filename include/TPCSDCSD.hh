// -*- C++ -*-

#ifndef TPC_SDC_SD_HH
#define TPC_SDC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCSDCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCSDCSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName( void );
  TPCSDCSD( const G4String& name );
  virtual ~TPCSDCSD( void );

private:
  G4THitsCollection<TPCSDCHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

//_____________________________________________________________________________
inline G4String
TPCSDCSD::ClassName( void )
{
  static G4String s_name("TPCSDCSD");
  return s_name;
}

#endif

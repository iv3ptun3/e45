// -*- C++ -*-

#ifndef TPC_SCH_SD_HH
#define TPC_SCH_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCSCHHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCSCHSD : public G4VSensitiveDetector
{
public:
  TPCSCHSD( const G4String& name );
  virtual ~TPCSCHSD( void );

private:
  G4THitsCollection<TPCSCHHit>* hitsCollection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

#endif

// -*- C++ -*-

#ifndef TPC_CH_SD_HH
#define TPC_CH_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCCHHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCCHSD : public G4VSensitiveDetector
{
public:
  TPCCHSD( const G4String& name );
  virtual ~TPCCHSD( void );

private:
  G4THitsCollection<TPCCHHit>* hitsCollection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

#endif

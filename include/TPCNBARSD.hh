// -*- C++ -*-

#ifndef TPC_NBAR_SD_HH
#define TPC_NBAR_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCNBARHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCNBARSD : public G4VSensitiveDetector
{
public:
  TPCNBARSD( const G4String& name );
  virtual ~TPCNBARSD( void );

private:
  G4THitsCollection<TPCNBARHit>* hitsCollection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );

};

#endif

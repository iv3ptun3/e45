// -*- C++ -*-

#ifndef TPC_AC_SD_HH
#define TPC_AC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCACHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCACSD : public G4VSensitiveDetector
{
public:
  TPCACSD( const G4String& name );
  virtual ~TPCACSD( void );

private:
  G4THitsCollection<TPCACHit>* hitsCollection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );

};

#endif

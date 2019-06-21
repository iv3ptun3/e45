// -*- C++ -*-

#ifndef TPC_SCINT_SD_HH
#define TPC_SCINT_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCScintHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCScintSD : public G4VSensitiveDetector
{
public:
  TPCScintSD( const G4String& name );
  virtual ~TPCScintSD( void );

private:
  G4THitsCollection<TPCScintHit>* hitsCollection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

#endif

// -*- C++ -*-

#ifndef TPC_FTOF_SD_HH
#define TPC_FTOF_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCFTOFHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCFTOFSD : public G4VSensitiveDetector
{
public:
  TPCFTOFSD( const G4String& name );
  virtual ~TPCFTOFSD( void );

private:
  G4THitsCollection<TPCFTOFHit>* hitsCollection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );

};

#endif

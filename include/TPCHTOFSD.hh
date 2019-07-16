// -*- C++ -*-

#ifndef TPC_HTOF_SD_HH
#define TPC_HTOF_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCHTOFHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCHTOFSD : public G4VSensitiveDetector
{
public:
  TPCHTOFSD( const G4String& name );
  virtual ~TPCHTOFSD( void );

private:
  G4THitsCollection<TPCHTOFHit>* hitsCollection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

#endif

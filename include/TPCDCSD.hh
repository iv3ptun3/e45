// -*- C++ -*-

#ifndef TPC_DC_SD_HH
#define TPC_DC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCDCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCDCSD : public G4VSensitiveDetector
{
public:
  TPCDCSD( const G4String& name );
  virtual ~TPCDCSD( void );

private:
  G4THitsCollection<TPCDCHit>* hitsCollection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );

};

#endif

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
  static G4String ClassName( void );
  TPCFTOFSD( const G4String& name );
  virtual ~TPCFTOFSD( void );

private:
  G4THitsCollection<TPCFTOFHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

//_____________________________________________________________________________
inline G4String
TPCFTOFSD::ClassName( void )
{
  static G4String s_name("TPCFTOFSD");
  return s_name;
}

#endif

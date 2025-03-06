// -*- C++ -*-

#ifndef TPC_VC1_SD_HH
#define TPC_VC1_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCVC1Hit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCVC1SD : public G4VSensitiveDetector
{
public:
  static G4String ClassName( void );
  TPCVC1SD( const G4String& name );
  virtual ~TPCVC1SD( void );

private:
  G4THitsCollection<TPCVC1Hit>* m_hits_collection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

//_____________________________________________________________________________
inline G4String
TPCVC1SD::ClassName( void )
{
  static G4String s_name("TPCVC1SD");
  return s_name;
}

#endif

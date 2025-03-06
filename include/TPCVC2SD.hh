// -*- C++ -*-

#ifndef TPC_VC2_SD_HH
#define TPC_VC2_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCVC2Hit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCVC2SD : public G4VSensitiveDetector
{
public:
  static G4String ClassName( void );
  TPCVC2SD( const G4String& name );
  virtual ~TPCVC2SD( void );

private:
  G4THitsCollection<TPCVC2Hit>* m_hits_collection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

//_____________________________________________________________________________
inline G4String
TPCVC2SD::ClassName( void )
{
  static G4String s_name("TPCVC2SD");
  return s_name;
}

#endif

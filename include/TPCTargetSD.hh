// -*- C++ -*-

#ifndef TPC_TARGET_SD_HH
#define TPC_TARGET_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCTargetHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCTargetSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName( void );
  TPCTargetSD( const G4String& name );
  virtual ~TPCTargetSD( void );

private:
  G4THitsCollection<TPCTargetHit>* m_hits_collection;

public:
  G4int ntrk;
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );

};

//_____________________________________________________________________________
inline G4String
TPCTargetSD::ClassName( void )
{
  static G4String s_name("TPCTargetSD");
  return s_name;
}

#endif

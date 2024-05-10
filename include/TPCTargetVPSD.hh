// -*- C++ -*-

#ifndef TPC_TGT_VP_SD_HH
#define TPC_TGT_VP_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCTargetVPHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCTargetVPSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName( void );
  TPCTargetVPSD( const G4String& name );
  virtual ~TPCTargetVPSD( void );

private:
  G4THitsCollection<TPCTargetVPHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

//_____________________________________________________________________________
inline G4String
TPCTargetVPSD::ClassName( void )
{
  static G4String s_name("TPCTargetVPSD");
  return s_name;
}

#endif

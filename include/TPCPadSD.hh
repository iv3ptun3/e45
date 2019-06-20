// -*- C++ -*-

#ifndef TPC_PAD_SD_HH
#define TPC_PAD_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCPadHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCPadSD : public G4VSensitiveDetector
{
public:
  TPCPadSD( const G4String& name );
  virtual ~TPCPadSD( void );

private:
  G4THitsCollection<TPCPadHit>* hitsCollection;
  //  G4int pidtr[50];
  //  G4int nparticle;
  G4double select_plane;
  G4int num_plane;
  G4double deadarea;
  G4int num_deadarea;
  G4int select_dead;

  G4int m_gem_discharge;
  G4int m_gem_fix_dead;
  G4int m_gem_dead_plane;
  G4int m_gem_dead_plane_division;
  G4double m_dead_area;

public:
  G4int ntrk;
  // virtual methods
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void Initialize( G4HCofThisEvent* HCTE );
  virtual void EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void DrawAll( void );
  virtual void PrintAll( void );
};

#endif

// ====================================================================
//   TPCScintSD.hh
//
// ====================================================================
#ifndef TPC_SCINT_SD_H
#define TPC_SCINT_SD_H
 
#include "G4VSensitiveDetector.hh"
#include "TPCScintHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class TPCScintSD : public G4VSensitiveDetector {
private:
  G4THitsCollection<TPCScintHit>* hitsCollection;

public:
  TPCScintSD(const G4String& name);
  virtual ~TPCScintSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};

#endif

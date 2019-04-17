// ====================================================================
//   TPCTargetSD.hh
//
// ====================================================================
#ifndef TPC_TARGET_SD_H
#define TPC_TARGET_SD_H
 
#include "G4VSensitiveDetector.hh"
#include "TPCTargetHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class TPCTargetSD : public G4VSensitiveDetector {
private:
  G4THitsCollection<TPCTargetHit>* hitsCollection;
  //  G4int pidtr[50];
  //  G4int nparticle;
public:
  G4int ntrk;
  TPCTargetSD(const G4String& name);
  virtual ~TPCTargetSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};

#endif

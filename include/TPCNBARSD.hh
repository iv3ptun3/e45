// ====================================================================
//   TPCNBARSD.hh
//
// ====================================================================
#ifndef TPC_NBAR_SD_H
#define TPC_NBAR_SD_H
 
#include "G4VSensitiveDetector.hh"
#include "TPCNBARHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class TPCNBARSD : public G4VSensitiveDetector {
private:
  G4THitsCollection<TPCNBARHit>* hitsCollection;

public:
  TPCNBARSD(const G4String& name);
  virtual ~TPCNBARSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};

#endif

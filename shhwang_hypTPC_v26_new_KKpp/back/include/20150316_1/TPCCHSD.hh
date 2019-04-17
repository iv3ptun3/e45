// ====================================================================
//   TPCCHSD.hh
//
// ====================================================================
#ifndef TPC_CH_SD_H
#define TPC_CH_SD_H
 
#include "G4VSensitiveDetector.hh"
#include "TPCCHHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class TPCCHSD : public G4VSensitiveDetector {
private:
  G4THitsCollection<TPCCHHit>* hitsCollection;

public:
  TPCCHSD(const G4String& name);
  virtual ~TPCCHSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};

#endif

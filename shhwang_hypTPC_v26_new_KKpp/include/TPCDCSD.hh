// ====================================================================
//   TPCDCSD.hh
//
// ====================================================================
#ifndef TPC_DC_SD_H
#define TPC_DC_SD_H
 
#include "G4VSensitiveDetector.hh"
#include "TPCDCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class TPCDCSD : public G4VSensitiveDetector {
private:
  G4THitsCollection<TPCDCHit>* hitsCollection;

public:
  TPCDCSD(const G4String& name);
  virtual ~TPCDCSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};

#endif

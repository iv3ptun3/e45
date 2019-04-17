// ====================================================================
//   TPCFTOFSD.hh
//
// ====================================================================
#ifndef TPC_FTOF_SD_H
#define TPC_FTOF_SD_H
 
#include "G4VSensitiveDetector.hh"
#include "TPCFTOFHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class TPCFTOFSD : public G4VSensitiveDetector {
private:
  G4THitsCollection<TPCFTOFHit>* hitsCollection;

public:
  TPCFTOFSD(const G4String& name);
  virtual ~TPCFTOFSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};

#endif

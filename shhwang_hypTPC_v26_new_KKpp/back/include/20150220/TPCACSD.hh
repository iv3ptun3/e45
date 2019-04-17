// ====================================================================
//   TPCACSD.hh
//
// ====================================================================
#ifndef TPC_AC_SD_H
#define TPC_AC_SD_H
 
#include "G4VSensitiveDetector.hh"
#include "TPCACHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class TPCACSD : public G4VSensitiveDetector {
private:
  G4THitsCollection<TPCACHit>* hitsCollection;

public:
  TPCACSD(const G4String& name);
  virtual ~TPCACSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};

#endif

// ====================================================================
//   TPCEventAction.hh
//
// ====================================================================
#ifndef TPC_EVENT_ACTION_H
#define TPC_EVENT_ACTION_H 

#include "G4UserEventAction.hh"

class G4Event;
class TPCAnaManager;

class TPCEventAction : public G4UserEventAction {
public:
  TPCEventAction(TPCAnaManager* ana);
  virtual ~TPCEventAction();

  virtual void BeginOfEventAction(const G4Event* anEvent);
  virtual void EndOfEventAction(const G4Event* anEvent);
private:
  TPCAnaManager* AnaManager;
  G4int ac_use;
  G4int n_bar_use;
  G4int experiment_num;
  G4int with_kurama;
  //  G4int pidtr[50];
  //  G4int nparticle;

};

#endif

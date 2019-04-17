// ====================================================================
//   TPCRunAction.hh
//
// ====================================================================
#ifndef TPC_RUN_ACTION_H
#define TPC_RUN_ACTION_H

#include "G4UserRunAction.hh"
#include "TPCAnaManager.hh"

class G4Run;
class TPCAnaManager;

class TPCRunAction : public G4UserRunAction {
public:
  //  TPCRunAction()
  TPCRunAction(TPCAnaManager* ana);
  virtual ~TPCRunAction();

  virtual void BeginOfRunAction(const G4Run* aRun);
  virtual void EndOfRunAction(const G4Run* aRun);
private:
  TPCAnaManager* AnaManager;
};

#endif

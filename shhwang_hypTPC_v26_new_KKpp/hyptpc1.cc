// ====================================================================
//    shhwang h-dibaryon code from Feb. 2012.
// ====================================================================
#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "TPCDetectorConstruction.hh"
#include "TPCPhysicsList.hh"
#include "TPCPrimaryGeneratorAction.hh"
#include "TPCRunAction.hh"
#include "TPCEventAction.hh"
#include "TPCSteppingAction.hh"

#include "TPCAnaManager.hh"
#include "getenv.hh"
#ifdef G4VIS_USE
#include "TPCVisManager.hh"
//#include "G4VisExecutive.hh"
#endif

#include <fstream>
std::ofstream ofs;

// ====================================================================
//     main
// ====================================================================
int main(int argc, char** argv) 
{
  // run manager
  G4RunManager* runManager= new G4RunManager;  G4cout << G4endl;
  // set mandatory user initialization classes...
  // detector setup
  runManager-> SetUserInitialization(new TPCDetectorConstruction);
  // particles and physics processes
  runManager-> SetUserInitialization(new TPCPhysicsList);

  getenv_read();

  TPCAnaManager* AnaManager = new TPCAnaManager();

  TPCPrimaryGeneratorAction* PrimaryGenerator =
    new TPCPrimaryGeneratorAction(AnaManager);

  TPCRunAction* RunAct = new TPCRunAction(AnaManager);

  TPCEventAction* EventAct = new TPCEventAction(AnaManager);

  TPCSteppingAction* StepAct = new TPCSteppingAction;


  // primary generator
  runManager-> SetUserAction(PrimaryGenerator);
  // user action classes... (optional)
  runManager-> SetUserAction(RunAct);
  runManager-> SetUserAction(StepAct);
  runManager-> SetUserAction(EventAct);

#ifdef G4VIS_USE
  // initialize visualization package
  G4VisManager* visManager= new TPCVisManager;
  //  G4VisManager* visManager= new G4VisExecutive;
  visManager-> Initialize();
  //  G4cout << G4endl;
#endif

  // user session...
  runManager-> Initialize();


  
    
  if(argc==1) { // interactive session, if no arguments given
    // tcsh-like
    G4UItcsh* tcsh= new G4UItcsh("hypTPC(%s)[%/][%h]:");
    G4UIterminal* session= new G4UIterminal(tcsh);
    tcsh-> SetLsColor(GREEN, CYAN);
    session-> SessionStart();
    delete session;

  } else { // batch mode
    G4UImanager* UImanager= G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager-> ApplyCommand(command+fileName);
  }
  


  // terminating...
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;  G4cout << G4endl;
  return 0;
}


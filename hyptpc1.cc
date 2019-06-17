// -*- C++ -*-

#include <G4RunManager.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>
#include <QGSP_BERT.hh>

#include <TFile.h>

#include "ConfMan.hh"
#include "TPCDetectorConstruction.hh"
#include "TPCPhysicsList.hh"
#include "TPCPrimaryGeneratorAction.hh"
#include "TPCRunAction.hh"
#include "TPCEventAction.hh"
#include "TPCSteppingAction.hh"

#ifdef G4VIS_USE
#include "TPCVisManager.hh"
//#include "G4VisExecutive.hh"
#endif

enum EArgv { kProcess, kConfFile, kOutputName, kG4Macro, kArgc };

//_____________________________________________________________________________
int
main( int argc, char** argv )
{
  if( argc != kArgc-1 && argc != kArgc ){
    G4cout << "Usage: " << argv[kProcess]
	   << " [ConfFile] [OutputName] (G4Macro)" << G4endl;
    return EXIT_SUCCESS;
  }

  ConfMan& gConf = ConfMan::GetInstance();
  if( !gConf.Initialize( argv[kConfFile] ) ||
      !gConf.InitializeParameterFiles() ){
    return EXIT_FAILURE;
  }

  new TFile( argv[kOutputName], "RECREATE" );

  auto runManager= new G4RunManager;
  runManager->SetUserInitialization( new TPCDetectorConstruction );
  // runManager-> SetUserInitialization( new QGSP_BERT );
  runManager->SetUserInitialization( new TPCPhysicsList );
  runManager->SetUserAction( new TPCPrimaryGeneratorAction );
  runManager->SetUserAction( new TPCRunAction );
  runManager->SetUserAction( new TPCSteppingAction );
  runManager->SetUserAction( new TPCEventAction );

#ifdef G4VIS_USE
  auto visManager = new TPCVisManager;
  // auto visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  runManager->Initialize();

  // interactive session, if no arguments given
  if( argc == kArgc-1 ) {
    auto tcsh = new G4UItcsh( "HypTPC(%s)[%/][%h]: " );
    auto session = new G4UIterminal( tcsh );
    tcsh->SetLsColor( GREEN, CYAN );
    session->SessionStart();
    delete session;
  }
  // batch mode
  else if( argc == kArgc ) {
    auto uiManager= G4UImanager::GetUIpointer();
    G4String command("/control/execute ");
    G4String fileName( argv[kG4Macro] );
    uiManager->ApplyCommand( command + fileName );
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  G4cout << "finished" << G4endl;
  return EXIT_SUCCESS;
}

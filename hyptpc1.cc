// -*- C++ -*-

#include <G4RunManager.hh>
#include <G4TrajectoryDrawByCharge.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>
#include <QGSP_BERT.hh>
#include <G4VisExecutive.hh>
#include <G4UIExecutive.hh>

#include <TFile.h>

#include "ConfMan.hh"
#include "TPCDetectorConstruction.hh"
#include "TPCPhysicsList.hh"
#include "TPCPrimaryGeneratorAction.hh"
#include "TPCRunAction.hh"
#include "TPCEventAction.hh"
#include "TPCSteppingAction.hh"

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
  if( !gConf.Initialize( argv[kConfFile] ) ){
    return EXIT_FAILURE;
  }

  new TFile( argv[kOutputName], "RECREATE" );

  auto runManager = new G4RunManager;
  runManager->SetUserInitialization( new TPCDetectorConstruction );
  if( gConf.Get<G4String>( "Physics" ) == "QGSP_BERT" )
    runManager->SetUserInitialization( new QGSP_BERT );
  else if( gConf.Get<G4String>( "Physics" ) == "USER" )
    runManager->SetUserInitialization( new TPCPhysicsList );
  runManager->SetUserAction( new TPCPrimaryGeneratorAction );
  runManager->SetUserAction( new TPCRunAction );
  runManager->SetUserAction( new TPCSteppingAction );
  runManager->SetUserAction( new TPCEventAction );

  // added by hjb
  auto visManager = new G4VisExecutive;
  visManager->SetVerboseLevel(0);
  visManager->Initialize();

  runManager->Initialize();
  auto uiManager = G4UImanager::GetUIpointer();

  if( gConf.Get<G4bool>( "EVDISP" ) )
    uiManager->ApplyCommand("/control/execute g4macro/vis.mac");

  // interactive session, if no arguments given
  if( argc == kArgc-1 ) {
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    uiManager->ApplyCommand("/control/execute g4macro/init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
  // batch mode
  else if( argc == kArgc ) {
    G4String command("/control/execute ");
    G4String fileName( argv[kG4Macro] );
    uiManager->ApplyCommand( command + fileName );
  }

  gFile->Write();
  gFile->Close();


  delete visManager;
  delete runManager;
  
  G4cout << "finished" << G4endl;
  return EXIT_SUCCESS;
}

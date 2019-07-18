// -*- C++ -*-

#include "TPCRunAction.hh"

#include <fstream>

#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4RunMessenger.hh>
#include <G4StateManager.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>

#include "FuncName.hh"
#include "GetNumberFromKernelEntropyPool.hh"
#include "TPCAnaManager.hh"

namespace
{
  TPCAnaManager& gAnaMan = TPCAnaManager::GetInstance();
}

//_____________________________________________________________________________
TPCRunAction::TPCRunAction( void )
  : G4UserRunAction()
{
}

//_____________________________________________________________________________
TPCRunAction::~TPCRunAction( void )
{
}

//_____________________________________________________________________________
void
TPCRunAction::BeginOfRunAction( const G4Run* aRun )
{
  G4cout << FUNC_NAME << G4endl
	 << "   Run# = " << aRun->GetRunID() << G4endl;
  gAnaMan.BeginOfRunAction(aRun->GetRunID());
  // currentEvent = GenerateEvent(i_event);
  // eventManager->ProcessOneEvent(currentEvent);

  //  if (G4VVisManager::GetConcreteInstance())
  //    {
  //      G4UImanager* UI = G4UImanager::GetUIpointer();
  //      UI->ApplyCommand("/vis/scene/notifyHandlers");
  //    }

  // Set the seed for randam function //Hwang-san
  // the_time = time((time_t *)0);
  // CLHEP::HepRandom::setTheSeed(the_time);

  // int initSeed = GetIntFromKernelEntropyPool()&0x7FFFFFFF;
  // G4Random::setTheSeed(initSeed);
  G4Random::setTheSeed( std::time( nullptr ) );
#ifdef DEBUG
  G4cout << "   Initial Seed = " << G4Random::getTheSeed() << G4endl;
  G4Random::showEngineStatus();
#endif
  // G4StateManager* stateManager = G4StateManager::GetStateManager();
  // G4RunManager* RunManager = G4RunManager::GetRunManager();
  // RunManager->DoEventLoop(5, "run_tmp.mac", 1);
}

//_____________________________________________________________________________
void
TPCRunAction::EndOfRunAction( const G4Run* aRun )
{
  gAnaMan.EndOfRunAction();
  G4cout << FUNC_NAME << G4endl
	 << "   Generated event# = " << aRun->GetNumberOfEvent() << G4endl;
}

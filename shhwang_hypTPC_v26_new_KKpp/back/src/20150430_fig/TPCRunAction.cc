// ====================================================================
//   TPCRunAction.cc
//
// ====================================================================
#include "G4Run.hh"
#include "TPCRunAction.hh"
#include "TPCAnaManager.hh"
#include <fstream>
#include "GetNumberFromKernelEntropyPool.hh"

extern std::ofstream ofs;

////////////////////////////
TPCRunAction::TPCRunAction(TPCAnaManager* ana)
  : AnaManager(ana)
////////////////////////////
{
}

/////////////////////////////
TPCRunAction::~TPCRunAction()
/////////////////////////////
{
}

//////////////////////////////////////////////////////
void TPCRunAction::BeginOfRunAction(const G4Run* aRun)
//////////////////////////////////////////////////////
{

  time_t the_time;

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  AnaManager->BeginOfRunAction(aRun->GetRunID());
  //  if (G4VVisManager::GetConcreteInstance())
  //    {
  //      G4UImanager* UI = G4UImanager::GetUIpointer();
  //      UI->ApplyCommand("/vis/scene/notifyHandlers");
  //    }

  //Set the seed for randam function //Hwang-san
  // the_time = time((time_t *)0);
  // CLHEP::HepRandom::setTheSeed(the_time);

  int initSeed=GetIntFromKernelEntropyPool()&0x7FFFFFFF;
  CLHEP::HepRandom::setTheSeed(initSeed);
  int startSeed=CLHEP::HepRandom::getTheSeed();
  G4cout << "*** Initial Seed = " << startSeed << G4endl;
  CLHEP::HepRandom::showEngineStatus();


  // ofs.open("a.dat", std::ios::out);
  // if(! ofs.good()) {
  //   G4String errorMessage= "*** fail to open a file (a.out).";
  //   G4Exception(errorMessage);
  // }
}

////////////////////////////////////////////////////
void TPCRunAction::EndOfRunAction(const G4Run* aRun)
////////////////////////////////////////////////////
{
  AnaManager->EndOfRunAction();
  //  ofs.close();
  G4cout << ">>> #events generated= " << aRun-> GetNumberOfEvent() << G4endl;
}


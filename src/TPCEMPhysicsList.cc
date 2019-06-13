// ====================================================================
//   TPCEMPhysicsList.cc
//
// ====================================================================
#include "TPCEMPhysics.hh"
#include "G4ProcessManager.hh"

////////////////////////////////////////////////
TPCEMPhysics::TPCEMPhysics(const G4String& name)
  :  G4VPhysicsConstructor(name)
////////////////////////////////////////////////
{
}

/////////////////////////////
TPCEMPhysics::~TPCEMPhysics()
/////////////////////////////
{
}

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

//////////////////////////////////////
void TPCEMPhysics::ConstructParticle()
//////////////////////////////////////
{
  // gamma
  G4Gamma::GammaDefinition();

  // electron
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();

  // nutrino
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

#include "G4GammaConversion.hh"
#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

/////////////////////////////////////
void TPCEMPhysics::ConstructProcess()
/////////////////////////////////////
{
  G4ProcessManager* pManager;

  // ===============================================================
  // gamma physics
  // ===============================================================
  pManager= G4Gamma::Gamma()-> GetProcessManager();

  pManager-> AddDiscreteProcess(new G4GammaConversion());
  pManager-> AddDiscreteProcess(new G4ComptonScattering());
  pManager-> AddDiscreteProcess(new G4PhotoElectricEffect());
  
  // ===============================================================
  // elecron physics
  // ===============================================================
  pManager= G4Electron::Electron()-> GetProcessManager();

  pManager-> AddProcess(new G4eMultipleScattering(), ordInActive, 1, 1);
  pManager-> AddProcess(new G4eIonisation(),        ordInActive, 2, 2);
  pManager-> AddDiscreteProcess(new G4eBremsstrahlung());
  
  // ===============================================================
  // positron physics
  // ===============================================================
  pManager= G4Positron::Positron()-> GetProcessManager();

  pManager-> AddProcess(new G4eMultipleScattering(), ordInActive, 1, 1);
  pManager-> AddProcess(new G4eIonisation(),        ordInActive, 2, 2);
  pManager-> AddDiscreteProcess(new G4eBremsstrahlung());

  G4eplusAnnihilation* eplusAnnihilation= new G4eplusAnnihilation();
  pManager-> AddDiscreteProcess(eplusAnnihilation);
  pManager-> AddRestProcess(eplusAnnihilation);
}


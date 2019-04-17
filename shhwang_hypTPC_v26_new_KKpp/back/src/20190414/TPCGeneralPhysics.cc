// ====================================================================
//   TPCGeneralPhysicsList.cc
//
// ====================================================================
#include "TPCGeneralPhysics.hh"
#include "G4ProcessManager.hh"

//////////////////////////////////////////////////////////
TPCGeneralPhysics::TPCGeneralPhysics(const G4String& name)
  :  G4VPhysicsConstructor(name)
//////////////////////////////////////////////////////////
{
}

///////////////////////////////////////
TPCGeneralPhysics::~TPCGeneralPhysics()
///////////////////////////////////////
{
}

#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"

///////////////////////////////////////////
void TPCGeneralPhysics::ConstructParticle()
///////////////////////////////////////////
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
}

#include "G4Decay.hh"

//////////////////////////////////////////
void TPCGeneralPhysics::ConstructProcess()
//////////////////////////////////////////
{
  G4Decay* decayProcess= new G4Decay();

  // add decay process...
  theParticleIterator-> reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle= theParticleIterator-> value();
    G4ProcessManager* pManager= particle-> GetProcessManager();
    if (decayProcess-> IsApplicable(*particle)) {
      pManager-> AddProcess(decayProcess);
      pManager-> SetProcessOrdering(decayProcess, idxPostStep);
      pManager-> SetProcessOrdering(decayProcess, idxAtRest);
    }
  }
}


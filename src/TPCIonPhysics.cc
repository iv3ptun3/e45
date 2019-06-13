// ====================================================================
//   TPCIonPhysicsList.cc
//
// ====================================================================
#include "TPCIonPhysics.hh"
#include "G4ProcessManager.hh"

//////////////////////////////////////////////////
TPCIonPhysics::TPCIonPhysics(const G4String& name)
  :  G4VPhysicsConstructor(name)
//////////////////////////////////////////////////
{
}

///////////////////////////////
TPCIonPhysics::~TPCIonPhysics()
///////////////////////////////
{
}


#include "G4IonConstructor.hh"

///////////////////////////////////////
void TPCIonPhysics::ConstructParticle()
///////////////////////////////////////
{
  // construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}


#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"


#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"

#include "G4DeuteronInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"

#include "G4TritonInelasticProcess.hh"
#include "G4LETritonInelastic.hh"

#include "G4AlphaInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"

//////////////////////////////////////
void TPCIonPhysics::ConstructProcess()
//////////////////////////////////////
{
  G4ProcessManager* pManager;

  // elastic model
  G4LElastic* elasticModel= new G4LElastic();  

  // ===============================================================
  // generic Ion
  // ===============================================================
  pManager= G4GenericIon::GenericIon()-> GetProcessManager();

  pManager-> AddProcess(new G4hMultipleScattering(), ordInActive, 1, 1);
  pManager-> AddProcess(new G4hIonisation(),        ordInActive, 2, 2);

  G4HadronElasticProcess* genericIonElasticProc= new G4HadronElasticProcess();
  genericIonElasticProc-> RegisterMe(elasticModel);
  pManager-> AddDiscreteProcess(genericIonElasticProc);

  // ===============================================================
  // deuteron
  // ===============================================================
  pManager= G4Deuteron::Deuteron()-> GetProcessManager();

  pManager-> AddProcess(new G4hMultipleScattering(), ordInActive, 1, 1);
  pManager-> AddProcess(new G4hIonisation(),        ordInActive, 2, 2);

  G4HadronElasticProcess* deuteronElasticProc= new G4HadronElasticProcess();
  deuteronElasticProc-> RegisterMe(elasticModel);
  pManager-> AddDiscreteProcess(deuteronElasticProc);

  G4DeuteronInelasticProcess* deuteronInelasticProc= 
    new G4DeuteronInelasticProcess();
  deuteronInelasticProc-> RegisterMe(new G4LEDeuteronInelastic());
  pManager-> AddDiscreteProcess(deuteronInelasticProc);

  // ===============================================================
  // triton
  // ===============================================================
  pManager= G4Triton::Triton()-> GetProcessManager();

  pManager-> AddProcess(new G4hMultipleScattering(), ordInActive, 1, 1);
  pManager-> AddProcess(new G4hIonisation(),        ordInActive, 2, 2);

  G4HadronElasticProcess* tritonElasticProc= new G4HadronElasticProcess();
  tritonElasticProc-> RegisterMe(elasticModel);
  pManager-> AddDiscreteProcess(tritonElasticProc);

  G4TritonInelasticProcess* tritonInelasticProc= 
    new G4TritonInelasticProcess();
  tritonInelasticProc-> RegisterMe(new G4LETritonInelastic());
  pManager-> AddDiscreteProcess(tritonInelasticProc);

  // ===============================================================
  // alpha
  // ===============================================================
  pManager= G4Alpha::Alpha()-> GetProcessManager();

  pManager-> AddProcess(new G4hMultipleScattering(), ordInActive, 1, 1);
  pManager-> AddProcess(new G4hIonisation(),        ordInActive, 2, 2);

  G4HadronElasticProcess* alphaElasticProc= new G4HadronElasticProcess();
  alphaElasticProc-> RegisterMe(elasticModel);
  pManager-> AddDiscreteProcess(alphaElasticProc);

  G4AlphaInelasticProcess* alphaInelasticProc= new G4AlphaInelasticProcess();
  alphaInelasticProc-> RegisterMe(new G4LEAlphaInelastic());
  pManager-> AddDiscreteProcess(alphaInelasticProc);

  // ===============================================================
  // He3
  // ===============================================================
  pManager= G4He3::He3()-> GetProcessManager();

  pManager-> AddProcess(new G4hMultipleScattering(), ordInActive, 1, 1);
  pManager-> AddProcess(new G4hIonisation(),        ordInActive, 2, 2);

  G4HadronElasticProcess* he3ElasticProc= new G4HadronElasticProcess();
  he3ElasticProc-> RegisterMe(elasticModel);
  pManager-> AddDiscreteProcess(he3ElasticProc);

}


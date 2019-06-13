// ====================================================================
//   TPCPhysicsList.hh
//
// ====================================================================
#ifndef TPC_PHYSICS_LIST_H
#define TPC_PHYSICS_LIST_H

#include "G4VModularPhysicsList.hh"

#include "G4VUserPhysicsList.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronInelasticDataSet.hh"

#include "G4LElastic.hh"  

#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4HEPionPlusInelastic.hh"

#include "G4PionMinusInelasticProcess.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"

#include "G4KaonPlusInelasticProcess.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"

#include "G4KaonMinusInelasticProcess.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"

#include "G4KaonZeroSInelasticProcess.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"

#include "G4KaonZeroLInelasticProcess.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"

#include "G4ProtonInelasticProcess.hh"
#include "G4LEProtonInelastic.hh"
#include "G4HEProtonInelastic.hh"

#include "G4NeutronInelasticProcess.hh"
#include "G4LENeutronInelastic.hh"
#include "G4HENeutronInelastic.hh"

#include "G4AntiProtonInelasticProcess.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"

#include "G4AntiNeutronInelasticProcess.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"

#include "G4LambdaInelasticProcess.hh"
#include "G4LELambdaInelastic.hh"
#include "G4HELambdaInelastic.hh"

#include "G4AntiLambdaInelasticProcess.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"

#include "G4SigmaPlusInelasticProcess.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"

#include "G4SigmaMinusInelasticProcess.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"

#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"

#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"

#include "G4hMultipleScattering.hh"
#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hIonisation.hh"

#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4IonConstructor.hh"

#include "globals.hh"

class TPCPhysicsList: public G4VModularPhysicsList {
public:

  TPCPhysicsList();
  virtual ~TPCPhysicsList();
  protected:

  // construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
  
  void SetCuts();
  
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructShortLived();
  void ConstructStableHyperons();
  void ConstructIons();

  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();
  void ConstructHadron();


protected:

  // Elastic Process
   G4HadronElasticProcess theElasticProcess;
   G4LElastic*            theElasticModel;

  // pi + 
  G4PionPlusInelasticProcess thePionPlusInelastic;
  G4LEPionPlusInelastic* theLEPionPlusModel;

  // pi - 
  G4PionMinusInelasticProcess thePionMinusInelastic;
  G4LEPionMinusInelastic* theLEPionMinusModel;

  // K + 
  G4KaonPlusInelasticProcess theKaonPlusInelastic;
  G4LEKaonPlusInelastic* theLEKaonPlusModel;
  G4HEKaonPlusInelastic* theHEKaonPlusModel;

  // K - 
  G4KaonMinusInelasticProcess theKaonMinusInelastic;
  G4LEKaonMinusInelastic* theLEKaonMinusModel;
  G4HEKaonMinusInelastic* theHEKaonMinusModel;

  // K0L
  G4KaonZeroLInelasticProcess theKaonZeroLInelastic;
  G4LEKaonZeroLInelastic* theLEKaonZeroLModel;
  G4HEKaonZeroInelastic* theHEKaonZeroLModel;

  // K0S
  G4KaonZeroSInelasticProcess theKaonZeroSInelastic;
  G4LEKaonZeroSInelastic* theLEKaonZeroSModel;
  G4HEKaonZeroInelastic* theHEKaonZeroSModel;

  // Proton
  G4ProtonInelasticProcess theProtonInelastic;
  G4LEProtonInelastic* theLEProtonModel;
  G4HEProtonInelastic* theHEProtonModel;

  // Neutron
  G4NeutronInelasticProcess theNeutronInelastic;
  G4LENeutronInelastic* theLENeutronModel;
  G4HENeutronInelastic* theHENeutronModel;
  
  // Anti-Proton
  G4AntiProtonInelasticProcess theAntiProtonInelastic;
  G4LEAntiProtonInelastic* theLEAntiProtonModel;
  G4HEAntiProtonInelastic* theHEAntiProtonModel;

  // anti-neutron
  G4AntiNeutronInelasticProcess  theAntiNeutronInelastic;
  G4LEAntiNeutronInelastic* theLEAntiNeutronModel;
  G4HEAntiNeutronInelastic* theHEAntiNeutronModel;

  // Lambda
  G4LambdaInelasticProcess  theLambdaInelastic;
  G4LELambdaInelastic*  theLELambdaModel;
  G4HELambdaInelastic*  theHELambdaModel;

  //  G4DecayTable* table;
  //  G4VDecayChannel* mode;

  
  // AntiLambda
  G4AntiLambdaInelasticProcess  theAntiLambdaInelastic;
  G4LEAntiLambdaInelastic*  theLEAntiLambdaModel;
  G4HEAntiLambdaInelastic*  theHEAntiLambdaModel;
  
  // SigmaMinus
  G4SigmaMinusInelasticProcess  theSigmaMinusInelastic;
  G4LESigmaMinusInelastic*  theLESigmaMinusModel;
  G4HESigmaMinusInelastic*  theHESigmaMinusModel;
  
  // AntiSigmaMinus
  G4AntiSigmaMinusInelasticProcess  theAntiSigmaMinusInelastic;
  G4LEAntiSigmaMinusInelastic*  theLEAntiSigmaMinusModel;
  G4HEAntiSigmaMinusInelastic*  theHEAntiSigmaMinusModel;
   
  // SigmaPlus
  G4SigmaPlusInelasticProcess  theSigmaPlusInelastic;
  G4LESigmaPlusInelastic*  theLESigmaPlusModel;
  G4HESigmaPlusInelastic*  theHESigmaPlusModel;
  
  // AntiSigmaPlus
  G4AntiSigmaPlusInelasticProcess  theAntiSigmaPlusInelastic;
  G4LEAntiSigmaPlusInelastic*  theLEAntiSigmaPlusModel;
  G4HEAntiSigmaPlusInelastic*  theHEAntiSigmaPlusModel;

};

#endif









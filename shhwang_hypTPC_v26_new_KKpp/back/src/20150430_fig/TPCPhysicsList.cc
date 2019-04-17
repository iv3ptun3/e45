// ====================================================================
//   TPCPhysicsList.cc
//
// ====================================================================
#include "TPCPhysicsList.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include "TPCEMPhysics.hh"
#include "TPCIonPhysics.hh"
#include "TPCGeneralPhysics.hh"
#include "G4BaryonConstructor.hh"
#include "G4ProcessManager.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "common.hh"

////////////////////////////////
TPCPhysicsList::TPCPhysicsList(): G4VUserPhysicsList()
////////////////////////////////
{
  defaultCutValue = 0.001*mm;
  SetVerboseLevel(1);
}
/////////////////////////////////
//TPCPhysicsList::TPCPhysicsList(): G4VUserPhysicsList()
//{

//}
/////////////////////////////////
TPCPhysicsList::~TPCPhysicsList()
/////////////////////////////////
{;}

void TPCPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructShortLived();
  ConstructStableHyperons();
  ConstructIons();
}

//////////////////////////////////////
void TPCPhysicsList::ConstructBosons()
//////////////////////////////////////
{
  // pseudo-particles
  //  G4Geantino::GeantinoDefinition();
  //  G4ChargedGeantino::ChargedGeantinoDefinition();
  // gamma
  G4Gamma::GammaDefinition();
}

///////////////////////////////////////
void TPCPhysicsList::ConstructLeptons()
///////////////////////////////////////
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}


//////////////////////////////////////
void TPCPhysicsList::ConstructMesons()
//////////////////////////////////////
{
  //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}

///////////////////////////////////////
void TPCPhysicsList::ConstructBaryons()
///////////////////////////////////////
{
//  barions
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}

#include "G4ShortLivedConstructor.hh"

///////////////////////////////////////
void TPCPhysicsList::ConstructShortLived()
///////////////////////////////////////
{
//  short lived particles
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

void TPCPhysicsList::ConstructProcess()
{
  // Define transportation process

  AddTransportation();
  //ConstructEM();
  ConstructGeneral(); //for test
  ConstructHadron();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"




///////////////////////////////////////
void TPCPhysicsList::ConstructIons()
///////////////////////////////////////
{
//  short lived particles
  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();  


  G4DecayTable* decayTable;
  G4VDecayChannel* mode;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Carbon12;
  G4int Z = 6, A = 12;
  G4double ionCharge   = 0.*eplus;
  G4double excitEnergy = 0.*keV;
  Carbon12
    = G4ParticleTable::GetParticleTable()->GetIon(Z,A,excitEnergy);
  //  Carbon12->DumpTable();
  //  C12->DumpTable();
  //  G4double rmC12=G4IonTable::GetNucleusMass( 6,12,0);
  //  G4ParticleDefinition* C12 = G4ParticleTable::G4IonTable::GetLightIon(6,12);
  //  G4ParticleDefinition* C12=G4IonTable::GetParticleTable()->GetIon(6,12,0.);
  //  G4ParticleDefinition* Be10=G4ParticleTable::FindIon(4,10,0.,0);
  //  G4double rmBe10= *G4IonTable::GetIonMass(4,10)/GeV;


  G4ParticleDefinition* kaonMinus=G4ParticleTable::GetParticleTable()->FindParticle("kaon-");
  //  G4ParticleDefinition* C12;
  //  C12 = G4ParticleTable::GetParticleTable()->GetIon(6, 12,0.*keV);
  //  G4ParticleDefinition* C12=G4ParticleTable::GetParticleTable()->FindParticle("C12[0,0]");
  G4double rmkn=kaonMinus->GetPDGMass()/GeV;
  G4double rmC12= 11.1749*GeV;
  //  G4double rmC12= C12->GetPDGMass()/GeV;

  G4String char_beam_mom = getenv("Beam_mom");

  G4ParticleDefinition* particle;
  G4double pbeam=atof(char_beam_mom.c_str());
  G4double Ebeam = sqrt(pow(pbeam,2)+pow(rmkn,2));
  G4double W=sqrt(pow(Ebeam+rmC12/GeV,2)-pow(pbeam,2));
  G4cout<<"---------------------------------"<<G4endl;
  G4cout<<"---------------------------------"<<G4endl;
  G4cout<<"---------------------------------"<<G4endl;
  G4cout<<"W:"<<W<<G4endl;
  G4cout<<"Ebeam:"<<Ebeam<<G4endl;
  G4cout<<"pbeam:"<<pbeam<<G4endl;
  G4cout<<"rmC12:"<<rmC12/GeV<<G4endl;
  G4cout<<"---------------------------------"<<G4endl;
  G4cout<<"---------------------------------"<<G4endl;
  G4cout<<"---------------------------------"<<G4endl;

  particle = new G4ParticleDefinition("phaseLL", W*GeV, 0.*GeV, 0,
				      0,              +0,             0,
				      0,              +0,             0,
				      "ion",        0,            +0,        101060120,
				      false,  0.*ns,          NULL);

  decayTable =  new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("phaseLL", 1.0,4,"lambda","lambda","Li8[0.0]","kaon+");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);





}


//////////////////////////////////
void TPCPhysicsList::ConstructEM()
//////////////////////////////////
{
  theParticleIterator-> reset();
  
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    if (particleName == "gamma") {
      // gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
      
    } 
    else if (particleName == "e-") {
      //electron
      G4VProcess* theeminusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeminusIonisation         = new G4eIonisation();
      G4VProcess* theeminusBremsstrahlung     = new G4eBremsstrahlung();
      //
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);
      //      
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, 
				   idxAlongStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         
					 idxAlongStep,2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung,     idxPostStep,3);
      
    }
    else if (particleName == "e+") {
      //positron
      G4VProcess* theeplusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeplusIonisation         = new G4eIonisation();
      G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
      G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();
      //
      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);
      //
      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,         idxAlongStep,2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(theeplusAnnihilation,       idxPostStep,4);
  
    }
    else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      G4VProcess* aMultipleScattering = new G4MuMultipleScattering();
      G4VProcess* aBremsstrahlung     = new G4MuBremsstrahlung();
      G4VProcess* aPairProduction     = new G4MuPairProduction();
      G4VProcess* anIonisation        = new G4MuIonisation();
      //
      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      pmanager->AddProcess(aBremsstrahlung);
      pmanager->AddProcess(aPairProduction);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
      pmanager->SetProcessOrdering(aBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(aPairProduction,     idxPostStep,4);

     }
    
    else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino     
     G4VProcess* aMultipleScattering = new G4hMultipleScattering();
     G4VProcess* anIonisation        = new G4hIonisation();
     //
     // add processes
     pmanager->AddProcess(anIonisation);
     pmanager->AddProcess(aMultipleScattering);
     //
     // set ordering for AlongStepDoIt
     // pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
     // pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
     pmanager->SetProcessOrderingToLast(anIonisation, idxAlongStep);
     //
     // set ordering for PostStepDoIt
     // pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
     // pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
     pmanager->SetProcessOrderingToLast(anIonisation, idxPostStep);
    }

    
  }

  
}

#include "G4Decay.hh"

///////////////////////////////////////
void TPCPhysicsList::ConstructGeneral()
///////////////////////////////////////
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();


  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    G4String particleName = particle->GetParticleName();
    if(particleName != "kaon+"){
      if (theDecayProcess->IsApplicable(*particle)) { 
	pmanager ->AddProcess(theDecayProcess);
	// set ordering for PostStepDoIt and AtRestDoIt
	pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
	pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
      }
    }
    else if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }

  ////shhwang; include decay process
    //comment out ichikawa
    // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    // G4VDecayChannel* mode;
    // G4DecayTable* Table = new G4DecayTable();
    // particle=particleTable->FindParticle("lambda");
    // mode = new G4PhaseSpaceDecayChannel("lambda",1.0000,2,"proton","pi-");
    // Table->Insert(mode);
    // particle->SetDecayTable(Table);
  //////


    ////k0
    //comment out ichikawa
    // G4VDecayChannel* mode1;
    // G4DecayTable* Table1 = new G4DecayTable();
    // particle=particleTable->FindParticle("kaon0S");
    // mode1 = new G4PhaseSpaceDecayChannel("kaon0S",1.0000,2,"pi+","pi-");
    // Table1->Insert(mode1);
    // particle->SetDecayTable(Table1);

    /*
    ////Decay mode control of lambda 1405. L1405 --> Lambda + gamma
    //    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    particle=particleTable->FindParticle("lambda(1405)");
    //    G4VDecayChannel* mode;
    //    G4DecayTable* Table = new G4DecayTable();
    mode = new G4PhaseSpaceDecayChannel("lambda(1405)",1.0000,2,"lambda","gamma");
    Table->Insert(mode);
    particle->SetDecayTable(Table);

    ////Decay mode control of sigma 1385. S1385 --> Lambda + gamma
    //    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    particle=particleTable->FindParticle("sigma(1385)");
    //    G4VDecayChannel* mode;
    //    G4DecayTable* Table = new G4DecayTable();
    mode = new G4PhaseSpaceDecayChannel("sigma(1385)",1.0000,2,"lambda","gamma");
    Table->Insert(mode);
    particle->SetDecayTable(Table);

  ////shhwang
  //  particleTable = G4ParticleTable::GetParticleTable();
  //  particleDef = particleTable->FindParticle("lambda");
  //  mode = new G4PhaseSpaceDecayChannel("lambda",1.00000,2,"proton","pi-");
  //  table = new G4DecayTable();
  //  table->Insert(mode);
  ///
  */

  }
}


void TPCPhysicsList::ConstructHadron()
{
  G4ProcessManager * pManager = 0;
  
  // Elastic Process
  theElasticModel = new G4LElastic();
  theElasticProcess.RegisterMe(theElasticModel);

  // PionPlus
  pManager = G4PionPlus::PionPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEPionPlusModel = new G4LEPionPlusInelastic();
  thePionPlusInelastic.RegisterMe(theLEPionPlusModel);
  pManager->AddDiscreteProcess(&thePionPlusInelastic);

  // PionMinus
  pManager = G4PionMinus::PionMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEPionMinusModel = new G4LEPionMinusInelastic();
  thePionMinusInelastic.RegisterMe(theLEPionMinusModel);
  pManager->AddDiscreteProcess(&thePionMinusInelastic);

  //KaonPlus
   pManager = G4KaonPlus::KaonPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonPlusModel = new G4LEKaonPlusInelastic();
  theHEKaonPlusModel = new G4HEKaonPlusInelastic();
  theKaonPlusInelastic.RegisterMe(theLEKaonPlusModel);
  theKaonPlusInelastic.RegisterMe(theHEKaonPlusModel);
  pManager->AddDiscreteProcess(&theKaonPlusInelastic);

  //  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //  G4ParticleDefinition* particle = particleTable->FindParticle("kaon+");

  //KaonMinus
  pManager = G4KaonMinus::KaonMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonMinusModel = new G4LEKaonMinusInelastic();
  theHEKaonMinusModel = new G4HEKaonMinusInelastic();
  theKaonMinusInelastic.RegisterMe(theLEKaonMinusModel);
  theKaonMinusInelastic.RegisterMe(theHEKaonMinusModel);
  pManager->AddDiscreteProcess(&theKaonMinusInelastic);

  // KaonZeroL
  pManager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonZeroLModel = new G4LEKaonZeroLInelastic();
  theHEKaonZeroLModel = new G4HEKaonZeroInelastic();
  theKaonZeroLInelastic.RegisterMe(theLEKaonZeroLModel);
  theKaonZeroLInelastic.RegisterMe(theHEKaonZeroLModel);
  pManager->AddDiscreteProcess(&theKaonZeroLInelastic);
 
  // KaonZeroS
  pManager = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonZeroSModel = new G4LEKaonZeroSInelastic();
  theHEKaonZeroSModel = new G4HEKaonZeroInelastic();
  theKaonZeroSInelastic.RegisterMe(theLEKaonZeroSModel);
  theKaonZeroSInelastic.RegisterMe(theHEKaonZeroSModel);
  pManager->AddDiscreteProcess(&theKaonZeroSInelastic);

  // Proton
  pManager = G4Proton::Proton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEProtonModel = new G4LEProtonInelastic();
  theHEProtonModel = new G4HEProtonInelastic();
  theProtonInelastic.RegisterMe(theLEProtonModel);
  theProtonInelastic.RegisterMe(theHEProtonModel);
  pManager->AddDiscreteProcess(&theProtonInelastic);

  // Neutron
  pManager = G4Neutron::Neutron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLENeutronModel = new G4LENeutronInelastic();
  theHENeutronModel = new G4HENeutronInelastic();
  theNeutronInelastic.RegisterMe(theLENeutronModel);
  theNeutronInelastic.RegisterMe(theHENeutronModel);
  pManager->AddDiscreteProcess(&theNeutronInelastic);

  // anti-Proton
  pManager = G4AntiProton::AntiProton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiProtonModel = new G4LEAntiProtonInelastic();
  theHEAntiProtonModel = new G4HEAntiProtonInelastic();
  theAntiProtonInelastic.RegisterMe(theLEAntiProtonModel);
  theAntiProtonInelastic.RegisterMe(theHEAntiProtonModel);
  pManager->AddDiscreteProcess(&theAntiProtonInelastic);

  // AntiNeutron
  pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiNeutronModel = new G4LEAntiNeutronInelastic();
  theHEAntiNeutronModel = new G4HEAntiNeutronInelastic();
  theAntiNeutronInelastic.RegisterMe(theLEAntiNeutronModel);
  theAntiNeutronInelastic.RegisterMe(theHEAntiNeutronModel);
  pManager->AddDiscreteProcess(&theAntiNeutronInelastic);

  // Lambda
  pManager = G4Lambda::Lambda()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);


  theLELambdaModel = new G4LELambdaInelastic();
  theHELambdaModel = new G4HELambdaInelastic();
  theLambdaInelastic.RegisterMe(theLELambdaModel);
  theLambdaInelastic.RegisterMe(theHELambdaModel);
  pManager->AddDiscreteProcess(&theLambdaInelastic);


  
  // AntiLambda
  pManager = G4AntiLambda::AntiLambda()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiLambdaModel = new G4LEAntiLambdaInelastic();
  theHEAntiLambdaModel = new G4HEAntiLambdaInelastic();
  theAntiLambdaInelastic.RegisterMe(theLEAntiLambdaModel);
  theAntiLambdaInelastic.RegisterMe(theHEAntiLambdaModel);
  pManager->AddDiscreteProcess(&theAntiLambdaInelastic);
    
  // SigmaMinus
  pManager = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLESigmaMinusModel = new G4LESigmaMinusInelastic();
  theHESigmaMinusModel = new G4HESigmaMinusInelastic();
  theSigmaMinusInelastic.RegisterMe(theLESigmaMinusModel);
  theSigmaMinusInelastic.RegisterMe(theHESigmaMinusModel);
  pManager->AddDiscreteProcess(&theSigmaMinusInelastic);

  // anti-SigmaMinus
  pManager = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiSigmaMinusModel = new G4LEAntiSigmaMinusInelastic();
  theHEAntiSigmaMinusModel = new G4HEAntiSigmaMinusInelastic();
  theAntiSigmaMinusInelastic.RegisterMe(theLEAntiSigmaMinusModel);
  theAntiSigmaMinusInelastic.RegisterMe(theHEAntiSigmaMinusModel);
  pManager->AddDiscreteProcess(&theAntiSigmaMinusInelastic);

  // SigmaPlus
  pManager = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLESigmaPlusModel = new G4LESigmaPlusInelastic();
  theHESigmaPlusModel = new G4HESigmaPlusInelastic();
  theSigmaPlusInelastic.RegisterMe(theLESigmaPlusModel);
  theSigmaPlusInelastic.RegisterMe(theHESigmaPlusModel);
  pManager->AddDiscreteProcess(&theSigmaPlusInelastic);

  // anti-SigmaPlus
  pManager = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiSigmaPlusModel = new G4LEAntiSigmaPlusInelastic();
  theHEAntiSigmaPlusModel = new G4HEAntiSigmaPlusInelastic();
  theAntiSigmaPlusInelastic.RegisterMe(theLEAntiSigmaPlusModel);
  theAntiSigmaPlusInelastic.RegisterMe(theHEAntiSigmaPlusModel);
  pManager->AddDiscreteProcess(&theAntiSigmaPlusInelastic);

}

void TPCPhysicsList::ConstructStableHyperons()
{
  G4DecayTable* decayTable;
  G4VDecayChannel* mode;
  G4ParticleDefinition* particle;

  // ssigma+ non-decay sigma+
  particle = new G4ParticleDefinition(
           "ssigma+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   true,                0,          NULL);


  // sigma1+  decay only to sigma+ -> pi+ neutron channel
  particle = new G4ParticleDefinition(
           "sigma1+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0.0799*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  // sigma+ -> neutron + pi+
  mode = new G4PhaseSpaceDecayChannel("sigma1+",1.0,2,"neutron","pi+");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // sigma2+  decay only to sigma+ -> pi0 proton channel
  particle = new G4ParticleDefinition(
           "sigma2+",    1.18937*GeV, 8.209e-12*MeV,       eplus,
                    1,              +1,             0,
                    2,              +2,             0,
             "baryon",               0,            +1,        3222,
	   false,                0.0799*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  // sigma+ -> proton + pi0
  mode = new G4PhaseSpaceDecayChannel("sigma2+", 1.0,2,"proton","pi0");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  /// Lambda1405 radioactive decay  
  particle = new G4ParticleDefinition(
           "lambda1405r",    1.4051*GeV, 50.*MeV,       0,
                    1,              -1,             0,
                    0,              +0,             0,
             "baryon",               0,            +1,        13122,
	   false,                0.*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("lambda1405r", 1.0,2,"lambda","gamma");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  /// Sigma1385 radioactive decay  
  particle = new G4ParticleDefinition(
           "sigma1385r",    1.3837*GeV, 36.*MeV,       0,
                    3,              +1,             0,
                    1,              +0,             0,
             "baryon",               0,            +1,        3214,
	   false,                0.*ns,          NULL);
  
  decayTable =  new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("sigma1385r", 1.0,2,"lambda","gamma");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);
  //		G4double mass, G4double	width,G4double charge,
  //		G4int iSpin, G4int iParity,G4int iConjugation,
  //		G4int iIsospin, G4int	iIsospinZ, G4int gParity,
  //		const G4String &pType, G4int lepton, G4int baryon, G4int encoding,
  //		G4bool 	stable,	G4double lifetime

  //, G4DecayTable *decaytable,
  //		G4bool  	shortlived = false,
  //		const G4String &  	subType = "",
  //		G4int  	anti_encoding = 0,
  //		G4double  	magneticMoment = 0.0	 
  /*
G4ParticleDefinition::G4ParticleDefinition 	( 	const G4String &  	aName,
		G4double  	mass,
		G4double  	width,
		G4double  	charge,
		G4int  	iSpin,
		G4int  	iParity,
		G4int  	iConjugation,
		G4int  	iIsospin,
		G4int  	iIsospinZ,
		G4int  	gParity,
		const G4String &  	pType,
		G4int  	lepton,
		G4int  	baryon,
		G4int  	encoding,
		G4bool  	stable,
		G4double  	lifetime,
		G4DecayTable *  	decaytable,
		G4bool  	shortlived = false,
		const G4String &  	subType = "",
		G4int  	anti_encoding = 0,
		G4double  	magneticMoment = 0.0	 
	) 	
*/


  //hybrid baryon mode
  //  particle = new G4ParticleDefinition(
  //				      "hybridb",  1.22*GeV, 0.*MeV, eplus,
  //				      0,              +0,             0,
  //				      0,              +0,             0,
  //				      "baryon",        0,            +1,        9223,
  //				      false,      0.0*ns,          NULL);
  //  decayTable =  new G4DecayTable();


  //  mode = new G4PhaseSpaceDecayChannel("hybridb", 1.0,3,"neutron","pi+","pi-"); //pi- p --> n pi+ pi-
  //  mode = new G4PhaseSpaceDecayChannel("hybridb", 1.0,3,"neutron","pi+","pi-"); //pi- p --> p pi0 pi-
  //  mode = new G4PhaseSpaceDecayChannel("hybridb", 1.0,3,"neutron","pi+","pi-"); //pi+ p --> n pi+ pi+
  //  mode = new G4PhaseSpaceDecayChannel("hybridb", 1.0,3,"proton","pi+","pi+"); //pi+ p --> p pi0 pi+
  //  mode = new G4PhaseSpaceDecayChannel("hybridb", 1.0,3,"proton","pi+","pi-"); //pi+ p --> p pi- pi+
  //  decayTable->Insert(mode);
  //  particle->SetDecayTable(decayTable);


  G4String h_decaytime = getenv("H_decaytime");
  G4String h_mass = getenv("H_mass");
  G4String h_width = getenv("H_width");
  //  atof(weak_decay_time.c_str())                                          

  particle = new G4ParticleDefinition(
				      "hdibaryon",  atof(h_mass.c_str())*GeV, atof(h_width.c_str())*GeV, 0,
				      0,              +0,             0,
				      0,              +0,             0,
				      "baryon",        0,            +2,        9223,
				      false,  atof(h_decaytime.c_str())*ns,          NULL);
  decayTable =  new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("hdibaryon", 1.0,3,"lambda","proton","pi-"); //pi+ p --> p pi- pi+
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);


  particle = new G4ParticleDefinition(
				      "hdibaryonS",  atof(h_mass.c_str())*GeV, atof(h_width.c_str())*GeV, 0,
				      0,              +0,             0,
				      0,              +0,             0,
				      "baryon",        0,            +2,        9224,
				      false,  atof(h_decaytime.c_str())*ns,          NULL);
  decayTable =  new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("hdibaryonS", 1.0,2,"sigma-","proton"); //pi+ p --> p pi- pi+
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  particle = new G4ParticleDefinition(
				      "hdibaryonLL", atof(h_mass.c_str())*GeV, atof(h_width.c_str())*GeV, 0,
				      0,              +0,             0,
				      0,              +0,             0,
				      "baryon",        0,            +2,        9225,
				      false,  atof(h_decaytime.c_str())*ns,          NULL);
  decayTable =  new G4DecayTable();
  mode = new G4PhaseSpaceDecayChannel("hdibaryonLL", 1.0,2,"lambda","lambda"); 
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);





}     



//////////////////////////////
void TPCPhysicsList::SetCuts()
//////////////////////////////
{
  // G4VUserPhysicsList::SetCutsWithDefault" method sets
  // the default cut value for all particle types
  SetCutsWithDefault();
}

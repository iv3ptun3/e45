// -*- C++ -*-

#include "TPCPhysicsList.hh"
#include "TPCDecayChannel.hh"
#include <G4ComptonScattering.hh>
#include <G4Decay.hh>
#include <G4DecayTable.hh>
#include <G4eBremsstrahlung.hh>
#include <G4eIonisation.hh>
#include <G4eMultipleScattering.hh>
#include <G4eplusAnnihilation.hh>
#include <G4GammaConversion.hh>
#include <G4hIonisation.hh>
#include <G4hMultipleScattering.hh>
#include <G4IonTable.hh>
#include <G4MuBremsstrahlung.hh>
#include <G4MuIonisation.hh>
#include <G4MuMultipleScattering.hh>
#include <G4MuPairProduction.hh>
#include <G4ParticleTable.hh>
#include <G4PhaseSpaceDecayChannel.hh>
#include <G4PhotoElectricEffect.hh>
#include <G4ProcessManager.hh>
// Particle Constructor
#include <G4BaryonConstructor.hh>
#include <G4BosonConstructor.hh>
#include <G4IonConstructor.hh>
#include <G4LeptonConstructor.hh>
#include <G4MesonConstructor.hh>
#include <G4ShortLivedConstructor.hh>
// Phycics Process
#include <G4EmStandardPhysics.hh>
#include <G4EmExtraPhysics.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4StoppingPhysics.hh>
#include <G4IonPhysics.hh>
#include <G4NeutronTrackingCut.hh>
#include <G4HadronPhysicsQGSP_BERT.hh>

#include "ConfMan.hh"

namespace
{
  using CLHEP::eplus;
  using CLHEP::GeV;
  using CLHEP::keV;
  using CLHEP::MeV;
  using CLHEP::ns;
  const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
TPCPhysicsList::TPCPhysicsList( void )
  : G4VUserPhysicsList(),
    m_em_physics_list(),
    m_hadron_physics_list()
{
  SetDefaultCutValue( 2.*CLHEP::mm );
  m_em_physics_list = new G4EmStandardPhysics( verboseLevel );
}

//_____________________________________________________________________________
TPCPhysicsList::~TPCPhysicsList( void )
{
  delete m_em_physics_list;
  for( auto&& phys : m_hadron_physics_list ){
    delete phys;
  }
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructParticle( void )
{
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructShortLived();
  ConstructStableHyperons();
  ConstructIons();
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructBosons( void )
{
  G4BosonConstructor pBosonConstructor;
  pBosonConstructor.ConstructParticle();
  // G4Geantino::GeantinoDefinition();
  // G4ChargedGeantino::ChargedGeantinoDefinition();
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructLeptons( void )
{
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructMesons( void )
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructBaryons( void )
{
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructShortLived( void )
{
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructProcess( void )
{
  AddTransportation();
  if( gConf.Get<G4bool>( "EM" ) )
    ConstructEM();
  if( gConf.Get<G4bool>( "Decay" ) )
    ConstructGeneral(); // for test
  if( gConf.Get<G4bool>( "Hadron" ) )
    ConstructHadron();
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructIons( void )
{
  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // G4int Z = 6, A = 12;
  // G4double ionCharge   = 0.*eplus;
  // G4double excitEnergy = 0.*keV;
  // auto Carbon12 = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
  //  Carbon12->DumpTable();
  //  C12->DumpTable();
  //  G4double rmC12=G4IonTable::GetNucleusMass( 6,12,0);
  //  G4ParticleDefinition* C12 = G4ParticleTable::G4IonTable::GetLightIon(6,12);
  //  G4ParticleDefinition* C12=G4IonTable::GetParticleTable()->GetIon(6,12,0.);
  //  G4ParticleDefinition* Be10=G4ParticleTable::FindIon(4,10,0.,0);
  //  G4double rmBe10= *G4IonTable::GetIonMass(4,10)/GeV;

  auto kaonMinus = G4ParticleTable::GetParticleTable()->FindParticle("kaon-");
  // auto C12 = G4ParticleTable::GetParticleTable()->GetIon(6, 12,0.*keV);
  // auto C12 = G4ParticleTable::GetParticleTable()->FindParticle("C12[0,0]");
  G4double rmkn = kaonMinus->GetPDGMass()/GeV;
  G4double rmC12 = 11.1749*GeV;
  //  G4double rmC12= C12->GetPDGMass()/GeV;
  G4double pbeam = gConf.Get<G4double>( "BeamMom" );
  G4double Ebeam = std::sqrt( std::pow( pbeam, 2 ) + std::pow( rmkn, 2 ) );
  G4double W = std::sqrt( std::pow( Ebeam+rmC12/GeV, 2 ) -
			  std::pow( pbeam, 2 ) );
  G4cout << "---------------------------------" << G4endl;
  G4cout << "---------------------------------" << G4endl;
  G4cout << "---------------------------------" << G4endl;
  G4cout << "W:" << W << G4endl;
  G4cout << "Ebeam:" << Ebeam << G4endl;
  G4cout << "pbeam:" << pbeam << G4endl;
  G4cout << "rmC12:" << rmC12/GeV << G4endl;
  G4cout << "---------------------------------" << G4endl;
  G4cout << "---------------------------------" << G4endl;
  G4cout << "---------------------------------" << G4endl;
  auto particle = new G4ParticleDefinition( "phaseLL", W*GeV, 0.*GeV, 0,
					    0, +0, 0,
					    0, +0, 0,
					    "ion", 0, +0, 101060120,
					    false, 0.*ns, nullptr );
  auto decayTable = new G4DecayTable;
  auto mode = new G4PhaseSpaceDecayChannel( "phaseLL", 1.0, 4, "lambda",
					    "lambda", "Li8[0.0]", "kaon+" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructEM( void )
{
  m_em_physics_list->ConstructProcess();
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructGeneral( void )
{
	G4cout<<"Ennabling Decays"<<G4endl;
  auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();
  auto theDecayProcess = new G4Decay;
  while( (*theParticleIterator)() ){
    auto particle = theParticleIterator->value();
    auto pmanager = particle->GetProcessManager();
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
    G4int Lambda_decay = gConf.Get<G4int>("LambdaDecay");
		G4int PolarizedDecay = gConf.Get<G4int>("PolarizedDecay");
		double XiPolarization = gConf.Get<G4double>("XiPolarization");
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    if(1){
      G4DecayTable* Table = new G4DecayTable();
      particle=particleTable->FindParticle("lambda");
			if(PolarizedDecay){
				TPCPolarizedDecayChannel* LdDecay;
				double DecayPar[2] = {0.732,-6.5};//Alpha,phi(in degree)
				LdDecay = new TPCPolarizedDecayChannel("lambda",1,DecayPar,1,"proton","pi-");
				Table->Insert(LdDecay);
			}
			else{
				G4VDecayChannel* mode;
				mode = new G4PhaseSpaceDecayChannel("lambda",1.0000,2,"proton","pi-");
				Table->Insert(mode);
			}
      particle->SetDecayTable(Table);
    }
    //////Xi-////
		if(1){
			G4DecayTable* Table = new G4DecayTable();
			particle = particleTable ->FindParticle("xi-");
			if(PolarizedDecay){
				TPCPolarizedDecayChannel* XiDecay;
				double DecayPar[2] = {-0.401,-2.1};//Alpha,phi(in degree)
				XiDecay = new TPCPolarizedDecayChannel("xi-",1,DecayPar,0,"lambda","pi-");
				XiDecay->SetPolarization(XiPolarization);
				Table->Insert(XiDecay);
			}
			else{
				G4VDecayChannel* mode;
				mode = new G4PhaseSpaceDecayChannel("xi-",1.0000,2,"lambda","pi-");
				Table->Insert(mode);
			}
      particle->SetDecayTable(Table);
		}
		////Xi0////
		if(1){
			G4VDecayChannel* mode;
			G4DecayTable* Table = new G4DecayTable();
			particle = particleTable ->FindParticle("xi0");
      mode = new G4PhaseSpaceDecayChannel("xi0",1.0000,2,"lambda","pi0");
      Table->Insert(mode);
      particle->SetDecayTable(Table);
		}
		////Xi1530////
    G4int XiStarToXiPi0 = gConf.Get<G4int>("XiStarDecayMode");
		if(XiStarToXiPi0==0){
			G4VDecayChannel* mode;
			G4DecayTable* Table = new G4DecayTable();
			particle = particleTable ->FindParticle("xi(1530)-");
      mode = new G4PhaseSpaceDecayChannel("xi(1530)-",1.0000,2,"xi-","pi0");
      Table->Insert(mode);
      particle->SetDecayTable(Table);
		}
		else if (XiStarToXiPi0==1){
			G4VDecayChannel* mode;
			G4DecayTable* Table = new G4DecayTable();
			particle = particleTable ->FindParticle("xi(1530)-");
      mode = new G4PhaseSpaceDecayChannel("xi(1530)-",1.0000,2,"xi0","pi-");
      Table->Insert(mode);
      particle->SetDecayTable(Table);
		}
		else{
			G4VDecayChannel* mode;
			G4DecayTable* Table = new G4DecayTable();
			particle = particleTable ->FindParticle("xi(1530)-");
      mode = new G4PhaseSpaceDecayChannel("xi(1530)-",1.0000,2,"xi-","gamma");
      Table->Insert(mode);
      particle->SetDecayTable(Table);
		}
		////
    ////k0
    G4int Ks_decay = gConf.Get<G4int>( "KsDecay" );
    if(Ks_decay==1){
      G4VDecayChannel* mode1;
      G4DecayTable* Table1 = new G4DecayTable();
      particle=particleTable->FindParticle("kaon0S");
      mode1 = new G4PhaseSpaceDecayChannel("kaon0S",1.0000,2,"pi+","pi-");
      Table1->Insert(mode1);
      particle->SetDecayTable(Table1);
    }
		//phi meson
		{
			particle = particleTable->FindParticle("phi");
      G4VDecayChannel* mode;
      G4DecayTable* Table = new G4DecayTable();
			mode = new G4PhaseSpaceDecayChannel("phi",1,2,"kaon+","kaon-");	
      Table->Insert(mode);
			particle -> SetDecayTable(Table);
		}
		// H-dibaryon
    {
      G4VDecayChannel* mode;
      G4DecayTable* Table = new G4DecayTable();
      particle=particleTable->FindParticle("hdibaryonLL");
      mode = new G4PhaseSpaceDecayChannel("hdibaryonLL",1.0000,2,"lambda","lambda");
      Table->Insert(mode);
      particle->SetDecayTable(Table);
    }
    {
      G4VDecayChannel* mode;
      G4DecayTable* Table = new G4DecayTable();
      particle=particleTable->FindParticle("hdibaryon");
      mode = new G4PhaseSpaceDecayChannel( "hdibaryon", 1.0, 3, "lambda",
	  "proton", "pi-" );
      Table->Insert( mode );
      particle->SetDecayTable( Table );
    }
    {
      G4VDecayChannel* mode;
      G4DecayTable* Table = new G4DecayTable();
      particle=particleTable->FindParticle("hdibaryonS");
       mode = new G4PhaseSpaceDecayChannel( "hdibaryonS", 1.0, 2,
				       "sigma-", "proton" );
      Table->Insert( mode );
      particle->SetDecayTable( Table );
    }
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

//_____________________________________________________________________________
void TPCPhysicsList::ConstructHadron( void )
{
  m_hadron_physics_list.push_back( new G4EmExtraPhysics );
  m_hadron_physics_list.push_back( new G4HadronElasticPhysics );
  m_hadron_physics_list.push_back( new G4StoppingPhysics );
  m_hadron_physics_list.push_back( new G4IonPhysics );
  m_hadron_physics_list.push_back( new G4NeutronTrackingCut );
  m_hadron_physics_list.push_back( new G4HadronPhysicsQGSP_BERT );
  for( auto&& phys : m_hadron_physics_list ){
    phys->ConstructProcess();
  }
}

//_____________________________________________________________________________
void
TPCPhysicsList::ConstructStableHyperons( void )
{
  G4ParticleDefinition* particle;
  G4DecayTable* decayTable;
  G4VDecayChannel* mode;

  // G4ParticleDefinition::G4ParticleDefinition(
  // const G4String& aName, G4double mass, G4double width,
  // G4double charge, G4int iSpin, G4int iParity, G4int iConjugation,
  // G4int iIsospin, G4int  iIsospinZ, G4int gParity,
  // const G4String& pType, G4int lepton, G4int baryon, G4int encoding,
  // G4bool stable, G4double lifetime, G4DecayTable* decaytable,
  // G4bool shortlived = false,
  // const G4String& subType = "",
  // G4int anti_encoding = 0,

  // G4double magneticMoment = 0.0
  // )
  // G4DecayTable::G4DecayTable(
  // G4bool shortlived = false,
  // const G4String& subType = "",
  // G4int anti_encoding = 0,
  // G4double magneticMoment = 0.0
  // )
  
	// skaon- non-decay K-
  particle = new G4ParticleDefinition( "skaon-", 0.493677*GeV, 5.315e-14*MeV,
				       -eplus, 0, -1, 0, 1, 1, 0,
				       "meson", 0, 0, 321,
				       true, 0, nullptr );

  // skaon+ non-decay K+
  particle = new G4ParticleDefinition( "skaon+", 0.493677*GeV, 5.315e-14*MeV,
				       eplus, 0, -1, 0, 1, 1, 0,
				       "meson", 0, 0, 321,
				       true, 0, nullptr );

  // ssigma+ non-decay sigma+
  particle = new G4ParticleDefinition( "ssigma+", 1.18937*GeV, 8.209e-12*MeV,
				       eplus, 1, +1, 0, 2, +2, 0,
				       "baryon", 0, +1, 3222,
				       true, 0, nullptr );
  // sigma1+  decay only to sigma+ -> pi+ neutron channel
  particle = new G4ParticleDefinition( "sigma1+", 1.18937*GeV, 8.209e-12*MeV,
				       eplus, 1, +1, 0, 2, +2, 0,
				       "baryon", 0, +1, 3222,
				       false, 0.0799*ns, nullptr );
  decayTable =  new G4DecayTable;
  // sigma+ -> neutron + pi+
  mode = new G4PhaseSpaceDecayChannel( "sigma1+", 1.0, 2, "neutron", "pi+" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );
  // sigma2+  decay only to sigma+ -> pi0 proton channel
  particle = new G4ParticleDefinition( "sigma2+", 1.18937*GeV, 8.209e-12*MeV,
				       eplus, 1, +1, 0, 2, +2, 0,
				       "baryon", 0, +1, 3222,
				       false, 0.0799*ns, nullptr );
  decayTable =  new G4DecayTable;
  // sigma+ -> proton + pi0
  mode = new G4PhaseSpaceDecayChannel( "sigma2+", 1.0, 2, "proton", "pi0" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );
  /// Lambda1405 radioactive decay
  particle = new G4ParticleDefinition( "lambda1405r", 1.4051*GeV, 50.*MeV,
				       0, 1, -1, 0, 0, +0, 0,
				       "baryon", 0, +1, 13122,
				       false, 0.*ns, nullptr );
  decayTable =  new G4DecayTable;
  mode = new G4PhaseSpaceDecayChannel( "lambda1405r", 1.0, 2,
				       "lambda", "gamma" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );
  /// Sigma1385 radioactive decay
  particle = new G4ParticleDefinition( "sigma1385r", 1.3837*GeV, 36.*MeV,
				       0, 3, +1, 0, 1, +0, 0,
				       "baryon", 0, +1, 3214,
				       false, 0.*ns, nullptr );
  decayTable =  new G4DecayTable;
  mode = new G4PhaseSpaceDecayChannel( "sigma1385r", 1.0, 2, "lambda", "gamma" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );
  // hybrid baryon mode
  // particle = new G4ParticleDefinition( "hybridb", 1.22*GeV, 0.*MeV,
  // 				       eplus, 0, +0, 0, 0, +0, 0,
  // 				       "baryon", 0, +1, 9223,
  // 				       false, 0.*ns, nullptr );
  // decayTable = new G4DecayTable;
  // mode = new G4PhaseSpaceDecayChannel( "hybridb", 1.0, 3, "neutron",
  // 				       "pi+", "pi-" ); // pi- p --> n pi+ pi-
  // mode = new G4PhaseSpaceDecayChannel( "hybridb", 1.0, 3, "neutron",
  // 				       "pi+", "pi-" ); // pi- p --> p pi0 pi-
  // mode = new G4PhaseSpaceDecayChannel( "hybridb", 1.0, 3, "neutron",
  // 				       "pi+", "pi-" ); // pi+ p --> n pi+ pi+
  // mode = new G4PhaseSpaceDecayChannel( "hybridb", 1.0, 3, "proton",
  // 				       "pi+", "pi+" ); // pi+ p --> p pi0 pi+
  // mode = new G4PhaseSpaceDecayChannel( "hybridb", 1.0, 3, "proton",
  // 				       "pi+", "pi-" ); // pi+ p --> p pi- pi+
  // decayTable->Insert( mode );
  // particle->SetDecayTable( decayTable );

  G4double h_lifetime = gConf.Get<G4double>("HdibaryonLifetime")* ns;
  G4double h_mass = gConf.Get<G4double>("HdibaryonMass") * GeV;
  G4double h_width = gConf.Get<G4double>("HdibaryonWidth") * keV;
  particle = new G4ParticleDefinition( "hdibaryon", h_mass, h_width,
				       0, 0, +0, 0, 0, +0, 0,
				       "baryon", 0, +2, 9223,
				       false, h_lifetime, nullptr );
  particle = new G4ParticleDefinition( "hdibaryonS", h_mass, h_width,
				       0, 0, +0, 0, 0, +0, 0,
				       "baryon", 0, +2, 9224,
				       false, h_lifetime, nullptr );
  particle = new G4ParticleDefinition( "hdibaryonLL", h_mass, h_width,
				       0, 0, +0, 0, 0, +0, 0,
				       "baryon", 0, +2, 9225,
				       false, 0.000000001*ns, nullptr );
}

//_____________________________________________________________________________
void
TPCPhysicsList::SetCuts( void )
{
  // G4VUserPhysicsList::SetCutsWithDefault" method sets
  // the default cut value for all particle types
  SetCutsWithDefault();
}

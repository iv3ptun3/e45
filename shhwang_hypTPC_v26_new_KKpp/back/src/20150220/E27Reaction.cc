//-----------------------------------------------------------
// E27Reaction.cc 
// for the EventGeneration for the E27 experiment
//-----------------------------------------------------------
#include "E27Reaction.hh"
#include "TPCPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "Kinema3Resonance.hh"
#include "KinemaHResonance.hh"
#include "Kinema3Body.hh"
#include "Kinema4Body.hh"
#include "KinemaHybrid.hh"
#include "KinemaHweak.hh"
#include "KinemaFermi.hh"
#include "KinemaKstar.hh"
#include "common.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4IonConstructor.hh"

void E27Reaction::E27_beamthrough(G4Event* anEvent){
  //  G4double  momk[3], mom[3],momkn[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  //  G4double pbm[4];
  G4double Energy_pip,  mom_pip_x, mom_pip_y, mom_pip_z;
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pionPlus;
  //  kaonMinus = particleTable->FindParticle("kaon-");
  //kaonMinus = particleTable->FindParticle("pi-");
  pionPlus = particleTable->FindParticle("pi+");
  pbeam=CLHEP::RandGauss::shoot(pGen->Get_env_Beam_mom(),0.01294*pGen->Get_env_Beam_mom());
  //  pbeam=CLHEP::RandGauss::shoot(env_Beam_mom,env_Beam_mom*3.3*0.0001/2.3548);
  //  pbeam=1.8;
  mom_pip_x=0;
  mom_pip_y=0;
  mom_pip_z=pbeam;
  Energy_pip=sqrt(pionPlus->GetPDGMass()/GeV*pionPlus->GetPDGMass()/GeV+pbeam*pbeam);

  pGen->anaManager->SetPrimaryBeam(0,0,pbeam);

  Ebeam = sqrt(pbeam*pbeam+pionPlus->GetPDGMass()/GeV*pionPlus->GetPDGMass()/GeV);   


  //  G4double vtx = CLHEP::RandFlat::shoot(-15.,15.)*mm;
  //  G4double vty = CLHEP::RandFlat::shoot(-5.,5.)*mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(env_target_pos_z-env_target_width/2,env_target_pos_z+env_target_width/2)*mm;

  G4double vtx = CLHEP::RandGauss::shoot(0,10.)*mm;
  G4double vty = CLHEP::RandFlat::shoot(0.,3.2)*mm;
  G4double vtz= CLHEP::RandFlat::shoot(pGen->Get_env_target_pos_z()-pGen->Get_env_target_width()/2,pGen->Get_env_target_pos_z()+pGen->Get_env_target_width()/2)*mm-250.*mm;
  std::cout<<"pbeam = "<<pbeam<<std::endl;
  //getchar();
  //beam pi+
  pGen->particleGun->SetParticleDefinition(pionPlus);
  pGen->particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pip_x,mom_pip_y,mom_pip_z));
  pGen->particleGun->SetParticleEnergy((Energy_pip - pionPlus->GetPDGMass()/GeV)*GeV);
  pGen->particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  //  anaManager->SetNumberOfPrimaryParticle(1);
  //  anaManager->SetPrimaryParticle(0,mom_pip_x,mom_pip_y,mom_pip_z,pionPlus->GetPDGMass()/GeV);
  //  anaManager->SetPrimaryVertex(0,vtx,vty,vtz);
}


void E27Reaction::E27_Kptest(G4Event* anEvent){
  //  G4double  momk[3], mom[3],momkn[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  //  G4double pbm[4];
  G4double Energy_Kp,  mom_Kp_x, mom_Kp_y, mom_Kp_z;
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* KaonPlus;
  //  kaonMinus = particleTable->FindParticle("kaon-");
  //kaonMinus = particleTable->FindParticle("pi-");
  KaonPlus = particleTable->FindParticle("kaon+");
  pbeam=CLHEP::RandGauss::shoot(pGen->Get_env_Beam_mom(),0.01294*pGen->Get_env_Beam_mom());
  //  pbeam=CLHEP::RandGauss::shoot(env_Beam_mom,env_Beam_mom*3.3*0.0001/2.3548);
  //  pbeam=1.8;
  mom_Kp_x=0;
  mom_Kp_y=0;
  mom_Kp_z=pbeam;
  Energy_Kp=sqrt(KaonPlus->GetPDGMass()/GeV*KaonPlus->GetPDGMass()/GeV+pbeam*pbeam);

  pGen->anaManager->SetPrimaryBeam(0,0,pbeam);

  Ebeam = sqrt(pbeam*pbeam+KaonPlus->GetPDGMass()/GeV*KaonPlus->GetPDGMass()/GeV);   


  //  G4double vtx = CLHEP::RandFlat::shoot(-15.,15.)*mm;
  //  G4double vty = CLHEP::RandFlat::shoot(-5.,5.)*mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(env_target_pos_z-env_target_width/2,env_target_pos_z+env_target_width/2)*mm;

  // G4double vtx = CLHEP::RandGauss::shoot(0,10.)*mm;
  // G4double vty = CLHEP::RandFlat::shoot(0.,3.2)*mm;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz= CLHEP::RandFlat::shoot(pGen->Get_env_target_pos_z()-pGen->Get_env_target_width()/2,pGen->Get_env_target_pos_z()+pGen->Get_env_target_width()/2)*mm-250.*mm;
  std::cout<<"pbeam = "<<pbeam<<std::endl;
  // getchar();
  //scat K+
  pGen->particleGun->SetParticleDefinition(KaonPlus);
  pGen->particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_Kp_x,mom_Kp_y,mom_Kp_z));
  pGen->particleGun->SetParticleEnergy((Energy_Kp - KaonPlus->GetPDGMass()/GeV)*GeV);
  pGen->particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  //  anaManager->SetNumberOfPrimaryParticle(1);
  //  anaManager->SetPrimaryParticle(0,mom_Kp_x,mom_Kp_y,mom_Kp_z,pionPlus->GetPDGMass()/GeV);
  //  anaManager->SetPrimaryVertex(0,vtx,vty,vtz);
}


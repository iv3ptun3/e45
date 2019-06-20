// -*- C++ -*-

/**
 * E27Reaction.cc
 * for the EventGeneration for the E27 experiment
 */

#include "E27Reaction.hh"

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
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
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4ParticleTypes.hh"
#include "G4IonConstructor.hh"
#include "GeneratorHelper.hh"
#include "AngDisGenerator.hh"
#include "TPCAnaManager.hh"

namespace
{
  using CLHEP::mm;
  using CLHEP::GeV;
  TPCAnaManager& gAnaMan = TPCAnaManager::GetInstance();
  const int MaxTry = 1000;
  const double AtomicMassUnit = 0.9314932;
  const auto& gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
}

//_____________________________________________________________________________
// reactio No #2701 pi+ beam through
void
E27Reaction::E27_beamthrough(G4Event* anEvent)
{
  //  G4double  momk[3], mom[3],momkn[3];
  //  G4double rmk=0.493677;
  //  G4double pbm[4];
  G4double Energy_pip,  mom_pip_x, mom_pip_y, mom_pip_z;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pionPlus;
  //  kaonMinus = particleTable->FindParticle("kaon-");
  //kaonMinus = particleTable->FindParticle("pi-");
  pionPlus = particleTable->FindParticle("pi+");
  G4double pbeam = CLHEP::RandGauss::shoot( gConf.Get<G4double>("BeamMom"),
					    0.01294*gConf.Get<G4double>("BeamMom") );
  //  pbeam=CLHEP::RandGauss::shoot(env_Beam_mom,env_Beam_mom*3.3*0.0001/2.3548);
  //  pbeam=1.8;
  mom_pip_x=0;
  mom_pip_y=0;
  mom_pip_z=pbeam;
  Energy_pip=sqrt(pionPlus->GetPDGMass()/GeV*pionPlus->GetPDGMass()/GeV+pbeam*pbeam);

  gAnaMan.SetPrimaryBeam(0,0,pbeam);

  // G4double Ebeam = sqrt(pbeam*pbeam+pionPlus->GetPDGMass()/GeV*pionPlus->GetPDGMass()/GeV);

  //  G4double vtx = CLHEP::RandFlat::shoot(-15.,15.)*mm;
  //  G4double vty = CLHEP::RandFlat::shoot(-5.,5.)*mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(env_target_pos_z-env_target_width/2,env_target_pos_z+env_target_width/2)*mm;

  G4double vtx = CLHEP::RandGauss::shoot(0,10.)*mm;
  G4double vty = CLHEP::RandFlat::shoot(0.,3.2)*mm;
  G4double vtz = gSize.Get("Target", ThreeVector::Z );
  //G4double vtz= CLHEP::RandFlat::shoot(pGen->Get_env_target_pos_z()-gSize.Get( "Target", ThreeVector::Z )/2,pGen->Get_env_target_pos_z()+gSize.Get( "Target", ThreeVector::Z )/2)*mm-250.*mm;
  std::cout<<"pbeam = "<<pbeam<<std::endl;
  //getchar();
  //beam pi+
  pGen->particleGun->SetParticleDefinition(pionPlus);
  pGen->particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pip_x,mom_pip_y,mom_pip_z));
  pGen->particleGun->SetParticleEnergy((Energy_pip - pionPlus->GetPDGMass()/GeV)*GeV);
  pGen->particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  //  gAnaMan->SetNumberOfPrimaryParticle(1);
  //  gAnaMan->SetPrimaryParticle(0,mom_pip_x,mom_pip_y,mom_pip_z,pionPlus->GetPDGMass()/GeV);
  //  gAnaMan->SetPrimaryVertex(0,vtx,vty,vtz);
}

//_____________________________________________________________________________
//reactio No #2702 K+ gun for test
void
E27Reaction::E27_Kptest(G4Event* anEvent)
{
  //  G4double  momk[3], mom[3],momkn[3];
  //  G4double rmk=0.493677;
  //  G4double pbm[4];
  G4double Energy_Kp,  mom_Kp_x, mom_Kp_y, mom_Kp_z;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* KaonPlus;
  //  kaonMinus = particleTable->FindParticle("kaon-");
  //kaonMinus = particleTable->FindParticle("pi-");
  KaonPlus = particleTable->FindParticle("kaon+");
  G4double pbeam=CLHEP::RandGauss::shoot(gConf.Get<G4double>( "BeamMom" ),0.01294*gConf.Get<G4double>( "BeamMom" ));
  //  pbeam=CLHEP::RandGauss::shoot(env_Beam_mom,env_Beam_mom*3.3*0.0001/2.3548);
  //  pbeam=1.8;
  mom_Kp_x=0;
  mom_Kp_y=0;
  mom_Kp_z=pbeam;
  Energy_Kp=sqrt(KaonPlus->GetPDGMass()/GeV*KaonPlus->GetPDGMass()/GeV+pbeam*pbeam);

  gAnaMan.SetPrimaryBeam(0,0,pbeam);

  // G4double Ebeam = sqrt(pbeam*pbeam+KaonPlus->GetPDGMass()/GeV*KaonPlus->GetPDGMass()/GeV);

  //  G4double vtx = CLHEP::RandFlat::shoot(-15.,15.)*mm;
  //  G4double vty = CLHEP::RandFlat::shoot(-5.,5.)*mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(env_target_pos_z-env_target_width/2,env_target_pos_z+env_target_width/2)*mm;

  // G4double vtx = CLHEP::RandGauss::shoot(0,10.)*mm;
  // G4double vty = CLHEP::RandFlat::shoot(0.,3.2)*mm;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(gGeom.GetGlobalPosition( "Target" ).z()-gSize.Get( "Target", ThreeVector::Z )/2,gGeom.GetGlobalPosition( "Target" ).z()+gSize.Get( "Target", ThreeVector::Z )/2)*mm-250.*mm;
  G4double vtz= gGeom.GetGlobalPosition( "Target" ).z();
  std::cout<<"pbeam = "<<pbeam<<std::endl;
  // getchar();
  //scat K+
  pGen->particleGun->SetParticleDefinition(KaonPlus);
  pGen->particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_Kp_x,mom_Kp_y,mom_Kp_z));
  pGen->particleGun->SetParticleEnergy((Energy_Kp - KaonPlus->GetPDGMass()/GeV)*GeV);
  pGen->particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  //  gAnaMan->SetNumberOfPrimaryParticle(1);
  //  gAnaMan->SetPrimaryParticle(0,mom_Kp_x,mom_Kp_y,mom_Kp_z,pionPlus->GetPDGMass()/GeV);
  //  gAnaMan->SetPrimaryVertex(0,vtx,vty,vtz);
}

//_____________________________________________________________________________
//reactio No #2703 pi+ d -> K+ K-pp, K-pp -> Lambda p reaction
void
E27Reaction::E27_Kpp_F_LambdaP(G4Event* anEvent)
{
  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4Lambda::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg( gConf.Get<G4double>( "BeamX0" ),
					     gConf.Get<G4double>( "BeamY0" ),
					     gGeom.GetGlobalPosition( "Target" ).z(),
					     gConf.Get<G4double>( "BeamDX" ),
					     gConf.Get<G4double>( "BeamDY" ),
					     gSize.Get( "Target", ThreeVector::X ),
					     gSize.Get( "Target", ThreeVector::Z ));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV( gConf.Get<G4double>( "BeamU0" ),
						gConf.Get<G4double>( "BeamV0" ),
						gConf.Get<G4double>( "BeamDU" ),
						gConf.Get<G4double>( "BeamDV" ));

  G4double pb = gConf.Get<G4double>( "BeamMom" )*GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>( "BeamMom" )!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
  pb += dpb;

  G4double Mm1 = BreitWigner(2.275, 0.162)*GeV; //E27

  G4double thetaK, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception("E27Reaction::E27_Kpp_F_LambdaP",
		  "Production under threshold",
		  RunMustBeAborted,
		  "E27Reaction::Production under Threshold!!");
     }

    status=Scattering2Body_theta( Mi1, Mi2, Mf1, Mm1,
				  pb*LBeamDir,LPini2,
				  LPf1, LPm1,theta_CM, gen1 );
    theta_CM = theta_CM*(180./(acos(-1.)));

     if(status ==true){
       thetaK = LPf1.theta()*(180./(acos(-1.)));
       G4cout<<"thetaK= " <<thetaK <<G4endl;
       if(thetaK<20.){
       status2=Decay2Body( Mm1, Mf2, Mf3, LPm1, LPf2, LPf3, gen2 );
       if(status2 == true)
	 break;
       else
	 std::cout<<"Mm1="<<Mm1<<std::endl;
       }
     }
     pb = gConf.Get<G4double>( "BeamMom" )*GeV;
     dpb = 0.;
     if(gConf.Get<G4double>( "BeamMom" )!=0.)
       dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
     pb += dpb;
     Mm1 = BreitWigner(2.275, 0.162)*GeV; //E27
  }



  double theta_scat;
  double cos_theta_scat;
  G4ThreeVector beam_mom = pb*LBeamDir;

  cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  theta_scat = acos(cos_theta_scat)*(180./acos(-1.));


  //  std::cout<<"theta_scat="<<theta_scat<<std::endl;

  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPf1, Mf1);

  double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();

  gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* KaonPlus;
  KaonPlus = particleTable->FindParticle("kaon+");
  pGen->particleGun->SetParticleDefinition(KaonPlus);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf1);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Lambda;
  Lambda= particleTable->FindParticle("lambda");
  pGen->particleGun->SetParticleDefinition(Lambda);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf2);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Proton;
  Proton= particleTable->FindParticle("proton");
  pGen->particleGun->SetParticleDefinition(Proton);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf3);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());

}

//_____________________________________________________________________________
// reactio No #2704 pi+ d -> K+ K-pp, K-pp -> SigmaZ p reaction
void
E27Reaction::E27_Kpp_F_SigmaZP(G4Event* anEvent)
{
  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4SigmaZero::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg( gConf.Get<G4double>( "BeamX0" ),
					     gConf.Get<G4double>( "BeamY0" ),
					     gGeom.GetGlobalPosition( "Target" ).z(),
					     gConf.Get<G4double>( "BeamDX" ),
					     gConf.Get<G4double>( "BeamDY" ),
					     gSize.Get( "Target", ThreeVector::X ),
					     gSize.Get( "Target", ThreeVector::Z ));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV( gConf.Get<G4double>( "BeamU0" ),
						gConf.Get<G4double>( "BeamV0" ),
						gConf.Get<G4double>( "BeamDU" ),
						gConf.Get<G4double>( "BeamDV" ));

  G4double pb = gConf.Get<G4double>( "BeamMom" )*GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>( "BeamMom" )!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
  pb += dpb;

  G4double Mm1 = BreitWigner(2.275, 0.162)*GeV; //E27

  G4double thetaK, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception("E27Reaction::E27_Kpp_F_LambdaP",
		  "Production under threshold",
		  RunMustBeAborted,
		  "E27Reaction::Production under Threshold!!");
     }

    status=Scattering2Body_theta( Mi1, Mi2, Mf1, Mm1,
				  pb*LBeamDir,LPini2,
				  LPf1, LPm1,theta_CM, gen1 );
    theta_CM = theta_CM*(180./(acos(-1.)));

     if(status ==true){
       thetaK = LPf1.theta()*(180./(acos(-1.)));
       G4cout<<"thetaK= " <<thetaK <<G4endl;
       if(thetaK<20.){

	 status2=Decay2Body( Mm1, Mf2, Mf3, LPm1, LPf2, LPf3, gen2 );
	 if(status2 == true)
	   break;
	 else
	   std::cout<<"Mm1="<<Mm1<<std::endl;
       }
     }
     pb = gConf.Get<G4double>( "BeamMom" )*GeV;
     dpb = 0.;
     if(gConf.Get<G4double>( "BeamMom" )!=0.)
       dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
     pb += dpb;
     Mm1 = BreitWigner(2.275, 0.162)*GeV; //E27
  }

  double theta_scat;
  double cos_theta_scat;
  G4ThreeVector beam_mom = pb*LBeamDir;

  cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  theta_scat = acos(cos_theta_scat)*(180./acos(-1.));

  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPf1, Mf1);

  double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();

  gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* KaonPlus;
  KaonPlus = particleTable->FindParticle("kaon+");
  pGen->particleGun->SetParticleDefinition(KaonPlus);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf1);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* SigmaZ;
  SigmaZ= particleTable->FindParticle("sigma0");
  pGen->particleGun->SetParticleDefinition(SigmaZ);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf2);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Proton;
  Proton= particleTable->FindParticle("proton");
  pGen->particleGun->SetParticleDefinition(Proton);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf3);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),SigmaZ->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());

}

//_____________________________________________________________________________
// reactio No #2705 pi+ d -> K+ K-pp, K-pp -> Lambda piz p reaction
void
E27Reaction::E27_Kpp_F_LambdaPizP(G4Event* anEvent)
{
  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4Lambda::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();
  G4double Mf4=G4PionZero::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg( gConf.Get<G4double>( "BeamX0" ),
					     gConf.Get<G4double>( "BeamY0" ),
					     gGeom.GetGlobalPosition( "Target" ).z(),
					     gConf.Get<G4double>( "BeamDX" ),
					     gConf.Get<G4double>( "BeamDY" ),
					     gSize.Get( "Target", ThreeVector::X ),
					     gSize.Get( "Target", ThreeVector::Z ));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV( gConf.Get<G4double>( "BeamU0" ),
						gConf.Get<G4double>( "BeamV0" ),
						gConf.Get<G4double>( "BeamDU" ),
						gConf.Get<G4double>( "BeamDV" ));

  G4double pb = gConf.Get<G4double>( "BeamMom" )*GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>( "BeamMom" )!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
  pb += dpb;

  G4double Mm1 = BreitWigner(2.275, 0.162)*GeV; //E27

  G4double thetaK, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3, LPf4;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception("E27Reaction::E27_Kpp_F_LambdaP",
		  "Production under threshold",
		  RunMustBeAborted,
		  "E27Reaction::Production under Threshold!!");
     }

    status=Scattering2Body_theta( Mi1, Mi2, Mf1, Mm1,
				  pb*LBeamDir,LPini2,
				  LPf1, LPm1,theta_CM, gen1 );
    theta_CM = theta_CM*(180./(acos(-1.)));

     if(status ==true){
       thetaK = LPf1.theta()*(180./(acos(-1.)));
       G4cout<<"thetaK= " <<thetaK <<G4endl;
       if(thetaK<20.){

	 status2=Decay3BodyPhaseSpace( Mm1, Mf2, Mf3, Mf4, LPm1, LPf2, LPf3, LPf4);
	 if(status2 == true)
	   break;
	 else
	   std::cout<<"Mm1="<<Mm1<<std::endl;
       }
     }
     pb = gConf.Get<G4double>( "BeamMom" )*GeV;
     dpb = 0.;
     if(gConf.Get<G4double>( "BeamMom" )!=0.)
       dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
     pb += dpb;
     Mm1 = BreitWigner(2.275, 0.162)*GeV; //E27
  }

  double theta_scat;
  double cos_theta_scat;
  G4ThreeVector beam_mom = pb*LBeamDir;

  cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  theta_scat = acos(cos_theta_scat)*(180./acos(-1.));

  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPf1, Mf1);

  double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();

  gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* KaonPlus;
  KaonPlus = particleTable->FindParticle("kaon+");
  pGen->particleGun->SetParticleDefinition(KaonPlus);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf1);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Lambda;
  Lambda= particleTable->FindParticle("lambda");
  pGen->particleGun->SetParticleDefinition(Lambda);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf2);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Proton;
  Proton= particleTable->FindParticle("proton");
  pGen->particleGun->SetParticleDefinition(Proton);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf3);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* PiZero;
  PiZero= particleTable->FindParticle("pi0");
  pGen->particleGun->SetParticleDefinition(PiZero);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf4);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(4);

  gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(3,LPf4.x(),LPf4.y(),LPf4.z(),PiZero->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
// reactio No #2706 pi+ d -> K+ K-pp, K-pp -> SigmaZ piz p reaction
void
E27Reaction::E27_Kpp_F_SigmaZPizP(G4Event* anEvent)
{
  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4SigmaZero::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();
  G4double Mf4=G4PionZero::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg( gConf.Get<G4double>( "BeamX0" ),
					     gConf.Get<G4double>( "BeamY0" ),
					     gGeom.GetGlobalPosition( "Target" ).z(),
					     gConf.Get<G4double>( "BeamDX" ),
					     gConf.Get<G4double>( "BeamDY" ),
					     gSize.Get( "Target", ThreeVector::X ),
					     gSize.Get( "Target", ThreeVector::Z ));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV( gConf.Get<G4double>( "BeamU0" ),
						gConf.Get<G4double>( "BeamV0" ),
						gConf.Get<G4double>( "BeamDU" ),
						gConf.Get<G4double>( "BeamDV" ));

  G4double pb = gConf.Get<G4double>( "BeamMom" )*GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>( "BeamMom" )!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
  pb += dpb;

  G4double Mm1 = BreitWigner(2.275, 0.162)*GeV; //E27

  G4double thetaK, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3, LPf4;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception("E27Reaction::E27_Kpp_F_LambdaP",
		  "Production under threshold",
		  RunMustBeAborted,
		  "E27Reaction::Production under Threshold!!");
     }

    status=Scattering2Body_theta( Mi1, Mi2, Mf1, Mm1,
				  pb*LBeamDir,LPini2,
				  LPf1, LPm1,theta_CM, gen1 );
    theta_CM = theta_CM*(180./(acos(-1.)));

     if(status ==true){
       thetaK = LPf1.theta()*(180./(acos(-1.)));
       G4cout<<"thetaK= " <<thetaK <<G4endl;
       if(thetaK<20.){

	 status2=Decay3BodyPhaseSpace( Mm1, Mf2, Mf3, Mf4, LPm1, LPf2, LPf3, LPf4);
	 if(status2 == true)
	   break;
	 else
	   std::cout<<"Mm1="<<Mm1<<std::endl;
       }
     }
     pb = gConf.Get<G4double>( "BeamMom" )*GeV;
     dpb = 0.;
     if(gConf.Get<G4double>( "BeamMom" )!=0.)
       dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
     pb += dpb;
     Mm1 = BreitWigner(2.275, 0.162)*GeV; //E27
  }

  double theta_scat;
  double cos_theta_scat;
  G4ThreeVector beam_mom = pb*LBeamDir;

  cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  theta_scat = acos(cos_theta_scat)*(180./acos(-1.));

  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPf1, Mf1);

  double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();

  gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* KaonPlus;
  KaonPlus = particleTable->FindParticle("kaon+");
  pGen->particleGun->SetParticleDefinition(KaonPlus);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf1);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* SigmaZero;
  SigmaZero= particleTable->FindParticle("sigma0");
  pGen->particleGun->SetParticleDefinition(SigmaZero);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf2);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Proton;
  Proton= particleTable->FindParticle("proton");
  pGen->particleGun->SetParticleDefinition(Proton);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf3);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* PiZero;
  PiZero= particleTable->FindParticle("pi0");
  pGen->particleGun->SetParticleDefinition(PiZero);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf4);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(4);

  gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),SigmaZero->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(3,LPf4.x(),LPf4.y(),LPf4.z(),PiZero->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
// reactio No #2707 pi+ d -> K+ K-pp, K-pp -> SigmaP pim p reaction
void
E27Reaction::E27_Kpp_F_SigmaPPimP(G4Event* anEvent)
{
  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4SigmaPlus::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();
  G4double Mf4=G4PionMinus::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg( gConf.Get<G4double>( "BeamX0" ),
					     gConf.Get<G4double>( "BeamY0" ),
					     gGeom.GetGlobalPosition( "Target" ).z(),
					     gConf.Get<G4double>( "BeamDX" ),
					     gConf.Get<G4double>( "BeamDY" ),
					     gSize.Get( "Target", ThreeVector::X ),
					     gSize.Get( "Target", ThreeVector::Z ));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV( gConf.Get<G4double>( "BeamU0" ),
						gConf.Get<G4double>( "BeamV0" ),
						gConf.Get<G4double>( "BeamDU" ),
						gConf.Get<G4double>( "BeamDV" ));

  G4double pb = gConf.Get<G4double>( "BeamMom" )*GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>( "BeamMom" )!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
  pb += dpb;

  G4double Mm1 = BreitWigner(2.275, 0.162)*GeV; //E27

  G4double thetaK, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3, LPf4;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception("E27Reaction::E27_Kpp_F_LambdaP",
		  "Production under threshold",
		  RunMustBeAborted,
		  "E27Reaction::Production under Threshold!!");
     }

    status=Scattering2Body_theta( Mi1, Mi2, Mf1, Mm1,
				  pb*LBeamDir,LPini2,
				  LPf1, LPm1,theta_CM, gen1 );
    theta_CM = theta_CM*(180./(acos(-1.)));

     if(status ==true){
       thetaK = LPf1.theta()*(180./(acos(-1.)));
       G4cout<<"thetaK= " <<thetaK <<G4endl;
       if(thetaK<20.){

	 status2=Decay3BodyPhaseSpace( Mm1, Mf2, Mf3, Mf4, LPm1, LPf2, LPf3, LPf4);
	 if(status2 == true)
	   break;
	 else
	   std::cout<<"Mm1="<<Mm1<<std::endl;
       }
     }
     pb = gConf.Get<G4double>( "BeamMom" )*GeV;
     dpb = 0.;
     if(gConf.Get<G4double>( "BeamMom" )!=0.)
       dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
     pb += dpb;
     Mm1 = BreitWigner(2.275, 0.162)*GeV; //E27
  }

  double theta_scat;
  double cos_theta_scat;
  G4ThreeVector beam_mom = pb*LBeamDir;

  cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  theta_scat = acos(cos_theta_scat)*(180./acos(-1.));

  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPf1, Mf1);

  double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();

  gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* KaonPlus;
  KaonPlus = particleTable->FindParticle("kaon+");
  pGen->particleGun->SetParticleDefinition(KaonPlus);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf1);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* SigmaPlus;
  SigmaPlus= particleTable->FindParticle("sigma+");
  pGen->particleGun->SetParticleDefinition(SigmaPlus);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf2);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Proton;
  Proton= particleTable->FindParticle("proton");
  pGen->particleGun->SetParticleDefinition(Proton);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf3);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* PiMinus;
  PiMinus= particleTable->FindParticle("pi-");
  pGen->particleGun->SetParticleDefinition(PiMinus);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf4);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(4);

  gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),SigmaPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(3,LPf4.x(),LPf4.y(),LPf4.z(),PiMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
// reactio No #2708 K- 12C -> p 11KB, 11KB -> Lambda 10Be reaction
void
E27Reaction::E27_K11B_Lambda10Be(G4Event* anEvent)
{
  G4double Mi1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mi2=12.*AtomicMassUnit*GeV;//12C
  G4double Mf1=G4Proton::Definition()->GetPDGMass();
  G4double Mf2=G4Lambda::Definition()->GetPDGMass();
  G4double Mf3=10.0135338*AtomicMassUnit*GeV;//10Be

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg( gConf.Get<G4double>( "BeamX0" ),
					     gConf.Get<G4double>( "BeamY0" ),
					     gGeom.GetGlobalPosition( "Target" ).z(),
					     gConf.Get<G4double>( "BeamDX" ),
					     gConf.Get<G4double>( "BeamDY" ),
					     gSize.Get( "Target", ThreeVector::X ),
					     gSize.Get( "Target", ThreeVector::Z ));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV( gConf.Get<G4double>( "BeamU0" ),
						gConf.Get<G4double>( "BeamV0" ),
						gConf.Get<G4double>( "BeamDU" ),
						gConf.Get<G4double>( "BeamDV" ));

  G4double pb = gConf.Get<G4double>( "BeamMom" )*GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>( "BeamMom" )!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>( "BeamWidth" ))*GeV;
  pb += dpb;

  double Boron11Mass = 11.0093054 * AtomicMassUnit*GeV;
  double MK = G4KaonPlus::Definition()->GetPDGMass();
  G4double Mm1 = BreitWigner(Boron11Mass + MK - 132.5
			     , 183.); //J. Yamagata et al.,

  G4double thetap, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception("E27Reaction::E27_Kpp_F_LambdaP",
		  "Production under threshold",
		  RunMustBeAborted,
		  "E27Reaction::Production under Threshold!!");
     }

    status=Scattering2Body_theta( Mi1, Mi2, Mf1, Mm1,
				  pb*LBeamDir,LPini2,
				  LPf1, LPm1,theta_CM, gen1 );
    theta_CM = theta_CM*(180./(acos(-1.)));

     if(status ==true){
       thetap = LPf1.theta()*(180./(acos(-1.)));
       G4cout<<"theta p= " <<thetap <<G4endl;
       if(thetap<20.){
       status2=Decay2Body( Mm1, Mf2, Mf3, LPm1, LPf2, LPf3, gen2 );
       if(status2 == true)
	 break;
       else
	 std::cout<<"Mm1="<<Mm1<<std::endl;
       }
     }
     pb = gConf.Get<G4double>( "BeamMom" )*GeV;
     dpb = 0.;


     if( gConf.Get<G4double>( "BeamWidth" ) != 0. ){
       dpb = G4RandGauss::shoot( 0., gConf.Get<G4double>( "BeamWidth" ) )*GeV ;
     }
     pb += dpb;

     Mm1 = BreitWigner(Boron11Mass + MK - 132.5
		       , 183.); //J. Yamagata et al.,
  }



  double theta_scat;
  double cos_theta_scat;
  G4ThreeVector beam_mom = pb*LBeamDir;

  cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  theta_scat = acos(cos_theta_scat)*(180./acos(-1.));


  //  std::cout<<"theta_scat="<<theta_scat<<std::endl;

  G4LorentzVector Lv_beam, Lv_targ, Lv_targ_p, Lv_p;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_p.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_p.setVectM(LPf1, Mf1);

  double mm = (Lv_beam + Lv_targ + (-1.)*Lv_p).mag();
  double mm_p = (Lv_beam + Lv_targ_p + (-1.)*Lv_p).mag();

  gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Proton;
  Proton = particleTable->FindParticle("proton");
  pGen->particleGun->SetParticleDefinition(Proton);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf1);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Lambda;
  Lambda= particleTable->FindParticle("lambda");
  pGen->particleGun->SetParticleDefinition(Lambda);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf2);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);


  G4IonTable* ionTable = G4IonTable::GetIonTable();
  G4int Z = 4;
  G4int A = 10;
  G4ParticleDefinition *Be10 = ionTable->GetIon(Z, A, 0.);
  if (0 == Be10) {
    G4cerr << " ERROR: Can't create nucleus with Z=" << Z << " A=" << A
	   << G4endl;
    ::exit(1);
  }
  // G4ParticleDefinition* ;
  // Proton= particleTable->FindParticle("proton");
  pGen->particleGun->SetParticleDefinition(Be10);
  pGen->particleGun->SetParticlePosition(LPos);
  pGen->particleGun->SetParticleMomentum(LPf3);
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryInfo(mm, mm_p, thetap, theta_scat, theta_CM);
  gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Be10->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());

}

//_____________________________________________________________________________
// reactio No #2709 K+ gun for test
void
E27Reaction::E27_Kptest2(G4Event* anEvent)
{
  int nev = anEvent->GetEventID();
  G4double pbeam;

  for(int i=0; i<14; ++i){
    int fac = nev%14;
    pbeam = (0.1+0.1*(double)fac);
  }

  // G4double Maxang = 30.;
  // G4double cost = cos(Maxang*G4UniformRand()*acos(-1.)/180.);
  G4double cost = cos(0.*acos(-1.)/180.);
  G4double sint =sqrt( 1.0 - cost *cost);
  //G4double phi = 360.*(0.5-G4UniformRand())*degree;
  G4double phi = 0.;


  G4double Energy_Kp,  mom_Kp_x, mom_Kp_y, mom_Kp_z;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* KaonPlus;
  KaonPlus = particleTable->FindParticle("kaon+");
  //KaonPlus = particleTable->FindParticle("proton");

  mom_Kp_x=pbeam * sint*cos(phi);
  mom_Kp_y=pbeam * sint*sin(phi);
  mom_Kp_z=pbeam * cost;
  Energy_Kp=sqrt(KaonPlus->GetPDGMass()/GeV*KaonPlus->GetPDGMass()/GeV+pbeam*pbeam);

  double theta = acos(cost)*180./acos(-1.);
  std::cout<<"theta ="<<theta<<", phi(rad)="<<phi<<std::endl;
  std::cout<<"Mom =("<<mom_Kp_x<<", "<<mom_Kp_y<<", "<<mom_Kp_z<<std::endl;


  gAnaMan.SetPrimaryBeam(mom_Kp_x,mom_Kp_y,mom_Kp_z);

  // G4double Ebeam = sqrt(pbeam*pbeam+KaonPlus->GetPDGMass()/GeV*KaonPlus->GetPDGMass()/GeV);


  //  G4double vtx = CLHEP::RandFlat::shoot(-15.,15.)*mm;
  //  G4double vty = CLHEP::RandFlat::shoot(-5.,5.)*mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(env_target_pos_z-env_target_width/2,env_target_pos_z+env_target_width/2)*mm;

  // G4double vtx = CLHEP::RandGauss::shoot(0,10.)*mm;
  // G4double vty = CLHEP::RandFlat::shoot(0.,3.2)*mm;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(gGeom.GetGlobalPosition( "Target" ).z()-gSize.Get( "Target", ThreeVector::Z )/2,gGeom.GetGlobalPosition( "Target" ).z()+gSize.Get( "Target", ThreeVector::Z )/2)*mm-250.*mm;
  G4double vtz= gGeom.GetGlobalPosition( "Target" ).z();
  std::cout<<"pbeam = "<<pbeam<<std::endl;
  // getchar();
  //scat K+
  pGen->particleGun->SetParticleDefinition(KaonPlus);
  pGen->particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_Kp_x,mom_Kp_y,mom_Kp_z));
  pGen->particleGun->SetParticleEnergy((Energy_Kp - KaonPlus->GetPDGMass()/GeV)*GeV);
  pGen->particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  pGen->particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_Kp_x,mom_Kp_y,mom_Kp_z,KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
}

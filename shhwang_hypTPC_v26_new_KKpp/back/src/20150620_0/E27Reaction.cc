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
#include "G4ParticleTypes.hh"
#include "G4IonConstructor.hh"
#include "GeneratorHelper.hh"
#include "AngDisGenerator.hh"

const int MaxTry=1000;

//reactio No #2701 pi+ beam through
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
  G4double vtz= pGen->Get_env_target_pos_z();
  //G4double vtz= CLHEP::RandFlat::shoot(pGen->Get_env_target_pos_z()-pGen->Get_env_target_width()/2,pGen->Get_env_target_pos_z()+pGen->Get_env_target_width()/2)*mm-250.*mm;
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

//reactio No #2702 K+ gun for test
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
  //  G4double vtz= CLHEP::RandFlat::shoot(pGen->Get_env_target_pos_z()-pGen->Get_env_target_width()/2,pGen->Get_env_target_pos_z()+pGen->Get_env_target_width()/2)*mm-250.*mm;
  G4double vtz= pGen->Get_env_target_pos_z();
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


//reactio No #2703 pi+ d -> K+ K-pp, K-pp -> Lambda p reaction
void E27Reaction::E27_Kpp_F_LambdaP(G4Event* anEvent){

  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4Lambda::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  G4double Mn=G4Neutron::Definition()->GetPDGMass();
  G4double Mpiz=G4PionZero::Definition()->GetPDGMass();
  

  G4ThreeVector LPos = GaussPosition_LqTarg( pGen->Get_env_Beam_x0(),
					     pGen->Get_env_Beam_y0(), 
					     pGen->Get_env_target_pos_z(),
					     pGen->Get_env_Beam_dx(),
					     pGen->Get_env_Beam_dy(), 
					     pGen->Get_env_target_size_x(),
					     pGen->Get_env_target_width());
  //Note!! env_target_width = Target_Size_z (height of target)
					     
  G4ThreeVector LBeamDir =  GaussDirectionInUV( pGen->Get_env_Beam_u0(), 
						pGen->Get_env_Beam_v0(), 
						pGen->Get_env_Beam_du(), 
						pGen->Get_env_Beam_dv()); 
  
  G4double pb = pGen->Get_env_Beam_mom()*GeV;
  G4double dpb = 0.;
  if(pGen->Get_env_Beam_mom()!=0.)
    dpb = G4RandGauss::shoot(0.,pGen->Get_env_Beam_width())*GeV;
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
     pb = pGen->Get_env_Beam_mom()*GeV;
     dpb = 0.;
     if(pGen->Get_env_Beam_mom()!=0.)
       dpb = G4RandGauss::shoot(0.,pGen->Get_env_Beam_width())*GeV;
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

  pGen->anaManager->SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
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

  pGen->anaManager->SetNumberOfPrimaryParticle(3);

  pGen->anaManager->SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  pGen->anaManager->SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),Lambda->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());

}


//reactio No #2704 pi+ d -> K+ K-pp, K-pp -> SigmaZ p reaction
void E27Reaction::E27_Kpp_F_SigmaZP(G4Event* anEvent){

  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4SigmaZero::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  G4double Mn=G4Neutron::Definition()->GetPDGMass();
  G4double Mpiz=G4PionZero::Definition()->GetPDGMass();
  

  G4ThreeVector LPos = GaussPosition_LqTarg( pGen->Get_env_Beam_x0(),
					     pGen->Get_env_Beam_y0(), 
					     pGen->Get_env_target_pos_z(),
					     pGen->Get_env_Beam_dx(),
					     pGen->Get_env_Beam_dy(), 
					     pGen->Get_env_target_size_x(),
					     pGen->Get_env_target_width());
  //Note!! env_target_width = Target_Size_z (height of target)
					     
  G4ThreeVector LBeamDir =  GaussDirectionInUV( pGen->Get_env_Beam_u0(), 
						pGen->Get_env_Beam_v0(), 
						pGen->Get_env_Beam_du(), 
						pGen->Get_env_Beam_dv()); 
  
  G4double pb = pGen->Get_env_Beam_mom()*GeV;
  G4double dpb = 0.;
  if(pGen->Get_env_Beam_mom()!=0.)
    dpb = G4RandGauss::shoot(0.,pGen->Get_env_Beam_width())*GeV;
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
     pb = pGen->Get_env_Beam_mom()*GeV;
     dpb = 0.;
     if(pGen->Get_env_Beam_mom()!=0.)
       dpb = G4RandGauss::shoot(0.,pGen->Get_env_Beam_width())*GeV;
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

  pGen->anaManager->SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
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

  pGen->anaManager->SetNumberOfPrimaryParticle(3);

  pGen->anaManager->SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  pGen->anaManager->SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),SigmaZ->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());

}


//reactio No #2705 pi+ d -> K+ K-pp, K-pp -> Lambda piz p reaction
void E27Reaction::E27_Kpp_F_LambdaPizP(G4Event* anEvent){

  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4Lambda::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();
  G4double Mf4=G4PionZero::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  G4double Mn=G4Neutron::Definition()->GetPDGMass();
  G4double Mpiz=G4PionZero::Definition()->GetPDGMass();
  

  G4ThreeVector LPos = GaussPosition_LqTarg( pGen->Get_env_Beam_x0(),
					     pGen->Get_env_Beam_y0(), 
					     pGen->Get_env_target_pos_z(),
					     pGen->Get_env_Beam_dx(),
					     pGen->Get_env_Beam_dy(), 
					     pGen->Get_env_target_size_x(),
					     pGen->Get_env_target_width());
  //Note!! env_target_width = Target_Size_z (height of target)
					     
  G4ThreeVector LBeamDir =  GaussDirectionInUV( pGen->Get_env_Beam_u0(), 
						pGen->Get_env_Beam_v0(), 
						pGen->Get_env_Beam_du(), 
						pGen->Get_env_Beam_dv()); 
  
  G4double pb = pGen->Get_env_Beam_mom()*GeV;
  G4double dpb = 0.;
  if(pGen->Get_env_Beam_mom()!=0.)
    dpb = G4RandGauss::shoot(0.,pGen->Get_env_Beam_width())*GeV;
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
     pb = pGen->Get_env_Beam_mom()*GeV;
     dpb = 0.;
     if(pGen->Get_env_Beam_mom()!=0.)
       dpb = G4RandGauss::shoot(0.,pGen->Get_env_Beam_width())*GeV;
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

  pGen->anaManager->SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
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

  pGen->anaManager->SetNumberOfPrimaryParticle(4);

  pGen->anaManager->SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  pGen->anaManager->SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),Lambda->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(3,LPf4.x(),LPf4.y(),LPf4.z(),PiZero->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}


//reactio No #2706 pi+ d -> K+ K-pp, K-pp -> SigmaZ piz p reaction
void E27Reaction::E27_Kpp_F_SigmaZPizP(G4Event* anEvent){

  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4SigmaZero::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();
  G4double Mf4=G4PionZero::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  G4double Mn=G4Neutron::Definition()->GetPDGMass();
  G4double Mpiz=G4PionZero::Definition()->GetPDGMass();
  

  G4ThreeVector LPos = GaussPosition_LqTarg( pGen->Get_env_Beam_x0(),
					     pGen->Get_env_Beam_y0(), 
					     pGen->Get_env_target_pos_z(),
					     pGen->Get_env_Beam_dx(),
					     pGen->Get_env_Beam_dy(), 
					     pGen->Get_env_target_size_x(),
					     pGen->Get_env_target_width());
  //Note!! env_target_width = Target_Size_z (height of target)
					     
  G4ThreeVector LBeamDir =  GaussDirectionInUV( pGen->Get_env_Beam_u0(), 
						pGen->Get_env_Beam_v0(), 
						pGen->Get_env_Beam_du(), 
						pGen->Get_env_Beam_dv()); 
  
  G4double pb = pGen->Get_env_Beam_mom()*GeV;
  G4double dpb = 0.;
  if(pGen->Get_env_Beam_mom()!=0.)
    dpb = G4RandGauss::shoot(0.,pGen->Get_env_Beam_width())*GeV;
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
     pb = pGen->Get_env_Beam_mom()*GeV;
     dpb = 0.;
     if(pGen->Get_env_Beam_mom()!=0.)
       dpb = G4RandGauss::shoot(0.,pGen->Get_env_Beam_width())*GeV;
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

  pGen->anaManager->SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
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

  pGen->anaManager->SetNumberOfPrimaryParticle(4);

  pGen->anaManager->SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  pGen->anaManager->SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),SigmaZero->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(3,LPf4.x(),LPf4.y(),LPf4.z(),PiZero->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}


//reactio No #2707 pi+ d -> K+ K-pp, K-pp -> SigmaP pim p reaction
void E27Reaction::E27_Kpp_F_SigmaPPimP(G4Event* anEvent){

  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4SigmaPlus::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();
  G4double Mf4=G4PionMinus::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  G4double Mn=G4Neutron::Definition()->GetPDGMass();
  G4double Mpiz=G4PionZero::Definition()->GetPDGMass();
  

  G4ThreeVector LPos = GaussPosition_LqTarg( pGen->Get_env_Beam_x0(),
					     pGen->Get_env_Beam_y0(), 
					     pGen->Get_env_target_pos_z(),
					     pGen->Get_env_Beam_dx(),
					     pGen->Get_env_Beam_dy(), 
					     pGen->Get_env_target_size_x(),
					     pGen->Get_env_target_width());
  //Note!! env_target_width = Target_Size_z (height of target)
					     
  G4ThreeVector LBeamDir =  GaussDirectionInUV( pGen->Get_env_Beam_u0(), 
						pGen->Get_env_Beam_v0(), 
						pGen->Get_env_Beam_du(), 
						pGen->Get_env_Beam_dv()); 
  
  G4double pb = pGen->Get_env_Beam_mom()*GeV;
  G4double dpb = 0.;
  if(pGen->Get_env_Beam_mom()!=0.)
    dpb = G4RandGauss::shoot(0.,pGen->Get_env_Beam_width())*GeV;
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
     pb = pGen->Get_env_Beam_mom()*GeV;
     dpb = 0.;
     if(pGen->Get_env_Beam_mom()!=0.)
       dpb = G4RandGauss::shoot(0.,pGen->Get_env_Beam_width())*GeV;
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

  pGen->anaManager->SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
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

  pGen->anaManager->SetNumberOfPrimaryParticle(4);

  pGen->anaManager->SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  pGen->anaManager->SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),SigmaPlus->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryParticle(3,LPf4.x(),LPf4.y(),LPf4.z(),PiMinus->GetPDGMass()/GeV);
  pGen->anaManager->SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  pGen->anaManager->SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}


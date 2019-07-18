// -*- C++ -*-

#include "TPCPrimaryGeneratorAction.hh"

#include <G4Event.hh>
#include <G4IonTable.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4UImanager.hh>
#include <G4IonConstructor.hh>
#include <Randomize.hh>

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TMath.h>

#include "BeamMan.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "Kinema3Resonance.hh"
#include "E27Reaction.hh"
#include "KKppReaction.hh"
#include "KinemaHResonance.hh"
#include "Kinema3Body.hh"
#include "Kinema4Body.hh"
#include "KinemaHybrid.hh"
#include "KinemaHweak.hh"
#include "KinemaFermi.hh"
#include "KinemaKstar.hh"
#include "TPCAnaManager.hh"
#include "ThreeVector.hh"

namespace
{
  using CLHEP::eplus;
  using CLHEP::GeV;
  using CLHEP::keV;
  using CLHEP::mm;
  auto& gAnaMan = TPCAnaManager::GetInstance();
  const auto& gBeam = BeamMan::GetInstance();
  const auto& gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
}

//_____________________________________________________________________________
TPCPrimaryGeneratorAction::TPCPrimaryGeneratorAction( void )
  : G4VUserPrimaryGeneratorAction(),
    particleGun( new G4ParticleGun ),
    m_target_pos( gGeom.GetGlobalPosition( "SHSTarget" )*mm ),
    m_target_size( gSize.GetSize( "Target" )*mm ),
    m_beam( new BeamInfo ),
    m_beam_p0( gConf.Get<G4double>( "BeamMom" ) /* 1.80*GeV */ )
{
}

//_____________________________________________________________________________
TPCPrimaryGeneratorAction::~TPCPrimaryGeneratorAction( void )
{
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent )
{
  E27Reaction E27Generator(this);
  KKppReaction KKppGenerator(this);

  G4int generator = gConf.Get<G4int>( "Generator" );

  static G4bool first = true;
  if( first ){
    G4cout << "   Generator : " << generator << G4endl;
    first = false;
  }

  TTree *t1 = nullptr;
  if( generator == 3101 ||
      generator == 3103 ||
      generator == 3104 ){
    TFile *fin = new TFile( gConf.Get<G4String>("Input_JAM_File_Name"), "READ");
    t1 = (TTree*)fin->Get("tree");
  }

  gAnaMan.SetGeneratorID( generator );

  *m_beam = gBeam.Get();
#ifdef DEBUG
   m_beam->Print();
#endif

  switch( generator  ){
  case 0:
    // no generation
    break;
  case 1:
    Generate_hanul(anEvent); // shhwang
    break;
  case 2:
    Generate_hdibaryon1(anEvent); // h-dibaryon --> LL
    break;
  case 3:
    Generate_hdibaryon2(anEvent); // h-dibaryon --> LL K+ 15deg
    break;
  case 4:
    Generate_test(anEvent); // test sako-san's code
    break;
  case 5:
    Generate_test2(anEvent); // test sako-san's code
    break;
  case 6:
    Generate_hdibaryon_PHSG(anEvent); // h weak. H->Lppi-
    break;
  case 7:
    Generate_hdibaryon_PHSG_S(anEvent); // h weak. H->S-p
    break;
  case 8:
    Generate_hdibaryon_non_reso(anEvent); // test non resonance LL??
    break;
  case 9:
    Generate_hdibaryon_PHSG_LL(anEvent); // H gen by using phsg
    break;
  case 10:
    // Generate_Kp_Kn(anEvent); // beam study
    GenerateBeam( anEvent ); // beam study
    break;
  case 11:
    Generate_pip_KsL(anEvent); // study pi-p --> KsL
    break;
  case 12:
    Generate_pip_KsS(anEvent); // study pi-p --> KsS
    break;
  case 13:
    Generate_pip_KstarL(anEvent); // study on pi-p --> KsS by using LL gen
    break;
  case 14:
    Generate_pip_KstarS(anEvent); // study on pi-p --> KsS by using LL gen
    break;
  case 21:
    Generate_hybrid3body_mode1(anEvent);
    break;
  case 22:
    Generate_hybrid3body_mode2(anEvent);
    break;
  case 23:
    Generate_hybrid3body_mode3(anEvent);
    break;
  case 24:
    Generate_hybrid3body_mode4(anEvent);
    break;
  case 25:
    Generate_E45_elastic_pip(anEvent);
    break;
  case 26:
    Generate_E45_elastic_pin(anEvent);
    break;
  case 30:
    Generate_PhaseSpace(anEvent); // phase space generator 12C + Kn --> 10Be L L Kp
    // 12C_amu = 12u --> 12*931.494061 MeV
    // 10C_amu = 10.012938u --> 10.012938*931.494061 MeV
   break;
  case 99:
    Generate_dedx_single(anEvent); // study on dedx. only single particle will be generated
    break;
  case 98:
    Generate_all(anEvent); // study on dedx. All particles will be generated
    break;
  case 70:
    Generate_E07_study(anEvent); // study on E07, generate beam and K+ from INC data
    break;
  case 71:
    Generate_E07_study_all(anEvent); // study on E07, generate K+ pi+ from INC data
    break;
  case 72:
    Generate_E07_study_knp(anEvent); // study on E07, generate K+ pi+ from INC data
    break;
  case 73:
    Generate_E07_study_kp(anEvent); // study on E07, generate K+ pi+ from INC data
    break;
  case 74:
    Generate_E07_study_knp_beam(anEvent); // study on E07, generate K+ pi+ from INC data
    break;
  case 75:
    Generate_E07_study_kp_beam(anEvent); // study on E07, generate K+ pi+ from INC data, beam size 1x3 cm^2
    break;
  case 76:
    Generate_E07_study_kpxi_beam(anEvent); // study on E07, generate K+ from isotropic, beam size 1x3 cm^2
  case 77:
    Generate_E07_study_kpxi_beam_only_kp(anEvent); // study on E07, generate K+ from isotropic, beam size 1x3 cm^2
    break;
  case 78:
    Generate_E07_study_pro_08_20(anEvent); // generate proton
    break;
  case 79:
    Generate_E07_study_kp_04_15(anEvent); // generate proton
    break;
  case 80:
    Generate_E07_study_kpxi1530(anEvent); // study on E07, Xi(1530)- generator, beam size 1x3 cm^2
    break;
  case 81:
    Generate_E07_study_Takahashi(anEvent); // study on E07, Xi-kp, by using Takahashi-san's code
    break;
    ///Study of lambda 1405
  case 60:
    Generate_Lambda1405_rad1(anEvent); // normal decay
    break;
  case 61:
    Generate_Lambda1405_rad2(anEvent); // radioactive decay
    break;
  case 62:
    Generate_Sigma1385(anEvent); // Sigma 1385 normal dcay
    break;
  case 63:
    Generate_Sigma1385_rad(anEvent); // Sigma 1385 radioactive decay
    break;
  case 64:
    // S+pi-,S0pi0, S-pi+ --> Lambda from threshold to 1.5 GeV
    Generate_Lambda1405_reso(anEvent); // pi+pi- --> K0
    break;
  case 2701:
    E27Generator.E27_beamthrough(anEvent);
    break;
  case 2702:
    E27Generator.E27_Kptest(anEvent);
    break;
  case 2703:
    E27Generator.E27_Kpp_F_LambdaP(anEvent);
    break;
  case 2704:
    E27Generator.E27_Kpp_F_SigmaZP(anEvent);
    break;
  case 2705:
    E27Generator.E27_Kpp_F_LambdaPizP(anEvent);
    break;
  case 2706:
    E27Generator.E27_Kpp_F_SigmaZPizP(anEvent);
    break;
  case 2707:
    E27Generator.E27_Kpp_F_SigmaPPimP(anEvent);
    break;
  case 2708:
    E27Generator.E27_K11B_Lambda10Be(anEvent);
    break;
  case 2709:
    E27Generator.E27_Kptest2(anEvent);
    break;
  case 3001:
    KKppGenerator.KKpp_LL1(anEvent);
    break;
  case 3002:
    KKppGenerator.KKpp_LL2(anEvent);
    break;
  case 3003:
    KKppGenerator.KKpp_LSmPip(anEvent);
    break;
  case 3004:
    KKppGenerator.KKpp_LSpPim(anEvent);
    break;
  case 3101:
    KKppGenerator.JAMInput(anEvent,t1);
    break;
  case 3102:
    KKppGenerator.KKpp_BeamThrough1(anEvent);
    break;
  case 3103:
    KKppGenerator.JAMInput_K0(anEvent,t1);
    break;
  case 3104:
    KKppGenerator.JAMInput_K0bar(anEvent,t1);
    break;
  default:
    G4cout << "Generator number error :" << generator << G4endl;
    break;
  }

  // G4UImanager* UImanager= G4UImanager::GetUIpointer();
  // UImanager-> ApplyCommand("/run/beamOn 3");
  // G4cout<<"hoge"<<G4endl;
  // getchar();


  //  Generate_hanulcut(anEvent);            // shhwang
  //  Generate_hybridPHSG(anEvent);            // isotropic h-dibaryon, K+ < 15*deg --> LL
  //    Generate_hybrid(anEvent);            // isotropic h-dibaryon, K+ < 15*deg --> LL
  //  Generate_hybrid3body(anEvent);            // isotropic h-dibaryon, K+ < 15*deg --> LL
  //  Generate_hybrid3body_mode1(anEvent);            // isotropic h-dibaryon, K+ < 15*deg --> LL
  //  Generate_hybrid3body_mode2(anEvent);            // isotropic h-dibaryon, K+ < 15*deg --> LL
  //  Generate_hybrid3body_mode3(anEvent);            // isotropic h-dibaryon, K+ < 15*deg --> LL
  //  Generate_hybrid3body_mode4(anEvent);            // isotropic h-dibaryon, K+ < 15*deg --> LL
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hdibaryon2( G4Event* anEvent )
{
  G4double mass_hdibaryon, width_hdibaryon;
  //  G4double Energy_h, momentum_h, mom_h_x, mom_h_y, mom_h_z;

  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  // G4double rmk=0.493677;
  auto particleTable = G4ParticleTable::GetParticleTable();
  auto proton = particleTable->FindParticle("proton");
  auto Lambda1 = particleTable->FindParticle("lambda");
  auto Lambda2 = particleTable->FindParticle("lambda");
  auto kaonPlus = particleTable->FindParticle("kaon+");
  auto kaonMinus = particleTable->FindParticle("kaon-");

  //  Ebeam = 1.9;       // Incident gamma energy (GeV)
  pg_x = 0.0;
  pg_y = 0.0;
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548)*GeV;
  //  pg_z = m_beam_p0;
  pg_z = G4RandGauss::shoot( m_beam_p0, 0.01294*m_beam_p0 );
  //  G4cout<<"Ebeam:"<<Ebeam<<G4endl;
  // G4double pbeam=sqrt(pow(pg_x,2)+pow(pg_y,2)+pow(pg_z,2));
  // G4double Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);       // Incident gamma energy (GeV)

  //  G4cout<<"env_beam_mom:"<<m_beam_p0<<G4endl;
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);
  mass_hdibaryon = gConf.Get<G4double>( "HdibaryonMass" );
  width_hdibaryon = gConf.Get<G4double>( "HdibaryonWidth" );

 up:
  Kinema3Resonance Hdibaryon( kaonMinus->GetPDGMass()/GeV,
			      ((proton->GetPDGMass()/GeV)*2),
			      Lambda1->GetPDGMass()/GeV,
			      Lambda2->GetPDGMass()/GeV,
			      kaonPlus->GetPDGMass()/GeV,
			      mass_hdibaryon, width_hdibaryon, pg_z, 0.0);

  Energy_kp = Hdibaryon.GetEnergy(5);
  // G4double momentum_kp = Hdibaryon.GetMomentum(5);
  Hdibaryon.GetMomentum(5,mom);

  //  if(atan((mom[1])/mom[2])*180/3.141592654 > 15. ) goto up;
  if( fabs(atan2(mom[1],mom[2])*180/3.141592654) > 15. || fabs(atan2(mom[0],mom[2])*180/3.141592654) > 20. ) goto up;
  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;

  /* L1 */
  Energy_L1 = Hdibaryon.GetEnergy(3);
  // G4double momentum_L1 = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];

  // G4double ThetaL1 = Hdibaryon.GetTheta(3);
  // G4double PhiL1 = Hdibaryon.GetPhi(3);

  /* L2 */
  Energy_L2 = Hdibaryon.GetEnergy(4);
  // G4double momentum_L2 = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double ThetaL2 = Hdibaryon.GetTheta(4);
  // G4double PhiL2 = Hdibaryon.GetPhi(4);

  /* Kp */
  Energy_kp = Hdibaryon.GetEnergy(5);
  // G4double momentum_kp = Hdibaryon.GetMomentum(5);
  Hdibaryon.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  //  Thetakp = Hdibaryon.GetThetaCM(5);
  // G4double Thetakp = Hdibaryon.GetThetaCM(1);
  //  Phikp = Hdibaryon.GetPhi(5);
  //  G4cout<<Phikp<<G4endl;
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4cout<<"phi test: "<<Phikp<<":"<<shphi<<G4endl;
  //  G4cout<<"test: "<<Hdibaryon.GetPhiCM(1)<<":"<<Hdibaryon.GetPhi(5)<<G4endl;
  // Vertex (Reaction) point
  //  G4cout<<Energy_kp<<G4endl;
  //  //  G4double beta= sqrt(1.8*1.8-0.493677*0.493677)/(2*0.93827203+Energy_kp);
  //  G4cout<<(proton->GetPDGMass()/GeV)<<G4endl;
  // G4double beta= pbeam/(2.*0.938272013+Ebeam);
  //  G4cout<<"Pri:"<<beta<<G4endl;
  //  G4double beta= sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2))/(2*0.93827203+1.8);
  // G4double test;
  // G4double momk[4]={0.};
  // G4double momkpp;
  // G4double momcmk[4]={0.};
  // G4double momcmkpp;

  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  // momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // test=lorentz(momk,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  cmk->Fill(Thetakp);
  //  G4cout<<"Theta LAB :"<<":"<<acos(momk[2]/momkpp)*180/3.141592654<<G4endl;
  //  G4cout<<"Theta CM gen:"<<Thetakp<<G4endl;
  //  G4cout<<"Theta CM cal:"<<acos(momcmk[2]/momcmkpp)*180/3.141592654<<G4endl;
  // G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

//  G4cout<<"sh phi:"<<shphi<<G4endl;
//  G4cout<<"cm phi:"<<cmphik<<G4endl;

  //coscmk->Fill((momcmk[2]/momcmkpp));

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);


  ////check angle in the CM frame
  G4double hlab[4];
  hlab[0]=mom_L1_x+mom_L2_x;
  hlab[1]=mom_L1_y+mom_L2_y;
  hlab[2]=mom_L1_z+mom_L2_z;
  //  hlab[3]=sqrt(pow(hlab[0],2)+pow(hlab[1],2)+pow(hlab[2],2)+pow(Lambda1->GetPDGMass()/GeV,2)*2);
  //  mass_hdibaryon = 2.250;    // Mass for H
  hlab[3]=sqrt(hlab[0]*hlab[0]+hlab[1]*hlab[1]+hlab[2]*hlab[2]+mass_hdibaryon*mass_hdibaryon);
  // G4double hcm[4];

  // test=lorentz(hlab,beta,hcm);
  // G4double hcmpp=sqrt(pow(hcm[0],2)+pow(hcm[1],2)+pow(hcm[2],2));
  // cmh->Fill(acos(hcm[2]/hcmpp)/3.141592654*180);

  // coscmh->Fill(hcm[2]/hcmpp);

  //  G4double checkphi;
  // G4double hphi=atan2(hcm[1],hcm[0])/3.141592654*180;
  //  G4double hphilab=atan2(hlab[1],hlab[0])/3.141592654*180;

  // phih->Fill(hphi);
  // thetadiff->Fill(acos(hcm[2]/hcmpp)/3.141592654*180+acos(momcmk[2]/momcmkpp)/3.141592654*180);

  // phidiff->Fill(hphi-cmphik);
  //  G4cout<<"theta diff H:"<<acos(hcm[2]/hcmpp)/3.141592654*180<<G4endl;
  //  G4cout<<"theta diff K+:"<<acos(momcmk[2]/momcmkpp)/3.141592654*180<<G4endl;
  // G4double momsum[3]={0};
  // momsum[0]=mom_kp_x+mom_L1_x+mom_L2_x;
  // momsum[1]=mom_kp_y+mom_L1_y+mom_L2_y;
  // momsum[2]=mom_kp_z+mom_L1_z+mom_L2_z;
  //  G4cout<<"mom sum:"<<momsum[0]<<":"<<momsum[1]<<":"<<momsum[2]<<G4endl;
  // G4double momH[3]={0};
  // momH[0]=mom_L1_x+mom_L2_x;
  // momH[1]=mom_L1_y+mom_L2_y;
  // momH[2]=mom_L1_z+mom_L2_z;
  //  G4cout<<"mom sum X H:kp:"<<momH[0]<<":"<<mom_kp_x<<G4endl;
  //  G4cout<<"mom sum Y H:kp:"<<momH[1]<<":"<<mom_kp_y<<G4endl;
  //  G4cout<<"mom sum Z H:kp:"<<momH[2]<<":"<<mom_kp_z<<G4endl;

  // G4double misskp = miss1(pbeam,rmk,momk);//pbeam, mass, momk
  //  G4cout<<rmk<<G4endl;
  //  G4cout<<"MM(K+):"<<misskp<<G4endl;
  //  G4cout<<"pp:"<<pbeam<<G4endl;
  //  G4cout<<"ppk:"<<sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z)<<G4endl;
  //missk->Fill(misskp);

  // G4double e1,e2,etot,invm2,ptot,invm;
  // invm=-1.;

  // ptot=pow(mom_L1_x+mom_L2_x,2)+pow(mom_L1_y+mom_L2_y,2)+pow(mom_L1_z+mom_L2_z,2);
  // e1=Energy_L1;
  // e2=Energy_L2;
  // etot=pow(e1+e2,2);
  // invm2=etot-ptot;
  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);

  //  G4cout<<"IM(LL)"<<invm<<G4endl;
  //  G4cout<<"-----------end-------------"<<G4endl;
  //  G4double vtr = 20.0*mm*((double) G4RandFlat::shoot());
  //  G4double vtphi = 2.0*pi*((double) G4RandFlat::shoot());
  //  G4double vtx = vtr*cos(vtphi);
  //  G4double vty = vtr*sin(vtphi);

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,
				       env_target_pos_z+m_target_size.z()/2)*mm;

  ///////kp PG
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////L1 PG
  particleGun->SetParticleDefinition(Lambda1);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  particleGun->SetParticleEnergy((Energy_L1 - Lambda1->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////L2 PG
  particleGun->SetParticleDefinition(Lambda2);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  particleGun->SetParticleEnergy((Energy_L2 - Lambda2->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,Lambda1->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,Lambda2->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hanul( G4Event* anEvent )
{

  double pbm[4]={-9999.9999}, pka[4]={-9999.9999}, vtx[3]={-9999.9999};
  double pL1[4]={-9999.9999}, pL2[4]={-9999.9999};
				//,xL1[3]={-9999.9999}, xL2[3]={-9999.9999};
  //  double pp1[4]={-9999.9999}, ppi1[4]={-9999.9999}, pp2[4]={-9999.9999}, ppi2[4]={-9999.9999};

  char fnam[30];

  auto particleTable = G4ParticleTable::GetParticleTable();
  // auto kaonMinus = particleTable->FindParticle("kaon-");
  auto kaonPlus = particleTable->FindParticle("kaon+");
  // auto proton = particleTable->FindParticle("proton");
  // auto piMinus = particleTable->FindParticle("pi-");
  auto Lambda = particleTable->FindParticle("lambda");
  auto lambda1 = particleTable->FindParticle("lambda");
  auto lambda2 = particleTable->FindParticle("lambda");

  FILE *fp;
  sprintf(fnam,"inc64cu.dat");

  if ((fp = fopen(fnam,"r")) == NULL){
    fprintf(stderr, "inc64cu.dat:: Cannot open file: %s\n", fnam);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check=0.;
  // 22478 lines --> remove NaN
  int ran = G4RandFlat::shoot(1.,22478.);

  while(1){
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
	   ,&data[0],&data[1],&data[2],&data[3],&data[4]
	   ,&data[5],&data[6],&data[7],&data[8],&data[9]
	   ,&data[10],&data[11],&data[12],&data[13],&data[14]
	   ,&data[15],&data[16],&data[17],&data[18],&data[19]
	   ,&data[20],&data[21],&data[22],&data[23],&data[24]
	   ,&data[25],&data[26],&data[27],&data[28],&data[29]
	   ,&data[30],&data[31],&data[32],&data[33],&data[34]
	   ,&data[35],&data[36],&data[37],&data[38],&data[39]
	   ,&data[40]);
    check=check+1.;
    if(check==ran) break;
    //    if(res==EOF) break;
  }

  for(int ii=0;ii<4;ii++){
    pbm[ii]=data[ii]/1000.;
    pka[ii]=data[ii+4]/1000.;
    pL1[ii]=data[ii+11]/1000.;
    pL2[ii]=data[ii+15]/1000.;
  }

  for(int ii=0;ii<3;ii++){
    vtx[ii]=data[ii+8];
  }

  fclose(fp);

  G4double Energy_L1;
  G4double Energy_L2;
  G4double Energy_ka;

  vtx[2] = G4RandFlat::shoot( -150. - m_target_size.z()/2,
				   -150. + m_target_size.z()/2 );

  // G4double Energy_beam = pbm[3];
  G4ThreeVector momentumBeam(-pbm[0],-pbm[1],-pbm[2]);
  //  G4ThreeVector vertexPos(vtx[0]*mm,vtx[1]*mm, vtx[2]*mm);
  G4ThreeVector vertexPos(0.*mm,0.*mm, vtx[2]*mm);

  Energy_ka = pka[3]; //total energy sqrt(pp^2+rmk^2)??
  Energy_L1 = pL1[3]; //lambda momenta ?? --> 4 momomenta sqrt(p^2+m^2)
  Energy_L2 = pL2[3]; //lambda momenta ?? --> 4 momomenta sqrt(p^2+m^2)

  // G4double beta= pbm[2]/(2*0.93827203+pbm[3]);
  // G4double test;
  // G4double momk[4]={0.};
  // G4double momkpp;
  // G4double momcmk[4]={0.};
  // momk[0]=pka[0];
  // momk[1]=pka[1];
  // momk[2]=pka[2];
  // momk[3]=sqrt(pka[0]*pka[0]+pka[1]*pka[1]+pka[2]*pka[2]+0.493677*0.493677);

  // momkpp=sqrt(pow(momk[0],2)+pow(momk[1],2)+pow(momk[2],2));

  // test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));

  // cmk->Fill(acos(momcmk[2]/momcmkpp)*180/3.141592654);
  // coscmk->Fill(momcmk[2]/momcmkpp);

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

  // G4double shphi=(atan2(momk[1],momk[0]))/3.141592654*180.;
  //  phik->Fill(shphi);


  // ---- K+ -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(momentumKaonPlus);
  particleGun->SetParticleEnergy((Energy_ka - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(vertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);

  // ---- Lambda 1 -------------
  G4ThreeVector momentumLambda1(pL1[0], pL1[1], pL1[2]);
  particleGun->SetParticleDefinition(lambda1);
  particleGun->SetParticleMomentumDirection(momentumLambda1);
  particleGun->SetParticleEnergy((Energy_L1 - Lambda->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(vertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);

  // ---- Lambda 2 -------------
  G4ThreeVector momentumLambda2(pL2[0], pL2[1], pL2[2]);
  particleGun->SetParticleDefinition(lambda2);
  particleGun->SetParticleMomentumDirection(momentumLambda2);
  particleGun->SetParticleEnergy((Energy_L2 - Lambda->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(vertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,pka[0],pka[1],pka[2],kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,pL1[0],pL1[1],pL1[2],Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,pL2[0],pL2[1],pL2[2],Lambda->GetPDGMass()/GeV);

  gAnaMan.SetPrimaryVertex(0,0.,0.,vtx[2]);
  gAnaMan.SetPrimaryVertex(1,0.,0.,vtx[2]);
  gAnaMan.SetPrimaryVertex(2,0.,0.,vtx[2]);


  // G4double e1=0.,e2=0.,etot=0.,invm2=0.,ptot=0.,invm=0.;
  // G4double rmk=0.493677;
  // G4double misskp = miss1(pbm,rmk,momk);//pbeam, mass, momk

  //  missk->Fill(misskp);

  // ptot=pow(pL1[0]+pL2[0],2)+pow(pL1[1]+pL2[1],2)+pow(pL1[2]+pL2[2],2);
  // e1=Energy_L1;
  // e2=Energy_L2;
  // etot=pow(e1+e2,2);
  // invm2=(etot-ptot);
  // if(invm2 > 0) invm=sqrt(invm2);
  //  gen_im->Fill(invm);
}

//_____________________________________________________________________________
// generator #30
void
TPCPrimaryGeneratorAction::Generate_PhaseSpace( G4Event* anEvent )
{

  // G4double Energy_Be10,momentum_Be10[4]={0.};
  // G4double Energy_L1, momentum_L1[4]={0.};
  // G4double Energy_L2, momentum_L2[4]={0.};
  // G4double Energy_kp, momentum_kp[4]={0.};

  // G4double ThetaBe10, PhiBe10;
  // G4double ThetaL1, PhiL1;
  // G4double ThetaL2, PhiL2;
  // G4double Thetakp1, Phikp1;

  // G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  auto particleTable = G4ParticleTable::GetParticleTable();
  // auto kaonPlus = particleTable->FindParticle("kaon+");
  // auto kaonMinus = particleTable->FindParticle("kaon-");
  // auto Lambda1 = particleTable->FindParticle("lambda");
  // auto Lambda2 = particleTable->FindParticle("lambda");
  auto LLphase = particleTable->FindParticle("phaseLL");
  // G4int Z, A;
  // G4double excitEnergy = 0.*keV;
  // G4double ionCharge   = 0.*eplus;
  // auto Carbon12 = G4IonTable::GetIonTable()->GetIon( Z=6, A=12, excitEnergy );
  // auto Beryllium10 = G4IonTable::GetIonTable()->GetIon( Z=4, A=10, excitEnergy );
  // kaonMinus->DumpTable();
  // Carbon12->DumpTable();
  // Beryllium10->DumpTable();
  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548)*GeV;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  // G4double Ebeam = sqrt(pow(pbeam/GeV,2)+pow(kaonMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  // G4double rmC12=Carbon12->GetPDGMass()/GeV;
  // G4double rmkp=kaonPlus->GetPDGMass()/GeV;
  // G4double rmBe10=Beryllium10->GetPDGMass()/GeV;
  // G4double rmL=Lambda1->GetPDGMass()/GeV;

  // G4double W=sqrt(pow(Ebeam+rmC12,2)-pow(pbeam/GeV,2));
  //  G4cout<<"C mass------->"<<Carbon12->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"Be mass------->"<<Beryllium10->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"W------->"<<W<<G4endl;
  //  G4cout<<"Ebeam------->"<<Ebeam<<G4endl;
  //  G4cout<<"pbeam------->"<<pbeam<<G4endl;
  G4double Energy_LLphase=sqrt(pow(LLphase->GetPDGMass()/GeV,2)+pow(pbeam/GeV,2));

  G4double vtx,vty,vtz;
  G4double rn_vtx,rn_vtz;
  while(1){
    rn_vtx = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    rn_vtz = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot( -m_target_size.z(), m_target_size.z() );
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;

  particleGun->SetParticleDefinition(LLphase);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy((Energy_LLphase - LLphase->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,pg_x,pg_y,pg_z, LLphase->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);


  /*
  Kinema4Body PhaseLL(W,0.,
		      Beryllium10->GetPDGMass()/GeV,
		      Lambda1->GetPDGMass()/GeV,
		      Lambda2->GetPDGMass()/GeV,
		      kaonPlus->GetPDGMass()/GeV,
		      pbeam, 0.0);
  */
  /*
//  Hdibaryon.Dump();
 // pro
  Energy_pro = Hybrid.GetEnergy(3);
  momentum_pro = Hybrid.GetMomentum(3);
  Hybrid.GetMomentum(3,mom);

  mom_pro_x = mom[0];
  mom_pro_y = mom[1];
  mom_pro_z = mom[2];

  Thetapro = Hybrid.GetTheta(3);
  Phipro = Hybrid.GetPhi(3);

  // L2
  Energy_pi1 = Hybrid.GetEnergy(4);
  momentum_pi1 = Hybrid.GetMomentum(4);
  Hybrid.GetMomentum(4,mom);

  mom_pi1_x = mom[0];
  mom_pi1_y = mom[1];
  mom_pi1_z = mom[2];

  Thetapi1 = Hybrid.GetTheta(4);
  Phipi1 = Hybrid.GetPhi(4);

  // pi2
  Energy_pi2 = Hybrid.GetEnergy(5);
  momentum_pi2 = Hybrid.GetMomentum(5);
  Hybrid.GetMomentum(5,mom);

  mom_pi2_x = mom[0];
  mom_pi2_y = mom[1];
  mom_pi2_z = mom[2];

  Thetapi2 = Hybrid.GetThetaCM(1);
  G4double shphi=(atan2(mom_pro_y,mom_pro_x))/3.141592654*180.;
  phik->Fill(shphi);

  G4double e1=0.,e2=0.,e3=0,etot=0.,invm2=0.,ptot=0.,invm=0.;

  ptot=pow(mom_pro_x+mom_pi1_x+mom_pi2_x,2)+pow(mom_pro_y+mom_pi1_y+mom_pi2_y,2)+pow(mom_pro_z+mom_pi1_z+mom_pi2_z,2);
  e1=Energy_pro;
  e2=Energy_pi1;
  e3=Energy_pi2;
  etot=pow(e1+e2+e3,2);
  invm2=etot-ptot;

  if(invm2 > 0) invm=sqrt(invm2);
  gen_im->Fill(invm);


  ///////neutron PG
  particleGun->SetParticleDefinition(neutron);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  particleGun->SetParticleEnergy((Energy_pro - neutron->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi1 PG
  particleGun->SetParticleDefinition(piMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi2 PG
  particleGun->SetParticleDefinition(piPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  particleGun->SetParticleEnergy((Energy_pi2 - piPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z, neutron->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z, piMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z, piPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

*/



  /*
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
 up:
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			(proton->GetPDGMass()/GeV)*2,
			hdibaryon->GetPDGMass()/GeV,
			kaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hdibaryon.GetEnergy(4);
  momentum_kp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];


  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hdibaryon.GetEnergy(3);
  momentum_h = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  Theta_h = Hdibaryon.GetTheta(3);
  Phi_h = Hdibaryon.GetPhi(3);

  ////check kinematics
  momk[0]=mom_kp_x;
  momk[1]=mom_kp_y;
  momk[2]=mom_kp_z;
  momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  coscmk->Fill(momcmk[2]/momcmkpp);


  //  G4double misskp = miss1(pbeam,rmk,momk);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;


  //  G4cout<<"env_target_width :"<<atof(env_env_target_width.c_str())<<G4endl;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-env_target_width/2,env_target_pos_z+env_target_width/2)*mm;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);

  //Kaon +
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //hdibaryon
  particleGun->SetParticleDefinition(hdibaryon);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - hdibaryon->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,hdibaryon->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  */
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hdibaryon_PHSG( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  auto particleTable = G4ParticleTable::GetParticleTable();
  auto kaonPlus = particleTable->FindParticle("kaon+");
  auto kaonMinus = particleTable->FindParticle("kaon-");
  auto proton = particleTable->FindParticle("proton");
  auto hdibaryon = particleTable->FindParticle("hdibaryon");
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			(proton->GetPDGMass()/GeV)*2,
			hdibaryon->GetPDGMass()/GeV,
			kaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hdibaryon.GetEnergy(4);
  // G4double momentum_kp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hdibaryon.GetEnergy(3);
  // G4double momentum_h = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hdibaryon.GetTheta(3);
  // G4double Phi_h = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);


  //  G4double misskp = miss1(pbeam,rmk,momk);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz = G4RandFlat::shoot( env_target_pos_z - m_target_size.z()/2,
					 env_target_pos_z + m_target_size.z()/2 )*mm;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);

  //Kaon +
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //hdibaryon
  particleGun->SetParticleDefinition(hdibaryon);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - hdibaryon->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,hdibaryon->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hdibaryon_PHSG_S( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  auto particleTable = G4ParticleTable::GetParticleTable();
  auto kaonPlus = particleTable->FindParticle("kaon+");
  auto kaonMinus = particleTable->FindParticle("kaon-");
  auto proton = particleTable->FindParticle("proton");
  auto hdibaryonS = particleTable->FindParticle("hdibaryonS");
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
 up:
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			(proton->GetPDGMass()/GeV)*2,
			hdibaryonS->GetPDGMass()/GeV,
			kaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hdibaryon.GetEnergy(4);
  // G4double momentum_kp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hdibaryon.GetEnergy(3);
  // G4double momentum_h = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hdibaryon.GetTheta(3);
  // G4double Phi_h = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);


  //  G4double misskp = miss1(pbeam,rmk,momk);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  //  G4cout<<"env_target_width :"<<atof(env_env_target_width.c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);
  G4double vtz= G4RandFlat::shoot( env_target_pos_z - m_target_size.z()/2,
					env_target_pos_z + m_target_size.z()/2 )*mm;

  //Kaon +
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //hdibaryon
  particleGun->SetParticleDefinition(hdibaryonS);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - hdibaryonS->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,hdibaryonS->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hdibaryon_PHSG_LL( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  auto particleTable = G4ParticleTable::GetParticleTable();
  auto kaonPlus = particleTable->FindParticle("kaon+");
  auto kaonMinus = particleTable->FindParticle("kaon-");
  auto proton = particleTable->FindParticle("proton");
  auto hdibaryonLL = particleTable->FindParticle("hdibaryonLL");
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
 up:
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			(proton->GetPDGMass()/GeV)*2,
			hdibaryonLL->GetPDGMass()/GeV,
			kaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hdibaryon.GetEnergy(4);
  // G4double momentum_kp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];


  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.  ) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hdibaryon.GetEnergy(3);
  // G4double momentum_h = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hdibaryon.GetTheta(3);
  // G4double Phi_h = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);


  //  G4double misskp = miss1(pbeam,rmk,momk);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz = G4RandFlat::shoot( env_target_pos_z - m_target_size.z()/2,
					 env_target_pos_z + m_target_size.z()/2 )*mm;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);

  //Kaon +
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //hdibaryon
  particleGun->SetParticleDefinition(hdibaryonLL);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - hdibaryonLL->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,hdibaryonLL->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_Kp_Kn( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double Energy_kn, mom_kn_x, mom_kn_y, mom_kn_z;
  auto particleTable = G4ParticleTable::GetParticleTable();
  auto kaonPlus = particleTable->FindParticle("kaon+");
  auto kaonMinus = particleTable->FindParticle("kaon-");
  auto proton = particleTable->FindParticle("proton");
  auto hdibaryonLL = particleTable->FindParticle("hdibaryonLL");
  pbeam=1.8;
  mom_kn_x=0;
  mom_kn_y=0;
  mom_kn_z=pbeam;
  Energy_kn=sqrt(kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV+pbeam*pbeam);

  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
 up:
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			(proton->GetPDGMass()/GeV)*2,
			hdibaryonLL->GetPDGMass()/GeV,
			kaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);


  Energy_kp=Hdibaryon.GetEnergy(4);
  // G4double momentum_kp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.  ) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  // G4double Energy_h = Hdibaryon.GetEnergy(3);
  // G4double momentum_h = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  // G4double mom_h_x = mom[0];
  // G4double mom_h_y = mom[1];
  // G4double mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hdibaryon.GetTheta(3);
  // G4double Phi_h = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);

  //  G4double misskp = miss1(pbeam,rmk,momk);

  //  G4double vtx = 0.*mm;
  //  G4double vty = 0.*mm;
  G4double vtx = G4RandFlat::shoot(-15.,15.)*mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*mm;
  G4double vtz = G4RandFlat::shoot( env_target_pos_z - m_target_size.z()/2,
					 env_target_pos_z + m_target_size.z()/2 )*mm;

  //Kaon -
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //beam K+
  particleGun->SetParticleDefinition(kaonMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kn_x,mom_kn_y,mom_kn_z));
  particleGun->SetParticleEnergy((Energy_kn - kaonMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz-150.*mm));
  particleGun->GeneratePrimaryVertex(anEvent);



  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_kn_x,mom_kn_y,mom_kn_z,kaonMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateBeam( G4Event* anEvent )
{
  static const auto particleTable = G4ParticleTable::GetParticleTable();
  static const auto kaonMinus = particleTable->FindParticle("kaon-");
  static const G4double mass = kaonMinus->GetPDGMass()/GeV;
  // pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  // pbeam=1.8;
  G4double energy = ( std::sqrt( mass*mass + m_beam->p.mag2() ) - mass )*GeV;
  gAnaMan.SetPrimaryBeam( m_beam->p );
  G4ThreeVector gen_pos( m_beam->x, m_beam->y, m_target_pos.z() );
  particleGun->SetParticleDefinition( kaonMinus );
  particleGun->SetParticleMomentumDirection( m_beam->p );
  particleGun->SetParticleEnergy( energy );
  particleGun->SetParticlePosition( gen_pos );
  particleGun->GeneratePrimaryVertex( anEvent );
  // gAnaMan.SetNumberOfPrimaryParticle( 1 );
  // gAnaMan.SetPrimaryParticle( 0, mom_kn_x, mom_kn_y, mom_kn_z, mass );
  // gAnaMan.SetPrimaryVertex( 0, vtx, vty, vtz );
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hdibaryon_non_reso( G4Event* anEvent )
{
  G4double  momk[4]={0.}, mom[3]={0.};
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  auto particleTable = G4ParticleTable::GetParticleTable();
  auto kaonPlus = particleTable->FindParticle("kaon+");
  auto kaonMinus = particleTable->FindParticle("kaon-");
  auto proton = particleTable->FindParticle("proton");
  auto hdibaryon = particleTable->FindParticle("hdibaryon");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
 up:
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			(proton->GetPDGMass()/GeV)*2,
			hdibaryon->GetPDGMass()/GeV,
			kaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hdibaryon.GetEnergy(4);
  // G4double momentum_kp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];


  if((acos(momk[2]/sqrt(pow(momk[0],2)+pow(momk[1],2)+pow(momk[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hdibaryon.GetEnergy(3);
  // G4double momentum_h = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hdibaryon.GetTheta(3);
  // G4double Phi_h = Hdibaryon.GetPhi(3);

  ////check kinematics
  momk[0]=mom_kp_x;
  momk[1]=mom_kp_y;
  momk[2]=mom_kp_z;
  momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);


  //  G4double misskp = miss1(pbeam,rmk,momk);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz = G4RandFlat::shoot( env_target_pos_z - m_target_size.z()/2,
					 env_target_pos_z + m_target_size.z()/2 )*mm;

  //Kaon +
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //hdibaryon
  particleGun->SetParticleDefinition(hdibaryon);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - hdibaryon->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,hdibaryon->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hybrid( G4Event* anEvent )
{
  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;
  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  // G4double rmpro = 0.93827203;
  // G4double rmpi = 0.13957018;

  auto particleTable = G4ParticleTable::GetParticleTable();
  auto proton = particleTable->FindParticle("proton");
  auto piPlus = particleTable->FindParticle("pi+");
  auto piMinus = particleTable->FindParticle("pi-");

  //  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  G4double pbeam=1.;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  // G4double Ebeam = sqrt(pow(pbeam,2)+pow(piMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  mass_hybrid = 1.22;    // Mass for H
  width_hybrid = 0.000;  // Width for H

  KinemaHybrid Hybrid(piMinus->GetPDGMass()/GeV,
		      proton->GetPDGMass()/GeV,
		      proton->GetPDGMass()/GeV,
		      piMinus->GetPDGMass()/GeV,//pi1
		      piPlus->GetPDGMass()/GeV, //pi2
		      mass_hybrid, width_hybrid, pg_z, 0.0);

  //  Hdibaryon.Dump();
  /* pro */
  Energy_pro = Hybrid.GetEnergy(3);
  // G4double momentum_pro = Hybrid.GetMomentum(3);
  Hybrid.GetMomentum(3,mom);

  mom_pro_x = mom[0];
  mom_pro_y = mom[1];
  mom_pro_z = mom[2];

  // G4double Thetapro = Hybrid.GetTheta(3);
  // G4double Phipro = Hybrid.GetPhi(3);

  /* L2 */
  Energy_pi1 = Hybrid.GetEnergy(4);
  // G4double momentum_pi1 = Hybrid.GetMomentum(4);
  Hybrid.GetMomentum(4,mom);

  mom_pi1_x = mom[0];
  mom_pi1_y = mom[1];
  mom_pi1_z = mom[2];

  // G4double Thetapi1 = Hybrid.GetTheta(4);
  // G4double Phipi1 = Hybrid.GetPhi(4);

  /* pi2 */
  Energy_pi2 = Hybrid.GetEnergy(5);
  // G4double momentum_pi2 = Hybrid.GetMomentum(5);
  Hybrid.GetMomentum(5,mom);

  mom_pi2_x = mom[0];
  mom_pi2_y = mom[1];
  mom_pi2_z = mom[2];

  // G4double Thetapi2 = Hybrid.GetThetaCM(1);
  // G4double shphi=(atan2(mom_pro_y,mom_pro_x))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4double beta= sqrt(pbeam*pbeam)/(0.938272013+Ebeam);
  /*  G4double test;
  G4double mom1[4];
  G4double mom2[4];
  G4double mom3[4];

  mom1[0]=mom_pro_x;
  mom1[1]=mom_pro_y;
  mom1[2]=mom_pro_z;
  mom1[3]=sqrt(mom_pro_x*mom_pro_x+mom_pro_y*mom_pro_y+mom_pro_z*mom_pro_z+rmpro*rmpro);
  momppro=sqrt(pow(mom_pro_x,2)+pow(mom_pro_y,2)+pow(mom_pro_z,2));
  */
  /*  test=lorentz(momk,beta,momcmk);
  momcmprop=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  cmk->Fill(acos(momcmk[2]/momcmprop)/3.141592654*180);

  coscmk->Fill((momcmk[2]/momcmprop));

  labk->Fill(acos(momk[2]/momprop)*180/3.141592654);
  coslabk->Fill(momk[2]/momprop);

  G4double misspro = miss1(pbeam,rmk,momk);//pbeam, mass, momk

  missk->Fill(misskp);
  */

  G4double e1=0.,e2=0.,e3=0,etot=0.,invm2=0.,ptot=0.,invm=0.;

  ptot=pow(mom_pro_x+mom_pi1_x+mom_pi2_x,2)+pow(mom_pro_y+mom_pi1_y+mom_pi2_y,2)+pow(mom_pro_z+mom_pi1_z+mom_pi2_z,2);
  e1=Energy_pro;
  e2=Energy_pi1;
  e3=Energy_pi2;
  etot=pow(e1+e2+e3,2);
  invm2=etot-ptot;

  if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);
  G4cout<<"invm:"<<invm<<G4endl;
  /////shhwang hybrid1

  double vtx = 0.*mm;
  double vty = 0.*mm;
  double vtz = -155.0+10.0*((double) G4RandFlat::shoot()); //--> 10 mm
  //  double vtz = -157.5+15.0*((double) G4RandFlat::shoot()); //--> 15 mm
  //  double vtz = -160.0+20.0*((double) G4RandFlat::shoot()); //--> 20 mm
  //  double vtz = -162.5+25.0*((double) G4RandFlat::shoot()); //--> 25 mm
  //    double vtz = -165.0+30.0*((double) G4RandFlat::shoot()); //--> 30 mm
  //  double vtz = -200.0+100.0*((double) G4RandFlat::shoot()); //--> test//

  ///////proton PG
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  particleGun->SetParticleEnergy((Energy_pro - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi1 PG
  particleGun->SetParticleDefinition(piMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi2 PG
  particleGun->SetParticleDefinition(piPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  particleGun->SetParticleEnergy((Energy_pi2 - piPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z,proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z,piMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z,piPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hybrid3body( G4Event* anEvent )
{
  //  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;
  G4double mom[3];
  G4double Ebeam,pg_x,pg_y,pg_z;
  G4double rmpro = 0.93827203;
  // G4double rmpi = 0.13957018;

  auto particleTable = G4ParticleTable::GetParticleTable();
  auto proton = particleTable->FindParticle("proton");
  auto piPlus = particleTable->FindParticle("pi+");
  auto piMinus = particleTable->FindParticle("pi-");

  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  G4double pbeam=0.7;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  Ebeam = sqrt(pow(pbeam,2)+pow(piMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  //  mass_hybrid = 1.22;    // Mass for H
  //  width_hybrid = 0.000;  // Width for H

  //  double temp= kin1.GetMomentumLab(3);
  //    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
  //  Kinema3Body kin2(kin3.M_res, 0.0, m3, m4, m5, temp, 0.0);
 G4double W=sqrt(pow(Ebeam+rmpro,2)-pow(pbeam,2));
 Kinema3Body Hybrid(W,0,
		    proton->GetPDGMass()/GeV,
		    piMinus->GetPDGMass()/GeV,
		    piPlus->GetPDGMass()/GeV,
		    pbeam, 0.0);
//  Hdibaryon.Dump();
  /* pro */
  Energy_pro = Hybrid.GetEnergy(3);
  // G4double momentum_pro = Hybrid.GetMomentum(3);
  Hybrid.GetMomentum(3,mom);

  mom_pro_x = mom[0];
  mom_pro_y = mom[1];
  mom_pro_z = mom[2];

  // G4double Thetapro = Hybrid.GetTheta(3);
  // G4double Phipro = Hybrid.GetPhi(3);

  /* L2 */
  Energy_pi1 = Hybrid.GetEnergy(4);
  // G4double momentum_pi1 = Hybrid.GetMomentum(4);
  Hybrid.GetMomentum(4,mom);

  mom_pi1_x = mom[0];
  mom_pi1_y = mom[1];
  mom_pi1_z = mom[2];

  // G4double Thetapi1 = Hybrid.GetTheta(4);
  // G4double Phipi1 = Hybrid.GetPhi(4);

  /* pi2 */
  Energy_pi2 = Hybrid.GetEnergy(5);
  // G4double momentum_pi2 = Hybrid.GetMomentum(5);
  Hybrid.GetMomentum(5,mom);

  mom_pi2_x = mom[0];
  mom_pi2_y = mom[1];
  mom_pi2_z = mom[2];

  // G4double Thetapi2 = Hybrid.GetThetaCM(1);
  // G4double shphi=(atan2(mom_pro_y,mom_pro_x))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4double beta= sqrt(pbeam*pbeam)/(0.938272013+Ebeam);
  /*  G4double test;
  G4double mom1[4];
  G4double mom2[4];
  G4double mom3[4];

  mom1[0]=mom_pro_x;
  mom1[1]=mom_pro_y;
  mom1[2]=mom_pro_z;
  mom1[3]=sqrt(mom_pro_x*mom_pro_x+mom_pro_y*mom_pro_y+mom_pro_z*mom_pro_z+rmpro*rmpro);
  momppro=sqrt(pow(mom_pro_x,2)+pow(mom_pro_y,2)+pow(mom_pro_z,2));
  */
  /*  test=lorentz(momk,beta,momcmk);
  momcmprop=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  cmk->Fill(acos(momcmk[2]/momcmprop)/3.141592654*180);

  coscmk->Fill((momcmk[2]/momcmprop));

  labk->Fill(acos(momk[2]/momprop)*180/3.141592654);
  coslabk->Fill(momk[2]/momprop);

  G4double misspro = miss1(pbeam,rmk,momk);//pbeam, mass, momk

  missk->Fill(misskp);
  */

  // G4double e1=0.,e2=0.,e3=0,etot=0.,invm2=0.,ptot=0.,invm=0.;

  // ptot=pow(mom_pro_x+mom_pi1_x+mom_pi2_x,2)+pow(mom_pro_y+mom_pi1_y+mom_pi2_y,2)+pow(mom_pro_z+mom_pi1_z+mom_pi2_z,2);
  // e1=Energy_pro;
  // e2=Energy_pi1;
  // e3=Energy_pi2;
  // etot=pow(e1+e2+e3,2);
  // invm2=etot-ptot;

  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);
  //  G4cout<<"invm:"<<invm<<G4endl;
  /////shhwang hybrid1

  double vtx = 0.*mm;
  double vty = 0.*mm;
  double vtz = -155.0+10.0*((double) G4RandFlat::shoot()); //--> 10 mm
  //  double vtz = -157.5+15.0*((double) G4RandFlat::shoot()); //--> 15 mm
  //  double vtz = -160.0+20.0*((double) G4RandFlat::shoot()); //--> 20 mm
  //  double vtz = -162.5+25.0*((double) G4RandFlat::shoot()); //--> 25 mm
  //    double vtz = -165.0+30.0*((double) G4RandFlat::shoot()); //--> 30 mm
  //  double vtz = -200.0+100.0*((double) G4RandFlat::shoot()); //--> test//

  ///////proton PG
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  particleGun->SetParticleEnergy((Energy_pro - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi1 PG
  particleGun->SetParticleDefinition(piMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi2 PG
  particleGun->SetParticleDefinition(piPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  particleGun->SetParticleEnergy((Energy_pi2 - piPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z,proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z,piMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z,piPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
// E45 elastic scattering pip
void
TPCPrimaryGeneratorAction::Generate_E45_elastic_pip( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam=0.;

  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  auto particleTable = G4ParticleTable::GetParticleTable();
  auto pip = particleTable->FindParticle("pi+");
  auto proton = particleTable->FindParticle("proton");

  G4double protonMass=proton->GetPDGMass()/GeV;//unit GeV
  G4double pipMass=pip->GetPDGMass()/GeV;//unit GeV

  /*
///first w/o beam
  ///beam optics from simulation file.
  char fname[100] ;
  FILE *fp;
  sprintf(fname,"./beam_simulation/profile_ve07-5.dat.txt");
  if ((fp = fopen(fname,"r")) == NULL){
    fprintf(stderr, ": Cannot open file: %s\n", fname);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check;
  int ran;
  int res;
 up1:
  check=0.;
  ran=G4RandFlat::shoot(1,17005);


  while(1){
    /// file structure : x, u, y, v, p, PID, ???
    /// unit           : cm, mrad, cm, mrad, gev, PID, ???
    res=fscanf(fp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf",&data[0],&data[1],&data[2],&data[3],&data[4],&data[5], &data[6]);
    check=check+1.;
    if(check==ran) break;
    //    if(res==EOF) break;
  }

  G4double dxdz,dydz,pp;
  dxdz=atan(data[1]*0.001);
  dydz=atan(data[3]*0.001);
  pp=data[4]/1.8*m_beam_p0;

  //  G4cout<<"pp:"<<pp<<G4endl;
  //  G4cout<<"data 4:"<<data[4]<<G4endl;


  */
  // G4double x0;
  // G4double y0;
  // G4double z0;
  // G4double mom_beam;
  G4String particleName;
  // G4ParticleDefinition* particle;

  G4double vtx;
  G4double vty;
  G4double vtz;
  ///for E45
  G4double rn_vtx,rn_vtz;
  while(1){
    rn_vtx = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    rn_vtz = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    //    G4cout<<rn_vtx<<G4endl;
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot( -m_target_size.z(), m_target_size.z() );
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;
  /*   ///E42
  vtz= G4RandFlat::shoot(env_Target_pos_z-env_Target_width/2.,env_Target_pos_z+env_Target_width/2.)*mm;
  vtx = (data[0]*10.+dxdz*(vtz-env_Target_pos_z))*mm;
  vty = (data[2]*10.+dydz*(vtz-env_Target_pos_z))*mm;

  //  G4cout<<"vtx and vty ::"<<fabs(vtx)<<", "<<fabs(vty)<<G4endl;
  //  G4cout<<"target_x, target_y::"<<env_Target_x<<", "<<env_Target_y<<G4endl;
  if(fabs(vtx)>env_Target_x/2. || fabs(vty)>env_Target_y/2.){
    goto up1;
  }
  */
  //  G4cout<<"passed-->vtx and vty ::"<<fabs(vtx)<<", "<<fabs(vty)<<G4endl;

  ///close file, if not it will make an error.
  //  fclose(fp);

  G4double pbeam_x;
  G4double pbeam_y;
  G4double pbeam_z;
  G4double pp;
  pp=m_beam_p0+G4RandGauss::shoot(0.,0.01294*m_beam_p0);
  pbeam_z=pp;
  pbeam_x=0.;
  pbeam_y=0.;

  gAnaMan.SetPrimaryBeam(pbeam_x,pbeam_y,pbeam_z);
  G4double cosx = G4RandFlat::shoot(-1.,1.);
  G4double p_proton[4]={0};
  //  gAnaMan.SetFermiMotion(p_proton);

  ///prepare real distribution
  Ebeam = sqrt(pbeam*pbeam+pipMass/GeV*pipMass/GeV);
  pbm[0]=pbeam_x;
  pbm[1]=pbeam_y;
  pbm[2]=pbeam_z;
  pbm[3]=Ebeam;
  ///first pip
  KinemaFermi Hdibaryon(pip->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			pip->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			pbm, p_proton,cosx);

  Energy_kp=Hdibaryon.GetEnergy(3);
  // G4double momentum_kp = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>30.) goto up_sh;
  //  G4cout<<"kp mom:"<<mom[2]<<G4endl;

  //  cross_section->Fill(cosx);
  //  gAnaMan.SetCrossSection(cosx);

  Energy_h = Hdibaryon.GetEnergy(4);
  // G4double momentum_h = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  cout<<mom_h_z<<endl;
  // G4double Theta_h = Hdibaryon.GetTheta(4);
  // G4double Phi_h = Hdibaryon.GetPhi(4);
  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+pipMass*pipMass);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0};
  //  G4double beta= pbeam/(proton->GetPDGMass()/GeV+Ebeam);
  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);

  //pi
  particleGun->SetParticleDefinition(pip);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - pipMass/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //proton
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - protonMass/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,pipMass/GeV);///pip
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,protonMass/GeV);///proton
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
// E45 elastic scattering pin
void
TPCPrimaryGeneratorAction::Generate_E45_elastic_pin( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam=0.;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pin;
  G4ParticleDefinition* proton;

  pin = particleTable->FindParticle("pi-");
  proton = particleTable->FindParticle("proton");

  G4double protonMass=proton->GetPDGMass()/GeV;//unit GeV
  G4double pinMass=pin->GetPDGMass()/GeV;//unit GeV

  /*
///first w/o beam
  ///beam optics from simulation file.
  char fname[100] ;
  FILE *fp;
  sprintf(fname,"./beam_simulation/profile_ve07-5.dat.txt");
  if ((fp = fopen(fname,"r")) == NULL){
    fprintf(stderr, ": Cannot open file: %s\n", fname);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check;
  int ran;
  int res;
 up1:
  check=0.;
  ran=G4RandFlat::shoot(1,17005);


  while(1){
    /// file structure : x, u, y, v, p, PID, ???
    /// unit           : cm, mrad, cm, mrad, gev, PID, ???
    res=fscanf(fp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf",&data[0],&data[1],&data[2],&data[3],&data[4],&data[5], &data[6]);
    check=check+1.;
    if(check==ran) break;
    //    if(res==EOF) break;
  }

  G4double dxdz,dydz,pp;
  dxdz=atan(data[1]*0.001);
  dydz=atan(data[3]*0.001);
  pp=data[4]/1.8*m_beam_p0;

  //  G4cout<<"pp:"<<pp<<G4endl;
  //  G4cout<<"data 4:"<<data[4]<<G4endl;


  */
  // G4double x0;
  // G4double y0;
  // G4double z0;
  // G4double mom_beam;
  G4String particleName;
  // G4ParticleDefinition* particle;

  G4double vtx;
  G4double vty;
  G4double vtz;
  ///for E45
  G4double rn_vtx,rn_vtz;
  while(1){
    rn_vtx = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    rn_vtz = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    //    G4cout<<rn_vtx<<G4endl;
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot( -m_target_size.z(), m_target_size.z() );
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;
  /*   ///E42
  vtz= G4RandFlat::shoot(env_Target_pos_z-env_Target_width/2.,env_Target_pos_z+env_Target_width/2.)*mm;
  vtx = (data[0]*10.+dxdz*(vtz-env_Target_pos_z))*mm;
  vty = (data[2]*10.+dydz*(vtz-env_Target_pos_z))*mm;

  //  G4cout<<"vtx and vty ::"<<fabs(vtx)<<", "<<fabs(vty)<<G4endl;
  //  G4cout<<"target_x, target_y::"<<env_Target_x<<", "<<env_Target_y<<G4endl;
  if(fabs(vtx)>env_Target_x/2. || fabs(vty)>env_Target_y/2.){
    goto up1;
  }
  */
  //  G4cout<<"passed-->vtx and vty ::"<<fabs(vtx)<<", "<<fabs(vty)<<G4endl;

  ///close file, if not it will make an error.
  //  fclose(fp);

  G4double pbeam_x;
  G4double pbeam_y;
  G4double pbeam_z;
  G4double pp;
  pp=m_beam_p0+G4RandGauss::shoot(0.,0.01294*m_beam_p0);
  pbeam_z=pp;
  pbeam_x=0.;
  pbeam_y=0.;

  gAnaMan.SetPrimaryBeam(pbeam_x,pbeam_y,pbeam_z);
  G4double cosx = G4RandFlat::shoot(-1.,1.);
  G4double p_proton[4]={0};
  //  gAnaMan.SetFermiMotion(p_proton);

  ///prepare real distribution
  Ebeam = sqrt(pbeam*pbeam+pinMass/GeV*pinMass/GeV);
  pbm[0]=pbeam_x;
  pbm[1]=pbeam_y;
  pbm[2]=pbeam_z;
  pbm[3]=Ebeam;
  ///first pin
  KinemaFermi Hdibaryon(pin->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			pin->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			pbm, p_proton,cosx);

  Energy_kp=Hdibaryon.GetEnergy(3);
  // G4double momentum_kp = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>30.) goto up_sh;
  //  G4cout<<"kp mom:"<<mom[2]<<G4endl;

  //  cross_section->Fill(cosx);
  //  gAnaMan.SetCrossSection(cosx);

  Energy_h = Hdibaryon.GetEnergy(4);
  // G4double   momentum_h = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  cout<<mom_h_z<<endl;
  // G4double   Theta_h = Hdibaryon.GetTheta(4);
  // G4double   Phi_h = Hdibaryon.GetPhi(4);
  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+pinMass*pinMass);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0};
  //  G4double beta= pbeam/(proton->GetPDGMass()/GeV+Ebeam);
  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);

  //pi
  particleGun->SetParticleDefinition(pin);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - pinMass/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //proton
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - protonMass/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,pinMass/GeV);///pin
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,protonMass/GeV);///proton
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hybrid3body_mode1( G4Event* anEvent )
{
  //  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;

  G4double mom[3];
  G4double Ebeam,pg_x,pg_y,pg_z;
  G4double rmpro=0.93827203;
  // G4double rmpi=0.13957018;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* neutron;
  G4ParticleDefinition* piMinus;
  G4ParticleDefinition* piPlus;

  neutron = particleTable->FindParticle("neutron");
  piPlus = particleTable->FindParticle("pi+");
  piMinus = particleTable->FindParticle("pi-");

  //  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*0.01294);
  //  G4double pbeam=0.7;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  Ebeam = sqrt(pow(pbeam,2)+pow(piMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

 G4double W=sqrt(pow(Ebeam+rmpro,2)-pow(pbeam,2));
 Kinema3Body Hybrid(W,0,
		    neutron->GetPDGMass()/GeV,
		    piMinus->GetPDGMass()/GeV,
		    piPlus->GetPDGMass()/GeV,
		    pbeam, 0.0);
//  Hdibaryon.Dump();
  /* pro */
  Energy_pro = Hybrid.GetEnergy(3);
  // G4double momentum_pro = Hybrid.GetMomentum(3);
  Hybrid.GetMomentum(3,mom);

  mom_pro_x = mom[0];
  mom_pro_y = mom[1];
  mom_pro_z = mom[2];

  // G4double   Thetapro = Hybrid.GetTheta(3);
  // G4double   Phipro = Hybrid.GetPhi(3);

  /* L2 */
  Energy_pi1 = Hybrid.GetEnergy(4);
  // G4double   momentum_pi1 = Hybrid.GetMomentum(4);
  Hybrid.GetMomentum(4,mom);

  mom_pi1_x = mom[0];
  mom_pi1_y = mom[1];
  mom_pi1_z = mom[2];

  // G4double   Thetapi1 = Hybrid.GetTheta(4);
  // G4double   Phipi1 = Hybrid.GetPhi(4);

  /* pi2 */
  Energy_pi2 = Hybrid.GetEnergy(5);
  // G4double   momentum_pi2 = Hybrid.GetMomentum(5);
  Hybrid.GetMomentum(5,mom);

  mom_pi2_x = mom[0];
  mom_pi2_y = mom[1];
  mom_pi2_z = mom[2];

  // G4double   Thetapi2 = Hybrid.GetThetaCM(1);
  // G4double shphi=(atan2(mom_pro_y,mom_pro_x))/3.141592654*180.;
  //phik->Fill(shphi);

  // G4double e1=0.,e2=0.,e3=0,etot=0.,invm2=0.,ptot=0.,invm=0.;

  // ptot=pow(mom_pro_x+mom_pi1_x+mom_pi2_x,2)+pow(mom_pro_y+mom_pi1_y+mom_pi2_y,2)+pow(mom_pro_z+mom_pi1_z+mom_pi2_z,2);
  // e1=Energy_pro;
  // e2=Energy_pi1;
  // e3=Energy_pi2;
  // etot=pow(e1+e2+e3,2);
  // invm2=etot-ptot;

  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);

  G4double vtx,vty,vtz;
  G4double rn_vtx,rn_vtz;
  while(1){
    rn_vtx = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    rn_vtz = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    //    G4cout<<rn_vtx<<G4endl;
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot( -m_target_size.z(), m_target_size.z() );
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;
  //env_target_pos_z

  ///////neutron PG
  particleGun->SetParticleDefinition(neutron);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  particleGun->SetParticleEnergy((Energy_pro - neutron->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi1 PG
  particleGun->SetParticleDefinition(piMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi2 PG
  particleGun->SetParticleDefinition(piPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  particleGun->SetParticleEnergy((Energy_pi2 - piPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z, neutron->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z, piMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z, piPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hybrid3body_mode2( G4Event* anEvent )
{
  //  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;

  G4double mom[3];
  G4double Ebeam,pg_x,pg_y,pg_z;
  G4double rmpro=0.93827203;
  // G4double rmpi=0.13957018;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* proton;
  G4ParticleDefinition* piMinus;
  G4ParticleDefinition* piPlus;

  proton = particleTable->FindParticle("proton");
  piPlus = particleTable->FindParticle("pi0");
  piMinus = particleTable->FindParticle("pi-");

  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*0.01294);
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  G4double pbeam=0.7;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  Ebeam = sqrt(pow(pbeam,2)+pow(piMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  //  mass_hybrid = 1.22;    // Mass for H
  //  width_hybrid = 0.000;  // Width for H

  //  double temp= kin1.GetMomentumLab(3);
  //    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
  //  Kinema3Body kin2(kin3.M_res, 0.0, m3, m4, m5, temp, 0.0);
 G4double W=sqrt(pow(Ebeam+rmpro,2)-pow(pbeam,2));
 Kinema3Body Hybrid(W,0,
		    proton->GetPDGMass()/GeV,
		    piMinus->GetPDGMass()/GeV,
		    piPlus->GetPDGMass()/GeV,
		    pbeam, 0.0);
//  Hdibaryon.Dump();
  /* pro */
  Energy_pro = Hybrid.GetEnergy(3);
  // G4double   momentum_pro = Hybrid.GetMomentum(3);
  Hybrid.GetMomentum(3,mom);

  mom_pro_x = mom[0];
  mom_pro_y = mom[1];
  mom_pro_z = mom[2];

  // G4double Thetapro = Hybrid.GetTheta(3);
  // G4double   Phipro = Hybrid.GetPhi(3);

  /* L2 */
  Energy_pi1 = Hybrid.GetEnergy(4);
  // G4double   momentum_pi1 = Hybrid.GetMomentum(4);
  Hybrid.GetMomentum(4,mom);

  mom_pi1_x = mom[0];
  mom_pi1_y = mom[1];
  mom_pi1_z = mom[2];

  // G4double   Thetapi1 = Hybrid.GetTheta(4);
  // G4double   Phipi1 = Hybrid.GetPhi(4);

  /* pi2 */
  Energy_pi2 = Hybrid.GetEnergy(5);
  // G4double   momentum_pi2 = Hybrid.GetMomentum(5);
  Hybrid.GetMomentum(5,mom);

  mom_pi2_x = mom[0];
  mom_pi2_y = mom[1];
  mom_pi2_z = mom[2];

  // G4double   Thetapi2 = Hybrid.GetThetaCM(1);
  // G4double shphi=(atan2(mom_pro_y,mom_pro_x))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4double beta= sqrt(pbeam*pbeam)/(0.938272013+Ebeam);
  /*  G4double test;
  G4double mom1[4];
  G4double mom2[4];
  G4double mom3[4];

  mom1[0]=mom_pro_x;
  mom1[1]=mom_pro_y;
  mom1[2]=mom_pro_z;
  mom1[3]=sqrt(mom_pro_x*mom_pro_x+mom_pro_y*mom_pro_y+mom_pro_z*mom_pro_z+rmpro*rmpro);
  momppro=sqrt(pow(mom_pro_x,2)+pow(mom_pro_y,2)+pow(mom_pro_z,2));
  */
  /*  test=lorentz(momk,beta,momcmk);
  momcmprop=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  cmk->Fill(acos(momcmk[2]/momcmprop)/3.141592654*180);

  coscmk->Fill((momcmk[2]/momcmprop));

  labk->Fill(acos(momk[2]/momprop)*180/3.141592654);
  coslabk->Fill(momk[2]/momprop);

  G4double misspro = miss1(pbeam,rmk,momk);//pbeam, mass, momk

  missk->Fill(misskp);
  */

  // G4double e1=0.,e2=0.,e3=0,etot=0.,invm2=0.,ptot=0.,invm=0.;

  // ptot=pow(mom_pro_x+mom_pi1_x+mom_pi2_x,2)+pow(mom_pro_y+mom_pi1_y+mom_pi2_y,2)+pow(mom_pro_z+mom_pi1_z+mom_pi2_z,2);
  // e1=Energy_pro;
  // e2=Energy_pi1;
  // e3=Energy_pi2;
  // etot=pow(e1+e2+e3,2);
  // invm2=etot-ptot;

  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);
  //  G4cout<<"invm:"<<invm<<G4endl;
  /////shhwang hybrid1
  G4double vtx,vty,vtz;
  G4double rn_vtx,rn_vtz;
  while(1){
    rn_vtx = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    rn_vtz = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    //    G4cout<<rn_vtx<<G4endl;
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot( -m_target_size.z(), m_target_size.z());
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;
  //  double vtz = -157.5+15.0*((double) G4RandFlat::shoot()); //--> 15 mm
  //  double vtz = -160.0+20.0*((double) G4RandFlat::shoot()); //--> 20 mm
  //  double vtz = -162.5+25.0*((double) G4RandFlat::shoot()); //--> 25 mm
  //    double vtz = -165.0+30.0*((double) G4RandFlat::shoot()); //--> 30 mm
  //  double vtz = -200.0+100.0*((double) G4RandFlat::shoot()); //--> test//

  ///////proton PG
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  particleGun->SetParticleEnergy((Energy_pro - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi1 PG
  particleGun->SetParticleDefinition(piMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi2 PG
  particleGun->SetParticleDefinition(piPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  particleGun->SetParticleEnergy((Energy_pi2 - piPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z,proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z,piMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z,piPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hybrid3body_mode3( G4Event* anEvent )
{
  //  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;

  G4double mom[3];
  G4double Ebeam,pg_x,pg_y,pg_z;
  G4double rmpro=0.93827203;
  // G4double rmpi=0.13957018;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* neutron;
  G4ParticleDefinition* piMinus;
  G4ParticleDefinition* piPlus;

  neutron = particleTable->FindParticle("neutron");
  piPlus = particleTable->FindParticle("pi+");
  piMinus = particleTable->FindParticle("pi+");
  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*0.01294);
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  G4double pbeam=0.7;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  Ebeam = sqrt(pow(pbeam,2)+pow(piMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  //  mass_hybrid = 1.22;    // Mass for H
  //  width_hybrid = 0.000;  // Width for H

  //  double temp= kin1.GetMomentumLab(3);
  //    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
  //  Kinema3Body kin2(kin3.M_res, 0.0, m3, m4, m5, temp, 0.0);
 G4double W=sqrt(pow(Ebeam+rmpro,2)-pow(pbeam,2));
 Kinema3Body Hybrid(W,0,
		    neutron->GetPDGMass()/GeV,
		    piMinus->GetPDGMass()/GeV,
		    piPlus->GetPDGMass()/GeV,
		    pbeam, 0.0);
//  Hdibaryon.Dump();
  /* pro */
  Energy_pro = Hybrid.GetEnergy(3);
  // G4double   momentum_pro = Hybrid.GetMomentum(3);
  Hybrid.GetMomentum(3,mom);

  mom_pro_x = mom[0];
  mom_pro_y = mom[1];
  mom_pro_z = mom[2];

  // G4double Thetapro = Hybrid.GetTheta(3);
  // G4double   Phipro = Hybrid.GetPhi(3);

  /* L2 */
  Energy_pi1 = Hybrid.GetEnergy(4);
  // G4double   momentum_pi1 = Hybrid.GetMomentum(4);
  Hybrid.GetMomentum(4,mom);

  mom_pi1_x = mom[0];
  mom_pi1_y = mom[1];
  mom_pi1_z = mom[2];

  // G4double   Thetapi1 = Hybrid.GetTheta(4);
  // G4double   Phipi1 = Hybrid.GetPhi(4);

  /* pi2 */
  Energy_pi2 = Hybrid.GetEnergy(5);
  // G4double   momentum_pi2 = Hybrid.GetMomentum(5);
  Hybrid.GetMomentum(5,mom);

  mom_pi2_x = mom[0];
  mom_pi2_y = mom[1];
  mom_pi2_z = mom[2];

  // G4double   Thetapi2 = Hybrid.GetThetaCM(1);
  // G4double shphi=(atan2(mom_pro_y,mom_pro_x))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4double beta= sqrt(pbeam*pbeam)/(0.938272013+Ebeam);
  /*  G4double test;
  G4double mom1[4];
  G4double mom2[4];
  G4double mom3[4];

  mom1[0]=mom_pro_x;
  mom1[1]=mom_pro_y;
  mom1[2]=mom_pro_z;
  mom1[3]=sqrt(mom_pro_x*mom_pro_x+mom_pro_y*mom_pro_y+mom_pro_z*mom_pro_z+rmpro*rmpro);
  momppro=sqrt(pow(mom_pro_x,2)+pow(mom_pro_y,2)+pow(mom_pro_z,2));
  */
  /*  test=lorentz(momk,beta,momcmk);
  momcmprop=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  cmk->Fill(acos(momcmk[2]/momcmprop)/3.141592654*180);

  coscmk->Fill((momcmk[2]/momcmprop));

  labk->Fill(acos(momk[2]/momprop)*180/3.141592654);
  coslabk->Fill(momk[2]/momprop);

  G4double misspro = miss1(pbeam,rmk,momk);//pbeam, mass, momk

  missk->Fill(misskp);
  */

  // G4double e1=0.,e2=0.,e3=0,etot=0.,invm2=0.,ptot=0.,invm=0.;
  // ptot=pow(mom_pro_x+mom_pi1_x+mom_pi2_x,2)+pow(mom_pro_y+mom_pi1_y+mom_pi2_y,2)+pow(mom_pro_z+mom_pi1_z+mom_pi2_z,2);
  // e1=Energy_pro;
  // e2=Energy_pi1;
  // e3=Energy_pi2;
  // etot=pow(e1+e2+e3,2);
  // invm2=etot-ptot;
  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);
  //  G4cout<<"invm:"<<invm<<G4endl;
  /////shhwang hybrid1
  G4double vtx,vty,vtz;
  G4double rn_vtx,rn_vtz;
  while(1){
    rn_vtx = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    rn_vtz = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    //    G4cout<<rn_vtx<<G4endl;
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot( -m_target_size.z(), m_target_size.z() );
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;
  //  double vtz = -157.5+15.0*((double) G4RandFlat::shoot()); //--> 15 mm
  //  double vtz = -160.0+20.0*((double) G4RandFlat::shoot()); //--> 20 mm
  //  double vtz = -162.5+25.0*((double) G4RandFlat::shoot()); //--> 25 mm
  //    double vtz = -165.0+30.0*((double) G4RandFlat::shoot()); //--> 30 mm
  //  double vtz = -200.0+100.0*((double) G4RandFlat::shoot()); //--> test//

  ///////neutron PG
  particleGun->SetParticleDefinition(neutron);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  particleGun->SetParticleEnergy((Energy_pro - neutron->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi1 PG
  particleGun->SetParticleDefinition(piMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi2 PG
  particleGun->SetParticleDefinition(piPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  particleGun->SetParticleEnergy((Energy_pi2 - piPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z,neutron->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z,piMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z,piPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hybrid3body_mode4( G4Event* anEvent )
{
  //  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;

  G4double mom[3];
  G4double Ebeam,pg_x,pg_y,pg_z;
  G4double rmpro = 0.93827203;
  // G4double rmpi = 0.13957018;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* proton;
  G4ParticleDefinition* piMinus;
  G4ParticleDefinition* piPlus;

  proton = particleTable->FindParticle("proton");
  piPlus = particleTable->FindParticle("pi+");
  piMinus = particleTable->FindParticle("pi0");
  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*0.01294);
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  G4double pbeam=0.7;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  Ebeam = sqrt(pow(pbeam,2)+pow(piMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  //  mass_hybrid = 1.22;    // Mass for H
  //  width_hybrid = 0.000;  // Width for H

  //  double temp= kin1.GetMomentumLab(3);
  //    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
  //  Kinema3Body kin2(kin3.M_res, 0.0, m3, m4, m5, temp, 0.0);
 G4double W=sqrt(pow(Ebeam+rmpro,2)-pow(pbeam,2));
 Kinema3Body Hybrid(W,0,
		    proton->GetPDGMass()/GeV,
		    piMinus->GetPDGMass()/GeV,
		    piPlus->GetPDGMass()/GeV,
		    pbeam, 0.0);
//  Hdibaryon.Dump();
  /* pro */
  Energy_pro = Hybrid.GetEnergy(3);
  // G4double momentum_pro = Hybrid.GetMomentum(3);
  Hybrid.GetMomentum(3,mom);

  mom_pro_x = mom[0];
  mom_pro_y = mom[1];
  mom_pro_z = mom[2];

  // G4double Thetapro = Hybrid.GetTheta(3);
  // G4double   Phipro = Hybrid.GetPhi(3);

  /* L2 */
  Energy_pi1 = Hybrid.GetEnergy(4);
  // G4double   momentum_pi1 = Hybrid.GetMomentum(4);
  Hybrid.GetMomentum(4,mom);

  mom_pi1_x = mom[0];
  mom_pi1_y = mom[1];
  mom_pi1_z = mom[2];

  // G4double   Thetapi1 = Hybrid.GetTheta(4);
  // G4double   Phipi1 = Hybrid.GetPhi(4);

  /* pi2 */
  Energy_pi2 = Hybrid.GetEnergy(5);
  // G4double   momentum_pi2 = Hybrid.GetMomentum(5);
  Hybrid.GetMomentum(5,mom);

  mom_pi2_x = mom[0];
  mom_pi2_y = mom[1];
  mom_pi2_z = mom[2];

  // G4double   Thetapi2 = Hybrid.GetThetaCM(1);
  // G4double shphi=(atan2(mom_pro_y,mom_pro_x))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4double beta= sqrt(pbeam*pbeam)/(0.938272013+Ebeam);
  /*  G4double test;
  G4double mom1[4];
  G4double mom2[4];
  G4double mom3[4];

  mom1[0]=mom_pro_x;
  mom1[1]=mom_pro_y;
  mom1[2]=mom_pro_z;
  mom1[3]=sqrt(mom_pro_x*mom_pro_x+mom_pro_y*mom_pro_y+mom_pro_z*mom_pro_z+rmpro*rmpro);
  momppro=sqrt(pow(mom_pro_x,2)+pow(mom_pro_y,2)+pow(mom_pro_z,2));
  */
  /*  test=lorentz(momk,beta,momcmk);
  momcmprop=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  cmk->Fill(acos(momcmk[2]/momcmprop)/3.141592654*180);

  coscmk->Fill((momcmk[2]/momcmprop));

  labk->Fill(acos(momk[2]/momprop)*180/3.141592654);
  coslabk->Fill(momk[2]/momprop);

  G4double misspro = miss1(pbeam,rmk,momk);//pbeam, mass, momk

  missk->Fill(misskp);
  */

  // G4double e1=0.,e2=0.,e3=0,etot=0.,invm2=0.,ptot=0.,invm=0.;

  // ptot=pow(mom_pro_x+mom_pi1_x+mom_pi2_x,2)+pow(mom_pro_y+mom_pi1_y+mom_pi2_y,2)+pow(mom_pro_z+mom_pi1_z+mom_pi2_z,2);
  // e1=Energy_pro;
  // e2=Energy_pi1;
  // e3=Energy_pi2;
  // etot=pow(e1+e2+e3,2);
  // invm2=etot-ptot;

  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);
  //  G4cout<<"invm:"<<invm<<G4endl;
  /////shhwang hybrid1
  G4double vtx,vty,vtz;
  G4double rn_vtx,rn_vtz;
  while(1){
    rn_vtx = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    rn_vtz = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    //    G4cout<<rn_vtx<<G4endl;
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot( -m_target_size.z(), m_target_size.z() );
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;
  //  double vtz = -157.5+15.0*((double) G4RandFlat::shoot()); //--> 15 mm
  //  double vtz = -160.0+20.0*((double) G4RandFlat::shoot()); //--> 20 mm
  //  double vtz = -162.5+25.0*((double) G4RandFlat::shoot()); //--> 25 mm
  //    double vtz = -165.0+30.0*((double) G4RandFlat::shoot()); //--> 30 mm
  //  double vtz = -200.0+100.0*((double) G4RandFlat::shoot()); //--> test//

  ///////proton PG
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  particleGun->SetParticleEnergy((Energy_pro - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi1 PG
  particleGun->SetParticleDefinition(piMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  particleGun->SetParticleEnergy((Energy_pi1 - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////pi2 PG
  particleGun->SetParticleDefinition(piPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  particleGun->SetParticleEnergy((Energy_pi2 - piPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z,proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z,piMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z,piPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hdibaryon1( G4Event* anEvent )
{
  G4double mass_hdibaryon, width_hdibaryon;

  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  // G4double rmk=0.493677;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonPlus;
  G4ParticleDefinition* kaonMinus;
  G4ParticleDefinition* Lambda1;
  G4ParticleDefinition* Lambda2;
  G4ParticleDefinition* proton;

  proton = particleTable->FindParticle("proton");
  Lambda1 = particleTable->FindParticle("lambda");
  Lambda2 = particleTable->FindParticle("lambda");
  kaonPlus = particleTable->FindParticle("kaon+");
  kaonMinus = particleTable->FindParticle("kaon-");

  // G4double Ebeam = sqrt(1.8*1.8+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);       // Incident gamma energy (GeV)

  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = 1.8;
  // G4double pbeam=1.8;

  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  mass_hdibaryon = gConf.Get<G4double>( "HdibaryonMass" );
  width_hdibaryon = gConf.Get<G4double>( "HdibaryonWidth" );

  Kinema3Resonance Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			     ((proton->GetPDGMass()/GeV)*2),
			     Lambda1->GetPDGMass()/GeV,
			     Lambda2->GetPDGMass()/GeV,
			     kaonPlus->GetPDGMass()/GeV,
			     mass_hdibaryon, width_hdibaryon, pg_z, 0.0);

//  Hdibaryon.Dump();

  Energy_kp = Hdibaryon.GetEnergy(5);
  // G4double   momentum_kp = Hdibaryon.GetMomentum(5);
  Hdibaryon.GetMomentum(5,mom);

  /* L1 */
  Energy_L1 = Hdibaryon.GetEnergy(3);
  // G4double   momentum_L1 = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];

  // G4double ThetaL1 = Hdibaryon.GetTheta(3);
  // G4double PhiL1 = Hdibaryon.GetPhi(3);

  /* L2 */
  Energy_L2 = Hdibaryon.GetEnergy(4);
  // G4double   momentum_L2 = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double   ThetaL2 = Hdibaryon.GetTheta(4);
  // G4double   PhiL2 = Hdibaryon.GetPhi(4);

  /* Kp */
  Energy_kp = Hdibaryon.GetEnergy(5);
  // G4double   momentum_kp = Hdibaryon.GetMomentum(5);
  Hdibaryon.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  // G4double Thetakp = Hdibaryon.GetThetaCM(1);
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  // G4double beta= sqrt(pbeam*pbeam)/(2.*0.938272013+Ebeam);
  // G4double momk[4]={};
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  // G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  // coscmk->Fill((momcmk[2]/momcmkpp));

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

  // G4double misskp = miss1(pbeam,rmk,momk);//pbeam, mass, momk
  //missk->Fill(misskp);

  //cal for invariant mass
  // G4double ptot=pow(mom_L1_x+mom_L2_x,2)+pow(mom_L1_y+mom_L2_y,2)+pow(mom_L1_z+mom_L2_z,2);
  // G4double e1=Energy_L1;
  // G4double e2=Energy_L2;
  // G4double etot=pow(e1+e2,2);
  // G4double invm2=etot-ptot;
  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);
  double vtx = 0.*mm;
  double vty = 0.*mm;
  /////shhwang hdibaryon1


  //  G4cout<<"env_target_width :"<<atof(env_env_target_width.c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);
  G4double vtz= G4RandFlat::shoot( env_target_pos_z - m_target_size.z()/2,
					env_target_pos_z+m_target_size.z()/2)*mm;

  ///////kp PG
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////L1 PG
  particleGun->SetParticleDefinition(Lambda1);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  particleGun->SetParticleEnergy((Energy_L1 - Lambda1->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////L2 PG
  particleGun->SetParticleDefinition(Lambda2);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  particleGun->SetParticleEnergy((Energy_L2 - Lambda2->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,Lambda1->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,Lambda2->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_hybridPHSG( G4Event* anEvent )
{
  //  G4double mom[3];
  G4double Ebeam;
  //  G4double rmk=0.493677;
  auto particleTable = G4ParticleTable::GetParticleTable();
  // auto piMinus = particleTable->FindParticle("pi-");
  // auto proton = particleTable->FindParticle("proton");
  auto hybridbaryon = particleTable->FindParticle("hybridb");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  G4double pimom=1.;
  gAnaMan.SetPrimaryBeam(0,0,pimom);
  //  Ebeam = sqrt(pimom*pimom+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pimom*pimom+0.13957018*0.13957018);

  // G4double pbm[4];
  // pbm[0]=0;
  // pbm[1]=0;
  // pbm[2]=pimom;
  // pbm[3]=Ebeam;

  // G4double beta= pimom/(0.93827203+pbm[3]);
  // G4double momhycm[4]={0};
  G4double rmp=0.93827203;
  G4double rmpi=0.13957018;
  G4double W=sqrt(pow(Ebeam+rmp,2)-pow(pimom,2));
  G4double W2=sqrt(pow(rmp,2)+pow(rmpi,2)+2*Ebeam*rmp);
  // momhycm[3]=hybridbaryon->GetPDGMass()/GeV;
  // momhycm[3]=W2;
  G4double momhylab[4]={0};
  // G4double test=lorentcmlab(momhycm,beta,momhylab);
  G4cout<<"momhylab x: "<<momhylab[0]<<G4endl;
  G4cout<<"momhylab y: "<<momhylab[1]<<G4endl;
  G4cout<<"momhylab z: "<<momhylab[2]<<G4endl;
  G4cout<<"momhylab energy"<<momhylab[3]<<G4endl;
  ///////gun
  G4cout<<"check energy"<<W<<G4endl;
  G4cout<<"w"<<W<<G4endl;
  G4cout<<"w2"<<W2<<G4endl;
  particleGun->SetParticleDefinition(hybridbaryon);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //  particleGun->SetParticleEnergy((Ebeam+proton->GetPDGMass()/GeV - hybridbaryon->GetPDGMass()/GeV)*GeV);
  // particleGun->SetParticleEnergy((momhylab[3] - hybridbaryon->GetPDGMass()/GeV)*GeV);
  // particleGun->SetParticleEnergy((W - hybridbaryon->GetPDGMass()/GeV)*GeV);
  // particleGun->SetParticleEnergy((momhylab[3])*GeV);
  // particleGun->SetParticleEnergy(0.5*GeV);
  //  particleGun->SetParticleEnergy((W - hybridbaryon->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticleMomentum(1.*GeV); //KE ~ 0.35746632293688GeV

  //  G4cout<<"sh test energy:"<<(Ebeam+proton->GetPDGMass()/GeV - hybridbaryon->GetPDGMass()/GeV)<<G4endl;
  //  G4cout<<"sh test energy:"<<momhylab[3]- hybridbaryon->GetPDGMass()/GeV<<G4endl;

  //  particleGun->SetParticleEnergy((W - hybridbaryon->GetPDGMass()/GeV)*GeV);

  particleGun->SetParticlePosition(G4ThreeVector(0.,0.,-150.*mm));
  particleGun->GeneratePrimaryVertex(anEvent);
  //  G4cout<<hybridbaryon->GetPDGMass()/GeV  <<G4endl;

  /*
  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,Lambda1->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,Lambda2->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
  */
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_test( G4Event* anEvent )
{
  //  G4double mass_hdibaryon, width_hdibaryon;
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //  G4double mom[3];
  //  G4double  pbeam;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* proton;
  G4ParticleDefinition* piminus;

  proton = particleTable->FindParticle("proton");
  piminus = particleTable->FindParticle("pi-");

  /* proton */
  G4double mompx=0.1;
  G4double mompz=0.0;
  G4double rand= G4RandFlat::shoot(-1.,1.)*3.141592654;
  mom_p_x = cos(rand)*mompx-sin(rand)*mompz;
  mom_p_y = 0.;
  mom_p_z = sin(rand)*mompx+cos(rand)*mompz;
  G4double mom_p=sqrt(pow(mom_p_x,2)+pow(mom_p_y,2)+pow(mom_p_z,2));
  //  Energy_p=pow(mom_p,2)+proton->GetPDGMass()/GeV;
  Energy_p=pow(mom_p,2)+piminus->GetPDGMass()/GeV;


  double vtz = -150+0.0*((double) G4RandFlat::shoot()); //--> 0 mm
  double vtx = 0.; //--> 0 mm
  double vty = -200.; //--> 0 mm

  ///////kp PG
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  //  particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticleEnergy((Energy_p - piminus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_dedx_single( G4Event* anEvent )
{
  //  G4double mass_hdibaryon, width_hdibaryon;
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //  G4double mom[3];
  //  G4double Ebeam,pg_x,pg_y,pg_z, pbeam;

  auto particleTable = G4ParticleTable::GetParticleTable();
  auto proton = particleTable->FindParticle("proton");
  // auto piplus = particleTable->FindParticle("pi+");
  // auto piminus = particleTable->FindParticle("pi-");

  /* proton */
  G4double phi=G4RandFlat::shoot(-1.,1.)*3.141592654;
  G4double mom_p=G4RandFlat::shoot(0.8,2.0);
  G4double theta=acos(G4RandFlat::shoot(-1.,1.));
  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);

  G4double rn_vtx=-1;  G4double rn_vtz=-1;
  G4double vtx=0;  G4double vty=0;   G4double vtz=0;
  if( gConf.Get<G4int>("Experiment") == 45. ){
    while(1){
      rn_vtx = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
      rn_vtz = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
      if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
    }
    vty = G4RandFlat::shoot(-m_target_size.z(),m_target_size.z());
    vtx=rn_vtx;
    vtz=rn_vtz+env_target_pos_z;
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-15.,15.)*mm;
    vty = G4RandFlat::shoot(-5.,5.)*mm;
    vtz = G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;
  }

  ///////kp PG
  Energy_p=sqrt(pow(mom_p,2)+pow(proton->GetPDGMass()/GeV,2));
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  /*
  Energy_p=sqrt(pow(mom_p,2)+pow(piplus->GetPDGMass()/GeV,2));
  particleGun->SetParticleDefinition(piplus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  particleGun->SetParticleEnergy((Energy_p - piplus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  */
  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  //  gAnaMan.SetPrimaryParticle(1,mom_p_x,mom_p_y,mom_p_z,piminus->GetPDGMass()/GeV);
  //  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_all( G4Event* anEvent )
{
  //  G4double mass_hdibaryon, width_hdibaryon;
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //  G4double mom[3];
  //  G4double Ebeam,pg_x,pg_y,pg_z, pbeam;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* proton;
  G4ParticleDefinition* piplus;
  G4ParticleDefinition* piminus;
  G4ParticleDefinition* kaonplus;
  G4ParticleDefinition* kaonminus;

  proton = particleTable->FindParticle("proton");
  piplus = particleTable->FindParticle("pi+");
  piminus = particleTable->FindParticle("pi-");
  kaonplus = particleTable->FindParticle("kaon+");
  kaonminus = particleTable->FindParticle("kaon-");

  //angle

  G4double phi=G4RandFlat::shoot(-1.,1.)*3.141592654;
  G4double mom_p=G4RandFlat::shoot(0.05,2.0);
  G4double theta=acos(G4RandFlat::shoot(0.,1.));
  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);
  Energy_p=sqrt(pow(mom_p,2)+pow(proton->GetPDGMass()/GeV,2))*GeV;

  G4double rn_vtx=-1;  G4double rn_vtz=-1;
  G4double vtx=0;  G4double vty=0;   G4double vtz=0;
  if( gConf.Get<G4int>("Experiment") == 45. ){
    while(1){
      rn_vtx = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
      rn_vtz = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
      if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
    }
    vty = G4RandFlat::shoot(-m_target_size.z(),m_target_size.z());
    vtx=rn_vtx;
    vtz=rn_vtz+env_target_pos_z;
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-15.,15.)*mm;
    vty = G4RandFlat::shoot(-5.,5.)*mm;
    vtz = G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;
  }


  /*
  ///////kp PG
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  */

  G4double ratio=G4RandFlat::shoot(0.,1.);
  //  ratio = 0.1;

  if(ratio >= 0.0 && ratio < 0.2 ){
  // proton//
    Energy_p=sqrt(pow(mom_p,2)+pow(proton->GetPDGMass()/GeV,2));
    particleGun->SetParticleDefinition(proton);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    particleGun->GeneratePrimaryVertex(anEvent);
    gAnaMan.SetModeID(1);
    gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,proton->GetPDGMass()/GeV);
    gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  } else  if(ratio >= 0.2 && ratio < 0.4 ){
    Energy_p=sqrt(pow(mom_p,2)+pow(kaonplus->GetPDGMass()/GeV,2));
    particleGun->SetParticleDefinition(kaonplus);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    particleGun->SetParticleEnergy((Energy_p - kaonplus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    particleGun->GeneratePrimaryVertex(anEvent);
    gAnaMan.SetModeID(2);
    gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,kaonplus->GetPDGMass()/GeV);
    gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  } else  if(ratio >= 0.4 && ratio< 0.6 ){
  // pi+ //
    Energy_p=sqrt(pow(mom_p,2)+ pow(piplus->GetPDGMass()/GeV,2));
    particleGun->SetParticleDefinition(piplus);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    particleGun->SetParticleEnergy((Energy_p - piplus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    particleGun->GeneratePrimaryVertex(anEvent);
    gAnaMan.SetModeID(3);
    gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,piplus->GetPDGMass()/GeV);
    gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  } else  if(ratio >= 0.6 && ratio< 0.8 ){
  // pi- //
    Energy_p=sqrt(pow(mom_p,2)+pow(piminus->GetPDGMass()/GeV,2));
    particleGun->SetParticleDefinition(piminus);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    particleGun->SetParticleEnergy((Energy_p - piminus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    particleGun->GeneratePrimaryVertex(anEvent);
    gAnaMan.SetModeID(4);
    gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,piminus->GetPDGMass()/GeV);
    gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  } else  if(ratio >= 0.8 && ratio <= 1.0 ){
  // k- //
    Energy_p=sqrt(pow(mom_p,2)+pow(kaonminus->GetPDGMass()/GeV,2));
    particleGun->SetParticleDefinition(kaonminus);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    particleGun->SetParticleEnergy((Energy_p - kaonminus->GetPDGMass()/GeV)*GeV);
    particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    particleGun->GeneratePrimaryVertex(anEvent);
    gAnaMan.SetModeID(5);
    gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,kaonminus->GetPDGMass()/GeV);
    gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  }

  ///////kp PG
}

//_____________________________________________________________________________
// generator 60
void
TPCPrimaryGeneratorAction::Generate_Lambda1405_rad1( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  G4double pbm[4];
  G4double Energy_L, mom_L_x, mom_L_y, mom_L_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Lambda1405;
  G4ParticleDefinition* proton;
  G4ParticleDefinition* k0;
  G4ParticleDefinition* piMinus;

  piMinus = particleTable->FindParticle("pi-");
  proton = particleTable->FindParticle("proton");
  k0 = particleTable->FindParticle("kaon0S");
  Lambda1405 = particleTable->FindParticle("lambda(1405)");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  pbeam=1.65;
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;

  KinemaKstar Hdibaryon(piMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			Lambda1405->GetPDGMass()/GeV,
			k0->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hdibaryon.GetEnergy(4);
  // G4double momentum_ksp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_L = Hdibaryon.GetEnergy(3);
  // G4double momentum_L = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_L_x = mom[0];
  mom_L_y = mom[1];
  mom_L_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_L = Hdibaryon.GetTheta(3);
  // G4double Phi_L = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+pow(k0->GetPDGMass()/GeV,2));

  //  G4cout<<"momks:"<<momks[0]<<":"<<momks[1]<<":"<<momks[2]<<":"<<momks[3]<<G4endl;
  //  G4cout<<"Energy_ksp:  "<<Energy_ksp<<G4endl;
  //  G4cout<<"mom_ks:"<<mom_ksp_x<<":"<<mom_ksp_y<<":"<<mom_ksp_z<<G4endl;
  // G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);
  //  G4double test=lorentz(momks,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;
  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_L<<G4endl;
  //coscmk->Fill(momcmk[2]/momcmkpp);


  // G4double misskp = missks(pbeam,piMinus->GetPDGMass()/GeV,k0->GetPDGMass()/GeV,momks);
  //  G4double misskp = missks(pbeam,kstar0->GetPDGMass()/GeV,momks);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz = 0*mm;
  G4double rn_vtx,rn_vtz;

  while(1){
    rn_vtx = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
    rn_vtz = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
    //    G4cout<<rn_vtx<<G4endl;
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot(-m_target_size.z(),m_target_size.z());
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;


  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);


  //Kstar0
  particleGun->SetParticleDefinition(k0);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  particleGun->SetParticleEnergy((Energy_ksp - k0->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //L
  particleGun->SetParticleDefinition(Lambda1405);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L_x,mom_L_y,mom_L_z));
  particleGun->SetParticleEnergy((Energy_L - Lambda1405->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,k0->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L_x,mom_L_y,mom_L_z,Lambda1405->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
// generator 61
void
TPCPrimaryGeneratorAction::Generate_Lambda1405_rad2( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  G4double pbm[4];
  G4double Energy_L, mom_L_x, mom_L_y, mom_L_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Lambda1405;
  G4ParticleDefinition* proton;
  G4ParticleDefinition* k0;
  G4ParticleDefinition* piMinus;

  piMinus = particleTable->FindParticle("pi-");
  proton = particleTable->FindParticle("proton");
  k0 = particleTable->FindParticle("kaon0S");
  Lambda1405 = particleTable->FindParticle("lambda1405r");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  pbeam=1.65;
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;

  KinemaKstar Hdibaryon(piMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			Lambda1405->GetPDGMass()/GeV,
			k0->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hdibaryon.GetEnergy(4);
  // G4double momentum_ksp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_L = Hdibaryon.GetEnergy(3);
  // G4double momentum_L = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_L_x = mom[0];
  mom_L_y = mom[1];
  mom_L_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_L = Hdibaryon.GetTheta(3);
  // G4double Phi_L = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+pow(k0->GetPDGMass()/GeV,2));

  //  G4cout<<"momks:"<<momks[0]<<":"<<momks[1]<<":"<<momks[2]<<":"<<momks[3]<<G4endl;
  //  G4cout<<"Energy_ksp:  "<<Energy_ksp<<G4endl;
  //  G4cout<<"mom_ks:"<<mom_ksp_x<<":"<<mom_ksp_y<<":"<<mom_ksp_z<<G4endl;
  // G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);
  //  G4double test=lorentz(momks,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;
  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_L<<G4endl;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  // G4double misskp = missks(pbeam,piMinus->GetPDGMass()/GeV,k0->GetPDGMass()/GeV,momks);
  //  G4double misskp = missks(pbeam,kstar0->GetPDGMass()/GeV,momks);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz = 0*mm;
  G4double rn_vtx,rn_vtz;

  while(1){
    rn_vtx = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
    rn_vtz = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
    //    G4cout<<rn_vtx<<G4endl;
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot(-m_target_size.z(),m_target_size.z());
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;


  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);


  //Kstar0
  particleGun->SetParticleDefinition(k0);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  particleGun->SetParticleEnergy((Energy_ksp - k0->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //L
  particleGun->SetParticleDefinition(Lambda1405);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L_x,mom_L_y,mom_L_z));
  particleGun->SetParticleEnergy((Energy_L - Lambda1405->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,k0->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L_x,mom_L_y,mom_L_z,Lambda1405->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_Lambda1405_reso( G4Event* anEvent )
{
  G4double mass_hdibaryon, width_hdibaryon;
  //  G4double Energy_h, momentum_h, mom_h_x, mom_h_y, mom_h_z;

  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  G4double mom[3];
  // G4double pg_x, pg_y, pg_z;
  G4double pbeam;
  // G4double rmk=0.493677;

  //  G4cout<<"check111"<<G4endl;
  auto particleTable = G4ParticleTable::GetParticleTable();

  // G4ParticleDefinition* sigmaPlus;
  // G4ParticleDefinition* sigmaMinus;
  // G4ParticleDefinition* sigmaZero;
  // G4ParticleDefinition* piPlus;
  // G4ParticleDefinition* piZero;
  G4ParticleDefinition* piMinus;
  G4ParticleDefinition* proton;
  G4ParticleDefinition* kaonZero;
  //  G4cout<<"check222"<<G4endl;
  G4ParticleDefinition* sigma = nullptr;
  G4ParticleDefinition* pion = nullptr;
  //  G4cout<<"check333"<<G4endl;

  //  G4cout<<"check444"<<G4endl;

  mass_hdibaryon = gConf.Get<G4double>( "HdibaryonMass" );
  width_hdibaryon = gConf.Get<G4double>( "HdibaryonWidth" );

  G4double check_decay_mode=G4RandFlat::shoot();
  G4int mode=0.;

  //  G4cout<<"check 1"<<G4endl;
  if(check_decay_mode<=1./3.){
    mode=1.; //S+pi-
    sigma = particleTable->FindParticle("sigma+");
    pion = particleTable->FindParticle("pi-");
  }else if(check_decay_mode>1./3. && check_decay_mode<=2./3.){
    mode=2.;     //S0pi0
    sigma = particleTable->FindParticle("sigma0");
    pion = particleTable->FindParticle("pi0");
  }else if(check_decay_mode>2./3. && check_decay_mode<=1.){
    mode=3.;     //S-pi+
    sigma = particleTable->FindParticle("sigma-");
    pion = particleTable->FindParticle("pi+");
  }else{
    G4cout<<"Break Lambda decay mode"<<G4endl;
  }
  //  G4cout<<"check 2"<<G4endl;
  gAnaMan.SetModeID(mode);

  proton = particleTable->FindParticle("proton");
  //  sigmaPlus = particleTable->FindParticle("sigma+");
  //  sigmaMinus = particleTable->FindParticle("sigma-");
  //  sigmaZero = particleTable->FindParticle("sigma0");

  // auto piPlus = particleTable->FindParticle("pi+");
  piMinus = particleTable->FindParticle("pi-");
  //  piZero = particleTable->FindParticle("pi0");

  kaonZero = particleTable->FindParticle("kaon0S");

  //  Ebeam = 1.9;       // Incident gamma energy (GeV)
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  // G4double Ebeam = sqrt(pow(pbeam,2)+pow(piMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(0,0,pbeam);// how about add mass information?

  //  G4cout<< piMinus->GetPDGMass()/GeV <<G4endl;
  //  G4cout<< proton->GetPDGMass()/GeV <<G4endl;
  //  G4cout<< pion->GetPDGMass()/GeV <<G4endl;
  //  G4cout<< sigma->GetPDGMass()/GeV <<G4endl;
  //  G4cout<< kaonZero->GetPDGMass()/GeV <<G4endl;
  //  G4cout<< mass_hdibaryon <<G4endl;
  //  G4cout<< width_hdibaryon <<G4endl;
  //  G4cout<< pbeam <<G4endl;
  Kinema3Resonance  Hdibaryon(piMinus->GetPDGMass()/GeV,
			     proton->GetPDGMass()/GeV,
			     pion->GetPDGMass()/GeV,
			     sigma->GetPDGMass()/GeV,
			     kaonZero->GetPDGMass()/GeV,
			     mass_hdibaryon, width_hdibaryon, pbeam, 0.0);

  //  G4cout<<"check 3"<<G4endl;
  Energy_kp = Hdibaryon.GetEnergy(5);
  // G4double momentum_kp = Hdibaryon.GetMomentum(5);
  Hdibaryon.GetMomentum(5,mom);

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;

  /* L1 */
  Energy_L1 = Hdibaryon.GetEnergy(3);
  // G4double momentum_L1 = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];

  // G4double ThetaL1 = Hdibaryon.GetTheta(3);
  // G4double PhiL1 = Hdibaryon.GetPhi(3);

  /* L2 */
  Energy_L2 = Hdibaryon.GetEnergy(4);
  // G4double momentum_L2 = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double ThetaL2 = Hdibaryon.GetTheta(4);
  // G4double PhiL2 = Hdibaryon.GetPhi(4);

  /* Kp */
  Energy_kp = Hdibaryon.GetEnergy(5);
  // G4double momentum_kp = Hdibaryon.GetMomentum(5);
  Hdibaryon.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  //  Thetakp = Hdibaryon.GetThetaCM(5);
  // G4double Thetakp = Hdibaryon.GetThetaCM(1);
  // G4double Phikp = Hdibaryon.GetPhi(5);
  //  G4cout<<Phikp<<G4endl;
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4cout<<"phi test: "<<Phikp<<":"<<shphi<<G4endl;
  //  G4cout<<"test: "<<Hdibaryon.GetPhiCM(1)<<":"<<Hdibaryon.GetPhi(5)<<G4endl;
  // Vertex (Reaction) point
  //  G4cout<<Energy_kp<<G4endl;
  //  //  G4double beta= sqrt(1.8*1.8-0.493677*0.493677)/(2*0.93827203+Energy_kp);
  //  G4cout<<(proton->GetPDGMass()/GeV)<<G4endl;

  // G4double beta= pbeam/(proton->GetPDGMass()/GeV+Ebeam);
  //  G4cout<<"Pri:"<<beta<<G4endl;
  //  G4double beta= sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2))/(2*0.93827203+1.8);

  // G4double momk[4]={};
  // G4double momcmk[4]={};
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+pow(piMinus->GetPDGMass()/GeV,2));
  // G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180.);


  // G4double cmphik=atan2(momk[1],momk[0])*180./3.141592654;


  // coscmk->Fill((momcmk[2]/momcmkpp));

  // labk->Fill(acos(momk[2]/momkpp)*180./3.141592654);
  // coslabk->Fill(momk[2]/momkpp);


  ////check angle in the CM frame
  G4double hlab[4];
  hlab[0]=mom_L1_x+mom_L2_x;
  hlab[1]=mom_L1_y+mom_L2_y;
  hlab[2]=mom_L1_z+mom_L2_z;

  hlab[3]=sqrt(hlab[0]*hlab[0]+hlab[1]*hlab[1]+hlab[2]*hlab[2]+mass_hdibaryon*mass_hdibaryon);
  // G4double hcm[4];

  // G4double test=lorentz(hlab,beta,hcm);
  // G4double hcmpp=sqrt(pow(hcm[0],2)+pow(hcm[1],2)+pow(hcm[2],2));
  // cmh->Fill(acos(hcm[2]/hcmpp)/3.141592654*180);

  // coscmh->Fill(hcm[2]/hcmpp);

  //  G4double checkphi;
  // G4double hphi=atan2(hcm[1],hcm[0])/3.141592654*180;

  // phih->Fill(hphi);
  // thetadiff->Fill(acos(hcm[2]/hcmpp)/3.141592654*180+acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //phidiff->Fill(hphi-cmphik);

  // G4double momsum[3]={};
  // momsum[0]=mom_kp_x+mom_L1_x+mom_L2_x;
  // momsum[1]=mom_kp_y+mom_L1_y+mom_L2_y;
  // momsum[2]=mom_kp_z+mom_L1_z+mom_L2_z;

  // G4double momH[3]={};
  // momH[0]=mom_L1_x+mom_L2_x;
  // momH[1]=mom_L1_y+mom_L2_y;
  // momH[2]=mom_L1_z+mom_L2_z;

  // G4double misskp = miss1(pbeam,rmk,momk);//pbeam, mass, momk

  //missk->Fill(misskp);

  // G4double e1=Energy_L1;
  // G4double e2=Energy_L2;
  // G4double invm2 = -9999.;
  // G4double invm = -9999.;
  // G4double ptot=pow(mom_L1_x+mom_L2_x,2)+pow(mom_L1_y+mom_L2_y,2)+pow(mom_L1_z+mom_L2_z,2);
  // G4double etot=pow(e1+e2,2);
  // invm2=etot-ptot;
  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);

  //  G4cout<<"IM(LL)"<<invm<<G4endl;
  //  G4cout<<"-----------end-------------"<<G4endl;
  //  G4double vtr = 20.0*mm*((double) G4RandFlat::shoot());
  //  G4double vtphi = 2.0*pi*((double) G4RandFlat::shoot());
  //  G4double vtx = vtr*cos(vtphi);
  //  G4double vty = vtr*sin(vtphi);

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  ///////kp PG
  particleGun->SetParticleDefinition(kaonZero);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonZero->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////L1 PG
  particleGun->SetParticleDefinition(pion);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  particleGun->SetParticleEnergy((Energy_L1 - pion->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////L2 PG
  particleGun->SetParticleDefinition(sigma);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  particleGun->SetParticleEnergy((Energy_L2 - sigma->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonZero->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,sigma->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,pion->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}

//_____________________________________________________________________________
// generator 63
void
TPCPrimaryGeneratorAction::Generate_Sigma1385_rad( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  G4double pbm[4];
  G4double Energy_L, mom_L_x, mom_L_y, mom_L_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Sigma1385;
  G4ParticleDefinition* proton;
  G4ParticleDefinition* k0;
  G4ParticleDefinition* piMinus;

  piMinus = particleTable->FindParticle("pi-");
  proton = particleTable->FindParticle("proton");
  k0 = particleTable->FindParticle("kaon0S");
  Sigma1385 = particleTable->FindParticle("sigma1385r");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  pbeam=1.65;
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;

  // up:
  KinemaKstar Hdibaryon(piMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			Sigma1385->GetPDGMass()/GeV,
			k0->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hdibaryon.GetEnergy(4);
  // G4double momentum_ksp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_L = Hdibaryon.GetEnergy(3);
  // G4double momentum_L = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_L_x = mom[0];
  mom_L_y = mom[1];
  mom_L_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_L = Hdibaryon.GetTheta(3);
  // G4double Phi_L = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+pow(k0->GetPDGMass()/GeV,2));

  //  G4cout<<"momks:"<<momks[0]<<":"<<momks[1]<<":"<<momks[2]<<":"<<momks[3]<<G4endl;
  //  G4cout<<"Energy_ksp:  "<<Energy_ksp<<G4endl;
  //  G4cout<<"mom_ks:"<<mom_ksp_x<<":"<<mom_ksp_y<<":"<<mom_ksp_z<<G4endl;
  // G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);
  //  G4double test=lorentz(momks,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;
  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_L<<G4endl;
  //coscmk->Fill(momcmk[2]/momcmkpp);


  // G4double misskp = missks(pbeam,piMinus->GetPDGMass()/GeV,k0->GetPDGMass()/GeV,momks);
  //  G4double misskp = missks(pbeam,kstar0->GetPDGMass()/GeV,momks);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz = 0*mm;
  G4double rn_vtx,rn_vtz;

  while(1){
    rn_vtx = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
    rn_vtz = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
    //    G4cout<<rn_vtx<<G4endl;
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot(-m_target_size.z(),m_target_size.z());
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;


  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);


  //Kstar0
  particleGun->SetParticleDefinition(k0);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  particleGun->SetParticleEnergy((Energy_ksp - k0->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //L
  particleGun->SetParticleDefinition(Sigma1385);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L_x,mom_L_y,mom_L_z));
  particleGun->SetParticleEnergy((Energy_L - Sigma1385->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,k0->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L_x,mom_L_y,mom_L_z,Sigma1385->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
// generator 62
void
TPCPrimaryGeneratorAction::Generate_Sigma1385( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  G4double pbm[4];
  G4double Energy_L, mom_L_x, mom_L_y, mom_L_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Sigma1385;
  G4ParticleDefinition* proton;
  G4ParticleDefinition* k0;
  G4ParticleDefinition* piMinus;

  piMinus = particleTable->FindParticle("pi-");
  proton = particleTable->FindParticle("proton");
  k0 = particleTable->FindParticle("kaon0S");
  Sigma1385 = particleTable->FindParticle("sigma(1385)0");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  pbeam=1.65;
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"check11"<<G4endl;
  KinemaKstar Hdibaryon(piMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			Sigma1385->GetPDGMass()/GeV,
			k0->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hdibaryon.GetEnergy(4);
  // G4double momentum_ksp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_L = Hdibaryon.GetEnergy(3);
  // G4double momentum_L = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_L_x = mom[0];
  mom_L_y = mom[1];
  mom_L_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_L = Hdibaryon.GetTheta(3);
  // G4double Phi_L = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4doulbe momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+pow(k0->GetPDGMass()/GeV,2));

  //  G4cout<<"momks:"<<momks[0]<<":"<<momks[1]<<":"<<momks[2]<<":"<<momks[3]<<G4endl;
  //  G4cout<<"Energy_ksp:  "<<Energy_ksp<<G4endl;
  //  G4cout<<"mom_ks:"<<mom_ksp_x<<":"<<mom_ksp_y<<":"<<mom_ksp_z<<G4endl;
  // G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);
  //  G4double test=lorentz(momks,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;
  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_L<<G4endl;
  //coscmk->Fill(momcmk[2]/momcmkpp);


  // G4double misskp = missks(pbeam,piMinus->GetPDGMass()/GeV,k0->GetPDGMass()/GeV,momks);
  //  G4double misskp = missks(pbeam,kstar0->GetPDGMass()/GeV,momks);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz = 0*mm;
  G4double rn_vtx,rn_vtz;

  while(1){
    rn_vtx = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
    rn_vtz = G4RandFlat::shoot(-m_target_size.x(),m_target_size.x());
    //    G4cout<<rn_vtx<<G4endl;
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot(-m_target_size.z(),m_target_size.z());
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;


  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);
  //  G4cout<<"check"<<G4endl;

  //Kstar0
  particleGun->SetParticleDefinition(k0);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  particleGun->SetParticleEnergy((Energy_ksp - k0->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //L
  particleGun->SetParticleDefinition(Sigma1385);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L_x,mom_L_y,mom_L_z));
  particleGun->SetParticleEnergy((Energy_L - Sigma1385->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,k0->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L_x,mom_L_y,mom_L_z,Sigma1385->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
// generator 11
void
TPCPrimaryGeneratorAction::Generate_pip_KsL( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  //  G4double rmks=0.896;
  G4double pbm[4];
  G4double Energy_L, mom_L_x, mom_L_y, mom_L_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Lambda;
  G4ParticleDefinition* kstar0;
  G4ParticleDefinition* piMinus;
  G4ParticleDefinition* proton;

  kstar0 = particleTable->FindParticle("k_star0");
  //  kstar0 = particleTable->FindParticle(313);
  piMinus = particleTable->FindParticle("pi-");
  proton = particleTable->FindParticle("proton");
  Lambda = particleTable->FindParticle("lambda");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"piMinus mass:"<<piMinus->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"KstarMinus mass:"<<kstar0->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"Lambda mass:"<<Lambda->GetPDGMass()/GeV<<G4endl;
  KinemaKstar Hdibaryon(piMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			Lambda->GetPDGMass()/GeV,
			kstar0->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hdibaryon.GetEnergy(4);
  // G4double momentum_ksp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_L = Hdibaryon.GetEnergy(3);
  // G4double momentum_L = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_L_x = mom[0];
  mom_L_y = mom[1];
  mom_L_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_L = Hdibaryon.GetTheta(3);
  // G4double Phi_L = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+0.896*0.896);

  //  G4cout<<"momks:"<<momks[0]<<":"<<momks[1]<<":"<<momks[2]<<":"<<momks[3]<<G4endl;
  //  G4cout<<"Energy_ksp:  "<<Energy_ksp<<G4endl;
  //  G4cout<<"mom_ks:"<<mom_ksp_x<<":"<<mom_ksp_y<<":"<<mom_ksp_z<<G4endl;
  //  G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);

  //  G4double test=lorentz(momks,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;
  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_L<<G4endl;
  //coscmk->Fill(momcmk[2]/momcmkpp);


  //  G4double misskp = missks(pbeam,piMinus->GetPDGMass()/GeV,kstar0->GetPDGMass()/GeV,momks);
  //  G4double misskp = missks(pbeam,kstar0->GetPDGMass()/GeV,momks);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  //Kstar0
  particleGun->SetParticleDefinition(kstar0);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  particleGun->SetParticleEnergy((Energy_ksp - kstar0->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //L
  particleGun->SetParticleDefinition(Lambda);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L_x,mom_L_y,mom_L_z));
  particleGun->SetParticleEnergy((Energy_L - Lambda->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,kstar0->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L_x,mom_L_y,mom_L_z,Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
// generator 12
void
TPCPrimaryGeneratorAction::Generate_pip_KsS( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  //  G4double rmks=0.896;
  G4double pbm[4];
  G4double Energy_S, mom_S_x, mom_S_y, mom_S_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Sigma0;
  G4ParticleDefinition* kstar0;
  G4ParticleDefinition* piMinus;
  G4ParticleDefinition* proton;

  //  kstar0 = particleTable->FindParticle("k_star0");
  kstar0 = particleTable->FindParticle(313);
  piMinus = particleTable->FindParticle("pi-");
  proton = particleTable->FindParticle("proton");
  Sigma0 = particleTable->FindParticle("sigma0");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.9;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  KinemaKstar Hdibaryon(piMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			Sigma0->GetPDGMass()/GeV,
			kstar0->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hdibaryon.GetEnergy(4);
  // G4double momentum_ksp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_S = Hdibaryon.GetEnergy(3);
  // G4double momentum_S = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_S_x = mom[0];
  mom_S_y = mom[1];
  mom_S_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_S = Hdibaryon.GetTheta(3);
  // G4double Phi_S = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+0.896*0.896);
  //  G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);

  //  G4double test=lorentz(momks,beta,momcmk);
  //  G4cout<<"test:"<<test<<G4endl;
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);

  //  G4double misskp = missks(pbeam,piMinus->GetPDGMass()/GeV,kstar0->GetPDGMass()/GeV,momks);
  //  G4double misskp = missks(pbeam,kstar0->GetPDGMass()/GeV,momks);
  //  G4cout<<"missing mass : "<<misskp<<G4endl;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  //Kstar0
  particleGun->SetParticleDefinition(kstar0);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  particleGun->SetParticleEnergy((Energy_ksp - kstar0->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //S
  particleGun->SetParticleDefinition(Sigma0);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_S_x,mom_S_y,mom_S_z));
  particleGun->SetParticleEnergy((Energy_S - Sigma0->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,kstar0->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_S_x,mom_S_y,mom_S_z,Sigma0->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_pip_KstarL( G4Event* anEvent )
{
  G4double mass_hdibaryon, width_hdibaryon;
  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  //  G4double rmk=0.493677;

  auto particleTable = G4ParticleTable::GetParticleTable();
  auto proton = particleTable->FindParticle("proton");
  auto Lambda = particleTable->FindParticle("lambda");
  auto kaonPlus = particleTable->FindParticle("kaon+");
  // auto kaonMinus = particleTable->FindParticle("kaon-");
  auto piMinus = particleTable->FindParticle("pi-");
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = 1.8;
  // G4double pbeam = 1.8;
  // G4double Ebeam = sqrt(pbeam*pbeam+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);       // Incident gamma energy (GeV)

  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);
  mass_hdibaryon = gConf.Get<G4double>( "HdibaryonMass" ); // 0.89594;
  width_hdibaryon = gConf.Get<G4double>( "HdibaryonWidth" ); // 0.0487;
  Kinema3Resonance Hdibaryon(piMinus->GetPDGMass()/GeV,
			     proton->GetPDGMass()/GeV,
			     kaonPlus->GetPDGMass()/GeV,
			     piMinus->GetPDGMass()/GeV,
			     Lambda->GetPDGMass()/GeV,
			     mass_hdibaryon, width_hdibaryon, pg_z, 0.0);
  //  Hdibaryon.Dump();
  Energy_kp = Hdibaryon.GetEnergy(5);
  // G4double momentum_kp = Hdibaryon.GetMomentum(5);
  Hdibaryon.GetMomentum(5,mom);
  /* L1 */
  Energy_L1 = Hdibaryon.GetEnergy(3);
  // G4double momentum_L1 = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);
  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];
  // G4double ThetaL1 = Hdibaryon.GetTheta(3);
  // G4double PhiL1 = Hdibaryon.GetPhi(3);
  /* L2 */
  Energy_L2 = Hdibaryon.GetEnergy(4);
  // G4double momentum_L2 = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double ThetaL2 = Hdibaryon.GetTheta(4);
  // G4double PhiL2 = Hdibaryon.GetPhi(4);

  /* Kp */
  // G4double Energy_kp = Hdibaryon.GetEnergy(5);
  // G4double momentum_kp = Hdibaryon.GetMomentum(5);
  Hdibaryon.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  // G4double Thetakp = Hdibaryon.GetThetaCM(1);
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  // G4double beta= sqrt(pbeam*pbeam)/(0.938272013+Ebeam);
  // G4double momcmk[4] = {};
  // G4double momk[4] = {};
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+Lambda->GetPDGMass()/GeV*Lambda->GetPDGMass()/GeV);
  // G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  // cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  // coscmk->Fill((momcmk[2]/momcmkpp));

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

  // G4double misskp = missks(pbeam,piMinus->GetPDGMass()/GeV,Lambda->GetPDGMass()/GeV,momk);//pbeam, mass, momk
  //missk->Fill(misskp);

  //cal for invariant mass
  // G4double ptot= pow(mom_L1_x+mom_L2_x,2) +
  //   pow(mom_L1_y+mom_L2_y,2) + pow(mom_L1_z+mom_L2_z,2);
  // G4double e1=Energy_L1;
  // G4double e2=Energy_L2;
  // G4double etot=pow(e1+e2,2);
  // G4double invm2=etot-ptot;
  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);
  double vtx = 0.*mm;
  double vty = 0.*mm;
  /////shhwang hdibaryon1

  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);

  ///////kp PG
  particleGun->SetParticleDefinition(Lambda);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - Lambda->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////L1 PG
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  particleGun->SetParticleEnergy((Energy_L1 - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////L2 PG
  particleGun->SetParticleDefinition(piMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  particleGun->SetParticleEnergy((Energy_L2 - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);



  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,piMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_pip_KstarS( G4Event* anEvent )
{
  G4double mass_hdibaryon, width_hdibaryon;

  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  //  G4double rmk=0.493677;

  auto particleTable = G4ParticleTable::GetParticleTable();
  auto proton = particleTable->FindParticle("proton");
  auto Lambda = particleTable->FindParticle("sigma0");
  auto kaonPlus = particleTable->FindParticle("kaon+");
  // auto kaonMinus = particleTable->FindParticle("kaon-");
  auto piMinus = particleTable->FindParticle("pi-");

  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = 1.8;
  // G4double pbeam=1.8;

  // G4double Ebeam = sqrt(pbeam*pbeam+piMinus->GetPDGMass()/GeV*piMinus->GetPDGMass()/GeV);       // Incident gamma energy (GeV)

  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);
  mass_hdibaryon = gConf.Get<G4double>( "HdibaryonMass" ); // 0.89594;
  width_hdibaryon = gConf.Get<G4double>( "HdibaryonWidth" ); // 0.0487;

  Kinema3Resonance Hdibaryon(piMinus->GetPDGMass()/GeV,
			     proton->GetPDGMass()/GeV,
			     kaonPlus->GetPDGMass()/GeV,
			     piMinus->GetPDGMass()/GeV,
			     Lambda->GetPDGMass()/GeV,
			     mass_hdibaryon, width_hdibaryon, pg_z, 0.0);

//  Hdibaryon.Dump();

  Energy_kp = Hdibaryon.GetEnergy(5);
  // G4double momentum_kp = Hdibaryon.GetMomentum(5);
  Hdibaryon.GetMomentum(5,mom);

  /* L1 */
  Energy_L1 = Hdibaryon.GetEnergy(3);
  // G4double momentum_L1 = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);

  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];

  // G4double ThetaL1 = Hdibaryon.GetTheta(3);
  // G4double PhiL1 = Hdibaryon.GetPhi(3);

  /* L2 */
  Energy_L2 = Hdibaryon.GetEnergy(4);
  // G4double momentum_L2 = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double ThetaL2 = Hdibaryon.GetTheta(4);
  // G4double PhiL2 = Hdibaryon.GetPhi(4);

  /* Kp */
  Energy_kp = Hdibaryon.GetEnergy(5);
  // G4double momentum_kp = Hdibaryon.GetMomentum(5);
  Hdibaryon.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  // G4double Thetakp = Hdibaryon.GetThetaCM(1);
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  // G4double beta= sqrt(pbeam*pbeam)/(0.938272013+Ebeam);
  // G4double momk[4] = {};
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+Lambda->GetPDGMass()/GeV*Lambda->GetPDGMass()/GeV);
  // G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  // cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  // coscmk->Fill((momcmk[2]/momcmkpp));

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

  // G4double misskp = missks(pbeam,piMinus->GetPDGMass()/GeV,Lambda->GetPDGMass()/GeV,
  // 			   momk);//pbeam, mass, momk
  //  missk->Fill(misskp);

  //cal for invariant mass
  // G4double ptot=pow(mom_L1_x+mom_L2_x,2)+pow(mom_L1_y+mom_L2_y,2)+pow(mom_L1_z+mom_L2_z,2);
  // G4double e1=Energy_L1;
  // G4double e2=Energy_L2;
  // G4double etot=pow(e1+e2,2);
  // G4double invm2=etot-ptot;
  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);

  double vtx = 0.*mm;
  double vty = 0.*mm;
  /////shhwang hdibaryon1
  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  ///////kp PG
  particleGun->SetParticleDefinition(Lambda);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - Lambda->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////L1 PG
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  particleGun->SetParticleEnergy((Energy_L1 - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  ///////L2 PG
  particleGun->SetParticleDefinition(piMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  particleGun->SetParticleEnergy((Energy_L2 - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,piMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study( G4Event* anEvent )
{

  double pbm[4]={-9999.9999}, pka[4]={-9999.9999}, vtx[3]={-9999.9999};
  // double pL1[4]={-9999.9999}, pL2[4]={-9999.9999};
  //,xL1[3]={-9999.9999}, xL2[3]={-9999.9999};
  //  double pp1[4]={-9999.9999}, ppi1[4]={-9999.9999}, pp2[4]={-9999.9999}, ppi2[4]={-9999.9999};

  char fnam[30];

  auto particleTable = G4ParticleTable::GetParticleTable();
  // auto kaonMinus = particleTable->FindParticle("kaon-");
  auto kaonPlus = particleTable->FindParticle("kaon+");
  // auto proton = particleTable->FindParticle("proton");
  // auto piMinus = particleTable->FindParticle("pi-");
  // auto Lambda = particleTable->FindParticle("lambda");
  // auto lambda1 = particleTable->FindParticle("lambda");
  // auto lambda2 = particleTable->FindParticle("lambda");

  FILE *fp;
  sprintf(fnam,"inc64cu.dat");

  if ((fp = fopen(fnam,"r")) == NULL){
    fprintf(stderr, "inc64cu.dat:: Cannot open file: %s\n", fnam);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check=0.;
  // 22478 lines --> remove NaN
  int ran = G4RandFlat::shoot(1.,22478.);

  while(1){
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
	   ,&data[0],&data[1],&data[2],&data[3],&data[4]
	   ,&data[5],&data[6],&data[7],&data[8],&data[9]
	   ,&data[10],&data[11],&data[12],&data[13],&data[14]
	   ,&data[15],&data[16],&data[17],&data[18],&data[19]
	   ,&data[20],&data[21],&data[22],&data[23],&data[24]
	   ,&data[25],&data[26],&data[27],&data[28],&data[29]
	   ,&data[30],&data[31],&data[32],&data[33],&data[34]
	   ,&data[35],&data[36],&data[37],&data[38],&data[39]
	   ,&data[40]);
    check=check+1.;
    if(check==ran) break;
    //    if(res==EOF) break;
  }

  for(int ii=0;ii<4;ii++){
    pbm[ii]=data[ii]/1000.;
    pka[ii]=data[ii+4]/1000.;
    // pL1[ii]=data[ii+11]/1000.;
    // pL2[ii]=data[ii+15]/1000.;
  }

  for(int ii=0;ii<3;ii++){
    vtx[ii]=data[ii+8];
  }

  fclose(fp);

  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  vtx[2]= G4RandFlat::shoot(100.-m_target_size.z()/2,100.+m_target_size.z()/2);

  // G4double Energy_beam = pbm[3];
  G4ThreeVector momentumBeam(pbm[0],pbm[1],pbm[2]);
  //  G4ThreeVector vertexPos(vtx[0]*mm,vtx[1]*mm, vtx[2]*mm);
  G4ThreeVector vertexPos(0.*mm,0.*mm, vtx[2]*mm);

  G4double Energy_ka = pka[3]; //total energy sqrt(pp^2+rmk^2)??
  // G4double Energy_L1 = pL1[3]; //lambda momenta ?? --> 4 momomenta sqrt(p^2+m^2)
  // G4double Energy_L2 = pL2[3]; //lambda momenta ?? --> 4 momomenta sqrt(p^2+m^2)

  // G4double beta= pbm[2]/(2*0.93827203+pbm[3]);
  // G4double momcmk[4] = {};
  // G4double momk[4];
  // momk[0]=pka[0];
  // momk[1]=pka[1];
  // momk[2]=pka[2];
  // momk[3]=sqrt(pka[0]*pka[0]+pka[1]*pka[1]+pka[2]*pka[2]+0.493677*0.493677);
  // G4double momkpp=sqrt(pow(momk[0],2)+pow(momk[1],2)+pow(momk[2],2));
  // G4double test = lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));

  // cmk->Fill(acos(momcmk[2]/momcmkpp)*180/3.141592654);
  // coscmk->Fill(momcmk[2]/momcmkpp);

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

  // G4double shphi=(atan2(momk[1],momk[0]))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4cout<<pka[1]<<G4endl;
  // ---- K+ -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(momentumKaonPlus);
  particleGun->SetParticleEnergy((Energy_ka - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(vertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);

  // ---- K-(beam)  -------------
  //  G4ThreeVector momentumLambda1(pbm[0], pbm[1], pbm[2]);
  //  particleGun->SetParticleDefinition(kaonMinus);
  //  particleGun->SetParticleMomentumDirection(momentumBeam);
  //  particleGun->SetParticleEnergy((Energy_beam - kaonMinus->GetPDGMass()/GeV)*GeV);
  //  particleGun->SetParticlePosition(vertexPos);
  //  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,pka[0],pka[1],pka[2],kaonPlus->GetPDGMass()/GeV);
  //  gAnaMan.SetPrimaryParticle(1,pbm[0],pbm[1],pbm[2],kaonMinus->GetPDGMass()/GeV);

  gAnaMan.SetPrimaryVertex(0,0.,0.,vtx[2]);
  //  gAnaMan.SetPrimaryVertex(1,0.,0.,vtx[2]);

  // G4double rmk=0.493677;
  // G4double misskp = miss1(pbm,rmk,momk);//pbeam, mass, momk

  //  missk->Fill(misskp);
  // G4double ptot=pow(pL1[0]+pL2[0],2)+pow(pL1[1]+pL2[1],2)+pow(pL1[2]+pL2[2],2);
  // G4double e1=Energy_L1;
  // G4double e2=Energy_L2;
  // G4double etot=pow(e1+e2,2);
  // G4double invm2=(etot-ptot);
  // if(invm2 > 0) invm=sqrt(invm2);
  // gen_im->Fill(invm);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_kp( G4Event* anEvent )
{

  double pbm[4]={-9999.9999}, pka[4]={-9999.9999}, vtx[3]={-9999.9999};
  // double pL1[4]={-9999.9999}, pL2[4]={-9999.9999};
  //,xL1[3]={-9999.9999}, xL2[3]={-9999.9999};
  //  double pp1[4]={-9999.9999}, ppi1[4]={-9999.9999}, pp2[4]={-9999.9999}, ppi2[4]={-9999.9999};

  char fnam[30];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonPlus;
  kaonPlus = particleTable->FindParticle("kaon+");

  FILE *fp;
  sprintf(fnam,"inc64cu.dat");

  if ((fp = fopen(fnam,"r")) == NULL){
    fprintf(stderr, "inc64cu.dat:: Cannot open file: %s\n", fnam);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check=0.;
  // 22478 lines --> remove NaN
  int ran = G4RandFlat::shoot(1.,22478.);

  while(1){
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
	   ,&data[0],&data[1],&data[2],&data[3],&data[4]
	   ,&data[5],&data[6],&data[7],&data[8],&data[9]
	   ,&data[10],&data[11],&data[12],&data[13],&data[14]
	   ,&data[15],&data[16],&data[17],&data[18],&data[19]
	   ,&data[20],&data[21],&data[22],&data[23],&data[24]
	   ,&data[25],&data[26],&data[27],&data[28],&data[29]
	   ,&data[30],&data[31],&data[32],&data[33],&data[34]
	   ,&data[35],&data[36],&data[37],&data[38],&data[39]
	   ,&data[40]);
    check=check+1.;
    if(check==ran) {
      //      G4cout<<check<<G4endl;
      break;
    }
    //    if(res==EOF) break;
  }

  for(int ii=0;ii<4;ii++){
    pbm[ii]=data[ii]/1000.;
    pka[ii]=data[ii+4]/1000.;
    // pL1[ii]=data[ii+11]/1000.;
    // pL2[ii]=data[ii+15]/1000.;
  }

  for(int ii=0;ii<3;ii++){
    vtx[ii]=data[ii+8];
  }

  fclose(fp);

  G4double Energy_ka;
  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  vtx[2]= G4RandFlat::shoot(100.-m_target_size.z()/2,100.+m_target_size.z()/2);
  // G4double Energy_beam = pbm[3];
  G4ThreeVector momentumBeam(pbm[0],pbm[1],pbm[2]);
  //  G4ThreeVector vertexPos(vtx[0]*mm,vtx[1]*mm, vtx[2]*mm);
  G4ThreeVector vertexPos(0.*mm,0.*mm, vtx[2]*mm);

  Energy_ka = pka[3]; //total energy sqrt(pp^2+rmk^2)??
  // ---- K+ -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(momentumKaonPlus);
  particleGun->SetParticleEnergy((Energy_ka - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(vertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,pka[0],pka[1],pka[2],kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,0.,0.,vtx[2]);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_kp_beam( G4Event* anEvent )
{

  double pbm[4]={-9999.9999}, pka[4]={-9999.9999}, vtx[3]={-9999.9999};
  // double pL1[4]={-9999.9999}, pL2[4]={-9999.9999};
  //,xL1[3]={-9999.9999}, xL2[3]={-9999.9999};
  //  double pp1[4]={-9999.9999}, ppi1[4]={-9999.9999}, pp2[4]={-9999.9999}, ppi2[4]={-9999.9999};

  char fnam[30];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonPlus;
  kaonPlus = particleTable->FindParticle("kaon+");

  FILE *fp;
  sprintf(fnam,"inc64cu.dat");

  if ((fp = fopen(fnam,"r")) == NULL){
    fprintf(stderr, "inc64cu.dat:: Cannot open file: %s\n", fnam);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check=0.;
  // 22478 lines --> remove NaN
  int ran = G4RandFlat::shoot(1.,22478.);

  while(1){
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
	   ,&data[0],&data[1],&data[2],&data[3],&data[4]
	   ,&data[5],&data[6],&data[7],&data[8],&data[9]
	   ,&data[10],&data[11],&data[12],&data[13],&data[14]
	   ,&data[15],&data[16],&data[17],&data[18],&data[19]
	   ,&data[20],&data[21],&data[22],&data[23],&data[24]
	   ,&data[25],&data[26],&data[27],&data[28],&data[29]
	   ,&data[30],&data[31],&data[32],&data[33],&data[34]
	   ,&data[35],&data[36],&data[37],&data[38],&data[39]
	   ,&data[40]);
    check=check+1.;
    if(check==ran) {
      //      G4cout<<check<<G4endl;
      break;
    }
    //    if(res==EOF) break;
  }

  for(int ii=0;ii<4;ii++){
    pbm[ii]=data[ii]/1000.;
    pka[ii]=data[ii+4]/1000.;
    // pL1[ii]=data[ii+11]/1000.;
    // pL2[ii]=data[ii+15]/1000.;
  }

  //  for(int ii=0;ii<3;ii++){
  //    vtx[ii]=data[ii+8];
  //  }
  fclose(fp);

  G4double Energy_ka;
  vtx[0] = G4RandFlat::shoot(-15.,15.)*mm;
  vtx[1] = G4RandFlat::shoot(-5.,5.)*mm;
  vtx[2]= G4RandFlat::shoot(100.-m_target_size.z()/2,100.+m_target_size.z()/2)*mm;

  // G4double Energy_beam = pbm[3];
  G4ThreeVector momentumBeam(pbm[0],pbm[1],pbm[2]);

  gAnaMan.SetPrimaryBeam(pbm[0],pbm[1],pbm[2]);

  G4ThreeVector vertexPos(vtx[0],vtx[1],vtx[2]);

  Energy_ka = pka[3]; //total energy sqrt(pp^2+rmk^2)??
  // ---- K+ -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(momentumKaonPlus);
  particleGun->SetParticleEnergy((Energy_ka - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(vertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,pka[0],pka[1],pka[2],kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx[0],vtx[1],vtx[2]);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_all( G4Event* anEvent )
{

  double pbm[4]={-9999.9999}, pka[4]={-9999.9999}, vtx[3]={-9999.9999};
  // double pL1[4]={-9999.9999}, pL2[4]={-9999.9999};
  //,xL1[3]={-9999.9999}, xL2[3]={-9999.9999};
  //  double pp1[4]={-9999.9999}, ppi1[4]={-9999.9999}, pp2[4]={-9999.9999}, ppi2[4]={-9999.9999};

  char fnam[30];

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //  G4ParticleDefinition* kaonMinus;
  G4ParticleDefinition* kaonPlus;
  //  G4ParticleDefinition* proton;
  //  G4ParticleDefinition* piMinus;
  //  G4ParticleDefinition* piPlus;

  //  G4ParticleDefinition* Lambda;
  //    G4ParticleDefinition* lambda1;
  //    G4ParticleDefinition* lambda2;

  //  kaonMinus = particleTable->FindParticle("kaon-");
  kaonPlus = particleTable->FindParticle("kaon+");
  //  proton = particleTable->FindParticle("proton");
  //  piMinus = particleTable->FindParticle("pi-");
  //  piPlus = particleTable->FindParticle("pi+");

  FILE *fp;
  sprintf(fnam,"inc64cu.dat");

  if ((fp = fopen(fnam,"r")) == NULL){
    fprintf(stderr, "inc64cu.dat:: Cannot open file: %s\n", fnam);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check=0.;
  // 22478 lines --> remove NaN
  int ran = G4RandFlat::shoot(1.,22478.);

  while(1){
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
	   ,&data[0],&data[1],&data[2],&data[3],&data[4]
	   ,&data[5],&data[6],&data[7],&data[8],&data[9]
	   ,&data[10],&data[11],&data[12],&data[13],&data[14]
	   ,&data[15],&data[16],&data[17],&data[18],&data[19]
	   ,&data[20],&data[21],&data[22],&data[23],&data[24]
	   ,&data[25],&data[26],&data[27],&data[28],&data[29]
	   ,&data[30],&data[31],&data[32],&data[33],&data[34]
	   ,&data[35],&data[36],&data[37],&data[38],&data[39]
	   ,&data[40]);
    check=check+1.;
    if(check==ran) break;
    //    if(res==EOF) break;
  }

  for(int ii=0;ii<4;ii++){
    pbm[ii]=data[ii]/1000.;
    pka[ii]=data[ii+4]/1000.;
    // pL1[ii]=data[ii+11]/1000.;
    // pL2[ii]=data[ii+15]/1000.;
  }

  for(int ii=0;ii<3;ii++){
    vtx[ii]=data[ii+8];
  }

  fclose(fp);

  //  G4double Energy_p;
  G4double Energy_kp;
  //  G4double Energy_pip;
  //  G4double Energy_pin;
  //  G4double Energy_kn;

  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  vtx[2]= G4RandFlat::shoot(100.-m_target_size.z()/2,100.+m_target_size.z()/2);

  // G4double Energy_beam = pbm[3];
  G4ThreeVector momentumBeam(pbm[0],pbm[1],pbm[2]);
  //  G4ThreeVector vertexPos(vtx[0]*mm,vtx[1]*mm, vtx[2]*mm);
  G4ThreeVector vertexPos(0.*mm,0.*mm, vtx[2]*mm);
  // G4double Energy_L1 = pL1[3]; //total energy sqrt(pp^2+rmk^2)??
  // G4double Energy_L2 = pL2[3]; //total energy sqrt(pp^2+rmk^2)??

  //  Energy_p = sqrt(pow(pka[0],2)+pow(pka[1],2)+pow(pka[2],2)+proton->GetPDGMass()/GeV*proton->GetPDGMass()/GeV); //total energy sqrt(pp^2+rmk^2)??
  Energy_kp = sqrt(pow(pka[0],2)+pow(pka[1],2)+pow(pka[2],2)+kaonPlus->GetPDGMass()/GeV*kaonPlus->GetPDGMass()/GeV); //total energy sqrt(pp^2+rmk^2)??
  //  Energy_pip = sqrt(pow(pka[0],2)+pow(pka[1],2)+pow(pka[2],2)+piPlus->GetPDGMass()/GeV*piPlus->GetPDGMass()/GeV); //total energy sqrt(pp^2+rmk^2)??
//  Energy_pin = sqrt(pow(pka[0],2)+pow(pka[1],2)+pow(pka[2],2)+proton->GetPDGMass()/GeV*proton->GetPDGMass()/GeV); //total energy sqrt(pp^2+rmk^2)??
//  Energy_kn = sqrt(pow(pka[0],2)+pow(pka[1],2)+pow(pka[2],2)+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV); //total energy sqrt(pp^2+rmk^2)??

  // G4double beta= pbm[2]/(2*0.93827203+pbm[3]);
  // G4double momcmk[4]={0.};
  // G4double momk[4]={0.};
  // momk[0]=pka[0];
  // momk[1]=pka[1];
  // momk[2]=pka[2];
  // momk[3]=sqrt(pka[0]*pka[0]+pka[1]*pka[1]+pka[2]*pka[2]+0.493677*0.493677);

  // G4double momkpp = sqrt(pow(momk[0],2)+pow(momk[1],2)+pow(momk[2],2));

  // G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));

  // cmk->Fill(acos(momcmk[2]/momcmkpp)*180/3.141592654);
  // coscmk->Fill(momcmk[2]/momcmkpp);

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

  // G4double shphi=(atan2(momk[1],momk[0]))/3.141592654*180.;
  //phik->Fill(shphi);

  // ---- proton -------------------
  //  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  //  particleGun->SetParticleDefinition(proton);
  //  particleGun->SetParticleMomentumDirection(momentumKaonPlus);
  //  particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
  //  particleGun->SetParticlePosition(vertexPos);
  //  particleGun->GeneratePrimaryVertex(anEvent);

  // ---- K+ -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(momentumKaonPlus);
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(vertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);

  // ---- pi+ -------------------
  //  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  //  particleGun->SetParticleDefinition(piPlus);
  //  particleGun->SetParticleMomentumDirection(momentumKaonPlus);
  //  particleGun->SetParticleEnergy((Energy_pip - piPlus->GetPDGMass()/GeV)*GeV);
  //  particleGun->SetParticlePosition(vertexPos);
  //  particleGun->GeneratePrimaryVertex(anEvent);
  /*
  // ---- pi- -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  particleGun->SetParticleDefinition(piMinus);
  particleGun->SetParticleMomentumDirection(momentumKaonPlus);
  particleGun->SetParticleEnergy((Energy_pin - piMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(vertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);

  // ---- K- -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  particleGun->SetParticleDefinition(kaonMinus);
  particleGun->SetParticleMomentumDirection(momentumKaonPlus);
  particleGun->SetParticleEnergy((Energy_kn - kaonMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(vertexPos);
  particleGun->GeneratePrimaryVertex(anEvent);
  */

  //  // ---- K-(beam)  -------------
  //  G4ThreeVector momentumLambda1(pbm[0], pbm[1], pbm[2]);
  //  particleGun->SetParticleDefinition(kaonMinus);
  //  particleGun->SetParticleMomentumDirection(momentumBeam);
  //  particleGun->SetParticleEnergy((Energy_beam - kaonMinus->GetPDGMass()/GeV)*GeV);
  //  particleGun->SetParticlePosition(vertexPos);
  //  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,pka[0],pka[1],pka[2],kaonPlus->GetPDGMass()/GeV);
  //  gAnaMan.SetPrimaryParticle(1,pka[0],pka[1],pka[2],kaonPlus->GetPDGMass()/GeV);
  //  gAnaMan.SetPrimaryParticle(2,pka[0],pka[1],pka[2],piPlus->GetPDGMass()/GeV);

  gAnaMan.SetPrimaryVertex(0,0.,0.,vtx[2]);
  //  gAnaMan.SetPrimaryVertex(1,0.,0.,vtx[2]);
  //  gAnaMan.SetPrimaryVertex(2,0.,0.,vtx[2]);

  // G4double rmk=0.493677;
  // G4double misskp = miss1(pbm,rmk,momk);//pbeam, mass, momk

  //missk->Fill(misskp);
  // G4double ptot=pow(pL1[0]+pL2[0],2)+pow(pL1[1]+pL2[1],2)+pow(pL1[2]+pL2[2],2);
  // G4double e1=Energy_L1;
  // G4double e2=Energy_L2;
  // G4double etot=pow(e1+e2,2);
  // G4double invm2=(etot-ptot);
  // if( invm2 > 0 )
  //   invm=sqrt(invm2);
  //gen_im->Fill(invm);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_knp( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonMinus;
  G4ParticleDefinition* proton;

  kaonMinus = particleTable->FindParticle("kaon-");
  proton = particleTable->FindParticle("proton");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			kaonMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hdibaryon.GetEnergy(3);
  // G4double momentum_kp = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hdibaryon.GetEnergy(4);
  // G4double momentum_h = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];
  //  cout<<mom_h_z<<endl;
  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.) goto up;

  // G4double Theta_h = Hdibaryon.GetTheta(4);
  // G4double Phi_h = Hdibaryon.GetPhi(4);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*mm;
  //  G4double vty = 0.*mm;

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  //Kaon -
  particleGun->SetParticleDefinition(kaonMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //proton
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_knp_beam( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonMinus;
  G4ParticleDefinition* proton;

  kaonMinus = particleTable->FindParticle("kaon-");
  proton = particleTable->FindParticle("proton");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			kaonMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hdibaryon.GetEnergy(3);
  // G4double momentum_kp = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hdibaryon.GetEnergy(4);
  // G4double momentum_h = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];
  //  cout<<mom_h_z<<endl;
  //    if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.) goto up;

  // G4double Theta_h = Hdibaryon.GetTheta(4);
  // G4double Phi_h = Hdibaryon.GetPhi(4);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*mm;
  //  G4double vty = 0.*mm;

  G4double vtx = G4RandFlat::shoot(-15.,15.)*mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  //Kaon -
  particleGun->SetParticleDefinition(kaonMinus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonMinus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //proton
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_kpxi_beam( G4Event* anEvent )
{
  G4double mom[4];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonMinus;
  G4ParticleDefinition* kaonPlus;
  G4ParticleDefinition* Xin;
  G4ParticleDefinition* proton;

  kaonMinus = particleTable->FindParticle("kaon-");
  proton = particleTable->FindParticle("proton");
  kaonPlus = particleTable->FindParticle("kaon+");
  Xin = particleTable->FindParticle("xi-");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  ////FWHM 3.3 x 10^-4, FWHM ~ 2.3548 sigma
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  pbeam=G4RandGauss::shoot(1.7,1.7*3.3*0.0001/2.3548);
  //  G4cout<<pbeam<<G4endl;

  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);
  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			Xin->GetPDGMass()/GeV,
			kaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hdibaryon.GetEnergy(4);
  // G4double momentum_kp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.)
    goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hdibaryon.GetEnergy(3);
  // G4double momentum_h = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];
  //  cout<<mom_h_z<<endl;

  // G4double Theta_h = Hdibaryon.GetTheta(3);
  // G4double Phi_h = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(proton->GetPDGMass()/GeV+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  // cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*mm;
  //  G4double vty = 0.*mm;

  G4double vtx = G4RandFlat::shoot(-15.,15.)*mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  //Kaon -
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //proton
  particleGun->SetParticleDefinition(Xin);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - Xin->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,Xin->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_kpxi_beam_only_kp( G4Event* anEvent )
{
  G4double mom[4];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonMinus;
  G4ParticleDefinition* kaonPlus;
  G4ParticleDefinition* Xin;
  G4ParticleDefinition* proton;

  kaonMinus = particleTable->FindParticle("kaon-");
  proton = particleTable->FindParticle("proton");
  kaonPlus = particleTable->FindParticle("kaon+");
  Xin = particleTable->FindParticle("xi-");

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  pbeam=1.7;
  //  pbeam=G4RandGauss::shoot(1.7,1.7*3.3*0.0001/2.3548);
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);
  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			Xin->GetPDGMass()/GeV,
			kaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hdibaryon.GetEnergy(4);
  // G4double momentum_kp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.)
    goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  // G4double Energy_h = Hdibaryon.GetEnergy(3);
  // G4double momentum_h = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);
  // G4double mom_h_x = mom[0];
  // G4double mom_h_y = mom[1];
  // G4double mom_h_z = mom[2];
  //  cout<<mom_h_z<<endl;

  // G4double Theta_h = Hdibaryon.GetTheta(3);
  // G4double Phi_h = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(proton->GetPDGMass()/GeV+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*mm;
  //  G4double vty = 0.*mm;

  G4double vtx = G4RandFlat::shoot(-15.,15.)*mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  //Kaon -
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //proton
  /*  particleGun->SetParticleDefinition(Xin);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - Xin->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);
  */

  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonMinus->GetPDGMass()/GeV);
  //  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,Xin->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  //  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_kpxi1530( G4Event* anEvent )
{
  G4double mom[4];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonMinus;
  G4ParticleDefinition* kaonPlus;
  G4ParticleDefinition* Xi1530;
  G4ParticleDefinition* proton;

  kaonMinus = particleTable->FindParticle("kaon-");
  proton = particleTable->FindParticle("proton");
  kaonPlus = particleTable->FindParticle("kaon+");
  //  G4cout<<"test111"<<G4endl;
  Xi1530 = particleTable->FindParticle(3314);
  //  G4cout<<"test21"<<G4endl;
  //  G4cout<<Xi1530->GetPDGMass()/GeV<<G4endl;
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  ////FWHM 3.3 x 10^-4, FWHM ~ 2.3548 sigma
  //  pbeam=G4RandGauss::shoot(1.7,1.7*3.3*0.0001/2.3548);
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  G4cout<<pbeam<<G4endl;

  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);
  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			proton->GetPDGMass()/GeV,
			Xi1530->GetPDGMass()/GeV,
			kaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hdibaryon.GetEnergy(4);
  // G4double momentum_kp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  if((atan2(mom_kp_y,mom_kp_z))/3.141592*180.<-18. ||
     (atan2(mom_kp_y,mom_kp_z))/3.141592*180.>20. )
    goto up;

  //    if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hdibaryon.GetEnergy(3);
  // G4double momentum_h = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];
  //  cout<<mom_h_z<<endl;

  // G4double Theta_h = Hdibaryon.GetTheta(3);
  // G4double Phi_h = Hdibaryon.GetPhi(3);
  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(proton->GetPDGMass()/GeV+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180.);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*mm;
  //  G4double vty = 0.*mm;

  G4double vtx = G4RandFlat::shoot(-15.,15.)*mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  //Kaon -
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //  //Xi1530-
  particleGun->SetParticleDefinition(Xi1530);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - Xi1530->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,Xi1530->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_Takahashi( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonMinus;
  G4ParticleDefinition* kaonPlus;
  G4ParticleDefinition* Xin;
  G4ParticleDefinition* proton;

  kaonMinus = particleTable->FindParticle("kaon-");
  proton = particleTable->FindParticle("proton");
  kaonPlus = particleTable->FindParticle("kaon+");
  Xin = particleTable->FindParticle("xi-");
  G4double p_proton[4]={0};
  G4int angular_mom = G4RandFlat::shoot() * 6. < 2. ? 0.:1.;
  //  G4cout<<"test111:"<<angular_mom<<G4endl;
  //  G4cout<<G4RandFlat::shoot()<<G4endl;
  HarmonicFermiMomentum( angular_mom, p_proton);
  G4double p_fermi = sqrt(pow(p_proton[0],2)+pow(p_proton[1],2)+pow(p_proton[2],2));
  p_proton[3] = sqrt(p_fermi * p_fermi + proton->GetPDGMass()/GeV * proton->GetPDGMass()/GeV);
  G4cout<<p_proton[3]<<G4endl;
  G4double Al[13]={1., -1.22, 1.55, -1.08, 0.37, -0.15, 0.16, -0.38, -0.18, 0.09, 0.05, -0.01, 0.20}; //from Dauber, et al., PR 179 5 1968 at K- beam energy 1.7 GeV/c

  G4double cross_section=0.;

  G4cout<<"p_proton:"<<p_proton[0]<<":"<<p_proton[1]<<":"<<p_proton[2]<<G4endl;
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  ////FWHM 3.3 x 10^-4, FWHM ~ 2.3548 sigma
    pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  pbeam=G4RandGauss::shoot(1.7,1.7*3.3*0.0001/2.3548);
  //  G4cout<<pbeam<<G4endl;

  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+kaonMinus->GetPDGMass()/GeV*kaonMinus->GetPDGMass()/GeV);
  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  // up:
  G4double cosx=G4RandFlat::shoot(-1.,1.);
  KinemaFermi Hdibaryon(kaonMinus->GetPDGMass()/GeV,
			p_proton[3],
			Xin->GetPDGMass()/GeV,
			kaonPlus->GetPDGMass()/GeV,
			pbm, p_proton,cosx);

  for(G4int le=0;le<13;le++){
    cross_section=cross_section+Al[le];
  }

  Energy_kp=Hdibaryon.GetEnergy(4);
  // G4double momentum_kp = Hdibaryon.GetMomentum(4);
  Hdibaryon.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  //    if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hdibaryon.GetEnergy(3);
  // G4double momentum_h = Hdibaryon.GetMomentum(3);
  Hdibaryon.GetMomentum(3,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  cout<<mom_h_z<<endl;

  // G4double Theta_h = Hdibaryon.GetTheta(3);
  // G4double Phi_h = Hdibaryon.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(proton->GetPDGMass()/GeV+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*mm;
  //  G4double vty = 0.*mm;
  G4double vtx = G4RandFlat::shoot(-15.,15.)*mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  //Kaon -
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  particleGun->SetParticleEnergy((Energy_kp - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  //proton
  particleGun->SetParticleDefinition(Xin);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  particleGun->SetParticleEnergy((Energy_h - Xin->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,kaonMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,Xin->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_pro_08_20( G4Event* anEvent )
{
  G4double Energy_p, mom_p_x, mom_p_y, mom_p_z;
  //  G4double mom[3];
  //  G4double Ebeam,pg_x,pg_y,pg_z, pbeam;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* proton;
  proton = particleTable->FindParticle("proton");

  /* proton */
 up:
  G4double phi=G4RandFlat::shoot(-1.,1.)*3.141592654;
  G4double mom_p=G4RandFlat::shoot(0.8,2.0);
  //  G4double theta=acos(G4RandFlat::shoot(-1.,1.));
  ///  G4double theta=acos(G4RandFlat::shoot(0.939692646,1.));
  //  G4double theta=acos(G4RandFlat::shoot(0.866025458,1.));//30 deg
  G4double theta=acos(G4RandFlat::shoot(0.80,1.));//

  //  G4double theta=acos(G4RandFlat::shoot(0.962646,0.97));
  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);

  if((atan2(mom_p_y,mom_p_z))/3.141592*180.<-18. || (atan2(mom_p_y,mom_p_z))/3.141592*180.>18. )
    goto up;

  //  G4cout<<atan2(mom_p_y,mom_p_z)*180./3.141592<<G4endl;

  // labk->Fill(theta*180/3.141592);
  // coslabk->Fill(cos(theta));
  // phik->Fill(phi/3.141592*180);

  //  G4cout<<mom_p_x<<":"<<mom_p_y<<":"<<mom_p_z<<G4endl;

  Energy_p=sqrt(pow(mom_p,2)+pow(proton->GetPDGMass()/GeV,2));
  /////vertex
  G4double vtx = G4RandFlat::shoot(-15.,15.)*mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*mm;
  G4double vtz = G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  ///////proton PG
  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_E07_study_kp_04_15( G4Event* anEvent )
{
  //  G4cout<<"e07_generator"<<G4endl;
  G4double Energy_p, mom_p_x, mom_p_y, mom_p_z;
  //  G4double mom[3];
  //  G4double Ebeam,pg_x,pg_y,pg_z, pbeam;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* kaonPlus;
  kaonPlus = particleTable->FindParticle("kaon+");

  /* kaon+ */
 up:
  G4double phi=G4RandFlat::shoot(-1.,1.)*3.141592654;
  G4double mom_p=G4RandFlat::shoot(0.4,1.5);
  //  G4double theta=acos(G4RandFlat::shoot(-1.,1.));
  //  G4double theta=acos(G4RandFlat::shoot(0.866025458,1.));//30 deg

  G4double theta=acos(G4RandFlat::shoot(0.8,1.));
  //  G4double theta=acos(G4RandFlat::shoot(0.939692646,1.));//25
  //  G4double theta=acos(G4RandFlat::shoot(0.962646,0.97));


  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);

  if((atan2(mom_p_y,mom_p_z))/3.141592*180.<-18. || (atan2(mom_p_y,mom_p_z))/3.141592*180.>18. )
    goto up;

  // labk->Fill(theta*180/3.141592);
  // coslabk->Fill(cos(theta));
  // phik->Fill(phi/3.141592*180);

  Energy_p=sqrt(pow(mom_p,2)+pow(kaonPlus->GetPDGMass()/GeV,2));
  //  Energy_p=pow(mom_p,2)+kaonPlus->GetPDGMass()/GeV;
  /////vertex
  G4double vtx = G4RandFlat::shoot(-15.,15.)*mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*mm;
  G4double vtz = G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;

  ///////proton PG
  particleGun->SetParticleDefinition(kaonPlus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  particleGun->SetParticleEnergy((Energy_p - kaonPlus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,kaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::Generate_test2( G4Event* anEvent )
{
  //  G4double mass_hdibaryon, width_hdibaryon;
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //  G4double mom[3];
  //  G4double Ebeam,pg_x,pg_y,pg_z, pbeam;

  auto particleTable = G4ParticleTable::GetParticleTable();
  // auto proton = particleTable->FindParticle("proton");
  auto piminus = particleTable->FindParticle("pi-");

  /* proton */
  //  G4double mompx=0.10;
  G4double mompx=0.10;
  G4double mompz=0.0;

  G4double rand= G4RandFlat::shoot(-1.,1.)*3.141592654;
  mom_p_x = cos(rand)*mompx-sin(rand)*mompz;
  mom_p_y = 0.;
  mom_p_z = sin(rand)*mompx+cos(rand)*mompz;
  G4double mom_p=sqrt(pow(mom_p_x,2)+pow(mom_p_y,2)+pow(mom_p_z,2));
  //  Energy_p=pow(mom_p,2)+proton->GetPDGMass()/GeV;
  Energy_p=pow(mom_p,2)+piminus->GetPDGMass()/GeV;

  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*mm;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);

  //  double vtz = -150+0.0*((double) G4RandFlat::shoot()); //--> 0 mm
  G4double vtx = 0.; //--> 0 mm
  G4double vty = 0.; //--> 0 mm

  ///////kp PG
  //  particleGun->SetParticleDefinition(proton);
  particleGun->SetParticleDefinition(piminus);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  //  particleGun->SetParticleEnergy((Energy_p - proton->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticleEnergy((Energy_p - piminus->GetPDGMass()/GeV)*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  particleGun->GeneratePrimaryVertex(anEvent);

}

//_____________________________________________________________________________
double
TPCPrimaryGeneratorAction::RandSin( void )
{
  int success=0;
  double x,fx;
  do {
    x = 180.0 * (double)G4RandFlat::shoot();
    fx = sin(TMath::DegToRad()*x);
    if (fx >= (double)G4RandFlat::shoot())
      success = 1;
  } while (success==0);
  return x;
}

//_____________________________________________________________________________
G4double
TPCPrimaryGeneratorAction::miss1( G4double pbeam, G4double m1,G4double *p1 )
{
  G4double miss = -9999.;
  G4double rmp = 0.938272013;;
  G4double px = p1[0];
  G4double py = p1[1];
  G4double pz = p1[2];
  G4double pp = sqrt(px*px+py*py+pz*pz);
  G4double ebeam = sqrt(pow(pbeam,2)+pow(m1,2));
  G4double pmiss = pow(px,2) +  pow(py,2) +  pow(pz-pbeam,2);
  G4double emiss = pow((ebeam+rmp*2-sqrt(pow(pp,2)+pow(m1,2))),2);
  if( emiss - pmiss > 0 ){
    miss = sqrt( emiss - pmiss );
  }
  return miss;
}

//_____________________________________________________________________________
G4double
TPCPrimaryGeneratorAction::missks( G4double pbeam, G4double mbeam,
				   G4double m1,G4double *p1 )
{
  G4double miss = -9999.;
  G4double rmp = 0.938272013;;
  G4double px = p1[0];
  G4double py = p1[1];
  G4double pz = p1[2];
  G4double pp = sqrt(px*px+py*py+pz*pz);
  G4double ebeam = sqrt(pow(pbeam,2)+pow(mbeam,2));
  G4double pmiss = pow(px,2) +  pow(py,2) +  pow(pz-pbeam,2);
  G4double emiss = pow((ebeam+rmp-sqrt(pow(pp,2)+pow(m1,2))),2);
  if( emiss - pmiss > 0 ){
    miss = sqrt( emiss - pmiss );
  }
  return miss;
}

//_____________________________________________________________________________
G4double
TPCPrimaryGeneratorAction::miss1( G4double *pbeam, G4double m1,G4double *p1 )
{
  G4double miss = -9999.;
  G4double rmp = 0.938272013;;
  G4double px = p1[0];
  G4double py = p1[1];
  G4double pz = p1[2];
  G4double pp = sqrt(px*px+py*py+pz*pz);
  G4double pbeamx=pbeam[0], pbeamy=pbeam[1], pbeamz=pbeam[2];
  G4double ebeam = sqrt(pow(pbeamx,2)+pow(pbeamy,2)+pow(pbeamz,2)+pow(m1,2));
  G4double pmiss = pow(px-pbeamx,2) +  pow(py-pbeamy,2) +  pow(pz-pbeamz,2);
  G4double emiss = pow((ebeam+rmp*2 - sqrt(pow(pp,2)+pow(m1,2)) ),2);
  if( emiss - pmiss > 0 ){
    miss = std::sqrt( emiss - pmiss );
  }
  return miss;
}

//_____________________________________________________________________________
G4double
TPCPrimaryGeneratorAction::lorentz( G4double *v1,G4double betaz,G4double *v2 )
{
  // boosts the vector v1 with beta
  //  G4double v1[4],betaz,v2[4];
  //  G4double v1[4],betaz,v2[4];
  G4double b2, gamma;
  //  lorentz=-1;
  b2= pow(betaz,2.);
  if(b2==0.){
    v2[0]=v1[0];
    v2[1]=v1[1];
    v2[2]=v1[2];
    v2[3]=v1[3];
    return -1;
  }else{
    gamma=1./sqrt(1.-b2);
    v2[0] = v1[0];
    v2[1] = v1[1];
    v2[2] = gamma*(v1[2] - betaz*v1[3]);
    v2[3] = gamma*(v1[3] - betaz*v1[2]);
    return 1;
  }
}

//_____________________________________________________________________________
G4double
TPCPrimaryGeneratorAction::lorentcmlab( G4double *v1, G4double betaz,
					G4double *v2 )
{
  // boosts the vector v1 with beta
  //  G4double v1[4],betaz,v2[4];
  //  G4double v1[4],betaz,v2[4];
  G4double b2, gamma;
  //  lorentz=-1;
  b2= pow(betaz,2.);
  if(b2==0.){
    v2[0]=v1[0];
    v2[1]=v1[1];
    v2[2]=v1[2];
    v2[3]=v1[3];
    return -1;
  }else{
    gamma=1./sqrt(1.-b2);
    v2[0] = v1[0];
    v2[1] = v1[1];
    v2[2] = gamma*(v1[2] + betaz*v1[3]);
    v2[3] = gamma*(v1[3] + betaz*v1[2]);
    return 1;
  }
}

//_____________________________________________________________________________
G4int
TPCPrimaryGeneratorAction::HarmonicFermiMomentum( G4int Angular_mom,
						  G4double *Kf )
{
  /////// Copy of H. Takahashi-san's code
  /////// revised to geant4 by S.Hwang
  G4double ymax, b;
  G4double x, y, yy, theta, phi;
  /*      THIS ROUTINE GENERATES FERMI MOMENTUM KF(3) BY HARMONIC */
  /*      OSCILATOR MODEL FOR 1S AND 1P */
  /*      INPUT; L=0 OR 1 */
  /*      OUTPUT; KF(3) */
  b = 1.62 / 0.197; /// unit GeV
  if (Angular_mom == 0) {
    ymax = exp(-1);
    do {
      x = G4RandFlat::shoot();
      y = x * x * exp(-b * b * x * x);
      yy = G4RandFlat::shoot() * ymax;
    } while (yy > y);
  } else {
    ymax = exp(-2) * 4 / (b * b);
    do {
      x = G4RandFlat::shoot();
      y = b * b * x * x * x * x * exp(-b * b * x * x);
      yy = G4RandFlat::shoot() * ymax;
    } while (yy > y);
  }
  theta=acos(G4RandFlat::shoot(-1.,1.));
  phi=(G4RandFlat::shoot(-1.,1.))*3.141592;
  //  IsotropicAngle(&theta, &phi);
  Kf[0] = x * sin(theta) * cos(phi);
  Kf[1] = x * sin(theta) * sin(phi);
  Kf[2] = x * cos(theta);
  //  G4cout<<x<<":"<<Kf[0]<<":"<<Kf[1]<<":"<<Kf[2]<<":"<<G4endl;
  return 0;
}

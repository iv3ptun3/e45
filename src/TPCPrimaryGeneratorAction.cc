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

#include "TMath.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "BeamMan.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "JamMan.hh"
#include "IncMan.hh"
#include "Kinema3Resonance.hh"
#include "KinemaHResonance.hh"
#include "Kinema3Body.hh"
#include "Kinema4Body.hh"
#include "KinemaHybrid.hh"
#include "KinemaHweak.hh"
#include "KinemaFermi.hh"
#include "KinemaKstar.hh"
#include "Kinematics.hh"
#include "PrintHelper.hh"
#include "TPCAnaManager.hh"
#include "TPCTrackBuffer.hh"

namespace
{
  using CLHEP::GeV;
  using CLHEP::keV;
  using CLHEP::mm;
  auto& gAnaMan = TPCAnaManager::GetInstance();
  const auto& gBeam = BeamMan::GetInstance();
  const auto& gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
  const auto& gJam  = JamMan::GetInstance();
  const auto& gInc  = IncMan::GetInstance();
  const auto particleTable = G4ParticleTable::GetParticleTable();
	auto& gTrackBuffer = TPCTrackBuffer::GetInstance();
}

//_____________________________________________________________________________
TPCPrimaryGeneratorAction::TPCPrimaryGeneratorAction( void )
  : G4VUserPrimaryGeneratorAction(),
    m_generator( gConf.Get<G4int>( "Generator" ) ),
    m_particle_gun( new G4ParticleGun ),
    m_target_pos( gGeom.GetGlobalPosition( "SHSTarget" )*mm ),
    m_target_size( gSize.GetSize( "Target" )*mm ),
    m_e45target_size( gSize.GetSize( "E45Target" )*mm ),
    m_beam( new BeamInfo ),
    m_mm_vert( new MMVertex),
    m_beam_p0( gConf.Get<G4double>( "BeamMom" ) /* 1.80*GeV */ ),
    m_jam(),
    m_inc(),
    m_Neutron( particleTable->FindParticle( "neutron" ) ),
    m_Proton( particleTable->FindParticle( "proton" ) ),
    m_Lambda( particleTable->FindParticle( "lambda" ) ),
    m_Lambda1405( particleTable->FindParticle( "lambda(1405)" ) ),
    m_Lambda1405R( particleTable->FindParticle( "lambda1405r" ) ),
    m_SigmaPlus( particleTable->FindParticle( "sigma+" ) ),
    m_SigmaMinus( particleTable->FindParticle( "sigma-" ) ),
    m_SigmaZero( particleTable->FindParticle( "sigma0" ) ),
    m_Sigma1385Zero( particleTable->FindParticle( "sigma(1385)0" ) ),
    m_Sigma1385R( particleTable->FindParticle( "sigma1385r" ) ),
    m_XiMinus( particleTable->FindParticle( "xi-" ) ),
    m_Xi1530Minus( particleTable->FindParticle( "xi(1530)-" ) ),
    m_PionPlus( particleTable->FindParticle( "pi+" ) ),
    m_PionMinus( particleTable->FindParticle( "pi-" ) ),
    m_PionZero( particleTable->FindParticle( "pi0" ) ),
    m_KaonPlus( particleTable->FindParticle( "kaon+" ) ),
    m_sKaonPlus( particleTable->FindParticle( "skaon+" ) ),
    m_KaonMinus( particleTable->FindParticle( "kaon-" ) ),
    m_KaonZeroS( particleTable->FindParticle( "kaon0S" ) ),
    m_KaonStarZero( particleTable->FindParticle( "k_star0" ) ),
    m_Hdibaryon( particleTable->FindParticle( "hdibaryon" ) ),
    m_HdibaryonS( particleTable->FindParticle( "hdibaryonS" ) ),
    m_HdibaryonLL( particleTable->FindParticle("hdibaryonLL") ),
    m_LLphase( particleTable->FindParticle( "phaseLL" ) ),
    m_HybridBaryon( particleTable->FindParticle( "hybridb" ) )
{
  G4cout << FUNC_NAME << G4endl
	 << "   Generator# = " << m_generator << G4endl;
  gAnaMan.SetGeneratorID( m_generator );
#ifdef DEBUG
  particleTable->DumpTable();
#endif
}

//_____________________________________________________________________________
TPCPrimaryGeneratorAction::~TPCPrimaryGeneratorAction( void )
{
  delete m_particle_gun;
  delete m_beam;
  delete m_mm_vert;
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent )
{
  if( gBeam.IsReady() ){
		if(gBeam.IsK18()){
    	G4int eventID = anEvent->GetEventID();
			*m_beam = gBeam.Get(eventID);
			if(m_beam->z == 0) return;
		}
		else if(gBeam.IsKurama()){
    	G4int eventID = anEvent->GetEventID();
			*m_beam = gBeam.Get(eventID);
		}
		else if (gBeam.IsMissMassXi()){
    	G4int eventID = anEvent->GetEventID();
			*m_mm_vert = gBeam.GetVertex(eventID);
		}
		else{
    	*m_beam = gBeam.Get();
		}
#ifdef DEBUG
    m_beam->Print();
#endif
  }

  if( gJam.IsReady() ){
    m_jam = gJam.Get();
#ifdef DEBUG
    m_jam->Print();
#endif
  }

  if( gInc.IsReady() ){
    m_inc = gInc.Get();
#ifdef DEBUG
    m_inc->Print();
#endif
  }

  switch( m_generator ){
  case  0: break; // no generation
  case  1: GenerateHanul( anEvent ); break; // shhwang
  case  2: GenerateHdibaryon1( anEvent ); break; // h-dibaryon --> LL
  case  3: GenerateHdibaryon2( anEvent ); break; // h-dibaryon --> LL K+ 15deg
  case  4: GenerateTest( anEvent ); break; // test sako-san's code
  case  5: GenerateTest2( anEvent ); break; // test sako-san's code
  case  6: GenerateHdibaryonPHSG( anEvent ); break; // h weak. H->Lppi-
  case  7: GenerateHdibaryonPHSGS( anEvent ); break; // h weak. H->S-p
  case  8: GenerateHdibaryonNonReso( anEvent ); break; // test non resonance LL??
  case  9: GenerateHdibaryonPHSGLL( anEvent ); break; // H gen by using phsg
  case 10: GenerateBeamVI( anEvent ); break;
  case 11: GenerateBeamVO( anEvent ); break;
  // case 11: GeneratePionPlusKsL( anEvent ); break; // Study pi-p --> KsL
  case 12: GeneratePionPlusKsS( anEvent ); break; // Study pi-p --> KsS
  case 13: GeneratePionPlusKstarL( anEvent ); break; // Study on pi-p --> KsS by using LL gen
  case 14: GeneratePionPlusKstarS( anEvent ); break; // Study on pi-p --> KsS by using LL gen
  case 15: GenerateKpXi2Body( anEvent ); break;
  case 150: GenerateKpXi2BodyUniform( anEvent ); break;
  case 151: GenerateKpXi1530Minus2BodyUniform( anEvent ); break;
  case 16: GenerateUniformProton( anEvent ); break;
  case 17: GenerateUniformPim( anEvent ); break;
  case 18: GenerateBeamProton( anEvent ); break;
  case 19: GenerateUniformProton_P( anEvent ); break;
  case 20: GenerateUniformProton_P_fixphi( anEvent ); break;
  case 21: GenerateHybrid3bodyMode1( anEvent ); break;
  case 22: GenerateHybrid3bodyMode2( anEvent ); break;
  case 23: GenerateHybrid3bodyMode3( anEvent ); break;
  case 24: GenerateHybrid3bodyMode4( anEvent ); break;
  case 25: GenerateUniformKaonPlus( anEvent ); break;
  case 30: GeneratePhaseSpace( anEvent ); break;
    /* phase space generator 12C + Kn --> 10Be L L Kp
       12C_amu = 12u --> 12*931.494061 MeV
       10C_amu = 10.012938u --> 10.012938*931.494061 MeV */
  case 31: GenerateJamInput( anEvent ); break;
  case 32: GenerateIncInput( anEvent ); break;
  case 33: GenerateJamInput_Randphi( anEvent ); break;
  case 34: GenerateLL_fromXiP( anEvent ); break;
  case 50: GenerateUniformProton_P_Multi( anEvent ); break;
  case 60: GenerateLambda1405Rad1( anEvent ); break; // normal decay
  case 61: GenerateLambda1405Rad2( anEvent ); break; // radioactive decay
  case 62: GenerateSigma1385( anEvent ); break; // Sigma 1385 normal dcay
  case 63: GenerateSigma1385Rad( anEvent ); break; // Sigma 1385 radioactive decay
  case 64: GenerateLambda1405Reso( anEvent ); break; // pi+pi- --> K0
    // S+pi-,S0pi0, S-pi+ --> Lambda from threshold to 1.5 GeV
  case 99: GenerateDedxSingle( anEvent ); break;
    // Study on dedx. only single particle will be generated
  case 98: GenerateAll( anEvent ); break; // Study on dedx. All particles will be generated
	case 100: GenerateKurama( anEvent); break;
	case 135: GeneratePionPlusBeamthrough( anEvent ); break;
	case -135: GeneratePionMinusBeamthrough( anEvent ); break;
	case 493: GenerateKaonPlusBeamthrough( anEvent ); break;
	case -493: GenerateKaonMinusBeamthrough( anEvent ); break;
	case 938: GenerateProtonBeamthrough( anEvent ); break;

	case 700: GenerateE07Study( anEvent ); break; // Study on E07, generate beam and K+ from INC data
  case 701: GenerateE07StudyAll( anEvent ); break; // Study on E07, generate K+ pi+ from INC data
  case 702: GenerateE07StudyKnP( anEvent ); break; // Study on E07, generate K+ pi+ from INC data
  case 703: GenerateE07StudyKp( anEvent ); break; // Study on E07, generate K+ pi+ from INC data
  case 704: GenerateE07StudyKnPBeam( anEvent ); break; // Study on E07, generate K+ pi+ from INC data
  case 705: GenerateE07StudyKpBeam( anEvent ); break; // Study on E07, generate K+ pi+ from INC data, beam size 1x3 cm^2
  case 706: GenerateE07StudyKpXiBeam( anEvent ); break;
    // Study on E07, generate K+ from isotropic, beam size 1x3 cm^2
  case 707: GenerateE07StudyKpXiBeamOnlyKp( anEvent ); break;
    // Study on E07, generate K+ from isotropic, beam size 1x3 cm^2
  case 708: GenerateE07StudyP08to20( anEvent ); break; // generate proton
  case 709: GenerateE07StudyKp04to15( anEvent ); break; // generate proton
  case 710: GenerateE07StudyKpxi1530( anEvent ); break; // Study on E07, Xi(1530)- generator, beam size 1x3 cm^2
  case 711: GenerateE07StudyTakahashi( anEvent ); break; // Study on E07, Xi-kp, by using Takahashi-san's code
  case 2701: GenerateE27BeamThrough( anEvent ); break;
  case 2702: GenerateE27Kptest( anEvent ); break;
  case 2703: GenerateE27KppFLambdaP( anEvent ); break;
  case 2704: GenerateE27KppFSigmaZP( anEvent ); break;
  case 2705: GenerateE27KppFLambdaPizP( anEvent ); break;
  case 2706: GenerateE27KppFSigmaZPizP( anEvent ); break;
  case 2707: GenerateE27KppFSigmaPPimP( anEvent ); break;
  case 2708: GenerateE27K11BLambda10Be( anEvent ); break;
  case 2709: GenerateE27Kptest2( anEvent ); break;
  case 3001: GenerateKKppLL1( anEvent ); break;
  case 3002: GenerateKKppLL2( anEvent ); break;
  case 3003: GenerateKKppLSmPip( anEvent ); break;
  case 3004: GenerateKKppLSpPim( anEvent ); break;
  // case 3101: GenerateKKppJAMInput(anEvent,t1); break;
  // case 3102: GenerateKKppKKpp_BeamThrough1( anEvent ); break;
  // case 3103: GenerateKKppJAMInputK0(anEvent,t1); break;
  // case 3104: GenerateKKppJAMInputK0bar(anEvent,t1); break;
  case 4501: GenerateE45ElasticPionPlus( anEvent ); break;
  case 4502: GenerateE45ElasticPionMinus( anEvent ); break;
  case -1321: GenerateXiMinus( anEvent ); break;
  case 1115: GenerateLambda( anEvent ); break;
 	case 181321:GenerateKuramaPKmKpXi( anEvent );	break; 
	default:
    G4cerr << " * Generator number error : " << m_generator << G4endl;
    break;
  }
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateKpXi2Body( G4Event* anEvent )
{
  static const G4double KaonMinusMass = m_KaonMinus->GetPDGMass()/CLHEP::GeV;
  static const G4double ProtonMass = m_Proton->GetPDGMass()/CLHEP::GeV;
  static const G4double KaonPlusMass = m_KaonPlus->GetPDGMass()/CLHEP::GeV;
  static const G4double XiMinusMass = m_XiMinus->GetPDGMass()/CLHEP::GeV;
  static const G4int KaonMinusID = m_KaonMinus->GetPDGEncoding();
  static const G4int ProtonID = m_Proton->GetPDGEncoding();
  static const G4int KaonPlusID = m_KaonPlus->GetPDGEncoding();
  static const G4int XiMinusID = m_XiMinus->GetPDGEncoding();

  G4LorentzVector KnLV, PLV, KpLV, XiLV;
  G4ThreeVector gen_pos( m_beam->GetX( -m_beam->z ),
			 m_beam->GetY( -m_beam->z ),
			 m_target_pos.z() );
  // K+ angular distribution
  static const G4int NLegendre = 13;
  static const G4double coeff[NLegendre] = { 1. , -1.22 ,1.55 ,-1.08 ,
					     0.37, -0.15, 0.16, -0.38,
					     -0.18, 0.09, 0.05, -0.01, 0.2 };
  G4double cosx = -100.;
  G4double comp = -100.;
  G4bool theta_flag = false;
  while( !theta_flag ){
    theta_flag = true;
    while( true ){
      cosx = G4RandFlat::shoot( 0., 1. );
      comp = G4RandFlat::shoot( 0., 1. );
      // G4cout << "test of legendre function" << G4endl;
      G4double cross = 0.;
      for( G4int jj=0; jj<NLegendre; ++jj ){
	// G4cout << Kinematics::Legendre( jj, -1. ) << G4endl;
	cross += coeff[jj]*Kinematics::Legendre( jj, cosx );
      }
      cross /= 5.9;
      if( comp <= cross ){
	// G4cout << "cross:comp --> " << cross << " : " << comp << G4endl;
	break;
      } else if( cross > 1. ){
	G4cout << "#W " << FUNC_NAME << " "
	       << "cross is larger than normalization" << G4endl;
      }
    }
    if( cosx < -2 )
      continue;
    ////////////// harmonic motion
    G4int type = G4RandFlat::shoot() * 6. < 2. ? 0 : 1 ;
    G4ThreeVector p_fermi = Kinematics::HarmonicFermiMomentum( type );
    gAnaMan.SetFermiMomentum( p_fermi );
    KinemaFermi kinema( KaonMinusMass, ProtonMass, KaonPlusMass, XiMinusMass,
			m_beam->p, p_fermi, cosx );
    // if( kinema.GetTheta( 3 )/CLHEP::deg > 60. ) //   theta_flag = false;
    KnLV = kinema.GetLorentzVector( 0 );
    PLV = kinema.GetLorentzVector( 1 );
    KpLV = kinema.GetLorentzVector( 2 );
    XiLV = kinema.GetLorentzVector( 3 );
		gTrackBuffer.SetLV(XiLV,0);
#ifdef DEBUG
    kinema.Print();
#endif
  }

  // G4cout << "#D " << FUNC_NAME << " k+ mom: " << mom[2] << G4endl;
  // m_cross_section->Fill(cosx);
  // gAnaMan.FillCrossSection(cosx);
  //  G4double momkpp=std::sqrt(std::pow(kp_mom_x,2)+std::pow(kp_mom_y,2)+std::pow(kp_mom_z,2));
  // G4double momcmk[4]={0};
  //  G4double beta= pbeam/(proton->GetPDGMass()/CLHEP::GeV+Ebeam);
  // G4double momcmkpp = std::sqrt( momcmk[0]*momcmk[0] +
  // 				 momcmk[1]*momcmk[1] +
  // 				 momcmk[2]*momcmk[2] );
  // m_cmk->Fill( std::acos(momcmk[2]/momcmkpp)/CLHEP::pi*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  // m_coscmk->Fill( momcmk[2]/momcmkpp );

  gAnaMan.SetPrimaryBeam( m_beam->p );
  gAnaMan.SetNumberOfPrimaryParticle( 3 );

  ///// KaonMinus (produced to upstream assumed with KaonPlus).
  // m_particle_gun->SetParticleDefinition( m_KaonPlus );
  // m_particle_gun->SetParticleMomentumDirection( -KnLV.v() );
  // m_particle_gun->SetParticleEnergy( ( KnLV.e() - KaonMinusMass )*CLHEP::GeV );
  // m_particle_gun->SetParticlePosition( gen_pos );
  // m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetPrimaryParticle( 0, -KnLV.v(), KaonMinusMass, KaonMinusID );
  gAnaMan.SetPrimaryVertex( 0, gen_pos );
  ///// Proton
  // m_particle_gun->SetParticleDefinition( m_Proton );
  // m_particle_gun->SetParticleMomentumDirection( PLV.v() );
  // m_particle_gun->SetParticleEnergy( ( PLV.e() - ProtonMass )*CLHEP::GeV );
  // m_particle_gun->SetParticlePosition( gen_pos );
  // m_particle_gun->GeneratePrimaryVertex( anEvent );
//  gAnaMan.SetPrimaryParticle( 1, PLV.v(), ProtonMass, ProtonID );
 // gAnaMan.SetPrimaryVertex( 1, gen_pos );
  ///// KaonPlus
  m_particle_gun->SetParticleDefinition( m_KaonPlus );
  m_particle_gun->SetParticleMomentumDirection( KpLV.v() );
  m_particle_gun->SetParticleEnergy( ( KpLV.e() - KaonPlusMass )*CLHEP::GeV );
  m_particle_gun->SetParticlePosition( gen_pos );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetPrimaryParticle( 1, KpLV.v(), KaonPlusMass, KaonPlusID );
  gAnaMan.SetPrimaryVertex( 1, gen_pos );
  ///// XiMinus
  m_particle_gun->SetParticleDefinition( m_XiMinus );
  m_particle_gun->SetParticleMomentumDirection( XiLV.v() );
  m_particle_gun->SetParticleEnergy( ( XiLV.e() - XiMinusMass )*CLHEP::GeV );
  m_particle_gun->SetParticlePosition( gen_pos );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetPrimaryParticle( 2, XiLV.v(), XiMinusMass, XiMinusID );
  gAnaMan.SetPrimaryVertex( 2, gen_pos );
}

//_____________________________________________________________________________
//case 150 Generate KpXi assuming uniform distribution on CM Frame
void
TPCPrimaryGeneratorAction::GenerateKpXi2BodyUniform( G4Event* anEvent )
{
  static const G4double KaonMinusMass = m_KaonMinus->GetPDGMass()/CLHEP::GeV;
  static const G4double ProtonMass = m_Proton->GetPDGMass()/CLHEP::GeV;
  static const G4double KaonPlusMass = m_KaonPlus->GetPDGMass()/CLHEP::GeV;
  static const G4double XiMinusMass = m_XiMinus->GetPDGMass()/CLHEP::GeV;
  static const G4int KaonMinusID = m_KaonMinus->GetPDGEncoding();
  static const G4int ProtonID = m_Proton->GetPDGEncoding();
  static const G4int KaonPlusID = m_KaonPlus->GetPDGEncoding();
  static const G4int XiMinusID = m_XiMinus->GetPDGEncoding();

  G4LorentzVector KnLV, PLV, KpLV, XiLV;
	double TargetDepth = G4RandFlat::shoot(-10.,10.); 
  G4ThreeVector gen_pos( m_beam->GetX( -m_beam->z ),
			 m_beam->GetY( -m_beam->z ),
			 m_target_pos.z() + TargetDepth);
	// K+ angular distribution
	G4double cosx = -100.;
	cosx = G4RandFlat::shoot( 0.5, 1. );//cos_CM= 0.83 -> theta_lab = 25 deg
	////////////// harmonic motion
	G4int type = G4RandFlat::shoot() * 6. < 2. ? 0 : 1 ;
	G4ThreeVector p_fermi = Kinematics::HarmonicFermiMomentum( type );
	p_fermi*=0;
	gAnaMan.SetFermiMomentum( p_fermi );
	KinemaFermi kinema( KaonMinusMass, ProtonMass, KaonPlusMass, XiMinusMass,
			m_beam->p, p_fermi, cosx );
	// if( kinema.GetTheta( 3 )/CLHEP::deg > 60. )
	//   theta_flag = false;
	KnLV = kinema.GetLorentzVector( 0 );
	PLV = kinema.GetLorentzVector( 1 );
	KpLV = kinema.GetLorentzVector( 2 );
	XiLV = kinema.GetLorentzVector( 3 );
#ifdef DEBUG
	kinema.Print();
#endif

  // G4cout << "#D " << FUNC_NAME << " k+ mom: " << mom[2] << G4endl;
  // m_cross_section->Fill(cosx);
  // gAnaMan.FillCrossSection(cosx);
  //  G4double momkpp=std::sqrt(std::pow(kp_mom_x,2)+std::pow(kp_mom_y,2)+std::pow(kp_mom_z,2));
  // G4double momcmk[4]={0};
  //  G4double beta= pbeam/(proton->GetPDGMass()/CLHEP::GeV+Ebeam);
  // G4double momcmkpp = std::sqrt( momcmk[0]*momcmk[0] +
  // 				 momcmk[1]*momcmk[1] +
  // 				 momcmk[2]*momcmk[2] );
  // m_cmk->Fill( std::acos(momcmk[2]/momcmkpp)/CLHEP::pi*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  // m_coscmk->Fill( momcmk[2]/momcmkpp );

  gAnaMan.SetPrimaryBeam( m_beam->p );
  gAnaMan.SetNumberOfPrimaryParticle( 3 );

  ///// KaonMinus (produced to upstream assumed with KaonPlus).
  m_particle_gun->SetParticleDefinition( m_KaonPlus );
	m_particle_gun->SetParticleMomentumDirection( -KnLV.v() );
	m_particle_gun->SetParticleEnergy( ( KnLV.e() - KaonMinusMass )*CLHEP::GeV );
	m_particle_gun->SetParticlePosition( gen_pos );
	m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetPrimaryParticle( 0, KnLV.v(), KaonMinusMass, KaonMinusID );
  gAnaMan.SetPrimaryVertex( 0, gen_pos );
  ///// Proton No Proton track
  ///// KaonPlus
  m_particle_gun->SetParticleDefinition( m_KaonPlus );
  m_particle_gun->SetParticleMomentumDirection( KpLV.v() );
  m_particle_gun->SetParticleEnergy( ( KpLV.e() - KaonPlusMass )*CLHEP::GeV );
  m_particle_gun->SetParticlePosition( gen_pos );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetPrimaryParticle( 1, KpLV.v(), KaonPlusMass, KaonPlusID );
  gAnaMan.SetPrimaryVertex( 1, gen_pos );
  ///// XiMinus
  m_particle_gun->SetParticleDefinition( m_XiMinus );
  m_particle_gun->SetParticleMomentumDirection( XiLV.v() );
  m_particle_gun->SetParticleEnergy( ( XiLV.e() - XiMinusMass )*CLHEP::GeV );
  m_particle_gun->SetParticlePosition( gen_pos );
	auto MomXi = XiLV.vect();
	auto MomKm = KnLV.vect();
	PLV = G4LorentzVector(PLV,ProtonMass);
	auto CMLV = KnLV + PLV;
//	G4cout<<CMLV.mag()<<G4endl;
	gTrackBuffer.SetCMLV(CMLV);
	auto SpinXi = MomXi.cross(MomKm);
	SpinXi = SpinXi*(1./SpinXi.mag());
	gTrackBuffer.SetPolarity(SpinXi,0);
	gTrackBuffer.SetMomentum(MomXi,0);
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetPrimaryParticle( 2, XiLV.v(), XiMinusMass, XiMinusID );
  gAnaMan.SetPrimaryVertex( 2, gen_pos );

}

//_____________________________________________________________________________
//case 151 Generate KpXi* assuming uniform distribution on CM Frame
void
TPCPrimaryGeneratorAction::GenerateKpXi1530Minus2BodyUniform( G4Event* anEvent )
{
  static const G4double KaonMinusMass = m_KaonMinus->GetPDGMass()/CLHEP::GeV;
  static const G4double ProtonMass = m_Proton->GetPDGMass()/CLHEP::GeV;
  static const G4double KaonPlusMass = m_KaonPlus->GetPDGMass()/CLHEP::GeV;
  static const G4double XiMinusMass = m_XiMinus->GetPDGMass()/CLHEP::GeV;
  static const G4double Xi1530MinusMass = m_Xi1530Minus->GetPDGMass()/CLHEP::GeV;
  static const G4int KaonMinusID = m_KaonMinus->GetPDGEncoding();
  static const G4int ProtonID = m_Proton->GetPDGEncoding();
  static const G4int KaonPlusID = m_KaonPlus->GetPDGEncoding();
  static const G4int Xi1530MinusID = m_Xi1530Minus->GetPDGEncoding();

  G4LorentzVector KnLV, PLV, KpLV, Xi1530MinusLV;
	double TargetDepth = G4RandFlat::shoot(-10.,10.); 
  G4ThreeVector gen_pos( m_beam->GetX( -m_beam->z ),
			 m_beam->GetY( -m_beam->z ),
			 m_target_pos.z() + TargetDepth);
	// K+ angular distribution
	G4double cosx = -100.;
	cosx = G4RandFlat::shoot( -0., 1. );//cos_CM= 0.83 -> theta_lab = 25 deg
	////////////// harmonic motion
	G4int type = G4RandFlat::shoot() * 6. < 2. ? 0 : 1 ;
	G4ThreeVector p_fermi = Kinematics::HarmonicFermiMomentum( type );
	p_fermi*=0;
	gAnaMan.SetFermiMomentum( p_fermi );
	KinemaFermi kinema( KaonMinusMass, ProtonMass, KaonPlusMass, Xi1530MinusMass,
			m_beam->p, p_fermi, cosx );
	// if( kinema.GetTheta( 3 )/CLHEP::deg > 60. )
	//   theta_flag = false;
	KnLV = kinema.GetLorentzVector( 0 );
	PLV = kinema.GetLorentzVector( 1 );
	KpLV = kinema.GetLorentzVector( 2 );
	Xi1530MinusLV = kinema.GetLorentzVector( 3 );
#ifdef DEBUG
	kinema.Print();
#endif

  // G4cout << "#D " << FUNC_NAME << " k+ mom: " << mom[2] << G4endl;
  // m_cross_section->Fill(cosx);
  // gAnaMan.FillCrossSection(cosx);
  //  G4double momkpp=std::sqrt(std::pow(kp_mom_x,2)+std::pow(kp_mom_y,2)+std::pow(kp_mom_z,2));
  // G4double momcmk[4]={0};
  //  G4double beta= pbeam/(proton->GetPDGMass()/CLHEP::GeV+Ebeam);
  // G4double momcmkpp = std::sqrt( momcmk[0]*momcmk[0] +
  // 				 momcmk[1]*momcmk[1] +
  // 				 momcmk[2]*momcmk[2] );
  // m_cmk->Fill( std::acos(momcmk[2]/momcmkpp)/CLHEP::pi*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  // m_coscmk->Fill( momcmk[2]/momcmkpp );

  gAnaMan.SetPrimaryBeam( m_beam->p );
  gAnaMan.SetNumberOfPrimaryParticle( 3 );

  ///// KaonMinus (produced to upstream assumed with KaonPlus).
  m_particle_gun->SetParticleDefinition( m_KaonPlus );
	m_particle_gun->SetParticleMomentumDirection( -KnLV.v() );
	m_particle_gun->SetParticleEnergy( ( KnLV.e() - KaonMinusMass )*CLHEP::GeV );
	m_particle_gun->SetParticlePosition( gen_pos );
	m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetPrimaryParticle( 0, KnLV.v(), KaonMinusMass, KaonMinusID );
  gAnaMan.SetPrimaryVertex( 0, gen_pos );
  ///// Proton No Proton track
  ///// KaonPlus
  m_particle_gun->SetParticleDefinition( m_KaonPlus );
  m_particle_gun->SetParticleMomentumDirection( KpLV.v() );
  m_particle_gun->SetParticleEnergy( ( KpLV.e() - KaonPlusMass )*CLHEP::GeV );
  m_particle_gun->SetParticlePosition( gen_pos );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetPrimaryParticle( 1, KpLV.v(), KaonPlusMass, KaonPlusID );
  gAnaMan.SetPrimaryVertex( 1, gen_pos );
  ///// Xi1530Minus
  m_particle_gun->SetParticleDefinition( m_Xi1530Minus );
  m_particle_gun->SetParticleMomentumDirection( Xi1530MinusLV.v() );
  m_particle_gun->SetParticleEnergy( ( Xi1530MinusLV.e() - Xi1530MinusMass )*CLHEP::GeV );
  m_particle_gun->SetParticlePosition( gen_pos );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetPrimaryParticle( 2, Xi1530MinusLV.v(), Xi1530MinusMass, Xi1530MinusID );
  gAnaMan.SetPrimaryVertex( 2, gen_pos );

}

//_____________________________________________________________________________
//case 16 based on GenerateAll
void
TPCPrimaryGeneratorAction::GenerateUniformProton( G4Event* anEvent )
{
  static const G4int ProtonID = m_Proton->GetPDGEncoding();
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //angle
  G4double phi=G4RandFlat::shoot(-1.,1.)*3.141592654;
  G4double mom_p=G4RandFlat::shoot(0.05,1.0);
  G4double theta=acos(G4RandFlat::shoot(0.,1.));
  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);
  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2))*GeV;

  G4double rn_vtx=-1;  G4double rn_vtz=-1;
  G4double vtx=0;  G4double vty=0;   G4double vtz=0;
  if( gConf.Get<G4int>("Experiment") == 45. ){
    while(1){
      rn_vtx = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
      rn_vtz = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
      if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_e45target_size.x()*m_e45target_size.x()) break;
    }
    vty = G4RandFlat::shoot(-m_e45target_size.z(),m_e45target_size.z());
    vtx=rn_vtx;
    vtz=rn_vtz+m_target_pos.z();
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-m_target_size.x()/2.,m_target_size.x()/2.)*mm;
    vty = G4RandFlat::shoot(-m_target_size.y()/2.,m_target_size.y()/2.)*mm;
    vtz = G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
  }

  //  G4double ratio=G4RandFlat::shoot(0.,1.);
  //  ratio = 0.1;
  // proton//
  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2));
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetModeID(1);
  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_Proton->GetPDGMass()/GeV,ProtonID);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);




}

//_____________________________________________________________________________
//case 17
void
TPCPrimaryGeneratorAction::GenerateUniformPim( G4Event* anEvent )
{
  G4double Energy_pi,  mom_pi_x, mom_pi_y, mom_pi_z;

  //angle
  G4double phi=G4RandFlat::shoot(-1.,1.)*3.141592654;
 // G4double mom_pi=G4RandFlat::shoot(0.05,1.0);
  G4double mom_pi=300;
  G4double theta=acos(G4RandFlat::shoot(0.,1.));
  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_pi_x = mom_pi*sin(theta)*cos(phi);
  mom_pi_y = mom_pi*sin(theta)*sin(phi);
  mom_pi_z = mom_pi*cos(theta);
  Energy_pi=sqrt(pow(mom_pi/GeV,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));


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
    vtz=rn_vtz+m_target_pos.z();
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-15.,15.)*mm;
    vty = G4RandFlat::shoot(-5.,5.)*mm;
    vtz = G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
  }



  //  G4double ratio=G4RandFlat::shoot(0.,1.);
  //  ratio = 0.1;

  // pi- //
  //Energy_p=sqrt(pow(mom_p,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi_x,mom_pi_y,mom_pi_z));
  m_particle_gun->SetParticleEnergy((Energy_pi - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetModeID(4);//?
  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_pi_x,mom_pi_y,mom_pi_z,m_PionMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);

}

//case 25
void
TPCPrimaryGeneratorAction::GenerateUniformKaonPlus( G4Event* anEvent )
{
  G4double Energy_pi,  mom_kp_x, mom_kp_y, mom_kp_z;

  //angle
  G4double phi=G4RandFlat::shoot(-1.,1.)*acos(-1);
  G4double mom_kp=G4RandFlat::shoot(0.05,1.5);
  G4double theta=acos(G4RandFlat::shoot(0.95,1.));
  mom_kp_x = mom_kp*sin(theta)*cos(phi);
  mom_kp_y = mom_kp*sin(theta)*sin(phi);
  mom_kp_z = mom_kp*cos(theta);
  Energy_pi=sqrt(pow(mom_kp,2)+pow(m_KaonPlus->GetPDGMass()/GeV,2))*GeV;

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
    vtz=rn_vtz+m_target_pos.z();
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-15.,15.)*mm;
    vty = G4RandFlat::shoot(-5.,5.)*mm;
    vtz = G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
  }

//	std::cout<<Form("KpMom:%g *(%g,%g,%g)",mom_kp,mom_kp_x/mom_kp,mom_kp_y/mom_kp,mom_kp_z/mom_kp)<<std::endl; 


  //  G4double ratio=G4RandFlat::shoot(0.,1.);
  //  ratio = 0.1;

  // pi- //
  //Energy_p=sqrt(pow(mom_p,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_pi - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetModeID(4);//?
  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);

}



//_____________________________________________________________________________
//case 18 based on GenerateAll
void
TPCPrimaryGeneratorAction::GenerateBeamProton( G4Event* anEvent )
{
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //angle
  G4double phi=G4RandFlat::shoot(-1.,1.)*acos(-1.);
  //  G4double mom_p=G4RandFlat::shoot(0.05,1.0);
  G4double mom_p=m_beam_p0;
  //G4double theta=acos(G4RandFlat::shoot(0.,1.));
  G4double theta=G4RandFlat::shoot(0.,3.)*acos(-1.)/180.;// 0-3 degree flat
  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);
  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2))*GeV;

  G4double rn_vtx=-1;  G4double rn_vtz=-1;
  G4double vtx=0;  G4double vty=0;   G4double vtz=0;
  if( gConf.Get<G4int>("Experiment") == 45. ){
    while(1){
      rn_vtx = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
      rn_vtz = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
      if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_e45target_size.x()*m_e45target_size.x()) break;
    }
    vty = G4RandFlat::shoot(-m_e45target_size.z(),m_e45target_size.z());
    vtx=rn_vtx;
    vtz=rn_vtz+m_target_pos.z();
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-m_target_size.x()/2.,m_target_size.x()/2.)*mm;
    vty = G4RandFlat::shoot(-m_target_size.y()/2.,m_target_size.y()/2.)*mm;
    vtz = G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
  }

  //  G4double ratio=G4RandFlat::shoot(0.,1.);
  //  ratio = 0.1;
  // proton//
  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2));
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetModeID(1);
  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);




}



//_____________________________________________________________________________
//case 19 based on GenerateAll
//beam mom is fixed
void
TPCPrimaryGeneratorAction::GenerateUniformProton_P( G4Event* anEvent )
{
  static const G4int ProtonID = m_Proton->GetPDGEncoding();
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //angle
  G4double phi=G4RandFlat::shoot(-1.,1.)*acos(-1.);
  //  G4double mom_p=G4RandFlat::shoot(0.05,1.0);
  G4double mom_p=m_beam_p0;
  //G4double theta=acos(G4RandFlat::shoot(0.,1.));
  //  G4double theta=G4RandFlat::shoot(0.,180.)*acos(-1.)/180.;// 0-180 degree flat
  G4double theta=acos(G4RandFlat::shoot(-1.,1.));// uniform generation in a sphere
  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);
  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2))*GeV;

  G4double rn_vtx=-1;  G4double rn_vtz=-1;
  G4double vtx=0;  G4double vty=0;   G4double vtz=0;
  if( gConf.Get<G4int>("Experiment") == 45. ){
    while(1){
      rn_vtx = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
      rn_vtz = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
      if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_e45target_size.x()*m_e45target_size.x()) break;
    }
    vty = G4RandFlat::shoot(-m_e45target_size.z(),m_e45target_size.z());
    vtx=rn_vtx;
    vtz=rn_vtz+m_target_pos.z();
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-m_target_size.x()/2.,m_target_size.x()/2.)*mm;
    vty = G4RandFlat::shoot(-m_target_size.y()/2.,m_target_size.y()/2.)*mm;
    vtz = G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
  }

  //  G4double ratio=G4RandFlat::shoot(0.,1.);
  //  ratio = 0.1;
  // proton//
  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2));
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetModeID(1);
  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_Proton->GetPDGMass()/GeV,ProtonID);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);




}

//_____________________________________________________________________________
//case 20 based on GenerateAll
//beam mom is fixed
//phi is fixed to be 0.
void
TPCPrimaryGeneratorAction::GenerateUniformProton_P_fixphi( G4Event* anEvent )
{
  static const G4int ProtonID = m_Proton->GetPDGEncoding();
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //angle
  //  G4double phi=G4RandFlat::shoot(-1.,1.)*acos(-1.);
  G4double phi=0.;
  //  G4double mom_p=G4RandFlat::shoot(0.05,1.0);
  G4double mom_p=m_beam_p0;
  //G4double theta=acos(G4RandFlat::shoot(0.,1.));
  G4double theta=G4RandFlat::shoot(0.,180.)*acos(-1.)/180.;// 0-180 degree flat
  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);
  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2))*GeV;

  G4double rn_vtx=-1;  G4double rn_vtz=-1;
  G4double vtx=0;  G4double vty=0;   G4double vtz=0;
  if( gConf.Get<G4int>("Experiment") == 45. ){
    while(1){
      rn_vtx = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
      rn_vtz = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
      if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_e45target_size.x()*m_e45target_size.x()) break;
    }
    vty = G4RandFlat::shoot(-m_e45target_size.z(),m_e45target_size.z());
    vtx=rn_vtx;
    vtz=rn_vtz+m_target_pos.z();
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-m_target_size.x()/2.,m_target_size.x()/2.)*mm;
    vty = G4RandFlat::shoot(-m_target_size.y()/2.,m_target_size.y()/2.)*mm;
    vtz = G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
  }

  //  G4double ratio=G4RandFlat::shoot(0.,1.);
  //  ratio = 0.1;
  // proton//
  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2));
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetModeID(1);
  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_Proton->GetPDGMass()/GeV,ProtonID);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);


}



//_____________________________________________________________________________
//case 50 based on GenerateUniformProton_P (case 19)
//beam mom is fixed
void
TPCPrimaryGeneratorAction::GenerateUniformProton_P_Multi( G4Event* anEvent )
{
  static const G4int ProtonID = m_Proton->GetPDGEncoding();
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  int n_event = 5;
  gAnaMan.SetNumberOfPrimaryParticle(n_event);
  for(int i=0; i<n_event; ++i){

    //angle
    G4double phi=G4RandFlat::shoot(-1.,1.)*acos(-1.);
    //  G4double mom_p=G4RandFlat::shoot(0.05,1.0);
    G4double mom_p=m_beam_p0;
    //G4double theta=acos(G4RandFlat::shoot(0.,1.));
    //  G4double theta=G4RandFlat::shoot(0.,180.)*acos(-1.)/180.;// 0-180 degree flat
    G4double theta=acos(G4RandFlat::shoot(-1.,1.));// uniform generation in a sphere
    //  G4cout<<theta*180/3.141592<<G4endl;
    mom_p_x = mom_p*sin(theta)*cos(phi);
    mom_p_y = mom_p*sin(theta)*sin(phi);
    mom_p_z = mom_p*cos(theta);
    Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2))*GeV;

    G4double rn_vtx=-1;  G4double rn_vtz=-1;
    G4double vtx=0;  G4double vty=0;   G4double vtz=0;
    if( gConf.Get<G4int>("Experiment") == 45. ){
      while(1){
	rn_vtx = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
	rn_vtz = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
	if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_e45target_size.x()*m_e45target_size.x()) break;
      }
      vty = G4RandFlat::shoot(-m_e45target_size.z(),m_e45target_size.z());
      vtx=rn_vtx;
      vtz=rn_vtz+m_target_pos.z();
    }else if( gConf.Get<G4int>("Experiment") == 42. ){
      vtx = G4RandFlat::shoot(-m_target_size.x()/2.,m_target_size.x()/2.)*mm;
      vty = G4RandFlat::shoot(-m_target_size.y()/2.,m_target_size.y()/2.)*mm;
      vtz = G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
    }

    //  G4double ratio=G4RandFlat::shoot(0.,1.);
    //  ratio = 0.1;
    // proton//
    Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2));
    m_particle_gun->SetParticleDefinition(m_Proton);
    m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
    m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    m_particle_gun->GeneratePrimaryVertex( anEvent );
    gAnaMan.SetModeID(1);
    //  gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(i,mom_p_x,mom_p_y,mom_p_z,m_Proton->GetPDGMass()/GeV,ProtonID);
    gAnaMan.SetPrimaryVertex(i,vtx,vty,vtz);

  }


}


//_____________________________________________________________________________
//case 34 LL production from K-"N" -> Xi K+, XiN->LL  
void
TPCPrimaryGeneratorAction::GenerateLL_fromXiP( G4Event* anEvent )
{
  static const G4double KaonMinusMass = m_KaonMinus->GetPDGMass()/CLHEP::GeV;
  static const G4double ProtonMass = m_Proton->GetPDGMass()/CLHEP::GeV;
  static const G4double KaonPlusMass = m_KaonPlus->GetPDGMass()/CLHEP::GeV;
  static const G4double XiMinusMass = m_XiMinus->GetPDGMass()/CLHEP::GeV;
  static const G4double LambdaMass = m_Lambda->GetPDGMass()/CLHEP::GeV;
  static const G4int KaonMinusID = m_KaonMinus->GetPDGEncoding();
  static const G4int ProtonID = m_Proton->GetPDGEncoding();
  static const G4int KaonPlusID = m_KaonPlus->GetPDGEncoding();
  static const G4int XiMinusID = m_XiMinus->GetPDGEncoding();
  static const G4int LambdaID = m_Lambda->GetPDGEncoding();

  G4double beam_p=m_beam_p0;
  double beam_E = sqrt(KaonMinusMass*KaonMinusMass + beam_p*beam_p);
  double masses1[2] = {KaonPlusMass, XiMinusMass};
  double masses2[2] = {LambdaMass, LambdaMass};
  TLorentzVector beam_L(0.0, 0.0, beam_p , beam_E );
  double proton_Fermi_p=0.;
  while(1){
    proton_Fermi_p= G4RandGauss::shoot(0.145,0.052291);//Tmp Fermimom of 12C
    if(proton_Fermi_p>0.)
      break;
  }
  G4double phi=G4RandFlat::shoot(-1.,1.)*acos(-1.);
  G4double theta=acos(G4RandFlat::shoot(-1.,1.));// uniform generation in a sphere
  
  G4double proton_Fermi_px = proton_Fermi_p*sin(theta)*cos(phi);
  G4double proton_Fermi_py = proton_Fermi_p*sin(theta)*sin(phi);
  G4double proton_Fermi_pz = proton_Fermi_p*cos(theta);
 
  G4double proton_Fermi_E = sqrt(proton_Fermi_p*proton_Fermi_p + ProtonMass);
  TLorentzVector target(proton_Fermi_px, proton_Fermi_py, 
			proton_Fermi_pz, proton_Fermi_E);

 
  double proton2_Fermi_p=0.;
  while(1){
    proton2_Fermi_p= G4RandGauss::shoot(0.145,0.052291);//Tmp Fermimom of 12C
    if(proton2_Fermi_p>0.)
      break;
  }
  G4double phi2=G4RandFlat::shoot(-1.,1.)*acos(-1.);
  G4double theta2=acos(G4RandFlat::shoot(-1.,1.));// uniform generation in a sphere
  
  G4double proton2_Fermi_px = proton2_Fermi_p*sin(theta)*cos(phi);
  G4double proton2_Fermi_py = proton2_Fermi_p*sin(theta)*sin(phi);
  G4double proton2_Fermi_pz = proton2_Fermi_p*cos(theta);
 
  G4double proton2_Fermi_E = sqrt(proton2_Fermi_p*proton2_Fermi_p + ProtonMass);
  TLorentzVector target2(proton2_Fermi_px, proton2_Fermi_py, 
			 proton2_Fermi_pz, proton2_Fermi_E);

  TLorentzVector W1 = beam_L + target;
  if(W1.Mag()<masses1[0]+masses1[1])
    G4cerr << " Error (PrimaryGeneratorAction 34): W1 < masses1"<< G4endl;
  

  TGenPhaseSpace event1;
  event1.SetDecay(W1, 2, masses1);

  TLorentzVector *kp_L = new TLorentzVector();
  TLorentzVector *XiM_L = new TLorentzVector();
  TLorentzVector *L1_L = new TLorentzVector();
  TLorentzVector *L2_L = new TLorentzVector();

  G4double weight1 = 0.;
  while(1){
    weight1 = event1.Generate();
    kp_L = event1.GetDecay(0);
    double kp_theta = kp_L->Theta()*180./acos(-1);
    if(kp_theta<20.)
      break;
  }
  XiM_L = event1.GetDecay(1);
  
  TLorentzVector W2 = *XiM_L + target2;
  if(W2.Mag()<masses2[0]+masses2[1])
    G4cerr << " Error (PrimaryGeneratorAction 34): W2 < masses2"<< G4endl;
  
  TGenPhaseSpace event2;
  event2.SetDecay(W2, 2, masses2);
  G4double weight2 = event2.Generate();
  L1_L = event2.GetDecay(0);
  L2_L = event2.GetDecay(1);

  

  G4double rn_vtx=-1;  G4double rn_vtz=-1;
  G4double vtx=0;  G4double vty=0;   G4double vtz=0;
  if( gConf.Get<G4int>("Experiment") == 45. ){
    while(1){
      rn_vtx = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
      rn_vtz = G4RandFlat::shoot(-m_e45target_size.x(),m_e45target_size.x());
      if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_e45target_size.x()*m_e45target_size.x()) break;
    }
    vty = G4RandFlat::shoot(-m_e45target_size.z(),m_e45target_size.z());
    vtx=rn_vtx;
    vtz=rn_vtz+m_target_pos.z();
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-m_target_size.x()/2.,m_target_size.x()/2.)*mm;
    vty = G4RandFlat::shoot(-m_target_size.y()/2.,m_target_size.y()/2.)*mm;
    vtz = G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
  }
  
  gAnaMan.SetNumberOfPrimaryParticle(3);

  //K+ gun, (skaon+ is injected)
  m_particle_gun->SetParticleDefinition(m_sKaonPlus);
  //  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  //  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
  std::cout<<"Kp Mom: "<<kp_L->P()<<std::endl;
  m_particle_gun->SetParticleMomentum(G4ThreeVector(kp_L->Px()*GeV,
						    kp_L->Py()*GeV,
						    kp_L->Pz()*GeV));
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetModeID(1);
  gAnaMan.SetPrimaryParticle(0,kp_L->Px(),kp_L->Py(),kp_L->Pz()
			     ,m_Proton->GetPDGMass()/GeV,KaonPlusID);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);

  
  //Lambda1 gun
  m_particle_gun->SetParticleDefinition(m_Lambda);
  std::cout<<"L1 Mom: "<<L1_L->P()<<std::endl;
  m_particle_gun->SetParticleMomentum(G4ThreeVector(L1_L->Px()*GeV,
						    L1_L->Py()*GeV,
						    L1_L->Pz()*GeV));
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetModeID(1);
  gAnaMan.SetPrimaryParticle(1,L1_L->Px(),L1_L->Py(),L1_L->Pz()
			     ,m_Lambda->GetPDGMass()/GeV,LambdaID);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);

    
  //Lambda2 gun
  m_particle_gun->SetParticleDefinition(m_Lambda);
  std::cout<<"L2 Mom: "<<L2_L->P()<<std::endl;
  m_particle_gun->SetParticleMomentum(G4ThreeVector(L2_L->Px()*GeV,
						    L2_L->Py()*GeV,
						    L2_L->Pz()*GeV));
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetModeID(1);
  gAnaMan.SetPrimaryParticle(2,L2_L->Px(),L2_L->Py(),L2_L->Pz()
			     ,m_Lambda->GetPDGMass()/GeV,LambdaID);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}




//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHdibaryon2( G4Event* anEvent )
{
  G4double mass_hdibaryon, width_hdibaryon;
  //  G4double Energy_h, momentum_h, mom_h_x, mom_h_y, mom_h_z;

  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  // G4double rmk=0.493677;
  //  Ebeam = 1.9;       // Incident gamma energy (GeV)
  pg_x = 0.0;
  pg_y = 0.0;
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548)*GeV;
  //  pg_z = m_beam_p0;
  pg_z = G4RandGauss::shoot( m_beam_p0, 0.01294*m_beam_p0 );
  //  G4cout<<"Ebeam:"<<Ebeam<<G4endl;
  // G4double pbeam=sqrt(pow(pg_x,2)+pow(pg_y,2)+pow(pg_z,2));
  // G4double Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/GeV*m_KaonMinus->GetPDGMass()/GeV);       // Incident gamma energy (GeV)

  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);
  mass_hdibaryon = gConf.Get<G4double>( "HdibaryonMass" );
  width_hdibaryon = gConf.Get<G4double>( "HdibaryonWidth" );

 up:
  Kinema3Resonance Hkinema( m_KaonMinus->GetPDGMass()/GeV,
			    ((m_Proton->GetPDGMass()/GeV)*2),
			    m_Lambda->GetPDGMass()/GeV,
			    m_Lambda->GetPDGMass()/GeV,
			    m_KaonPlus->GetPDGMass()/GeV,
			    mass_hdibaryon, width_hdibaryon, pg_z, 0.0);

  Energy_kp = Hkinema.GetEnergy(5);
  // G4double momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  //  if(atan((mom[1])/mom[2])*180/3.141592654 > 15. ) goto up;
  if( fabs(atan2(mom[1],mom[2])*180/3.141592654) > 15. || fabs(atan2(mom[0],mom[2])*180/3.141592654) > 20. ) goto up;
  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;

  /* L1 */
  Energy_L1 = Hkinema.GetEnergy(3);
  // G4double momentum_L1 = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];

  // G4double ThetaL1 = Hkinema.GetTheta(3);
  // G4double PhiL1 = Hkinema.GetPhi(3);

  /* L2 */
  Energy_L2 = Hkinema.GetEnergy(4);
  // G4double momentum_L2 = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double ThetaL2 = Hkinema.GetTheta(4);
  // G4double PhiL2 = Hkinema.GetPhi(4);

  /* Kp */
  Energy_kp = Hkinema.GetEnergy(5);
  // G4double momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  //  Thetakp = Hkinema.GetThetaCM(5);
  // G4double Thetakp = Hkinema.GetThetaCM(1);
  //  Phikp = Hkinema.GetPhi(5);
  //  G4cout<<Phikp<<G4endl;
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4cout<<"phi test: "<<Phikp<<":"<<shphi<<G4endl;
  //  G4cout<<"test: "<<Hkinema.GetPhiCM(1)<<":"<<Hkinema.GetPhi(5)<<G4endl;
  // Vertex (Reaction) point
  //  G4cout<<Energy_kp<<G4endl;
  //  //  G4double beta= sqrt(1.8*1.8-0.493677*0.493677)/(2*0.93827203+Energy_kp);
  //  G4cout<<(m_Proton->GetPDGMass()/GeV)<<G4endl;
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
  //  hlab[3]=sqrt(pow(hlab[0],2)+pow(hlab[1],2)+pow(hlab[2],2)+pow(m_Lambda->GetPDGMass()/GeV,2)*2);
  //  mass_hdibaryon = 2.250;    // Mass for H
  hlab[3]=sqrt(hlab[0]*hlab[0]+hlab[1]*hlab[1]+hlab[2]*hlab[2]+mass_hdibaryon*mass_hdibaryon);
  // G4double hcm[4];

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
  G4double vtz= G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,
				  m_target_pos.z()+m_target_size.z()/2)*mm;

  ///////kp PG
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////L1 PG
  m_particle_gun->SetParticleDefinition(m_Lambda);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  m_particle_gun->SetParticleEnergy((Energy_L1 - m_Lambda->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////L2 PG
  m_particle_gun->SetParticleDefinition(m_Lambda);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  m_particle_gun->SetParticleEnergy((Energy_L2 - m_Lambda->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,m_Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,m_Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHanul( G4Event* anEvent )
{

  double pbm[4]={-9999.9999}, pka[4]={-9999.9999}, vtx[3]={-9999.9999};
  double pL1[4]={-9999.9999}, pL2[4]={-9999.9999};
				//,xL1[3]={-9999.9999}, xL2[3]={-9999.9999};
  //  double pp1[4]={-9999.9999}, ppi1[4]={-9999.9999}, pp2[4]={-9999.9999}, ppi2[4]={-9999.9999};

  char fnam[30];

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

  // ---- K+ -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(momentumKaonPlus);
  m_particle_gun->SetParticleEnergy((Energy_ka - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(vertexPos);
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  // ---- Lambda 1 -------------
  G4ThreeVector momentumLambda1(pL1[0], pL1[1], pL1[2]);
  m_particle_gun->SetParticleDefinition( m_Lambda );
  m_particle_gun->SetParticleMomentumDirection(momentumLambda1);
  m_particle_gun->SetParticleEnergy((Energy_L1 - m_Lambda->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(vertexPos);
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  // ---- Lambda 2 -------------
  G4ThreeVector momentumLambda2(pL2[0], pL2[1], pL2[2]);
  m_particle_gun->SetParticleDefinition( m_Lambda );
  m_particle_gun->SetParticleMomentumDirection(momentumLambda2);
  m_particle_gun->SetParticleEnergy((Energy_L2 - m_Lambda->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(vertexPos);
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,pka[0],pka[1],pka[2],m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,pL1[0],pL1[1],pL1[2],m_Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,pL2[0],pL2[1],pL2[2],m_Lambda->GetPDGMass()/GeV);

  gAnaMan.SetPrimaryVertex(0,0.,0.,vtx[2]);
  gAnaMan.SetPrimaryVertex(1,0.,0.,vtx[2]);
  gAnaMan.SetPrimaryVertex(2,0.,0.,vtx[2]);
}

//_____________________________________________________________________________
// generator #30
void
TPCPrimaryGeneratorAction::GeneratePhaseSpace( G4Event* anEvent )
{
  G4double pg_x,pg_y,pg_z;
  G4double pbeam = G4RandGauss::shoot( m_beam_p0,
				       m_beam_p0*3.3*0.0001/2.3548 )*GeV;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;
  gAnaMan.SetPrimaryBeam( pg_x, pg_y, pg_z );

  // G4double Ebeam = sqrt(pow(pbeam/GeV,2)+pow(m_KaonMinus->GetPDGMass()/GeV,2));
  // G4double rmC12=Carbon12->GetPDGMass()/GeV;
  // G4double rmkp=m_KaonPlus->GetPDGMass()/GeV;
  // G4double rmBe10=Beryllium10->GetPDGMass()/GeV;
  // G4double rmL=m_Lambda->GetPDGMass()/GeV;

  // G4double W=sqrt(pow(Ebeam+rmC12,2)-pow(pbeam/GeV,2));
  //  G4cout<<"C mass------->"<<Carbon12->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"Be mass------->"<<Beryllium10->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"W------->"<<W<<G4endl;
  //  G4cout<<"Ebeam------->"<<Ebeam<<G4endl;
  //  G4cout<<"pbeam------->"<<pbeam<<G4endl;
  G4double Energy_LLphase = std::sqrt( std::pow( m_LLphase->GetPDGMass()/GeV, 2 ) +
				       std::pow( pbeam/GeV, 2 ) );
  G4double vtx,vty,vtz;
  G4double rn_vtx,rn_vtz;
  while( true ){
    rn_vtx = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    rn_vtz = G4RandFlat::shoot( -m_target_size.x(), m_target_size.x() );
    if( (rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x() )
      break;
  }
  vty = G4RandFlat::shoot( -m_target_size.z(), m_target_size.z() );
  vtx = rn_vtx;
  vtz = rn_vtz+m_target_pos.z();

  m_particle_gun->SetParticleDefinition(m_LLphase);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  m_particle_gun->SetParticleEnergy((Energy_LLphase - m_LLphase->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,pg_x,pg_y,pg_z, m_LLphase->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);

  /*
  Kinema4Body PhaseLL(W,0.,
		      Beryllium10->GetPDGMass()/GeV,
		      m_Lambda->GetPDGMass()/GeV,
		      m_Lambda->GetPDGMass()/GeV,
		      m_KaonPlus->GetPDGMass()/GeV,
		      pbeam, 0.0);
  */
  /*
//  Hkinema.Dump();
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


  ///////Neutron PG
  m_particle_gun->SetParticleDefinition(m_Neutron);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  m_particle_gun->SetParticleEnergy((Energy_pro - m_Neutron->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi1 PG
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  m_particle_gun->SetParticleEnergy((Energy_pi1 - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi2 PG
  m_particle_gun->SetParticleDefinition(m_PionPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  m_particle_gun->SetParticleEnergy((Energy_pi2 - m_PionPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z, m_Neutron->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z, m_PionMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z, m_PionPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

*/



  /*
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/GeV*m_KaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/GeV,
			(m_Proton->GetPDGMass()/GeV)*2,
			m_Hdibaryon->GetPDGMass()/GeV,
			m_KaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];


  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  Theta_h = Hkinema.GetTheta(3);
  Phi_h = Hkinema.GetPhi(3);

  ////check kinematics
  momk[0]=mom_kp_x;
  momk[1]=mom_kp_y;
  momk[2]=mom_kp_z;
  momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4cout<<"beta:"<<beta<<G4endl;
  G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  coscmk->Fill(momcmk[2]/momcmkpp);


  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;


  //  G4cout<<"env_target_width :"<<atof(env_env_target_width.c_str())<<G4endl;
  G4double vtz= G4RandFlat::shoot(m_target_pos.z()-env_target_width/2,m_target_pos.z()+env_target_width/2)*mm;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);

  //Kaon +
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //hdibaryon
  m_particle_gun->SetParticleDefinition(m_Hdibaryon);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_Hdibaryon->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_Hdibaryon->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  */
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHdibaryonPHSG( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/GeV*m_KaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/GeV,
			(m_Proton->GetPDGMass()/GeV)*2,
			m_Hdibaryon->GetPDGMass()/GeV,
			m_KaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hkinema.GetTheta(3);
  // G4double Phi_h = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);


  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz = G4RandFlat::shoot( m_target_pos.z() - m_target_size.z()/2,
					 m_target_pos.z() + m_target_size.z()/2 )*mm;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);

  //Kaon +
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //Hdibaryon
  m_particle_gun->SetParticleDefinition(m_Hdibaryon);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_Hdibaryon->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_Hdibaryon->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHdibaryonPHSGS( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/GeV*m_KaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"m_Hdibaryon mass"<<m_Hdibaryon->GetPDGMass()/GeV<<G4endl;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/GeV,
		      (m_Proton->GetPDGMass()/GeV)*2,
		      m_HdibaryonS->GetPDGMass()/GeV,
		      m_KaonPlus->GetPDGMass()/GeV,
		      pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hkinema.GetTheta(3);
  // G4double Phi_h = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);


  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  //  G4cout<<"env_target_width :"<<atof(env_env_target_width.c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);
  G4double vtz= G4RandFlat::shoot( m_target_pos.z() - m_target_size.z()/2,
					m_target_pos.z() + m_target_size.z()/2 )*mm;

  //Kaon +
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //Hdibaryon
  m_particle_gun->SetParticleDefinition(m_HdibaryonS);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_HdibaryonS->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_HdibaryonS->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHdibaryonPHSGLL( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/GeV*m_KaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/GeV,
			(m_Proton->GetPDGMass()/GeV)*2,
			m_HdibaryonLL->GetPDGMass()/GeV,
			m_KaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];


  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.  ) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hkinema.GetTheta(3);
  // G4double Phi_h = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz = G4RandFlat::shoot( m_target_pos.z() - m_target_size.z()/2,
					 m_target_pos.z() + m_target_size.z()/2 )*mm;

  //Kaon +
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //Hdibaryon
  m_particle_gun->SetParticleDefinition(m_HdibaryonLL);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_HdibaryonLL->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_HdibaryonLL->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateKpKn( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double Energy_kn, mom_kn_x, mom_kn_y, mom_kn_z;
  pbeam=1.8;
  mom_kn_x=0;
  mom_kn_y=0;
  mom_kn_z=pbeam;
  Energy_kn=sqrt(m_KaonMinus->GetPDGMass()/GeV*m_KaonMinus->GetPDGMass()/GeV+pbeam*pbeam);

  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/GeV*m_KaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/GeV,
			(m_Proton->GetPDGMass()/GeV)*2,
			m_HdibaryonLL->GetPDGMass()/GeV,
			m_KaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);


  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.  ) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  // G4double Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  // G4double mom_h_x = mom[0];
  // G4double mom_h_y = mom[1];
  // G4double mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hkinema.GetTheta(3);
  // G4double Phi_h = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);

  //  G4double vtx = 0.*mm;
  //  G4double vty = 0.*mm;
  G4double vtx = G4RandFlat::shoot(-15.,15.)*mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*mm;
  G4double vtz = G4RandFlat::shoot( m_target_pos.z() - m_target_size.z()/2,
					 m_target_pos.z() + m_target_size.z()/2 )*mm;

  //Kaon -
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //beam K+
  m_particle_gun->SetParticleDefinition(m_KaonMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kn_x,mom_kn_y,mom_kn_z));
  m_particle_gun->SetParticleEnergy((Energy_kn - m_KaonMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz-150.*mm));
  m_particle_gun->GeneratePrimaryVertex( anEvent );



  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_kn_x,mom_kn_y,mom_kn_z,m_KaonMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateBeamVI( G4Event* anEvent )
{
  static const G4double mass = m_KaonMinus->GetPDGMass()/GeV;
  static const G4int pid = m_KaonMinus->GetPDGEncoding();
  static const G4double D4BendAngle = gSize.Get( "D4BendAngle" )*CLHEP::deg;
  m_beam->p.rotateY( - D4BendAngle );
  G4double energy = ( std::sqrt( mass*mass + m_beam->p.mag2() ) - mass )*GeV;
  gAnaMan.SetPrimaryBeam( m_beam->p );
  G4ThreeVector gen_pos( m_beam->x, m_beam->y, m_beam->z );
  gen_pos.rotateY( - D4BendAngle );
  gen_pos += gBeam.GetVIPosition();
  m_particle_gun->SetParticleDefinition( m_KaonMinus );
  m_particle_gun->SetParticleMomentumDirection( m_beam->p );
  m_particle_gun->SetParticleEnergy( energy );
  m_particle_gun->SetParticlePosition( gen_pos );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetNumberOfPrimaryParticle( 1 );
  gAnaMan.SetPrimaryParticle( 0, m_beam->p, mass, pid );
  gAnaMan.SetPrimaryVertex( 0, gen_pos );
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateBeamVO( G4Event* anEvent )
{
  static const G4double mass = m_KaonMinus->GetPDGMass()/GeV;
  static const G4int pid = m_KaonMinus->GetPDGEncoding();
  G4double energy = ( std::sqrt( mass*mass + m_beam->p.mag2() ) - mass )*GeV;
  gAnaMan.SetPrimaryBeam( m_beam->p );
  //std::cout<<"beam mom="<<m_beam->p.mag()<<std::endl;
  G4ThreeVector gen_pos( m_beam->x, m_beam->y, m_target_pos.z() + m_beam->z );
  m_particle_gun->SetParticleDefinition( m_KaonMinus );
  m_particle_gun->SetParticleMomentumDirection( m_beam->p );
  m_particle_gun->SetParticleEnergy( energy );
  m_particle_gun->SetParticlePosition( gen_pos );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  gAnaMan.SetNumberOfPrimaryParticle( 1 );
  gAnaMan.SetPrimaryParticle( 0, m_beam->p, mass, pid );
  gAnaMan.SetPrimaryVertex( 0, gen_pos );
}

//_____________________________________________________________________________
// Generator# 31
void
TPCPrimaryGeneratorAction::GenerateJamInput( G4Event* anEvent )
{
  static const auto particleTable = G4ParticleTable::GetParticleTable();
  if( !m_jam )
    return;
  gAnaMan.SetNumberOfPrimaryParticle( m_jam->np );
  for( G4int i=0; i<m_jam->np; ++i ){
    auto particle = particleTable->FindParticle( m_jam->pid[i] );
    //G4ThreeVector x( 0., 0., gGeom.GetGlobalPosition( "Target" ).z() );
    G4ThreeVector x( m_beam->GetX( -m_beam->z ),
		     m_beam->GetY( -m_beam->z ),
		     G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2, m_target_pos.z()+m_target_size.z()/2)*mm);
    G4ThreeVector p( m_jam->px[i]*GeV, m_jam->py[i]*GeV, m_jam->pz[i]*GeV );
    double p_u = p.x()/p.z();
    double p_v = p.y()/p.z();
    double beam_u = m_beam->p.x()/m_beam->p.z();
    double beam_v = m_beam->p.y()/m_beam->p.z();
    double p_u_cor = p_u + beam_u;
    double p_v_cor = p_v + beam_v;
    double p_z_cor = p.mag()/std::sqrt( p_u_cor*p_u_cor + p_v_cor*p_v_cor + 1. );
    G4ThreeVector p_cor( p_z_cor*p_u_cor, p_z_cor*p_v_cor, p_z_cor);
    // std::cout<<"u="<<p_u<<", v="<<p_v
    // 	     <<", b_u="<<beam_u<<", b_v="<<beam_v
    // 	     <<", p="<<p.mag()<<std::endl;
    //    double theta = p_cor.theta()*180./acos(-1.);
    // if(theta>30.)
    //   continue;

    m_particle_gun->SetParticleDefinition( particle );
    //    m_particle_gun->SetParticleMomentumDirection( p );
    m_particle_gun->SetParticleMomentumDirection( p_cor );
    G4double m = particle->GetPDGMass()/GeV;
    //    G4double ke = std::sqrt( m*m + p.mag2() ) - m;
    G4double ke = std::sqrt( m*m + p_cor.mag2() ) - m;
    m_particle_gun->SetParticleEnergy( ke );
    m_particle_gun->SetParticlePosition( x );
    m_particle_gun->GeneratePrimaryVertex( anEvent );
    gAnaMan.SetPrimaryParticle( i, p_cor, m, m_jam->pid[i] );
    gAnaMan.SetPrimaryVertex( i, x );
  }
}

// Generator# 33
void
TPCPrimaryGeneratorAction::GenerateJamInput_Randphi( G4Event* anEvent )
{
  static const auto particleTable = G4ParticleTable::GetParticleTable();
  if( !m_jam )
    return;
  gAnaMan.SetNumberOfPrimaryParticle( m_jam->np );
  for( G4int i=0; i<m_jam->np; ++i ){
    auto particle = particleTable->FindParticle( m_jam->pid[i] );
    //G4ThreeVector x( 0., 0., gGeom.GetGlobalPosition( "Target" ).z() );
    G4ThreeVector x( m_beam->GetX( -m_beam->z ),
		     m_beam->GetY( -m_beam->z ),
		     G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2, m_target_pos.z()+m_target_size.z()/2)*mm);
    G4ThreeVector p( m_jam->px[i]*GeV, m_jam->py[i]*GeV, m_jam->pz[i]*GeV );
    double theta = p.theta()*180./acos(-1.);
    if(theta>30.)
      continue;
    double phi_rand = p.phi()+G4RandFlat::shoot(0., 2.*acos(-1.));
    double p_x_rand = p.mag()*sin(p.theta())*cos(phi_rand);
    double p_y_rand = p.mag()*sin(p.theta())*sin(phi_rand);

    double p_u = p_x_rand/p.z();
    double p_v = p_y_rand/p.z();
    double beam_u = m_beam->p.x()/m_beam->p.z();
    double beam_v = m_beam->p.y()/m_beam->p.z();
    double p_u_cor = p_u + beam_u;
    double p_v_cor = p_v + beam_v;
    double p_z_cor = p.mag()/std::sqrt( p_u_cor*p_u_cor + p_v_cor*p_v_cor + 1. );
    G4ThreeVector p_cor( p_z_cor*p_u_cor, p_z_cor*p_v_cor, p_z_cor);
    // std::cout<<"u="<<p_u<<", v="<<p_v
    // 	     <<", b_u="<<beam_u<<", b_v="<<beam_v
    // 	     <<", p="<<p.mag()<<std::endl;




    m_particle_gun->SetParticleDefinition( particle );
    //    m_particle_gun->SetParticleMomentumDirection( p );
    m_particle_gun->SetParticleMomentumDirection( p_cor );
    G4double m = particle->GetPDGMass()/GeV;
    //    G4double ke = std::sqrt( m*m + p.mag2() ) - m;
    G4double ke = std::sqrt( m*m + p_cor.mag2() ) - m;
    m_particle_gun->SetParticleEnergy( ke );
    m_particle_gun->SetParticlePosition( x );
    m_particle_gun->GeneratePrimaryVertex( anEvent );
    gAnaMan.SetPrimaryParticle( i, p_cor, m, m_jam->pid[i] );
    gAnaMan.SetPrimaryVertex( i, x );
  }
}


//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateIncInput( G4Event* anEvent )
{
  static const auto particleTable = G4ParticleTable::GetParticleTable();
  if( !m_inc )
    return;
  gAnaMan.SetIncID( m_inc->ich );
  gAnaMan.SetPrimaryBeam( m_inc->bpx, m_inc->bpy, m_inc->bpz );
  gAnaMan.SetNumberOfPrimaryParticle( m_inc->np );
  for( G4int i=0; i<m_inc->np; ++i ){
    auto particle = particleTable->FindParticle( m_inc->pid[i] );
    G4ThreeVector x( m_target_pos.x(), m_target_pos.y(), m_target_pos.z() );
    G4ThreeVector p( m_inc->px[i]*GeV, m_inc->py[i]*GeV, m_inc->pz[i]*GeV );
    m_particle_gun->SetParticleDefinition( particle );
    m_particle_gun->SetParticleMomentumDirection( p );
    G4double m = particle->GetPDGMass()/GeV;
    G4double ke = std::sqrt( m*m + p.mag2() ) - m;
    m_particle_gun->SetParticleEnergy( ke );
    m_particle_gun->SetParticlePosition( x );
    m_particle_gun->GeneratePrimaryVertex( anEvent );
    gAnaMan.SetPrimaryParticle( i, p, m, m_inc->pid[i] );
    gAnaMan.SetPrimaryVertex( i, x );
  }
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHdibaryonNonReso( G4Event* anEvent )
{
  G4double  momk[4]={0.}, mom[3]={0.};
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/GeV*m_KaonMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/GeV,
			(m_Proton->GetPDGMass()/GeV)*2,
			m_Hdibaryon->GetPDGMass()/GeV,
			m_KaonPlus->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];


  if((acos(momk[2]/sqrt(pow(momk[0],2)+pow(momk[1],2)+pow(momk[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hkinema.GetTheta(3);
  // G4double Phi_h = Hkinema.GetPhi(3);

  ////check kinematics
  momk[0]=mom_kp_x;
  momk[1]=mom_kp_y;
  momk[2]=mom_kp_z;
  momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz = G4RandFlat::shoot( m_target_pos.z() - m_target_size.z()/2,
				    m_target_pos.z() + m_target_size.z()/2 )*mm;

  //Kaon +
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //Hdibaryon
  m_particle_gun->SetParticleDefinition(m_Hdibaryon);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_Hdibaryon->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_Hdibaryon->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHybrid( G4Event* anEvent )
{
  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;
  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  // G4double rmpro = 0.93827203;
  // G4double rmpi = 0.13957018;

  //  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  G4double pbeam=1.;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  // G4double Ebeam = sqrt(pow(pbeam,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  mass_hybrid = 1.22;    // Mass for H
  width_hybrid = 0.000;  // Width for H

  KinemaHybrid Hybrid(m_PionMinus->GetPDGMass()/GeV,
		      m_Proton->GetPDGMass()/GeV,
		      m_Proton->GetPDGMass()/GeV,
		      m_PionMinus->GetPDGMass()/GeV,//pi1
		      m_PionPlus->GetPDGMass()/GeV, //pi2
		      mass_hybrid, width_hybrid, pg_z, 0.0);

  //  Hkinema.Dump();
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
  /*
  momcmprop=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  cmk->Fill(acos(momcmk[2]/momcmprop)/3.141592654*180);

  coscmk->Fill((momcmk[2]/momcmprop));

  labk->Fill(acos(momk[2]/momprop)*180/3.141592654);
  coslabk->Fill(momk[2]/momprop);
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
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  m_particle_gun->SetParticleEnergy((Energy_pro - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi1 PG
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  m_particle_gun->SetParticleEnergy((Energy_pi1 - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi2 PG
  m_particle_gun->SetParticleDefinition(m_PionPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  m_particle_gun->SetParticleEnergy((Energy_pi2 - m_PionPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z,m_Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z,m_PionMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z,m_PionPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHybrid3body( G4Event* anEvent )
{
  //  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;
  G4double mom[3];
  G4double Ebeam,pg_x,pg_y,pg_z;
  G4double rmpro = 0.93827203;
  // G4double rmpi = 0.13957018;

  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  G4double pbeam=0.7;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  Ebeam = sqrt(pow(pbeam,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  //  mass_hybrid = 1.22;    // Mass for H
  //  width_hybrid = 0.000;  // Width for H

  //  double temp= kin1.GetMomentumLab(3);
  //    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
  //  Kinema3Body kin2(kin3.M_res, 0.0, m3, m4, m5, temp, 0.0);
 G4double W=sqrt(pow(Ebeam+rmpro,2)-pow(pbeam,2));
 Kinema3Body Hybrid(W,0,
		    m_Proton->GetPDGMass()/GeV,
		    m_PionMinus->GetPDGMass()/GeV,
		    m_PionPlus->GetPDGMass()/GeV,
		    pbeam, 0.0);
//  Hkinema.Dump();
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
  /*
  momcmprop=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  cmk->Fill(acos(momcmk[2]/momcmprop)/3.141592654*180);

  coscmk->Fill((momcmk[2]/momcmprop));

  labk->Fill(acos(momk[2]/momprop)*180/3.141592654);
  coslabk->Fill(momk[2]/momprop);
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

  ///////Proton PG
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  m_particle_gun->SetParticleEnergy((Energy_pro - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi1 PG
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  m_particle_gun->SetParticleEnergy((Energy_pi1 - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi2 PG
  m_particle_gun->SetParticleDefinition(m_PionPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  m_particle_gun->SetParticleEnergy((Energy_pi2 - m_PionPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z,m_Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z,m_PionMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z,m_PionPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHybrid3bodyMode1( G4Event* anEvent )
{
  //  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;

  G4double mom[3];
  G4double Ebeam,pg_x,pg_y,pg_z;
  G4double rmpro=0.93827203;
  // G4double rmpi=0.13957018;

  //  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*0.01294);
  //  G4double pbeam=0.7;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  Ebeam = sqrt(pow(pbeam,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

 G4double W=sqrt(pow(Ebeam+rmpro,2)-pow(pbeam,2));
 Kinema3Body Hybrid(W,0,
		    m_Neutron->GetPDGMass()/GeV,
		    m_PionMinus->GetPDGMass()/GeV,
		    m_PionPlus->GetPDGMass()/GeV,
		    pbeam, 0.0);
//  Hkinema.Dump();
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
  vtz=rn_vtz+m_target_pos.z();
  //m_target_pos.z()

  ///////Neutron PG
  m_particle_gun->SetParticleDefinition(m_Neutron);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  m_particle_gun->SetParticleEnergy((Energy_pro - m_Neutron->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi1 PG
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  m_particle_gun->SetParticleEnergy((Energy_pi1 - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi2 PG
  m_particle_gun->SetParticleDefinition(m_PionPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  m_particle_gun->SetParticleEnergy((Energy_pi2 - m_PionPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z, m_Neutron->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z, m_PionMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z, m_PionPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHybrid3bodyMode2( G4Event* anEvent )
{
  //  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;

  G4double mom[3];
  G4double Ebeam,pg_x,pg_y,pg_z;
  G4double rmpro=0.93827203;
  // G4double rmpi=0.13957018;

  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*0.01294);
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  G4double pbeam=0.7;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  Ebeam = sqrt(pow(pbeam,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  //  mass_hybrid = 1.22;    // Mass for H
  //  width_hybrid = 0.000;  // Width for H

  //  double temp= kin1.GetMomentumLab(3);
  //    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
  //  Kinema3Body kin2(kin3.M_res, 0.0, m3, m4, m5, temp, 0.0);
 G4double W=sqrt(pow(Ebeam+rmpro,2)-pow(pbeam,2));
 Kinema3Body Hybrid(W,0,
		    m_Proton->GetPDGMass()/GeV,
		    m_PionMinus->GetPDGMass()/GeV,
		    m_PionPlus->GetPDGMass()/GeV,
		    pbeam, 0.0);
//  Hkinema.Dump();
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
  /*
  momcmprop=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  cmk->Fill(acos(momcmk[2]/momcmprop)/3.141592654*180);

  coscmk->Fill((momcmk[2]/momcmprop));

  labk->Fill(acos(momk[2]/momprop)*180/3.141592654);
  coslabk->Fill(momk[2]/momprop);
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
  vtz=rn_vtz+m_target_pos.z();
  //  double vtz = -157.5+15.0*((double) G4RandFlat::shoot()); //--> 15 mm
  //  double vtz = -160.0+20.0*((double) G4RandFlat::shoot()); //--> 20 mm
  //  double vtz = -162.5+25.0*((double) G4RandFlat::shoot()); //--> 25 mm
  //    double vtz = -165.0+30.0*((double) G4RandFlat::shoot()); //--> 30 mm
  //  double vtz = -200.0+100.0*((double) G4RandFlat::shoot()); //--> test//

  ///////proton PG
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  m_particle_gun->SetParticleEnergy((Energy_pro - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi1 PG
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  m_particle_gun->SetParticleEnergy((Energy_pi1 - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi2 PG
  m_particle_gun->SetParticleDefinition(m_PionPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  m_particle_gun->SetParticleEnergy((Energy_pi2 - m_PionPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z,m_Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z,m_PionMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z,m_PionPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHybrid3bodyMode3( G4Event* anEvent )
{
  //  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;

  G4double mom[3];
  G4double Ebeam,pg_x,pg_y,pg_z;
  G4double rmpro=0.93827203;
  // G4double rmpi=0.13957018;

  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*0.01294);
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  G4double pbeam=0.7;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  Ebeam = sqrt(pow(pbeam,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  //  mass_hybrid = 1.22;    // Mass for H
  //  width_hybrid = 0.000;  // Width for H

  //  double temp= kin1.GetMomentumLab(3);
  //    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
  //  Kinema3Body kin2(kin3.M_res, 0.0, m3, m4, m5, temp, 0.0);
 G4double W=sqrt(pow(Ebeam+rmpro,2)-pow(pbeam,2));
 Kinema3Body Hybrid(W,0,
		    m_Neutron->GetPDGMass()/GeV,
		    m_PionMinus->GetPDGMass()/GeV,
		    m_PionPlus->GetPDGMass()/GeV,
		    pbeam, 0.0);
//  Hkinema.Dump();
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
  /*
  momcmprop=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  cmk->Fill(acos(momcmk[2]/momcmprop)/3.141592654*180);

  coscmk->Fill((momcmk[2]/momcmprop));

  labk->Fill(acos(momk[2]/momprop)*180/3.141592654);
  coslabk->Fill(momk[2]/momprop);
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
  vtz=rn_vtz+m_target_pos.z();
  //  double vtz = -157.5+15.0*((double) G4RandFlat::shoot()); //--> 15 mm
  //  double vtz = -160.0+20.0*((double) G4RandFlat::shoot()); //--> 20 mm
  //  double vtz = -162.5+25.0*((double) G4RandFlat::shoot()); //--> 25 mm
  //    double vtz = -165.0+30.0*((double) G4RandFlat::shoot()); //--> 30 mm
  //  double vtz = -200.0+100.0*((double) G4RandFlat::shoot()); //--> test//

  ///////Neutron PG
  m_particle_gun->SetParticleDefinition(m_Neutron);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  m_particle_gun->SetParticleEnergy((Energy_pro - m_Neutron->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi1 PG
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  m_particle_gun->SetParticleEnergy((Energy_pi1 - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi2 PG
  m_particle_gun->SetParticleDefinition(m_PionPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  m_particle_gun->SetParticleEnergy((Energy_pi2 - m_PionPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z,m_Neutron->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z,m_PionMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z,m_PionPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHybrid3bodyMode4( G4Event* anEvent )
{
  //  G4double mass_hybrid, width_hybrid;

  G4double Energy_pro, mom_pro_x, mom_pro_y, mom_pro_z;
  G4double Energy_pi1, mom_pi1_x, mom_pi1_y, mom_pi1_z;
  G4double Energy_pi2, mom_pi2_x, mom_pi2_y, mom_pi2_z;

  G4double mom[3];
  G4double Ebeam,pg_x,pg_y,pg_z;
  G4double rmpro = 0.93827203;
  // G4double rmpi = 0.13957018;

  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*0.01294);
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  G4double pbeam=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  G4double pbeam=0.7;
  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = pbeam;

  Ebeam = sqrt(pow(pbeam,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  //  mass_hybrid = 1.22;    // Mass for H
  //  width_hybrid = 0.000;  // Width for H

  //  double temp= kin1.GetMomentumLab(3);
  //    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
  //  Kinema3Body kin2(kin3.M_res, 0.0, m3, m4, m5, temp, 0.0);
 G4double W=sqrt(pow(Ebeam+rmpro,2)-pow(pbeam,2));
 Kinema3Body Hybrid(W,0,
		    m_Proton->GetPDGMass()/GeV,
		    m_PionMinus->GetPDGMass()/GeV,
		    m_PionPlus->GetPDGMass()/GeV,
		    pbeam, 0.0);
//  Hkinema.Dump();
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
  /*
  momcmprop=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  cmk->Fill(acos(momcmk[2]/momcmprop)/3.141592654*180);

  coscmk->Fill((momcmk[2]/momcmprop));

  labk->Fill(acos(momk[2]/momprop)*180/3.141592654);
  coslabk->Fill(momk[2]/momprop);
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
  vtz=rn_vtz+m_target_pos.z();
  //  double vtz = -157.5+15.0*((double) G4RandFlat::shoot()); //--> 15 mm
  //  double vtz = -160.0+20.0*((double) G4RandFlat::shoot()); //--> 20 mm
  //  double vtz = -162.5+25.0*((double) G4RandFlat::shoot()); //--> 25 mm
  //    double vtz = -165.0+30.0*((double) G4RandFlat::shoot()); //--> 30 mm
  //  double vtz = -200.0+100.0*((double) G4RandFlat::shoot()); //--> test//

  ///////Proton PG
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pro_x,mom_pro_y,mom_pro_z));
  m_particle_gun->SetParticleEnergy((Energy_pro - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi1 PG
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi1_x,mom_pi1_y,mom_pi1_z));
  m_particle_gun->SetParticleEnergy((Energy_pi1 - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////pi2 PG
  m_particle_gun->SetParticleDefinition(m_PionPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pi2_x,mom_pi2_y,mom_pi2_z));
  m_particle_gun->SetParticleEnergy((Energy_pi2 - m_PionPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_pro_x,mom_pro_y,mom_pro_z,m_Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_pi1_x,mom_pi1_y,mom_pi1_z,m_PionMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_pi2_x,mom_pi2_y,mom_pi2_z,m_PionPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}



//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHdibaryon1( G4Event* anEvent )
{
  G4double mass_hdibaryon, width_hdibaryon;

  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  // G4double rmk=0.493677;

  // G4double Ebeam = sqrt(1.8*1.8+m_KaonMinus->GetPDGMass()/GeV*m_KaonMinus->GetPDGMass()/GeV);       // Incident gamma energy (GeV)

  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = 1.8;
  // G4double pbeam=1.8;

  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  mass_hdibaryon = gConf.Get<G4double>( "HdibaryonMass" );
  width_hdibaryon = gConf.Get<G4double>( "HdibaryonWidth" );

  Kinema3Resonance Hkinema(m_KaonMinus->GetPDGMass()/GeV,
			     ((m_Proton->GetPDGMass()/GeV)*2),
			     m_Lambda->GetPDGMass()/GeV,
			     m_Lambda->GetPDGMass()/GeV,
			     m_KaonPlus->GetPDGMass()/GeV,
			     mass_hdibaryon, width_hdibaryon, pg_z, 0.0);

//  Hkinema.Dump();

  Energy_kp = Hkinema.GetEnergy(5);
  // G4double   momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  /* L1 */
  Energy_L1 = Hkinema.GetEnergy(3);
  // G4double   momentum_L1 = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];

  // G4double ThetaL1 = Hkinema.GetTheta(3);
  // G4double PhiL1 = Hkinema.GetPhi(3);

  /* L2 */
  Energy_L2 = Hkinema.GetEnergy(4);
  // G4double   momentum_L2 = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double   ThetaL2 = Hkinema.GetTheta(4);
  // G4double   PhiL2 = Hkinema.GetPhi(4);

  /* Kp */
  Energy_kp = Hkinema.GetEnergy(5);
  // G4double   momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  // G4double Thetakp = Hkinema.GetThetaCM(1);
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  // G4double beta= sqrt(pbeam*pbeam)/(2.*0.938272013+Ebeam);
  // G4double momk[4]={};
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  // G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  // coscmk->Fill((momcmk[2]/momcmkpp));

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

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
  G4double vtz= G4RandFlat::shoot( m_target_pos.z() - m_target_size.z()/2,
					m_target_pos.z()+m_target_size.z()/2)*mm;

  ///////kp PG
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////L1 PG
  m_particle_gun->SetParticleDefinition( m_Lambda );
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  m_particle_gun->SetParticleEnergy((Energy_L1 - m_Lambda->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////L2 PG
  m_particle_gun->SetParticleDefinition( m_Lambda );
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  m_particle_gun->SetParticleEnergy((Energy_L2 - m_Lambda->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z, m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z, m_Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z, m_Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateHybridPHSG( G4Event* anEvent )
{
  //  G4double mom[3];
  G4double Ebeam;
  //  G4double rmk=0.493677;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  G4double pimom=1.;
  gAnaMan.SetPrimaryBeam(0,0,pimom);
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);
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
  // momhycm[3]=HybridBaryon->GetPDGMass()/GeV;
  // momhycm[3]=W2;
  G4double momhylab[4]={0};
  G4cout<<"momhylab x: "<<momhylab[0]<<G4endl;
  G4cout<<"momhylab y: "<<momhylab[1]<<G4endl;
  G4cout<<"momhylab z: "<<momhylab[2]<<G4endl;
  G4cout<<"momhylab energy"<<momhylab[3]<<G4endl;
  ///////gun
  G4cout<<"check energy"<<W<<G4endl;
  G4cout<<"w"<<W<<G4endl;
  G4cout<<"w2"<<W2<<G4endl;
  m_particle_gun->SetParticleDefinition( m_HybridBaryon );
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //  m_particle_gun->SetParticleEnergy((Ebeam+m_Proton->GetPDGMass()/GeV - HybridBaryon->GetPDGMass()/GeV)*GeV);
  // m_particle_gun->SetParticleEnergy((momhylab[3] - HybridBaryon->GetPDGMass()/GeV)*GeV);
  // m_particle_gun->SetParticleEnergy((W - HybridBaryon->GetPDGMass()/GeV)*GeV);
  // m_particle_gun->SetParticleEnergy((momhylab[3])*GeV);
  // m_particle_gun->SetParticleEnergy(0.5*GeV);
  //  m_particle_gun->SetParticleEnergy((W - HybridBaryon->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticleMomentum(1.*GeV); //KE ~ 0.35746632293688GeV

  //  G4cout<<"sh test energy:"<<(Ebeam+m_Proton->GetPDGMass()/GeV - HybridBaryon->GetPDGMass()/GeV)<<G4endl;
  //  G4cout<<"sh test energy:"<<momhylab[3]- HybridBaryon->GetPDGMass()/GeV<<G4endl;

  //  m_particle_gun->SetParticleEnergy((W - HybridBaryon->GetPDGMass()/GeV)*GeV);

  m_particle_gun->SetParticlePosition(G4ThreeVector(0.,0.,-150.*mm));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  //  G4cout<<HybridBaryon->GetPDGMass()/GeV  <<G4endl;

  /*
  gAnaMan.SetNumberOfPrimaryParticle(3);
  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,m_Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,m_Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
  */
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateTest( G4Event* anEvent )
{
  //  G4double mass_hdibaryon, width_hdibaryon;
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //  G4double mom[3];
  //  G4double  pbeam;

  /* proton */
  G4double mompx=0.1;
  G4double mompz=0.0;
  G4double rand= G4RandFlat::shoot(-1.,1.)*3.141592654;
  mom_p_x = cos(rand)*mompx-sin(rand)*mompz;
  mom_p_y = 0.;
  mom_p_z = sin(rand)*mompx+cos(rand)*mompz;
  G4double mom_p=sqrt(pow(mom_p_x,2)+pow(mom_p_y,2)+pow(mom_p_z,2));
  //  Energy_p=pow(mom_p,2)+m_Proton->GetPDGMass()/GeV;
  Energy_p=pow(mom_p,2)+m_PionMinus->GetPDGMass()/GeV;


  double vtz = -150+0.0*((double) G4RandFlat::shoot()); //--> 0 mm
  double vtx = 0.; //--> 0 mm
  double vty = -200.; //--> 0 mm

  ///////kp PG
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  //  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticleEnergy((Energy_p - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateDedxSingle( G4Event* anEvent )
{
  //  G4double mass_hdibaryon, width_hdibaryon;
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //  G4double mom[3];
  //  G4double Ebeam,pg_x,pg_y,pg_z, pbeam;

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
    vtz=rn_vtz+m_target_pos.z();
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-15.,15.)*mm;
    vty = G4RandFlat::shoot(-5.,5.)*mm;
    vtz = G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
  }

  ///////kp PG
  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2));
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  /*
  Energy_p=sqrt(pow(mom_p,2)+pow(piplus->GetPDGMass()/GeV,2));
  m_particle_gun->SetParticleDefinition(piplus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  m_particle_gun->SetParticleEnergy((Energy_p - piplus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  */
  gAnaMan.SetNumberOfPrimaryParticle(1);
  gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_Proton->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  //  gAnaMan.SetPrimaryParticle(1,mom_p_x,mom_p_y,mom_p_z,m_PionMinus->GetPDGMass()/GeV);
  //  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateAll( G4Event* anEvent )
{
  //  G4double mass_hdibaryon, width_hdibaryon;
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //  G4double mom[3];
  //  G4double Ebeam,pg_x,pg_y,pg_z, pbeam;

  //angle

  G4double phi=G4RandFlat::shoot(-1.,1.)*3.141592654;
  G4double mom_p=G4RandFlat::shoot(0.05,2.0);
  G4double theta=acos(G4RandFlat::shoot(0.,1.));
  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);
  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2))*GeV;

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
    vtz=rn_vtz+m_target_pos.z();
  }else if( gConf.Get<G4int>("Experiment") == 42. ){
    vtx = G4RandFlat::shoot(-15.,15.)*mm;
    vty = G4RandFlat::shoot(-5.,5.)*mm;
    vtz = G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
  }


  /*
  ///////kp PG
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  */

  G4double ratio=G4RandFlat::shoot(0.,1.);
  //  ratio = 0.1;

  if(ratio >= 0.0 && ratio < 0.2 ){
    // proton//
    Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/GeV,2));
    m_particle_gun->SetParticleDefinition(m_Proton);
    m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
    m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    m_particle_gun->GeneratePrimaryVertex( anEvent );
    gAnaMan.SetModeID(1);
    gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_Proton->GetPDGMass()/GeV);
    gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  } else  if(ratio >= 0.2 && ratio < 0.4 ){
    Energy_p=sqrt(pow(mom_p,2)+pow(m_KaonPlus->GetPDGMass()/GeV,2));
    m_particle_gun->SetParticleDefinition(m_KaonPlus);
    m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    m_particle_gun->SetParticleEnergy((Energy_p - m_KaonPlus->GetPDGMass()/GeV)*GeV);
    m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    m_particle_gun->GeneratePrimaryVertex( anEvent );
    gAnaMan.SetModeID(2);
    gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_KaonPlus->GetPDGMass()/GeV);
    gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  } else  if(ratio >= 0.4 && ratio< 0.6 ){
  // pi+ //
    Energy_p=sqrt(pow(mom_p,2)+ pow(m_PionPlus->GetPDGMass()/GeV,2));
    m_particle_gun->SetParticleDefinition(m_PionPlus);
    m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    m_particle_gun->SetParticleEnergy((Energy_p - m_PionPlus->GetPDGMass()/GeV)*GeV);
    m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    m_particle_gun->GeneratePrimaryVertex( anEvent );
    gAnaMan.SetModeID(3);
    gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_PionPlus->GetPDGMass()/GeV);
    gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  } else  if(ratio >= 0.6 && ratio< 0.8 ){
  // pi- //
    Energy_p=sqrt(pow(mom_p,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));
    m_particle_gun->SetParticleDefinition(m_PionMinus);
    m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    m_particle_gun->SetParticleEnergy((Energy_p - m_PionMinus->GetPDGMass()/GeV)*GeV);
    m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    m_particle_gun->GeneratePrimaryVertex( anEvent );
    gAnaMan.SetModeID(4);
    gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_PionMinus->GetPDGMass()/GeV);
    gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  } else  if(ratio >= 0.8 && ratio <= 1.0 ){
  // k- //
    Energy_p=sqrt(pow(mom_p,2)+pow(m_KaonMinus->GetPDGMass()/GeV,2));
    m_particle_gun->SetParticleDefinition(m_KaonMinus);
    m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
    m_particle_gun->SetParticleEnergy((Energy_p - m_KaonMinus->GetPDGMass()/GeV)*GeV);
    m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
    m_particle_gun->GeneratePrimaryVertex( anEvent );
    gAnaMan.SetModeID(5);
    gAnaMan.SetNumberOfPrimaryParticle(1);
    gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_KaonMinus->GetPDGMass()/GeV);
    gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  }

  ///////kp PG
}
void
TPCPrimaryGeneratorAction::GenerateKurama( G4Event* anEvent){
  static const G4double PionPlusMass = m_PionPlus->GetPDGMass()/CLHEP::GeV;
  static const G4double PionMinusMass = m_PionMinus->GetPDGMass()/CLHEP::GeV;
  static const G4double KaonPlusMass = m_KaonPlus->GetPDGMass()/CLHEP::GeV;
  static const G4double KaonMinusMass = m_KaonMinus->GetPDGMass()/CLHEP::GeV;
  static const G4double ProtonMass = m_Proton->GetPDGMass()/CLHEP::GeV;
  static const G4int PionPlusID = m_PionPlus->GetPDGEncoding();
  static const G4int PionMinusID = m_PionMinus->GetPDGEncoding();
  static const G4int KaonPlusID = m_KaonPlus->GetPDGEncoding();
  static const G4int KaonMinusID = m_KaonMinus->GetPDGEncoding();
  static const G4int ProtonID = m_Proton->GetPDGEncoding();

  G4ThreeVector gen_pos( m_beam->GetX( 0 ),
			 m_beam->GetY( 0 ),
			 m_beam->z );//Z = zTgt;
	G4int ntBeam = m_beam->ntBeam;
	if(ntBeam != 1 ) return;
	G4ThreeVector BeamMom = m_beam->p;
	G4double BeamMass2 = m_beam->m2;
	G4double BeamQ = m_beam->q;
	G4double BeamEnergy,BeamMass;
	if(isnan(BeamMass2) or isnan(gen_pos.x())) return;
//	G4cout<<Form("Ev %d, gen_pos(%f,%f)",m_beam->evnum,m_beam->GetX(0),m_beam->GetY(0))<<G4endl;
//	G4cout<<Form("BeamMom(%f,%f,%f) , Mass2 %f, Q %f",BeamMom.x(),BeamMom.y(),BeamMom.z(),BeamMass2,BeamQ)<<G4endl;
	G4int BeamID;
	if(BeamMass2 < 0.13 ){
		BeamEnergy = hypot(PionPlusMass,BeamMom.mag()) - PionPlusMass;
		BeamID = PionPlusID;
		BeamMass = PionPlusMass;
	}
	else if(BeamMass2 < 0.4 ){
		BeamEnergy = hypot(KaonPlusMass,BeamMom.mag()) - KaonPlusMass;
		BeamID = KaonPlusID;
		BeamMass = KaonPlusMass;
	}
	else if (BeamQ > 0){
		BeamEnergy = hypot(ProtonMass,BeamMom.mag()) - ProtonMass;
		BeamID = ProtonID;
		BeamMass = ProtonMass;
	}
	else{
		BeamEnergy = hypot(KaonMinusMass,BeamMom.mag()) - KaonMinusMass;
		BeamID = KaonMinusID;
		BeamMass = KaonMinusMass;
	}
	if(BeamQ < 0){
		BeamID = -BeamID;
	}
  gAnaMan.SetPrimaryVertex( 0, gen_pos );
	gAnaMan.SetPrimaryBeam( BeamMom );
  gAnaMan.SetNumberOfPrimaryParticle( 1 );
  gAnaMan.SetPrimaryParticle( 0, BeamMom, BeamMass, BeamID );
 	gAnaMan.SetRealBeamData(m_beam->runnum,m_beam->evnum,m_beam->trigpat);	
	if( BeamID == PionPlusID){	
			m_particle_gun->SetParticleDefinition( m_PionPlus );
	}else if( BeamID == PionMinusID){
			m_particle_gun->SetParticleDefinition( m_PionMinus );
	}else if( BeamID == KaonPlusID){
			m_particle_gun->SetParticleDefinition( m_KaonPlus );
	}else if( BeamID == KaonMinusID){
			m_particle_gun->SetParticleDefinition( m_KaonMinus );
	}else if( BeamID == ProtonID){
			m_particle_gun->SetParticleDefinition( m_Proton );
	}else {
			m_particle_gun->SetParticleDefinition( m_KaonMinus );
	}
	m_particle_gun->SetParticleEnergy(BeamEnergy * GeV);
	m_particle_gun->SetParticlePosition(gen_pos );
	m_particle_gun->SetParticleMomentumDirection(BeamMom );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
		
}

//Generator 135, PionPlusBeamthrough with RealData
	void
TPCPrimaryGeneratorAction::GeneratePionPlusBeamthrough( G4Event* anEvent )
{
  static const G4double PionPlusMass = m_PionPlus->GetPDGMass()/CLHEP::GeV;
  
  static const G4int PionPlusID = m_PionPlus->GetPDGEncoding();
  G4ThreeVector gen_pos( m_beam->GetX( 0 ),
			 m_beam->GetY( 0 ),
			 m_beam->z );//Z = V0
	G4ThreeVector BeamMom = m_beam->p;
	double BeamEnergy = hypot(PionPlusMass,BeamMom.mag()) - PionPlusMass;
	gAnaMan.SetPrimaryBeam( m_beam->p );
  gAnaMan.SetNumberOfPrimaryParticle( 1 );
  gAnaMan.SetPrimaryParticle( 0, BeamMom, PionPlusMass, PionPlusID );
 	gAnaMan.SetRealBeamData(m_beam->runnum,m_beam->evnum,m_beam->trigpat);	

	m_particle_gun->SetParticleDefinition( m_PionPlus );
	m_particle_gun->SetParticleEnergy(BeamEnergy * GeV);
	m_particle_gun->SetParticlePosition(gen_pos );
	m_particle_gun->SetParticleMomentumDirection(BeamMom );
  m_particle_gun->GeneratePrimaryVertex( anEvent );

}
	
//Generator -135, PionMinusBeamthrough with RealData
	void
TPCPrimaryGeneratorAction::GeneratePionMinusBeamthrough( G4Event* anEvent )
{
  static const G4double PionMinusMass = m_PionMinus->GetPDGMass()/CLHEP::GeV;
  
  static const G4int PionMinusID = m_PionMinus->GetPDGEncoding();

  G4ThreeVector gen_pos( m_beam->GetX( 0 ),
			 m_beam->GetY( 0 ),
			 m_beam->z );//Z = V0
	G4ThreeVector BeamMom = m_beam->p;
	double BeamEnergy = hypot(PionMinusMass,BeamMom.mag()) - PionMinusMass;
  gAnaMan.SetPrimaryBeam( m_beam->p );
  gAnaMan.SetNumberOfPrimaryParticle( 1 );
  gAnaMan.SetPrimaryParticle( 0, BeamMom, PionMinusMass, PionMinusID );
 	gAnaMan.SetRealBeamData(m_beam->runnum,m_beam->evnum,m_beam->trigpat);	
	m_particle_gun->SetParticleDefinition( m_PionMinus );
	m_particle_gun->SetParticleEnergy(BeamEnergy * GeV);
	m_particle_gun->SetParticlePosition(gen_pos );
	m_particle_gun->SetParticleMomentumDirection(BeamMom );
  m_particle_gun->GeneratePrimaryVertex( anEvent );

}

//Generator 493, KaonPlusBeamthrough with RealData
	void
TPCPrimaryGeneratorAction::GenerateKaonPlusBeamthrough( G4Event* anEvent )
{
  static const G4double KaonPlusMass = m_KaonPlus->GetPDGMass()/CLHEP::GeV;
  
  static const G4int KaonPlusID = m_KaonPlus->GetPDGEncoding();

  G4ThreeVector gen_pos( m_beam->GetX( 0 ),
			 m_beam->GetY( 0 ),
			 m_beam->z );//Z = V0
	G4ThreeVector BeamMom = m_beam->p;
	double BeamEnergy = hypot(KaonPlusMass,BeamMom.mag()) - KaonPlusMass;
  gAnaMan.SetPrimaryBeam( m_beam->p );
  gAnaMan.SetNumberOfPrimaryParticle( 1 );
  gAnaMan.SetPrimaryParticle( 0, BeamMom, KaonPlusMass, KaonPlusID );
 	gAnaMan.SetRealBeamData(m_beam->runnum,m_beam->evnum,m_beam->trigpat);	
  
	m_particle_gun->SetParticleDefinition( m_KaonPlus );
	m_particle_gun->SetParticleEnergy(BeamEnergy * GeV);
	m_particle_gun->SetParticlePosition(gen_pos );
	m_particle_gun->SetParticleMomentumDirection(BeamMom );
  m_particle_gun->GeneratePrimaryVertex( anEvent );

}
	
//Generator -493 KaonMinusBeamthrough with RealData
	void
TPCPrimaryGeneratorAction::GenerateKaonMinusBeamthrough( G4Event* anEvent )
{
  static const G4double KaonMinusMass = m_KaonMinus->GetPDGMass()/CLHEP::GeV;
  
  static const G4int KaonMinusID = m_KaonMinus->GetPDGEncoding();

  G4ThreeVector gen_pos( m_beam->GetX( 0 ),
			 m_beam->GetY( 0 ),
			 m_beam->z );//Z = V0
	G4ThreeVector BeamMom = m_beam->p;
	double BeamEnergy = hypot(KaonMinusMass,BeamMom.mag()) - KaonMinusMass;
  gAnaMan.SetPrimaryBeam( m_beam->p );
  gAnaMan.SetNumberOfPrimaryParticle( 1 );
  gAnaMan.SetPrimaryParticle( 0, BeamMom, KaonMinusMass, KaonMinusID );
 	gAnaMan.SetRealBeamData(m_beam->runnum,m_beam->evnum,m_beam->trigpat);	
  
	m_particle_gun->SetParticleDefinition( m_KaonMinus );
	m_particle_gun->SetParticleEnergy(BeamEnergy * GeV);
	m_particle_gun->SetParticlePosition(gen_pos );
	m_particle_gun->SetParticleMomentumDirection(BeamMom );
  m_particle_gun->GeneratePrimaryVertex( anEvent );

}

//Generator 938, ProtonBeamthrough with RealData
	void
TPCPrimaryGeneratorAction::GenerateProtonBeamthrough( G4Event* anEvent )
{
  static const G4double ProtonMass = m_Proton->GetPDGMass()/CLHEP::GeV;
  
  static const G4int ProtonID = m_Proton->GetPDGEncoding();

  G4ThreeVector gen_pos( m_beam->GetX( 0 ),
			 m_beam->GetY( 0 ),
			 m_beam->z );//Z = V0
	G4ThreeVector BeamMom = m_beam->p;
	double BeamEnergy = hypot(ProtonMass,BeamMom.mag()) - ProtonMass;
  gAnaMan.SetPrimaryBeam( m_beam->p );
  gAnaMan.SetNumberOfPrimaryParticle( 1 );
  gAnaMan.SetPrimaryParticle( 0, BeamMom, ProtonMass, ProtonID );
  
	m_particle_gun->SetParticleDefinition( m_Proton );
	m_particle_gun->SetParticleEnergy(BeamEnergy * GeV);
	m_particle_gun->SetParticlePosition(gen_pos );
	m_particle_gun->SetParticleMomentumDirection(BeamMom );
  m_particle_gun->GeneratePrimaryVertex( anEvent );

}

//_____________________________________________________________________________
// generator 60
void
TPCPrimaryGeneratorAction::GenerateLambda1405Rad1( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  G4double pbm[4];
  G4double Energy_L, mom_L_x, mom_L_y, mom_L_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  pbeam=1.65;
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;

  KinemaKstar Hkinema(m_PionMinus->GetPDGMass()/GeV,
			m_Proton->GetPDGMass()/GeV,
			m_Lambda1405->GetPDGMass()/GeV,
			m_KaonZeroS->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hkinema.GetEnergy(4);
  // G4double momentum_ksp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_L = Hkinema.GetEnergy(3);
  // G4double momentum_L = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_L_x = mom[0];
  mom_L_y = mom[1];
  mom_L_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_L = Hkinema.GetTheta(3);
  // G4double Phi_L = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+pow(m_KaonZeroS->GetPDGMass()/GeV,2));

  //  G4cout<<"momks:"<<momks[0]<<":"<<momks[1]<<":"<<momks[2]<<":"<<momks[3]<<G4endl;
  //  G4cout<<"Energy_ksp:  "<<Energy_ksp<<G4endl;
  //  G4cout<<"mom_ks:"<<mom_ksp_x<<":"<<mom_ksp_y<<":"<<mom_ksp_z<<G4endl;
  // G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;
  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_L<<G4endl;
  //coscmk->Fill(momcmk[2]/momcmkpp);

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
  vtz=rn_vtz+m_target_pos.z();


  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);


  //KaonStarZero
  m_particle_gun->SetParticleDefinition(m_KaonZeroS);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  m_particle_gun->SetParticleEnergy((Energy_ksp - m_KaonZeroS->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //L
  m_particle_gun->SetParticleDefinition(m_Lambda1405);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L_x,mom_L_y,mom_L_z));
  m_particle_gun->SetParticleEnergy((Energy_L - m_Lambda1405->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,m_KaonZeroS->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L_x,mom_L_y,mom_L_z,m_Lambda1405->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
// generator 61
void
TPCPrimaryGeneratorAction::GenerateLambda1405Rad2( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  G4double pbm[4];
  G4double Energy_L, mom_L_x, mom_L_y, mom_L_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  pbeam=1.65;
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;

  KinemaKstar Hkinema(m_PionMinus->GetPDGMass()/GeV,
			m_Proton->GetPDGMass()/GeV,
			m_Lambda1405R->GetPDGMass()/GeV,
			m_KaonZeroS->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hkinema.GetEnergy(4);
  // G4double momentum_ksp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_L = Hkinema.GetEnergy(3);
  // G4double momentum_L = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_L_x = mom[0];
  mom_L_y = mom[1];
  mom_L_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_L = Hkinema.GetTheta(3);
  // G4double Phi_L = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+pow(m_KaonZeroS->GetPDGMass()/GeV,2));

  //  G4cout<<"momks:"<<momks[0]<<":"<<momks[1]<<":"<<momks[2]<<":"<<momks[3]<<G4endl;
  //  G4cout<<"Energy_ksp:  "<<Energy_ksp<<G4endl;
  //  G4cout<<"mom_ks:"<<mom_ksp_x<<":"<<mom_ksp_y<<":"<<mom_ksp_z<<G4endl;
  // G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;
  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_L<<G4endl;
  //coscmk->Fill(momcmk[2]/momcmkpp);

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
  vtz=rn_vtz+m_target_pos.z();


  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);


  //KaonStarZero
  m_particle_gun->SetParticleDefinition(m_KaonZeroS);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  m_particle_gun->SetParticleEnergy((Energy_ksp - m_KaonZeroS->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //L
  m_particle_gun->SetParticleDefinition(m_Lambda1405R);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L_x,mom_L_y,mom_L_z));
  m_particle_gun->SetParticleEnergy((Energy_L - m_Lambda1405R->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,m_KaonZeroS->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L_x,mom_L_y,mom_L_z,m_Lambda1405R->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateLambda1405Reso( G4Event* anEvent )
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
    sigma = m_SigmaPlus;
    pion = m_PionMinus;
  }else if(check_decay_mode>1./3. && check_decay_mode<=2./3.){
    mode=2.;     //S0pi0
    sigma = m_SigmaZero;
    pion = m_PionZero;
  }else if(check_decay_mode>2./3. && check_decay_mode<=1.){
    mode=3.;     //S-pi+
    sigma = m_SigmaMinus;
    pion = m_PionPlus;
  }else{
    G4cout<<"Break Lambda decay mode"<<G4endl;
  }
  //  G4cout<<"check 2"<<G4endl;
  gAnaMan.SetModeID(mode);

  //  Ebeam = 1.9;       // Incident gamma energy (GeV)
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  // G4double Ebeam = sqrt(pow(pbeam,2)+pow(m_PionMinus->GetPDGMass()/GeV,2));
  gAnaMan.SetPrimaryBeam(0,0,pbeam);// how about add mass information?

  //  G4cout<< m_PionMinus->GetPDGMass()/GeV <<G4endl;
  //  G4cout<< m_Proton->GetPDGMass()/GeV <<G4endl;
  //  G4cout<< pion->GetPDGMass()/GeV <<G4endl;
  //  G4cout<< sigma->GetPDGMass()/GeV <<G4endl;
  //  G4cout<< m_KaonZeroS->GetPDGMass()/GeV <<G4endl;
  //  G4cout<< mass_hdibaryon <<G4endl;
  //  G4cout<< width_hdibaryon <<G4endl;
  //  G4cout<< pbeam <<G4endl;
  Kinema3Resonance  Hkinema(m_PionMinus->GetPDGMass()/GeV,
			     m_Proton->GetPDGMass()/GeV,
			     pion->GetPDGMass()/GeV,
			     sigma->GetPDGMass()/GeV,
			     m_KaonZeroS->GetPDGMass()/GeV,
			     mass_hdibaryon, width_hdibaryon, pbeam, 0.0);

  //  G4cout<<"check 3"<<G4endl;
  Energy_kp = Hkinema.GetEnergy(5);
  // G4double momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;

  /* L1 */
  Energy_L1 = Hkinema.GetEnergy(3);
  // G4double momentum_L1 = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];

  // G4double ThetaL1 = Hkinema.GetTheta(3);
  // G4double PhiL1 = Hkinema.GetPhi(3);

  /* L2 */
  Energy_L2 = Hkinema.GetEnergy(4);
  // G4double momentum_L2 = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double ThetaL2 = Hkinema.GetTheta(4);
  // G4double PhiL2 = Hkinema.GetPhi(4);

  /* Kp */
  Energy_kp = Hkinema.GetEnergy(5);
  // G4double momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  //  Thetakp = Hkinema.GetThetaCM(5);
  // G4double Thetakp = Hkinema.GetThetaCM(1);
  // G4double Phikp = Hkinema.GetPhi(5);
  //  G4cout<<Phikp<<G4endl;
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4cout<<"phi test: "<<Phikp<<":"<<shphi<<G4endl;
  //  G4cout<<"test: "<<Hkinema.GetPhiCM(1)<<":"<<Hkinema.GetPhi(5)<<G4endl;
  // Vertex (Reaction) point
  //  G4cout<<Energy_kp<<G4endl;
  //  //  G4double beta= sqrt(1.8*1.8-0.493677*0.493677)/(2*0.93827203+Energy_kp);
  //  G4cout<<(m_Proton->GetPDGMass()/GeV)<<G4endl;

  // G4double beta= pbeam/(m_Proton->GetPDGMass()/GeV+Ebeam);
  //  G4cout<<"Pri:"<<beta<<G4endl;
  //  G4double beta= sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2))/(2*0.93827203+1.8);

  // G4double momk[4]={};
  // G4double momcmk[4]={};
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+pow(m_PionMinus->GetPDGMass()/GeV,2));
  // G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
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
  G4double vtz= G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;

  ///////kp PG
  m_particle_gun->SetParticleDefinition(m_KaonZeroS);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonZeroS->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////L1 PG
  m_particle_gun->SetParticleDefinition(pion);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  m_particle_gun->SetParticleEnergy((Energy_L1 - pion->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////L2 PG
  m_particle_gun->SetParticleDefinition(sigma);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  m_particle_gun->SetParticleEnergy((Energy_L2 - sigma->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonZeroS->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,sigma->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,pion->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}

//_____________________________________________________________________________
// generator 63
void
TPCPrimaryGeneratorAction::GenerateSigma1385Rad( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  G4double pbm[4];
  G4double Energy_L, mom_L_x, mom_L_y, mom_L_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  pbeam=1.65;
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;

  // up:
  KinemaKstar Hkinema(m_PionMinus->GetPDGMass()/GeV,
			m_Proton->GetPDGMass()/GeV,
			m_Sigma1385R->GetPDGMass()/GeV,
			m_KaonZeroS->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hkinema.GetEnergy(4);
  // G4double momentum_ksp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_L = Hkinema.GetEnergy(3);
  // G4double momentum_L = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_L_x = mom[0];
  mom_L_y = mom[1];
  mom_L_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_L = Hkinema.GetTheta(3);
  // G4double Phi_L = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+pow(m_KaonZeroS->GetPDGMass()/GeV,2));

  //  G4cout<<"momks:"<<momks[0]<<":"<<momks[1]<<":"<<momks[2]<<":"<<momks[3]<<G4endl;
  //  G4cout<<"Energy_ksp:  "<<Energy_ksp<<G4endl;
  //  G4cout<<"mom_ks:"<<mom_ksp_x<<":"<<mom_ksp_y<<":"<<mom_ksp_z<<G4endl;
  // G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;
  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_L<<G4endl;
  //coscmk->Fill(momcmk[2]/momcmkpp);

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
  vtz=rn_vtz+m_target_pos.z();


  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);


  //KaonStarZero
  m_particle_gun->SetParticleDefinition(m_KaonZeroS);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  m_particle_gun->SetParticleEnergy((Energy_ksp - m_KaonZeroS->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //L
  m_particle_gun->SetParticleDefinition(m_Sigma1385R);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L_x,mom_L_y,mom_L_z));
  m_particle_gun->SetParticleEnergy((Energy_L - m_Sigma1385R->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,m_KaonZeroS->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L_x,mom_L_y,mom_L_z,m_Sigma1385R->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
// generator 62
void
TPCPrimaryGeneratorAction::GenerateSigma1385( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  G4double pbm[4];
  G4double Energy_L, mom_L_x, mom_L_y, mom_L_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  pbeam=1.65;
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"check11"<<G4endl;
  KinemaKstar Hkinema(m_PionMinus->GetPDGMass()/GeV,
			m_Proton->GetPDGMass()/GeV,
			m_Sigma1385Zero->GetPDGMass()/GeV,
			m_KaonZeroS->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hkinema.GetEnergy(4);
  // G4double momentum_ksp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_L = Hkinema.GetEnergy(3);
  // G4double momentum_L = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_L_x = mom[0];
  mom_L_y = mom[1];
  mom_L_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_L = Hkinema.GetTheta(3);
  // G4double Phi_L = Hkinema.GetPhi(3);

  ////check kinematics
  // G4doulbe momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+pow(m_KaonZeroS->GetPDGMass()/GeV,2));

  //  G4cout<<"momks:"<<momks[0]<<":"<<momks[1]<<":"<<momks[2]<<":"<<momks[3]<<G4endl;
  //  G4cout<<"Energy_ksp:  "<<Energy_ksp<<G4endl;
  //  G4cout<<"mom_ks:"<<mom_ksp_x<<":"<<mom_ksp_y<<":"<<mom_ksp_z<<G4endl;
  // G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;
  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_L<<G4endl;
  //coscmk->Fill(momcmk[2]/momcmkpp);

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
  vtz=rn_vtz+m_target_pos.z();


  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);
  //  G4cout<<"check"<<G4endl;

  //KaonStarZero
  m_particle_gun->SetParticleDefinition(m_KaonZeroS);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  m_particle_gun->SetParticleEnergy((Energy_ksp - m_KaonZeroS->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //L
  m_particle_gun->SetParticleDefinition(m_Sigma1385Zero);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L_x,mom_L_y,mom_L_z));
  m_particle_gun->SetParticleEnergy((Energy_L - m_Sigma1385Zero->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,m_KaonZeroS->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L_x,mom_L_y,mom_L_z,m_Sigma1385Zero->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
// generator 11
void
TPCPrimaryGeneratorAction::GeneratePionPlusKsL( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  //  G4double rmks=0.896;
  G4double pbm[4];
  G4double Energy_L, mom_L_x, mom_L_y, mom_L_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"m_PionMinus mass:"<<m_PionMinus->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"KstarMinus mass:"<<m_KaonStarZero->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"Lambda mass:"<<m_Lambda->GetPDGMass()/GeV<<G4endl;
  KinemaKstar Hkinema(m_PionMinus->GetPDGMass()/GeV,
			m_Proton->GetPDGMass()/GeV,
			m_Lambda->GetPDGMass()/GeV,
			m_KaonStarZero->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hkinema.GetEnergy(4);
  // G4double momentum_ksp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_L = Hkinema.GetEnergy(3);
  // G4double momentum_L = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_L_x = mom[0];
  mom_L_y = mom[1];
  mom_L_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_L = Hkinema.GetTheta(3);
  // G4double Phi_L = Hkinema.GetPhi(3);

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
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;
  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_L<<G4endl;
  //coscmk->Fill(momcmk[2]/momcmkpp);

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  G4double vtz= G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;

  //KaonStarZero
  m_particle_gun->SetParticleDefinition(m_KaonStarZero);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  m_particle_gun->SetParticleEnergy((Energy_ksp - m_KaonStarZero->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //L
  m_particle_gun->SetParticleDefinition(m_Lambda);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L_x,mom_L_y,mom_L_z));
  m_particle_gun->SetParticleEnergy((Energy_L - m_Lambda->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,m_KaonStarZero->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L_x,mom_L_y,mom_L_z,m_Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
// generator 12
void
TPCPrimaryGeneratorAction::GeneratePionPlusKsS( G4Event* anEvent )
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  //  G4double rmks=0.896;
  G4double pbm[4];
  G4double Energy_S, mom_S_x, mom_S_y, mom_S_z;
  G4double Energy_ksp, mom_ksp_x, mom_ksp_y, mom_ksp_z;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.9;
  gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);
  Ebeam = sqrt(pbeam*pbeam+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  //  G4cout<<"hdibaryon mass"<<hdibaryon->GetPDGMass()/GeV<<G4endl;
  KinemaKstar Hkinema(m_PionMinus->GetPDGMass()/GeV,
			m_Proton->GetPDGMass()/GeV,
			m_SigmaZero->GetPDGMass()/GeV,
			m_KaonStarZero->GetPDGMass()/GeV,
			pbm[2], 0.0);

  Energy_ksp=Hkinema.GetEnergy(4);
  // G4double momentum_ksp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_ksp_x = mom[0];
  mom_ksp_y = mom[1];
  mom_ksp_z = mom[2];

  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_S = Hkinema.GetEnergy(3);
  // G4double momentum_S = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_S_x = mom[0];
  mom_S_y = mom[1];
  mom_S_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_S = Hkinema.GetTheta(3);
  // G4double Phi_S = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momks[4];
  // momks[0]=mom_ksp_x;
  // momks[1]=mom_ksp_y;
  // momks[2]=mom_ksp_z;
  // momks[3]=sqrt(mom_ksp_x*mom_ksp_x+mom_ksp_y*mom_ksp_y+mom_ksp_z*mom_ksp_z+0.896*0.896);
  //  G4double momkpp=sqrt(pow(mom_ksp_x,2)+pow(mom_ksp_y,2)+pow(mom_ksp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);
  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momks[1],momks[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);

  G4double vtx = 0.*mm;
  G4double vty = 0.*mm;
  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);
  G4double vtz= G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;

  //KaonStarZero
  m_particle_gun->SetParticleDefinition(m_KaonStarZero);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_ksp_x,mom_ksp_y,mom_ksp_z));
  m_particle_gun->SetParticleEnergy((Energy_ksp - m_KaonStarZero->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  //S
  m_particle_gun->SetParticleDefinition(m_SigmaZero);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_S_x,mom_S_y,mom_S_z));
  m_particle_gun->SetParticleEnergy((Energy_S - m_SigmaZero->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );


  gAnaMan.SetNumberOfPrimaryParticle(2);
  gAnaMan.SetPrimaryParticle(0,mom_ksp_x,mom_ksp_y,mom_ksp_z,m_KaonStarZero->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_S_x,mom_S_y,mom_S_z,m_SigmaZero->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GeneratePionPlusKstarL( G4Event* anEvent )
{
  G4double mass_hdibaryon, width_hdibaryon;
  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  //  G4double rmk=0.493677;

  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = 1.8;
  // G4double pbeam = 1.8;
  // G4double Ebeam = sqrt(pbeam*pbeam+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);       // Incident gamma energy (GeV)

  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);
  mass_hdibaryon = gConf.Get<G4double>( "HdibaryonMass" ); // 0.89594;
  width_hdibaryon = gConf.Get<G4double>( "HdibaryonWidth" ); // 0.0487;
  Kinema3Resonance Hkinema(m_PionMinus->GetPDGMass()/GeV,
			     m_Proton->GetPDGMass()/GeV,
			     m_KaonPlus->GetPDGMass()/GeV,
			     m_PionMinus->GetPDGMass()/GeV,
			     m_Lambda->GetPDGMass()/GeV,
			     mass_hdibaryon, width_hdibaryon, pg_z, 0.0);
  //  Hkinema.Dump();
  Energy_kp = Hkinema.GetEnergy(5);
  // G4double momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);
  /* L1 */
  Energy_L1 = Hkinema.GetEnergy(3);
  // G4double momentum_L1 = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);
  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];
  // G4double ThetaL1 = Hkinema.GetTheta(3);
  // G4double PhiL1 = Hkinema.GetPhi(3);
  /* L2 */
  Energy_L2 = Hkinema.GetEnergy(4);
  // G4double momentum_L2 = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double ThetaL2 = Hkinema.GetTheta(4);
  // G4double PhiL2 = Hkinema.GetPhi(4);

  /* Kp */
  // G4double Energy_kp = Hkinema.GetEnergy(5);
  // G4double momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  // G4double Thetakp = Hkinema.GetThetaCM(1);
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  // G4double beta= sqrt(pbeam*pbeam)/(0.938272013+Ebeam);
  // G4double momcmk[4] = {};
  // G4double momk[4] = {};
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+m_Lambda->GetPDGMass()/GeV*m_Lambda->GetPDGMass()/GeV);
  // G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  // cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  // coscmk->Fill((momcmk[2]/momcmkpp));

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

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
  m_particle_gun->SetParticleDefinition(m_Lambda);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_Lambda->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////L1 PG
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  m_particle_gun->SetParticleEnergy((Energy_L1 - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////L2 PG
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  m_particle_gun->SetParticleEnergy((Energy_L2 - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);



  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,m_PionMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GeneratePionPlusKstarS( G4Event* anEvent )
{
  G4double mass_hdibaryon, width_hdibaryon;

  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  //  G4double rmk=0.493677;

  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = 1.8;
  // G4double pbeam=1.8;

  // G4double Ebeam = sqrt(pbeam*pbeam+m_PionMinus->GetPDGMass()/GeV*m_PionMinus->GetPDGMass()/GeV);       // Incident gamma energy (GeV)

  gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);
  mass_hdibaryon = gConf.Get<G4double>( "HdibaryonMass" ); // 0.89594;
  width_hdibaryon = gConf.Get<G4double>( "HdibaryonWidth" ); // 0.0487;

  Kinema3Resonance Hkinema(m_PionMinus->GetPDGMass()/GeV,
			     m_Proton->GetPDGMass()/GeV,
			     m_KaonPlus->GetPDGMass()/GeV,
			     m_PionMinus->GetPDGMass()/GeV,
			     m_Lambda->GetPDGMass()/GeV,
			     mass_hdibaryon, width_hdibaryon, pg_z, 0.0);

//  Hkinema.Dump();

  Energy_kp = Hkinema.GetEnergy(5);
  // G4double momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  /* L1 */
  Energy_L1 = Hkinema.GetEnergy(3);
  // G4double momentum_L1 = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];

  // G4double ThetaL1 = Hkinema.GetTheta(3);
  // G4double PhiL1 = Hkinema.GetPhi(3);

  /* L2 */
  Energy_L2 = Hkinema.GetEnergy(4);
  // G4double momentum_L2 = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double ThetaL2 = Hkinema.GetTheta(4);
  // G4double PhiL2 = Hkinema.GetPhi(4);

  /* Kp */
  Energy_kp = Hkinema.GetEnergy(5);
  // G4double momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  // G4double Thetakp = Hkinema.GetThetaCM(1);
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  // G4double beta= sqrt(pbeam*pbeam)/(0.938272013+Ebeam);
  // G4double momk[4] = {};
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+m_Lambda->GetPDGMass()/GeV*m_Lambda->GetPDGMass()/GeV);
  // G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  // cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  // coscmk->Fill((momcmk[2]/momcmkpp));

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

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
  G4double vtz= G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;

  ///////kp PG
  m_particle_gun->SetParticleDefinition(m_Lambda);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_Lambda->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////L1 PG
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  m_particle_gun->SetParticleEnergy((Energy_L1 - m_KaonPlus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );
  ///////L2 PG
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  m_particle_gun->SetParticleEnergy((Energy_L2 - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

  gAnaMan.SetNumberOfPrimaryParticle(3);

  gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_Lambda->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,m_KaonPlus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,m_PionMinus->GetPDGMass()/GeV);
  gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);

}

//_____________________________________________________________________________
void
TPCPrimaryGeneratorAction::GenerateTest2( G4Event* anEvent )
{
  //  G4double mass_hdibaryon, width_hdibaryon;
  G4double Energy_p,  mom_p_x, mom_p_y, mom_p_z;

  //  G4double mom[3];
  //  G4double Ebeam,pg_x,pg_y,pg_z, pbeam;

  /* proton */
  //  G4double mompx=0.10;
  G4double mompx=0.10;
  G4double mompz=0.0;

  G4double rand= G4RandFlat::shoot(-1.,1.)*3.141592654;
  mom_p_x = cos(rand)*mompx-sin(rand)*mompz;
  mom_p_y = 0.;
  mom_p_z = sin(rand)*mompx+cos(rand)*mompz;
  G4double mom_p=sqrt(pow(mom_p_x,2)+pow(mom_p_y,2)+pow(mom_p_z,2));
  //  Energy_p=pow(mom_p,2)+m_Proton->GetPDGMass()/GeV;
  Energy_p=pow(mom_p,2)+m_PionMinus->GetPDGMass()/GeV;

  G4double vtz= G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,m_target_pos.z()+m_target_size.z()/2)*mm;
  //  G4double vtz= G4RandFlat::shoot(-150.-m_target_size.z()/2,-150.+m_target_size.z()/2);

  //  double vtz = -150+0.0*((double) G4RandFlat::shoot()); //--> 0 mm
  G4double vtx = 0.; //--> 0 mm
  G4double vty = 0.; //--> 0 mm

  ///////kp PG
  //  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  //  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticleEnergy((Energy_p - m_PionMinus->GetPDGMass()/GeV)*GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex( anEvent );

}

void
TPCPrimaryGeneratorAction::GenerateLambda(G4Event* anEvent){
  static const G4double LambdaMass = m_Lambda->GetPDGMass()/CLHEP::GeV;
  m_particle_gun->SetParticleDefinition( m_Lambda );
  m_particle_gun->SetParticleMomentumDirection( G4ThreeVector(0,0,1) );
  m_particle_gun->SetParticleEnergy( (sqrt(2.)-1.)*LambdaMass*CLHEP::GeV );
  m_particle_gun->SetParticlePosition(G4ThreeVector(0,0,-143) );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
}
void
TPCPrimaryGeneratorAction::GenerateXiMinus(G4Event* anEvent){
  static const G4double XiMinusMass = m_XiMinus->GetPDGMass()/CLHEP::GeV;
  m_particle_gun->SetParticleDefinition( m_XiMinus );
  m_particle_gun->SetParticleMomentumDirection( G4ThreeVector(0,0,1) );
  m_particle_gun->SetParticleEnergy( (sqrt(2.)-1.)*XiMinusMass*CLHEP::GeV );
  m_particle_gun->SetParticlePosition(G4ThreeVector(0,0,-143) );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
}

//_____________________________________________________________________________
G4double
TPCPrimaryGeneratorAction::RandSin( void )
{
  G4int success = 0;
  G4double x,fx;
  while( success==0 ){
    x = 180.0 * G4RandFlat::shoot();
    //fx = sin(TMath::DegToRad()*x);
    fx = sin(x*TMath::Pi()/180.);
    if (fx >= G4RandFlat::shoot())
      success = 1;
  }
  return x;
}

void
TPCPrimaryGeneratorAction::GenerateKuramaPKmKpXi( G4Event* anEvent ){
	
	
  static const G4int KaonMinusID = m_KaonMinus->GetPDGEncoding();
  static const G4int KaonPlusID = m_KaonPlus->GetPDGEncoding();
  static const G4int XiMinusID = m_XiMinus->GetPDGEncoding();
  
	
	static const G4double KaonPlusMass = m_KaonPlus->GetPDGMass()/CLHEP::GeV;
  static const G4double XiMinusMass = m_XiMinus->GetPDGMass()/CLHEP::GeV;
	auto MomKm=	m_mm_vert->Moms[0];	
	auto MomKp=	m_mm_vert->Moms[1];	
	auto MomXi=	m_mm_vert->Moms[2];	
	double phi = G4RandFlat::shoot(-acos(-1),acos(-1));
	double KpEnergy = hypot(KaonPlusMass,MomKp.mag()) - KaonPlusMass;
	double XiEnergy = hypot(XiMinusMass,MomKp.mag()) - XiMinusMass;
	
	G4double dx = G4RandFlat::shoot(-15,15);
	G4double dy = G4RandFlat::shoot(-25,25);
	G4double dz = G4RandFlat::shoot(-153,-133);
	G4ThreeVector gen_pos(dx,dy,dz);
	MomKp.rotateZ(phi);
	MomKm.rotateZ(phi);
	MomXi.rotateZ(phi);

	auto SpinXi = MomXi.cross(MomKm);
	SpinXi = SpinXi*(1./SpinXi.mag());
	gTrackBuffer.SetPolarity(SpinXi,0);
	gTrackBuffer.SetMomentum(MomXi,0);
  gAnaMan.SetNumberOfPrimaryParticle( 3 );
	gAnaMan.SetPrimaryParticle(0,MomKm,KaonPlusMass,KaonMinusID);
	gAnaMan.SetPrimaryParticle(1,MomKp,KaonPlusMass,KaonPlusID);
	gAnaMan.SetPrimaryParticle(2,MomXi,XiMinusMass,XiMinusID);
	m_particle_gun->SetParticleDefinition( m_KaonPlus );
	m_particle_gun->SetParticleEnergy(KpEnergy * GeV);
	m_particle_gun->SetParticlePosition(gen_pos );
	m_particle_gun->SetParticleMomentumDirection(MomKp );
  m_particle_gun->GeneratePrimaryVertex( anEvent );
	
	m_particle_gun->SetParticleDefinition( m_XiMinus );
	m_particle_gun->SetParticleEnergy(XiEnergy * GeV);
	m_particle_gun->SetParticlePosition(gen_pos );
	m_particle_gun->SetParticleMomentumDirection(MomXi );
  m_particle_gun->GeneratePrimaryVertex( anEvent );




}

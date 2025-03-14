// -*- C++ -*-

#ifndef TPC_PRIMARY_GENERATOR_ACTION_HH
#define TPC_PRIMARY_GENERATOR_ACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ThreeVector.hh>
#include <G4Types.hh>

class G4ParticleGun;
class G4ParticleDefinition;
class TFile;
class TTree;

struct BeamInfo;
struct MMVertex;
struct KmKpL;
struct JamInfo;
struct IncInfo;

//_____________________________________________________________________________
class TPCPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  static G4String ClassName( void );
  TPCPrimaryGeneratorAction( void );
  ~TPCPrimaryGeneratorAction( void );

private:
  G4int                 m_generator;
  G4ParticleGun*        m_particle_gun;
  G4ThreeVector         m_target_pos;
  G4ThreeVector         m_target_size;
  G4ThreeVector         m_e45target_size;
  BeamInfo*             m_beam;
  MMVertex*             m_mm_vert;
  KmKpL*                m_kmkpl;
  KmKpL*                m_kmkpl2;
  G4double              m_beam_p0;
  JamInfo*              m_jam;
  IncInfo*              m_inc;
	std::vector<double>   rand_cont;
	G4ParticleDefinition* m_Neutron;
  G4ParticleDefinition* m_Proton;
  G4ParticleDefinition* m_AntiProton;
  G4ParticleDefinition* m_Lambda;
  G4ParticleDefinition* m_Lambda1405;
  G4ParticleDefinition* m_Lambda1405R;
  G4ParticleDefinition* m_SigmaPlus;
  G4ParticleDefinition* m_SigmaMinus;
  G4ParticleDefinition* m_SigmaZero;
  G4ParticleDefinition* m_Sigma1385Zero;
  G4ParticleDefinition* m_Sigma1385R;
  G4ParticleDefinition* m_XiMinus;
  G4ParticleDefinition* m_Xi1530Minus;
  G4ParticleDefinition* m_XiZero;
  G4ParticleDefinition* m_PionPlus;
  G4ParticleDefinition* m_PionMinus;
  G4ParticleDefinition* m_PionZero;
  G4ParticleDefinition* m_KaonPlus;
  G4ParticleDefinition* m_sKaonPlus;
  G4ParticleDefinition* m_KaonMinus;
  G4ParticleDefinition* m_sKaonMinus;
  G4ParticleDefinition* m_KaonZeroS;
  G4ParticleDefinition* m_KaonStarZero;
  G4ParticleDefinition* m_Phi;
  G4ParticleDefinition* m_Hdibaryon;
  G4ParticleDefinition* m_HdibaryonS;
  G4ParticleDefinition* m_HdibaryonLL;
  G4ParticleDefinition* m_LLphase;
  G4ParticleDefinition* m_HybridBaryon;
  G4double env_target_pos_z;

public:
  virtual void GeneratePrimaries( G4Event* anEvent );

public:
  void     GenerateKpXi2Body( G4Event* anEvent );
  void     GenerateKpXi2BodyUniform( G4Event* anEvent );
  void     GenerateKpXi1530Minus2BodyUniform( G4Event* anEvent );
	void     GeneratePPBar2Phi( G4Event* anEvent );
	void     GenerateHanul( G4Event* anEvent );
  void     GeneratePhaseSpace( G4Event* anEvent );
  void     GenerateHdibaryon1( G4Event* anEvent );
  void     GenerateHdibaryon2( G4Event* anEvent );
  void     GenerateHdibaryonPHSG( G4Event* anEvent );
  void     GenerateHdibaryonPHSGS( G4Event* anEvent );
  void     GenerateHdibaryonPHSGLL( G4Event* anEvent );
  void     GenerateKpKn( G4Event* anEvent );
  void     GenerateBeamVI( G4Event* anEvent );
  void     GenerateBeamVO( G4Event* anEvent );
  void     GenerateUniformProton( G4Event* anEvent );
  void     GenerateUniformPim( G4Event* anEvent );
  void     GenerateUniformKaonPlus( G4Event* anEvent );
  void     GenerateBeamProton( G4Event* anEvent );
  void     GenerateUniformProton_P( G4Event* anEvent );
  void     GenerateUniformProton_P_fixphi( G4Event* anEvent );
  void     GenerateUniformProton_P_Multi( G4Event* anEvent );
  void     GenerateJamInput( G4Event* anEvent );
  void     GenerateIncInput( G4Event* anEvent );
  void     GenerateIncInput_Vertex( G4Event* anEvent );
  void     GenerateJamInput_Randphi( G4Event* anEvent );
  void     GenerateLL_fromXiP( G4Event* anEvent );
  void     GenerateHdibaryonNonReso( G4Event* anEvent );
  void     GenerateTest( G4Event* anEvent );
  void     GenerateTest2( G4Event* anEvent );
  void     GenerateHybrid( G4Event* anEvent );
  void     GenerateHybridPHSG( G4Event* anEvent );
  void     GenerateHybrid3body( G4Event* anEvent );
  void     GenerateHybrid3bodyMode1( G4Event* anEvent );
  void     GenerateHybrid3bodyMode2( G4Event* anEvent );
  void     GenerateHybrid3bodyMode3( G4Event* anEvent );
  void     GenerateHybrid3bodyMode4( G4Event* anEvent );
  /// Study of Lambda 1405
  //  void GenerateLambda1405( G4Event* anEvent ); //#60
  void     GenerateLambda1405Rad1( G4Event* anEvent ); //#60
  void     GenerateLambda1405Rad2( G4Event* anEvent ); //#61  --> normal decay
  void     GenerateSigma1385( G4Event* anEvent ); //#62
  void     GenerateSigma1385Rad( G4Event* anEvent ); //#63
  void     GenerateLambda1405Reso( G4Event* anEvent ); //#64
  //Study of Kstar production
  void     GeneratePionPlusKsL( G4Event* anEvent ); //#11
  void     GeneratePionPlusKsS( G4Event* anEvent ); //#12
  void     GeneratePionPlusKstarL( G4Event* anEvent ); //#13
  void     GeneratePionPlusKstarS( G4Event* anEvent ); //#14
  void     GenerateAll( G4Event* anEvent ); //#98
  void     GenerateKurama( G4Event* anEvent ); //#100
  void     GeneratePionPlusBeamthrough( G4Event* anEvent ); //#135
  void     GeneratePionMinusBeamthrough( G4Event* anEvent ); //#-135
  void     GenerateKaonPlusBeamthrough( G4Event* anEvent ); //#493
  void     GenerateKaonMinusBeamthrough( G4Event* anEvent ); //#-493
  void     GenerateProtonBeamthrough( G4Event* anEvent ); //#938
  void     GenerateDedxSingle( G4Event* anEvent ); //#99
	void		 GenerateKuramaPKmKpXi(G4Event* anEvent);//181321;
	void		 GenerateKuramaPKmKpXi1530(G4Event* anEvent);//181530;
  void     GenerateTPCXiKmKp(G4Event* anEvent);//1001321

  // E07
  void     GenerateE07Study( G4Event* anEvent ); //#700
  void     GenerateE07StudyAll( G4Event* anEvent ); //#701
  void     GenerateE07StudyKnP( G4Event* anEvent ); //#702 knp --> knp, elastic
  void     GenerateE07StudyKnPBeam( G4Event* anEvent ); //#704 knp --> knp, elastic
  void     GenerateE07StudyKp( G4Event* anEvent ); //#703 knp --> kp, INC
  void     GenerateE07StudyKpBeam( G4Event* anEvent ); //#705 knp --> kp, INC, beam
  void     GenerateE07StudyKpXiBeam( G4Event* anEvent ); //#706 knp --> kpXi-, isotropic, beam
  void     GenerateE07StudyKpXiBeamOnlyKp( G4Event* anEvent ); //#707 knp -> kpXi-, isotropic, beam, only Kp
  void     GenerateE07StudyP08to20( G4Event* anEvent ); //#708 proton 0.4 GeV to 2 GeV
  void     GenerateE07StudyKp04to15( G4Event* anEvent ); //#709 k+ from 0.4 GeV to 1.5 GeV
  void     GenerateE07StudyKpxi1530( G4Event* anEvent ); //#710 knp --> kpxi1530-
  void     GenerateE07StudyTakahashi( G4Event* anEvent ); //#711 knp --> kpXi-, Takahashi-san's code
  // E27
  void     GenerateE27BeamThrough( G4Event* anEvent );//#2701
  void     GenerateE27Kptest( G4Event* anEvent );//#2702
  void     GenerateE27KppFLambdaP( G4Event* anEvent );//#2703
  void     GenerateE27KppFSigmaZP( G4Event* anEvent );//#2704
  void     GenerateE27KppFLambdaPizP( G4Event* anEvent );//#2705
  void     GenerateE27KppFSigmaZPizP( G4Event* anEvent );//#2706
  void     GenerateE27KppFSigmaPPimP( G4Event* anEvent );//#2707
  void     GenerateE27K11BLambda10Be( G4Event* anEvent );//#2708
  void     GenerateE27Kptest2( G4Event* anEvent );//#2709
  // KKpp
  void     GenerateKKppLL1( G4Event* anEvent );//#3001
  void     GenerateKKppLL2( G4Event* anEvent );//#3002
  void     GenerateKKppLSmPip( G4Event* anEvent );//#3003
  void     GenerateKKppLSpPim( G4Event* anEvent );//#3004
  void     GenerateJAMInput(G4Event* anEvent, TTree *t1);//#3101
  void     GenerateKKppBeamThrough1( G4Event* anEvent );//#3102
  void     GenerateJAMInputK0(G4Event* anEvent, TTree* t1);//#3103
  void     GenerateJAMInputK0bar(G4Event* anEvent, TTree* t1);//#3104
  // E45
  void     GenerateE45ElasticPionPlus( G4Event* anEvent ); //#4501
  void     GenerateE45ElasticPionMinus( G4Event* anEvent ); //#4502
  // Other methods
  void     GenerateXiMinus( G4Event* anEvent ); //#-1321
  void     GenerateLambda( G4Event* anEvent ); //#-1115
  void     GenerateKaonMinus( G4Event* anEvent ); //#-4930
  void     GenerateKmKpLL_BE( G4Event* anEvent ); //#1811161116
  void     GenerateTPCXi0nUniform(G4Event* anEvent);
  void     GenerateTPCLLUniform(G4Event* anEvent);

  void     GenerateAccidentals( int nAccidentals ,G4Event* anEvent); //Not intended for stand-alone use

	double   RandSin(void);
};

//_____________________________________________________________________________
inline G4String
TPCPrimaryGeneratorAction::ClassName( void )
{
  static G4String s_name("TPCPrimaryGeneratorAction");
  return s_name;
}

#endif

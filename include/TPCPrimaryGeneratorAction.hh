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
struct JamInfo;

//_____________________________________________________________________________
class TPCPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  static G4String ClassName( void );
  TPCPrimaryGeneratorAction( void );
  ~TPCPrimaryGeneratorAction( void );

private:
  G4int          m_generator;
  G4ParticleGun* m_particle_gun;
  G4ThreeVector  m_target_pos;
  G4ThreeVector  m_target_size;
  BeamInfo*      m_beam;
  G4double       m_beam_p0;
  JamInfo*       m_jam;
  G4ParticleDefinition* m_Neutron;
  G4ParticleDefinition* m_Proton;
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
  G4ParticleDefinition* m_PionPlus;
  G4ParticleDefinition* m_PionMinus;
  G4ParticleDefinition* m_PionZero;
  G4ParticleDefinition* m_KaonPlus;
  G4ParticleDefinition* m_KaonMinus;
  G4ParticleDefinition* m_KaonZeroS;
  G4ParticleDefinition* m_KaonStarZero;
  G4ParticleDefinition* m_Hdibaryon;
  G4ParticleDefinition* m_HdibaryonS;
  G4ParticleDefinition* m_HdibaryonLL;
  G4ParticleDefinition* m_LLphase;
  G4ParticleDefinition* m_HybridBaryon;

  G4double env_target_pos_z;
  G4double env_Beam_width;
  G4double env_Beam_x0;
  G4double env_Beam_y0;
  G4double env_Beam_dx;
  G4double env_Beam_dy;
  G4double env_Beam_u0;
  G4double env_Beam_v0;
  G4double env_Beam_du;
  G4double env_Beam_dv;

public:
  virtual void GeneratePrimaries( G4Event* anEvent );

public:
  G4double miss1(G4double pbeam, G4double m1,G4double *p1);
  G4double missks(G4double pbeam, G4double mbeam, G4double m1,G4double *p1);
  G4double miss1(G4double *pbeam, G4double m1,G4double *p1);
  G4double lorentz(G4double *v1,G4double betaz,G4double *v2);
  G4double lorentcmlab(G4double *v1,G4double betaz, G4double *v2);
  G4int HarmonicFermiMomentum(G4int Angular_mom, G4double *Kf);
  void Generate_hanul( G4Event* anEvent );
  void Generate_PhaseSpace( G4Event* anEvent );
  void Generate_hdibaryon1( G4Event* anEvent );
  void Generate_hdibaryon2( G4Event* anEvent );
  void Generate_hdibaryon_PHSG( G4Event* anEvent );
  void Generate_hdibaryon_PHSG_S( G4Event* anEvent );
  void Generate_hdibaryon_PHSG_LL( G4Event* anEvent );
  void GenerateKpKn( G4Event* anEvent );
  void GenerateBeamVI( G4Event* anEvent );
  void GenerateBeamVO( G4Event* anEvent );
  void GenerateJamInput( G4Event* anEvent );
  void Generate_hdibaryon_non_reso( G4Event* anEvent );
  void Generate_test( G4Event* anEvent );
  void Generate_test2( G4Event* anEvent );
  void Generate_hybrid( G4Event* anEvent );
  void Generate_hybridPHSG( G4Event* anEvent );
  void Generate_hybrid3body( G4Event* anEvent );
  void Generate_hybrid3body_mode1( G4Event* anEvent );
  void Generate_hybrid3body_mode2( G4Event* anEvent );
  void Generate_hybrid3body_mode3( G4Event* anEvent );
  void Generate_hybrid3body_mode4( G4Event* anEvent );
  /// Study of Lambda 1405
  //  void Generate_Lambda1405( G4Event* anEvent ); //#60
  void Generate_Lambda1405_rad1( G4Event* anEvent ); //#60
  void Generate_Lambda1405_rad2( G4Event* anEvent ); //#61  --> normal decay
  void Generate_Sigma1385( G4Event* anEvent ); //#62
  void Generate_Sigma1385_rad( G4Event* anEvent ); //#63
  void Generate_Lambda1405_reso( G4Event* anEvent ); //#64
  //Study of Kstar production
  void Generate_pip_KsL( G4Event* anEvent ); //#11
  void Generate_pip_KsS( G4Event* anEvent ); //#12
  void Generate_pip_KstarL( G4Event* anEvent ); //#13
  void Generate_pip_KstarS( G4Event* anEvent ); //#14
  //Study on Kstar production
  //  void Generate_dedx_all( G4Event* anEvent ); //#98
  void Generate_all( G4Event* anEvent ); //#98
  void Generate_dedx_single( G4Event* anEvent ); //#99
  // E07
  void GenerateE07Study( G4Event* anEvent ); //#700
  void GenerateE07StudyAll( G4Event* anEvent ); //#701
  void GenerateE07StudyKnP( G4Event* anEvent ); //#702 knp --> knp, elastic
  void GenerateE07StudyKnPBeam( G4Event* anEvent ); //#704 knp --> knp, elastic
  void GenerateE07StudyKp( G4Event* anEvent ); //#703 knp --> kp, INC
  void GenerateE07StudyKpBeam( G4Event* anEvent ); //#705 knp --> kp, INC, beam
  void GenerateE07StudyKpXiBeam( G4Event* anEvent ); //#706 knp --> kpXi-, isotropic, beam
  void GenerateE07StudyKpXiBeamOnlyKp( G4Event* anEvent ); //#707 knp -> kpXi-, isotropic, beam, only Kp
  void GenerateE07StudyP08to20( G4Event* anEvent ); //#708 proton 0.4 GeV to 2 GeV
  void GenerateE07StudyKp04to15( G4Event* anEvent ); //#709 k+ from 0.4 GeV to 1.5 GeV
  void GenerateE07StudyKpxi1530( G4Event* anEvent ); //#710 knp --> kpxi1530-
  void GenerateE07StudyTakahashi( G4Event* anEvent ); //#711 knp --> kpXi-, Takahashi-san's code
  // E27
  void GenerateE27BeamThrough( G4Event* anEvent );//#2701
  void GenerateE27Kptest( G4Event* anEvent );//#2702
  void GenerateE27KppFLambdaP( G4Event* anEvent );//#2703
  void GenerateE27KppFSigmaZP( G4Event* anEvent );//#2704
  void GenerateE27KppFLambdaPizP( G4Event* anEvent );//#2705
  void GenerateE27KppFSigmaZPizP( G4Event* anEvent );//#2706
  void GenerateE27KppFSigmaPPimP( G4Event* anEvent );//#2707
  void GenerateE27K11BLambda10Be( G4Event* anEvent );//#2708
  void GenerateE27Kptest2( G4Event* anEvent );//#2709
  // KKpp
  void GenerateKKppLL1( G4Event* anEvent );//#3001
  void GenerateKKppLL2( G4Event* anEvent );//#3002
  void GenerateKKppLSmPip( G4Event* anEvent );//#3003
  void GenerateKKppLSpPim( G4Event* anEvent );//#3004
  void GenerateJAMInput(G4Event* anEvent, TTree *t1);//#3101
  void GenerateKKppBeamThrough1( G4Event* anEvent );//#3102
  void GenerateJAMInputK0(G4Event* anEvent, TTree* t1);//#3103
  void GenerateJAMInputK0bar(G4Event* anEvent, TTree* t1);//#3104
  // E45
  void GenerateE45ElasticPionPlus( G4Event* anEvent ); //#4501
  void GenerateE45ElasticPionMinus( G4Event* anEvent ); //#4502

  double RandSin(void);
};

//_____________________________________________________________________________
inline G4String
TPCPrimaryGeneratorAction::ClassName( void )
{
  static G4String s_name("TPCPrimaryGeneratorAction");
  return s_name;
}

#endif

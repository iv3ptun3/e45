// -*- C++ -*-

#ifndef TPC_PRIMARY_GENERATOR_ACTION_HH
#define TPC_PRIMARY_GENERATOR_ACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ThreeVector.hh>
#include <G4Types.hh>

class G4ParticleGun;

struct BeamInfo;

//_____________________________________________________________________________
class TPCPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  static G4String ClassName( void );
  TPCPrimaryGeneratorAction( void );
  ~TPCPrimaryGeneratorAction( void );

private:
  G4ParticleGun* particleGun;   // particle gun provided by Geant4
  G4ThreeVector m_target_pos;
  G4ThreeVector m_target_size;
  BeamInfo*     m_beam;
  G4double      m_beam_p0;

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
  virtual void GeneratePrimaries(G4Event* anEvent);

public:
  G4double miss1(G4double pbeam, G4double m1,G4double *p1);
  G4double missks(G4double pbeam, G4double mbeam, G4double m1,G4double *p1);
  G4double miss1(G4double *pbeam, G4double m1,G4double *p1);
  G4double lorentz(G4double *v1,G4double betaz,G4double *v2);
  G4double lorentcmlab(G4double *v1,G4double betaz, G4double *v2);
  G4int HarmonicFermiMomentum(G4int Angular_mom, G4double *Kf);
  void Generate_hanul(G4Event* anEvent);
  void Generate_PhaseSpace(G4Event* anEvent); ///#30
  void Generate_hdibaryon1(G4Event* anEvent);
  void Generate_hdibaryon2(G4Event* anEvent);
  void Generate_hdibaryon_PHSG(G4Event* anEvent); //6
  void Generate_hdibaryon_PHSG_S(G4Event* anEvent); //7
  void Generate_hdibaryon_PHSG_LL(G4Event* anEvent); //9
  void Generate_Kp_Kn(G4Event* anEvent); //10
  void GenerateBeam( G4Event* anEvent ); //10
  void Generate_hdibaryon_non_reso(G4Event* anEvent);
  void Generate_test(G4Event* anEvent);
  void Generate_test2(G4Event* anEvent);
  void Generate_hybrid(G4Event* anEvent);
  void Generate_hybridPHSG(G4Event* anEvent);
  void Generate_hybrid3body(G4Event* anEvent);
  void Generate_E45_elastic_pip(G4Event* anEvent); //#25
  void Generate_E45_elastic_pin(G4Event* anEvent); //#26
  void Generate_hybrid3body_mode1(G4Event* anEvent);
  void Generate_hybrid3body_mode2(G4Event* anEvent);
  void Generate_hybrid3body_mode3(G4Event* anEvent);
  void Generate_hybrid3body_mode4(G4Event* anEvent);
  /// study of Lambda 1405
  //  void Generate_Lambda1405(G4Event* anEvent); //#60
  void Generate_Lambda1405_rad1(G4Event* anEvent); //#60
  void Generate_Lambda1405_rad2(G4Event* anEvent); //#61  --> normal decay
  void Generate_Sigma1385(G4Event* anEvent); //#62
  void Generate_Sigma1385_rad(G4Event* anEvent); //#63
  void Generate_Lambda1405_reso(G4Event* anEvent); //#64
  //study of Kstar production
  void Generate_pip_KsL(G4Event* anEvent); //#11
  void Generate_pip_KsS(G4Event* anEvent); //#12
  void Generate_pip_KstarL(G4Event* anEvent); //#13
  void Generate_pip_KstarS(G4Event* anEvent); //#14
  //study on Kstar production
  //  void Generate_dedx_all(G4Event* anEvent); //#98
  void Generate_all(G4Event* anEvent); //#98
  void Generate_dedx_single(G4Event* anEvent); //#99
  void Generate_E07_study(G4Event* anEvent); //#70
  void Generate_E07_study_all(G4Event* anEvent); //#71
  void Generate_E07_study_knp(G4Event* anEvent); //#72 knp --> knp, elastic
  void Generate_E07_study_knp_beam(G4Event* anEvent); //#74 knp --> knp, elastic
  void Generate_E07_study_kp(G4Event* anEvent); //#73 knp --> kp, INC
  void Generate_E07_study_kp_beam(G4Event* anEvent); //#75 knp --> kp, INC, beam
  void Generate_E07_study_kpxi_beam(G4Event* anEvent); //#76 knp --> kpXi-, isotropic, beam
  void Generate_E07_study_kpxi_beam_only_kp(G4Event* anEvent); //#77 knp --> kpXi-, isotropic, beam, only Kp
  void Generate_E07_study_pro_08_20(G4Event* anEvent); //#78 proton 0.4 GeV to 2 GeV
  void Generate_E07_study_kp_04_15(G4Event* anEvent); //#79 k+ from 0.4 GeV to 1.5 GeV
  void Generate_E07_study_kpxi1530(G4Event* anEvent); //#80 knp --> kpxi1530-
  void Generate_E07_study_Takahashi(G4Event* anEvent); //#81 knp --> kpXi-, Takahashi-san's code
  double RandSin(void);

  friend class E27Reaction;
  friend class KKppReaction;
};

//_____________________________________________________________________________
inline G4String
TPCPrimaryGeneratorAction::ClassName( void )
{
  static G4String s_name("TPCPrimaryGeneratorAction");
  return s_name;
}

#endif

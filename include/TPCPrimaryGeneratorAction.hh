// ====================================================================
//   TPCPrimaryGeneratorAction.hh
//
// ====================================================================
#ifndef TPC_PRIMARY_GENERATOR_ACTION_H
#define TPC_PRIMARY_GENERATOR_ACTION_H
 
#include "G4VUserPrimaryGeneratorAction.hh"
#include "TPCAnaManager.hh"
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>

class G4ParticleGun;
class TPCAnaManager;

class TPCPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:
  G4ParticleGun* particleGun;   // particle gun provided by Geant4
  TPCAnaManager* anaManager;

  G4double env_mass_hdibaryon;
  G4double env_width_hdibaryon;
  G4int env_Generator;
  //  G4double env_Target_z;
  G4double env_Kurama_gap;
  G4double env_target_width;
  G4double env_target_size_x;
  G4double env_target_size_y;

  G4double env_target_pos_z;
  G4double env_Beam_mom;
  G4double env_Beam_width;

  G4double env_Beam_x0;
  G4double env_Beam_y0;
  G4double env_Beam_dx;
  G4double env_Beam_dy;

  G4double env_Beam_u0;
  G4double env_Beam_v0;
  G4double env_Beam_du;
  G4double env_Beam_dv;




  G4int env_Experiment_num;

  TFile* file;
  TH1F* gen_im;
  TH1F* cmk;
  TH1F* phik;
  TH1F* phih;

  TH1F* phidiff;
  TH1F* thetadiff;

  TH1F* missk;
  TH1F* coscmk;
  TH1F* coscmh;

  TH1F* coslabk;
  TH1F* labk;

  TH1F* cmh;
public:
  TPCPrimaryGeneratorAction(TPCAnaManager* ana);
  ~TPCPrimaryGeneratorAction();


  G4double Get_env_Beam_mom() const {return env_Beam_mom;};
  G4double Get_env_Beam_width() const {return env_Beam_width;};
  G4double Get_env_target_width() const {return env_target_width;};
  G4double Get_env_target_size_x() const {return env_target_size_x;};
  G4double Get_env_target_size_y() const {return env_target_size_y;};

  G4double Get_env_target_pos_z() const {return env_target_pos_z;};

  G4double Get_env_Beam_x0() const {return env_Beam_x0;};
  G4double Get_env_Beam_y0() const {return env_Beam_y0;};
  G4double Get_env_Beam_dx() const {return env_Beam_dx;};
  G4double Get_env_Beam_dy() const {return env_Beam_dy;};
  G4double Get_env_Beam_u0() const {return env_Beam_u0;};
  G4double Get_env_Beam_v0() const {return env_Beam_v0;};
  G4double Get_env_Beam_du() const {return env_Beam_du;};
  G4double Get_env_Beam_dv() const {return env_Beam_dv;};



  G4double miss1(G4double pbeam, G4double m1,G4double *p1);
  G4double missks(G4double pbeam, G4double mbeam, G4double m1,G4double *p1);
  G4double miss1(G4double *pbeam, G4double m1,G4double *p1);
  G4double lorentz(G4double *v1,G4double betaz,G4double *v2);
  G4double lorentcmlab(G4double *v1,G4double betaz, G4double *v2);
  G4int HarmonicFermiMomentum(G4int Angular_mom, G4double *Kf);

  virtual void GeneratePrimaries(G4Event* anEvent);

  void Generate_hanul(G4Event* anEvent);

  void Generate_PhaseSpace(G4Event* anEvent); ///#30

  void Generate_hdibaryon1(G4Event* anEvent);
  void Generate_hdibaryon2(G4Event* anEvent); 
  void Generate_hdibaryon_PHSG(G4Event* anEvent); //6
  void Generate_hdibaryon_PHSG_S(G4Event* anEvent); //7
  void Generate_hdibaryon_PHSG_LL(G4Event* anEvent); //9
  void Generate_Kp_Kn(G4Event* anEvent); //10
  void Generate_beam(G4Event* anEvent); //10

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
  double deg2rad(double theta);
  double RandSin(void);


  friend class E27Reaction;
  friend class KKppReaction;


};


#endif

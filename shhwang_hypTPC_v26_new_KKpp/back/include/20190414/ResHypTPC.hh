// Author : Hiroyuki Sako
#ifndef RES_HYP_TPC_H
#define RES_HYP_TPC_H
#include <math.h>
#include "globals.hh"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "common.hh"

class TH1F;
class TF1;

class ResHypTPC {
private:
  static const bool debug = false;
  G4double pad_size;
  G4double pad_length;
  G4double threshold;

  //  static const G4double y_resolution = 0.5;
  //  static const G4double y_resolution;// = 0.5;
  //  static const G4double y_resolution;// = 0.5;
  G4double y_resolution;// = 0.5;

  //sigma/amp ratio of max amp dist (from Run 238)
  //static const G4double sigma_amp;// = 13.3/53.9;
  G4double sigma_amp;// = 13.3/53.9;

  //  static const G4double neff = 38;//proton T=400MeV (RCNP test)
  //  static const G4double neff;// = 26;//MIP
  G4double neff;// = 26;//MIP

  //  static const G4double neff_nmpv_correction;// = 0.7*1.2;
  G4double neff_nmpv_correction;// = 0.7*1.2;

  G4double diffusion_T;//mm/sqrt(cm)
  G4double const_smearing;
  //  static const G4double diffuse_GEM;//; = 0.1;
  G4double diffuse_GEM;//; = 0.1;

  G4double nmpv;
  G4double ncutoff;
  G4double nsigma;

  G4double n_sum;
  G4double x_mean;
  G4double x2_mean;

    TF1 *f_n_drift_electron;
  //    TF1 *f_diffusion;
    TF1 *f_GEM_avalanche;
  //    TF1 *f_const_smearing;

  void init() {
    //initialization  by sh
    // y_resolution = 0.5;
    //    sigma_amp = 13.3/53.9;
    //    neff = 26;//MIP
    //    neff_nmpv_correction;// = 0.7*1.2;
    //    diffuse_GEM = 0.1;
    /// init end
    nmpv = neff * neff_nmpv_correction;
    ncutoff = nmpv * 5;
    nsigma = nmpv * sigma_amp;

  };
  

public:
  //  ResHypTPC(G4double pad_size=2., G4double pad_length=10., G4double threshold=0.1, G4double diff_T=0.18, G4double smearing=0);
  ResHypTPC();
  ResHypTPC(G4double pad_size, G4double pad_length, G4double threshold, G4double diff_T, G4double smearing);
  ~ResHypTPC() {
    delete f_n_drift_electron;
    //    delete f_diffusion;
    delete f_GEM_avalanche;
    //    delete f_const_smearing;
};
G4double getXDeviation(G4int &n_electron, G4int &n_pad, G4double &x_rms, G4double x, G4double y, G4double dxdz, G4double dydz);
  G4double getYDeviation(G4double y);
  
  G4double getDiffusionX(G4double y) const {
    if (y>=0) {
      return diffusion_T*sqrt(y/10.);
    } else {
      return -1;
    }
  };



  G4double getNsum() {
    return n_sum;
  };

  G4double getNcutoff() {
    return ncutoff;
  }

};

G4double GetTransverseRes(G4double y);
//const G4double ResHypTPC::y_resolution = 0.5;
#endif


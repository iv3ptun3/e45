#ifndef __KINEMAFERMI_HH__
#define __KINEMAFERMI_HH__

#include "Kinema2Body.hh"
struct KINEMA_FERMI{
  double E_1_lab;
  double p_1_lab[4];
  double M_1;

  double E_2_lab;
  double p_2_lab[4];
  double M_2;

  double E_3_lab;
  double p_3_lab;
  double M_3;
  double P_3_lab[3];
  double theta3, phi3;

  double E_4_lab;
  double p_4_lab;
  double M_4;
  double P_4_lab[3];
  double theta4, phi4;

  double Theta1CM;
  double Phi1;
};

class KinemaFermi {
private:
  Kinema2Body kin1;
  KINEMA_FERMI kin3;

public:
  KinemaFermi(double m1, double m2, double m3, double m4, double *p1, double *p2, double cos_theta);
  double p2E(double p,double m);
  void CalcDistoribution(double unitx, double unity, double unitz, double *theta, double *phi);
  double deg2rad(double theta);
  double rag2deg(double rag);
  double RandSin(void);
  void Dump(void);
  double GetEnergy(int i);
  double GetMomentum(int i);
  void GetMomentum(int i, double *mom);
  double GetTheta(int i);
  double GetPhi(int i);
  double GetThetaCM(int i);
  double GetPhiCM(int i);
  void  RotateMom(int i, double deg, double *mom);
};

#endif

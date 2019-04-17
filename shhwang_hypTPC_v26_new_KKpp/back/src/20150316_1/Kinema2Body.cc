#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Kinema2Body.hh"
#define PI acos(-1.)
Kinema2Body::Kinema2Body(void)
{
}

Kinema2Body::Kinema2Body(double m1, double m2, double m3, double m4)
{
  kin.M_1 = m1;
  kin.M_2 = m2;
  kin.M_3 = m3;
  kin.M_4 = m4;
  //  G4cout<<"mass test 2body:"<<kin.M_1<<kin.M_2<<kin.M_3<<kin.M_4<<G4endl;
}

void Kinema2Body::SetMass(double *mass)
{
  kin.M_1 = mass[0];
  kin.M_2 = mass[1];
  kin.M_3 = mass[2];
  kin.M_4 = mass[3];
}
  
void Kinema2Body::SetMass(int i, double mass)
{
  switch (i) {
  case 1:
    kin.M_1 = mass;
    break;
  case 2:
    kin.M_2 = mass;
    break;
  case 3:
    kin.M_3 = mass;
    break;
  case 4:
    kin.M_4 = mass;
    break;
  default:
    fprintf(stderr,"Kinema2Body::SetMass No such particle:%d\n",i);
  }
  return;
}

void Kinema2Body::SetMomentum(int i,double mom)
{
  switch (i) {
  case 1:
    kin.p_1_lab = mom;
    break;
  case 2:
    kin.p_2_lab = mom;
    break;
  case 3:
    kin.p_3_lab = mom;
    break;
  case 4:
    kin.p_4_lab = mom;
    break;
  default:
    fprintf(stderr,"Kinema2Body::SetMomentum No such particle:%d\n",i);
  }
  return;
}


void Kinema2Body::SetTheta(int i,double theta)
{
  switch (i) {
  case 1:
    kin.p_1_theta = theta;
    break;
  case 2:
    kin.p_2_theta = theta;
    break;
  default:
    fprintf(stderr,"Kinema2Body::SetTheta No such particle:%d\n",i);
  }
  return;
}


void Kinema2Body::SetThetaCM(double theta_cm)
{
  kin.theta_cm = theta_cm;
  return;
}

int Kinema2Body::calc_kinema(void)
{
  
  kin.E_1_lab = p2E(kin.p_1_lab/cos(kin.p_1_theta),kin.M_1); //energy m1
  kin.E_2_lab = p2E(kin.p_2_lab/cos(kin.p_2_theta),kin.M_2); //energy m2
  /*
  G4cout<<"kin.p_1_theta:"<<kin.p_1_theta<<G4endl;
  G4cout<<"kin.p_2_theta:"<<kin.p_2_theta<<G4endl;
  G4cout<<"calc_kinema"<<G4endl;  
  G4cout<<"E1"<<G4endl;  
  G4cout<<kin.E_1_lab<<G4endl;
  G4cout<<kin.p_1_lab<<":"<<kin.M_1<<G4endl;
  G4cout<<"E2"<<G4endl;  
  G4cout<<kin.E_2_lab<<G4endl;
  G4cout<<kin.p_2_lab<<":"<<kin.M_2<<G4endl;
  */
  //  kin.beta_cm = pE2beta(kin.p_1_lab,kin.E_1_lab,kin.p_2_lab,kin.M_2);
  kin.beta_cm = pE2beta(kin.p_1_lab,kin.E_1_lab,kin.M_1, kin.p_2_lab,kin.E_2_lab,kin.M_2);
  //  G4cout<<"beta:"<<kin.beta_cm<<G4endl;
  //  G4cout<<"beta 2body:"<<kin.beta_cm<<G4endl;
  //double Kinema2Body::pE2beta(double p1,double E1,double p2, double E2)
  kin.gamma_cm = beta2gamma(kin.beta_cm);

  //calculate of Kinematics in CM system
  kin.E_1_cm = kin.gamma_cm * kin.E_1_lab - kin.beta_cm * kin.gamma_cm * kin.p_1_lab;
  kin.E_2_cm = kin.gamma_cm * kin.E_2_lab - kin.beta_cm * kin.gamma_cm * kin.p_2_lab;
  //  kin.E_2_cm = kin.gamma_cm * kin.M_2;
  //  kin.p_12_cm = E2p(kin.E_1_cm+kin.E_2_cm,kin.M_1+kin.M_2);//?? 12?? 1+2?
  kin.p_12_cm= -kin.gamma_cm*kin.beta_cm*kin.E_1_lab + kin.gamma_cm*kin.p_1_lab/cos(kin.p_1_theta);

  //  kin.E_12_cm = kin.gamma_cm * kin.E_2_lab - kin.beta_cm * kin.gamma_cm * kin.p_2_lab;
  /*
  G4cout<<"-----------------------------------------"<<G4endl;  
  G4cout<<"center of mass e1 e2"<<G4endl;  
  G4cout<<"p1:"<<sqrt(pow(kin.E_1_cm,2)-pow(kin.M_1,2))<<G4endl;
  G4cout<<"p2:"<<sqrt(pow(kin.E_2_cm,2)-pow(kin.M_2,2))<<G4endl;
  G4cout<<kin.E_1_cm<<G4endl;
  G4cout<<kin.E_2_cm<<G4endl;
  G4cout<<kin.p_12_cm<<G4endl;
  G4cout<<"E1+E2:"<<kin.E_1_cm+kin.E_2_cm<<G4endl;
  */

  //  G4cout<<kin.p_12_cm<<G4endl;
  //  kin.E_3_cm = (pow(kin.E_1_cm + kin.E_2_cm,2.0) + pow(kin.M_3,2.0)  - pow(kin.M_4,2.0))/(2.0*(kin.E_1_cm + kin.E_2_cm));
  //  kin.E_4_cm = (pow(kin.E_1_cm + kin.E_2_cm,2.0) + pow(kin.M_4,2.0)  - pow(kin.M_3,2.0))/(2.0*(kin.E_1_cm + kin.E_2_cm));

  //  double E12_cm=sqrt(pow(kin.M_1,2)+2.*kin.E_1_cm*kin.E_2_cm -2.*kin.p_1_lab*kin.p_2_lab + pow(kin.M_2,2));
  double E12_cm=sqrt( pow( sqrt(pow(kin.M_1,2)+pow(kin.p_1_lab/cos(kin.p_1_theta),2)) +sqrt(pow(kin.M_2,2)+pow(kin.p_2_lab/cos(kin.p_2_theta),2) ),2)-pow(kin.p_1_lab+kin.p_2_lab,2));
  kin.E_3_cm = (pow(E12_cm,2) +pow(kin.M_3,2)-pow(kin.M_4,2))/(2.*E12_cm);
  kin.E_4_cm = (pow(E12_cm,2) -pow(kin.M_3,2)+pow(kin.M_4,2))/(2.*E12_cm);

  //  kin.p_34_cm = sqrt((pow((kin.M_3 + kin.E_1_cm + kin.E_2_cm),2.0)-pow(kin.M_4,2.0))*(pow((kin.M_3 - (kin.E_1_cm + kin.E_2_cm)),2.0)-pow(kin.M_4,2.0))/(4.0*pow(kin.E_1_cm + kin.E_2_cm,2.0)));
  //  kin.p_34_cm = sqrt( pow(E12_cm,4)+pow( pow(kin.M_1,2)-pow(kin.M_2,2), 2) -2.*pow(E12_cm,2)*(pow(kin.M_1,2)+pow(kin.M_2,2)   )   )/(2.*E12_cm);
  //  kin.p_34_cm=sqrt(pow(E12_cm,2)-pow(kin.M_4+kin.M_3,2));


  kin.p_34_cm = sqrt( pow(E12_cm,4)+pow( pow(kin.M_3,2)-pow(kin.M_4,2), 2) -2.*pow(E12_cm,2)*(pow(kin.M_3,2)+pow(kin.M_4,2)   )   )/(2.*E12_cm);

  /*
  //  kin.p_34_cm=sqrt(pow(kin.E_3_cm,2)-pow(kin.M_3,2));
  G4cout<<"p3:"<<sqrt(pow(kin.E_3_cm,2)-pow(kin.M_3,2))<<G4endl;
  G4cout<<"p4:"<<sqrt(pow(kin.E_4_cm,2)-pow(kin.M_4,2))<<G4endl;
  G4cout<<"p_34:"<<kin.p_34_cm<<G4endl;
  G4cout<<"E3+E4:"<<kin.E_3_cm+kin.E_4_cm<<G4endl;
  */

  /*
    G4cout<<"------------------------------------------"<<G4endl;
    G4cout<<"2body inside, energy cm 3 4:"<<kin.E_3_cm<<":"<<kin.E_4_cm<<G4endl;
    G4cout<<"2body inside, mom 34:"<<kin.p_34_cm<<G4endl;
    G4cout<<"------------------------------------------"<<G4endl;
    G4cout<<"comp cm 12 and 34:"<<kin.E_1_cm+kin.E_2_cm<<", "<<kin.E_3_cm+kin.E_4_cm<<G4endl;
  //  G4cout<<kin.p_34_cm<<G4endl;
  */



  //calculate theta,phi of LAB system from CM system
  kin.theta_lab = theta_cm2theta_lab(kin.theta_cm,kin.p_34_cm,kin.E_3_cm,kin.gamma_cm,kin.beta_cm);
  kin.phi_lab   = theta_cm2phi_lab(kin.theta_cm,kin.p_34_cm,kin.E_4_cm,kin.gamma_cm,kin.beta_cm);

  //    G4cout<<"------------------------------------------"<<G4endl;
  //    G4cout<<"2body inside, theta_lab 1 2:"<<kin.theta_lab<<":"<<kin.phi_lab<<G4endl;
  //    G4cout<<"------------------------------------------"<<G4endl;

  //calculate pormentum ,energy of LAB system
  double p_lab_z,p_lab_xy;

  //  p_lab_z = kin.beta_cm * kin.gamma_cm * kin.E_3_cm + kin.gamma_cm * kin.p_34_cm * cos(deg2rad(kin.theta_cm));
  p_lab_z = kin.beta_cm * kin.gamma_cm * kin.E_3_cm + kin.gamma_cm * kin.p_34_cm * cos(deg2rad(kin.theta_cm));
  p_lab_xy = kin.p_34_cm * sin(deg2rad(kin.theta_cm));
  kin.p_3_lab = sqrt(pow(p_lab_z,2.0) + pow(p_lab_xy,2.0));

  
  p_lab_z = kin.beta_cm * kin.gamma_cm * kin.E_4_cm - kin.gamma_cm * kin.p_34_cm * cos(deg2rad(kin.theta_cm));
  kin.p_4_lab = sqrt(pow(p_lab_z,2.0) + pow(-p_lab_xy,2.0));

  kin.E_3_lab = p2E(kin.p_3_lab,kin.M_3);
  kin.E_4_lab = p2E(kin.p_4_lab,kin.M_4);

  /*
  G4cout<<"kin.p_3_lab:"<<kin.p_3_lab<<G4endl;
  G4cout<<"kin.p_4_lab:"<<kin.p_4_lab<<G4endl;
  G4cout<<"kin.E_3_lab:"<<kin.E_3_lab<<G4endl;
  G4cout<<"kin.E_4_lab:"<<kin.E_4_lab<<G4endl;
  */
  return 0;
}

double Kinema2Body::p2E(double p,double m)
{
  return sqrt(p*p + m*m);
}

double Kinema2Body::E2p(double E,double m)
{
  if (E*E - m*m<0.000000001)
    return 0.0;
  else
    return sqrt(E*E - m*m);
}

double Kinema2Body::pE2beta(double p,double E,double m_2)
{
  //target stops in rest
  return p/(E + m_2);
}

double Kinema2Body::pE2beta(double p1,double E1,double m1, double p2, double E2, double m2)
{
  //move --> it is correct
  //Ecm=sqrt((E1+E2)^2-(p1+p2)^2)
  //beta = p/E
  /*
  G4cout<<"-------------------------"<<G4endl;
  G4cout<<"in the pE2beta"<<G4endl;
  G4cout<<"p1:"<<p1<<G4endl;
  G4cout<<"E1:"<<E1<<G4endl;
  G4cout<<"m1:"<<m1<<G4endl;
  G4cout<<"p2:"<<p2<<G4endl;
  G4cout<<"E2:"<<E2<<G4endl;
  G4cout<<"m2:"<<m2<<G4endl;
  G4cout<<"-------------------------"<<G4endl;
  */
  //  return (p1+p2)/sqrt(pow(E1+E2,2)-pow(p1+p2,2));
  return (p1+p2)/(E1+E2);
}

double Kinema2Body::beta2gamma(double beta)
{
  return 1.0/sqrt(1.0 - pow(beta,2.0));
}

double Kinema2Body::theta2theta_cm(double theta,double p,double E,double kin_gamma,double beta)
{
  double value;

  value = atan(p*sin(deg2rad(theta))/(-kin_gamma*beta*E + kin_gamma*p*cos(deg2rad(theta))))*180.0/PI;

  if ( value < 0.0 )
    value = value + 180.0;

  return value;
}

double Kinema2Body::theta_cm2theta_lab(double theta,double p,double E,double gamma_cm,double beta_cm)
{
  double value;

  value = atan(p*sin(deg2rad(theta))/(beta_cm*gamma_cm*E + gamma_cm*p*cos(deg2rad(theta))))*180.0/PI;

  if (value < 0.0)
    value = value + 180.0;

  return value;
}

double Kinema2Body::theta_cm2phi_lab(double theta,double p,double E,double gamma_cm,double beta_cm)
{
  double value;

  value = atan(p*sin(deg2rad(theta))/(beta_cm*gamma_cm*E - gamma_cm*p*cos(deg2rad(theta))))*180.0/PI;

  if (value < 0.0)
    value = value + 180.0;

  return value;
}

double Kinema2Body::deg2rad(double theta) {
  return PI*theta/180.0;
}

double Kinema2Body::rag2deg(double rag)
{
  return 360.0 * rag/ (2.0 * PI);
}

void Kinema2Body::Dump(void)
{
  printf("===theta_cm:%f===\n",kin.theta_cm);
  printf("E_1_cm:%f p_1_cm:%f\n",kin.E_1_cm,kin.p_12_cm);
  printf("E_2_cm:%f p_2_cm:%f\n",kin.E_2_cm,kin.p_12_cm);
  printf("E_3_cm:%f p_3_cm:%f\n",kin.E_3_cm,kin.p_34_cm);
  printf("E_4_cm:%f p_4_cm:%f\n",kin.E_4_cm,kin.p_34_cm);
  printf("\n");
  printf("***LAB SYSTEM***\n");
  printf("E_1_lab:%f p_1_lab:%f\n",kin.E_1_lab,kin.p_1_lab);
  printf("E_2_lab:%f\n",kin.M_2);
  printf("E_3_lab:%f p_3_lab:%f\n",kin.E_3_lab,kin.p_3_lab);
  printf("E_4_lab:%f p_4_lab:%f\n",kin.E_4_lab,kin.p_4_lab);
  printf("theta_lab:%f phi_lab:%f",kin.theta_lab,kin.phi_lab);
  printf("\n\n");

}

double Kinema2Body::GetMomentumLab(int i)
{
  switch (i) {
  case 1:
    return kin.p_1_lab;
    break;
  case 2:
    return kin.p_2_lab;
    break;
  case 3:
    return kin.p_3_lab;
    break;
  case 4:
    return kin.p_4_lab;
    break;
  default:
    fprintf(stderr,"Kinema2Body::GetMomLab No such particle:%d\n",i);
    exit(1);
  }

}

double Kinema2Body::GetTheta(int i)
{
  switch (i) {
  case 1:
    return kin.p_1_theta;
    break;
  case 2:
    return kin.p_2_theta;
    break;
  default:
    fprintf(stderr,"Kinema2Body::GetTheta No such particle:%d\n",i);
    exit(1);
  }

}


double Kinema2Body::GetEnergyLab(int i)
{
  switch (i) {
  case 1:
    return kin.E_1_lab;
    break;
  case 2:
    return kin.E_2_lab;
    break;
  case 3:
    return kin.E_3_lab;
    break;
  case 4:
    return kin.E_4_lab;
    break;
  default:
    fprintf(stderr,"Kinema2Body::GetEnergyLab No such particle:%d\n",i);
    exit(1);
  }

}



double Kinema2Body::GetMomentumCM(int i)
{
  switch (i) {
  case 1:
    return kin.p_12_cm;
    break;
  case 2:
    return -kin.p_12_cm;
    break;
  case 3:
    return kin.p_34_cm;
    break;
  case 4:
    return -kin.p_34_cm;
    break;
  default:
    fprintf(stderr,"Kinema2Body::GetMomCm No such particle:%d\n",i);
    exit(1);
  }

}

double Kinema2Body::GetEnergyCM(int i)
{
  switch (i) {
  case 1:
    return kin.E_1_cm;
    break;
  case 2:
    return kin.E_2_cm;
    break;
  case 3:
    return kin.E_3_cm;
    break;
  case 4:
    return kin.E_4_cm;
    break;
  default:
    fprintf(stderr,"Kinema2Body::GetEnergyCM No such particle:%d\n",i);
    exit(1);
  }

}

double Kinema2Body::GetThetaLab(void)
{
  return kin.theta_lab;
}

double Kinema2Body::GetPhiLab(void)
{
  return kin.phi_lab;
}

double Kinema2Body::GetThetaCM(void)
{
  return kin.theta_cm;
}

double Kinema2Body::GetMass(int i) 
{
  switch (i) {
  case 1:
    return kin.M_1;
    break;
  case 2:
    return kin.M_2;
    break;
  case 3:
    return kin.M_3;
    break;
  case 4:
    return kin.M_4;
    break;
  default:
    fprintf(stderr,"Kinema2Body::GetMass No such particle:%d\n",i);
    exit(1);
  }

}

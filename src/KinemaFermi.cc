// -*- C++ -*-

#include "KinemaFermi.hh"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <G4ios.hh>
#include <Randomize.hh>

#include <CLHEP/Units/SystemOfUnits.h>

//_____________________________________________________________________________
KinemaFermi::KinemaFermi( double m1, double m2, double m3, double m4,
			  double *p1, double *p2,double cos_theta )
{
  double ECM;
  //  double vx0, vy0, vz0;            /* unit vector */
  double vx3, vy3, vz3;            /* unit vector */
  double vx4, vy4, vz4;            /* unit vector */
  double theta3, theta4;
  // double theta2,phi5;
  double phi4;
  double phi3;                     /* 2phi(CM system)*/
  double Theta3,Phi3;              /*heta,Ph*/
  double Theta4,Phi4;              /* eta,Ph*/
  // double P_Theta,P_Phi;              /* production angle*/
  double theta1, phi1;
  //
  /// Let's think about vector
  ///
  /*
    G4cout<<"m1 :"<<m1<<G4endl;
    G4cout<<"m2 :"<<m2<<G4endl;
    G4cout<<"m3 :"<<m3<<G4endl;
    G4cout<<"m4 :"<<m4<<G4endl;
    G4cout<<"shhwang test test 1"<<G4endl;
  */

  kin3.M_1 = m1; //--> beam
  kin3.M_2 = m2; //--> target
  kin3.M_3 = m3; //--> K+
  kin3.M_4 = m4; //--> Xi-

  kin3.p_1_lab[0] = p1[0];
  kin3.p_1_lab[1] = p1[1];
  kin3.p_1_lab[2] = p1[2];
  //  kin3.p_1_lab[3] = p1[3];//kinetic energy

  kin3.p_2_lab[0] = p2[0];
  kin3.p_2_lab[1] = p2[1];
  kin3.p_2_lab[2] = p2[2];

  //  kin3.p_2_lab[3] = p2[3];//kinetic energy
  ////////////shhwang ///////////////////////////////////
  /////production angle, total angle will be rotate
  double Plab1=sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);
  double Plab2=sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);

  ////// sum p1 p2
  double p12_sum[4];
  p12_sum[0]=kin3.p_1_lab[0]+kin3.p_2_lab[0];
  p12_sum[1]=kin3.p_1_lab[1]+kin3.p_2_lab[1];
  p12_sum[2]=kin3.p_1_lab[2]+kin3.p_2_lab[2];
  p12_sum[3]=sqrt(p12_sum[0]*p12_sum[0] + p12_sum[1]*p12_sum[1] + p12_sum[2]*p12_sum[2] );

  /*
  G4cout<<"p1:"<<p1[0]<<", "<<p1[1]<<", "<<p1[2]<<G4endl;
  G4cout<<"p2:"<<p2[0]<<", "<<p2[1]<<", "<<p2[2]<<G4endl;
  */

  //  theta1=atan(sqrt(p12_sum[0]*p12_sum[0]+p12_sum[1]*p12_sum[1]),p12_sum[2])*180./CLHEP::pi);
  //  G4cout<<"theta atan2:"<<theta1<<G4endl;


  //  G4cout<<"theta acos:"<<theta1<<G4endl;

  theta1=atan2(sqrt(pow(p12_sum[0],2)+pow(p12_sum[1],2)),p12_sum[2])*180./CLHEP::pi;
  phi1=atan2(p12_sum[1],p12_sum[0])*180./CLHEP::pi;
  // theta1=acos(fabs(p12_sum[2])/p12_sum[3])*180./CLHEP::pi;
  //  phi1=atan(p12_sum[1]/p12_sum[0])*180./CLHEP::pi;

  //  if(phi1<0) phi1+180.;
  //  G4cout<<"phi1:"<<phi1<<G4endl;

  //    vx0 = sin(deg2rad(P_Theta))*cos(deg2rad(P_Phi));
  //    vy0 = sin(deg2rad(P_Theta))*sin(deg2rad(P_Phi));
  //    vz0 = cos(deg2rad(P_Theta));
  //  vx0=vy0=vz0=1.;

  /*
  //////////////check//////////////
  G4cout<<"p1:"<<p1[0]<<", "<<p1[1]<<", "<<p1[2]<<G4endl;
  G4cout<<"p2:"<<p2[0]<<", "<<p2[1]<<", "<<p2[2]<<G4endl;

  G4cout<<"p12_sum:"<<p12_sum[0]<<", "<<p12_sum[1]<<", "<<p12_sum[2]<<", "<<p12_sum[3]<<G4endl;
  G4cout<<"theta1:"<<theta1<<G4endl;
  G4cout<<"phi1:"<<phi1<<G4endl;
  */
  //  G4cout<<"vx0, vy0, vz0:"<<vx0<<", "<<vy0<<", "<<vz0<<G4endl;
  //  CalcDistoribution(vx0, vy0, vz0, &P_Theta, &P_Phi);

  kin3.E_1_lab = p2E(sqrt(pow(kin3.p_1_lab[0],2)+pow(kin3.p_1_lab[1],2)+pow(kin3.p_1_lab[2],2)), kin3.M_1);///beam
  kin3.E_2_lab = p2E(sqrt(pow(kin3.p_2_lab[0],2)+pow(kin3.p_2_lab[1],2)+pow(kin3.p_2_lab[2],2)), kin3.M_2);///target

  ECM = sqrt(pow(kin3.E_1_lab+kin3.E_2_lab,2)
	     -(pow(kin3.p_1_lab[0]+kin3.p_2_lab[0],2)
	       +pow(kin3.p_1_lab[1]+kin3.p_2_lab[1],2)
	       +pow(kin3.p_1_lab[2]+kin3.p_2_lab[2],2)));

  //  do {
  //    kin3.M_res = CLHEP::RandBreitWigner::shoot(m_res, width);
  //  } while (kin3.M_3+kin3.M_4 > kin3.M_res || kin3.M_res > ECM-kin3.M_5);
  //  kin3.M_res = m_res;


  if(ECM<m3+m4){
    G4cout<<"Center of energy less than the mass sum for the produced particles"<<G4endl;
    return;
  }

  kin1 = Kinema2Body(m1, m2, m3, m4);
  //    kin1.SetMomentum(1,Plab1);
  //    kin1.SetMomentum(2,Plab2);
  /////sh.  2body also should change////
  //  cos(theta)=|ab|/

  double theta_p1_psum;
  double theta_p2_psum;
  // double theta_p1_p2;
  if(Plab2==0){
    theta_p1_psum=0.;
    theta_p2_psum=0.;
    // theta_p1_p2=0.;
  }else {
    theta_p1_psum=acos( (p1[0]*p12_sum[0]+p1[1]*p12_sum[1]+p1[2]*p12_sum[2]) /(Plab1*p12_sum[3]));//beam
    theta_p2_psum=acos(  ( p2[0]*p12_sum[0] + p2[1]*p12_sum[1] +p2[2]*p12_sum[2]) /(Plab2*p12_sum[3]));//proton
    // theta_p1_p2=acos((p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2])/(Plab1*Plab2));
  }
  /*
  G4cout<<"theta_p1_psum:"<<theta_p1_psum<<G4endl;
  G4cout<<"theta_p2_psum:"<<theta_p2_psum<<G4endl;
  G4cout<<"comp theta_p1_p2:"<<theta_p1_p2<<":"<<theta_p1_psum+theta_p2_psum<<G4endl;
  */
  /*
  G4cout<<"------------->"<<( p1[0]*p12_sum[0] + p1[1]*p12_sum[1] +p1[2]*p12_sum[2])<<G4endl;
  G4cout<<"------------->"<<Plab1*p12_sum[3]<<G4endl;
  G4cout<<"------------->"<<CLHEP::pi<<G4endl;

  G4cout<<"theta_p1_p2:"<<theta_p1_p2*180./CLHEP::pi<<G4endl;
  G4cout<<"cos(theta_p1_psum)*Plab1:"<<cos(theta_p1_psum)*Plab1<<G4endl;
  G4cout<<"cos(theta_p2_psum)*Plab2:"<<cos(theta_p2_psum)*Plab2<<G4endl;
  //    G4cout<<p12_sum[3]<<":"<<cos(theta_p1_psum)*Plab1+cos(theta_p2_psum)*Plab2<<G4endl;
  G4cout<<"plab1:"<<Plab1<<G4endl;
  G4cout<<"plab2:"<<Plab2<<G4endl;
  G4cout<<"cos(theta_p1_psum)*Plab1:"<<cos(theta_p1_psum)*Plab1<<G4endl;
  G4cout<<"cos(theta_p2_psum)*Plab2:"<<cos(theta_p2_psum)*Plab2<<G4endl;
  */
  kin1.SetMomentum(1,cos(theta_p1_psum)*Plab1);
  kin1.SetMomentum(2,cos(theta_p2_psum)*Plab2);

  kin1.SetTheta(1,theta_p1_psum);
  kin1.SetTheta(2,theta_p2_psum);

  //  kin1.SetMomentum(1,Plab1);
  //  kin1.SetMomentum(2,Plab2);


  //  G4cout<<"comp cos:"<<sin(theta_p1_psum)*Plab1<<", "<<sin(theta_p2_psum)*Plab2<<G4endl;



  /*
  G4cout<<"p12_sum:"<<p12_sum[0]<<", "<<p12_sum[1]<<", "<<p12_sum[2]<<", "<<p12_sum[3]<<G4endl;
  G4cout<<"Plab1:"<<Plab1<<G4endl;
  G4cout<<"cos(theta_p1_psum)*Plab1:"<<cos(theta_p1_psum)*Plab1<<G4endl;
  G4cout<<"input 2body mom:"<<cos(theta_p1_psum)*Plab1+cos(theta_p2_psum)*Plab2<<G4endl;
  G4cout<<"Plab2:"<<Plab2<<G4endl;
  G4cout<<"cos(theta_p2_psum)*Plab2:"<<cos(theta_p2_psum)*Plab2<<G4endl;
  */
  ///////////  kin1.SetThetaCM((double)RandSin()); //set to zero
  kin1.SetThetaCM(acos(cos_theta)*180./CLHEP::pi); //input cross section value

  //  kin1.SetThetaCM(170.); //input cross section value
  //  G4cout<<"rand sin:cos_theta-->"<<(double)RandSin()<<", "<<acos(cos_theta)*180./CLHEP::pi <<G4endl;
  kin1.CalcKinema();
  //  phi3 = (-180.+360.0*(double)CLHEP::RandFlat::shoot());
  phi3 = (360.0*(double)CLHEP::RandFlat::shoot());


  kin3.Theta1CM = kin1.GetThetaCM();
  kin3.Phi1     = phi3;
  /*
  G4cout<<"kin3.Theta1CM:"<<kin3.Theta1CM<<G4endl;
  G4cout<<"kin3.Phi1:"<<kin3.Phi1<<G4endl;
  */

  /* calculate m4 */
  //  theta4 = kin1.GetThetaLab();
  theta3 = kin1.GetThetaLab();
  //  theta3 = kin1.GetPhi();
  /*
    vx3 = sin(deg2rad(theta3))*cos(deg2rad(phi3));
    vy3 = sin(deg2rad(theta3))*sin(deg2rad(phi3));
    vz3 = cos(deg2rad(theta3));
  */

//////////new cal by shhwang
/*
  vx3=sin(deg2rad(theta1))*cos(deg2rad(phi1))*sin(deg2rad(theta3))*cos(deg2rad(phi3))
    +cos(deg2rad(phi1))*cos(deg2rad(theta1))*sin(deg2rad(theta3))*sin(deg2rad(phi3))
    -sin(deg2rad(phi1))*cos(deg2rad(theta3));

  vy3=sin(deg2rad(theta1))*sin(deg2rad(phi1))*sin(deg2rad(theta3))*cos(deg2rad(phi3))
    +sin(deg2rad(phi1))*cos(deg2rad(theta1))*sin(deg2rad(theta3))*sin(deg2rad(phi3))
    +cos(deg2rad(phi1))*cos(deg2rad(theta3));

  vz3=cos(deg2rad(theta1))*sin(deg2rad(theta3))*cos(deg2rad(phi3))
    -sin(deg2rad(theta1))*sin(deg2rad(theta3))*sin(deg2rad(phi3));
*/
/*
////////////test shhwang
  vx3=sin(deg2rad(theta3))*cos(deg2rad(phi3))*sin(deg2rad(theta1))*cos(deg2rad(phi1))
    +cos(deg2rad(phi3))*cos(deg2rad(theta3))*sin(deg2rad(theta1))*sin(deg2rad(phi1))
    -sin(deg2rad(phi3))*cos(deg2rad(theta1));

  vy3=sin(deg2rad(theta3))*sin(deg2rad(phi3))*sin(deg2rad(theta1))*cos(deg2rad(phi1))
    +sin(deg2rad(phi3))*cos(deg2rad(theta3))*sin(deg2rad(theta1))*sin(deg2rad(phi1))
    +cos(deg2rad(phi3))*cos(deg2rad(theta1));

  vz3=cos(deg2rad(theta3))*sin(deg2rad(theta1))*cos(deg2rad(phi1))
    -sin(deg2rad(theta3))*sin(deg2rad(theta1))*sin(deg2rad(phi1));
*/


  ///////////from 3body, eular rotation

  vx3 = (cos(deg2rad(phi1))*cos(deg2rad(theta3))*sin(deg2rad(theta1)) +
    cos(deg2rad(theta1))*cos(deg2rad(phi1))*cos(deg2rad(phi3))*sin(deg2rad(theta3)) -
    sin(deg2rad(phi1))*sin(deg2rad(phi3))*sin(deg2rad(theta3)));
  vy3 = (sin(deg2rad(phi1))*cos(deg2rad(theta3))*sin(deg2rad(theta1)) +
    cos(deg2rad(theta1))*sin(deg2rad(phi1))*cos(deg2rad(phi3))*sin(deg2rad(theta3))+
    cos(deg2rad(phi1))*sin(deg2rad(phi3))*sin(deg2rad(theta3)));
  vz3 = cos(deg2rad(theta3))*cos(deg2rad(theta1)) -
    sin(deg2rad(theta1))*cos(deg2rad(phi3))*sin(deg2rad(theta3));


  //  G4cout<<"|v3|:"<<sqrt(vx3*vx3+vy3*vy3+vz3*vz3)<<G4endl;
  /*
  vx3 = sin(deg2rad(theta3))*cos(deg2rad(phi3));
  vy3 = sin(deg2rad(theta3))*sin(deg2rad(phi3));
  vz3 = cos(deg2rad(theta3));
  */
  /*
  G4cout<<"theta1:phi1:"<<theta1<<" : "<<phi1 << G4endl;
  G4cout<<"vx3, vy3, vz3:"<<vx3<<", "<<vy3<<", "<<vz3<<G4endl;
  */
  CalcDistoribution(vx3, vy3, vz3, &Theta3, &Phi3);
  //  G4cout<<"vx3, vy3, vz3:"<<vx3<<", "<<vy3<<", "<<vz3<<G4endl;
  kin3.E_3_lab = kin1.GetEnergyLab(3);
  kin3.p_3_lab = kin1.GetMomentumLab(3);

  kin3.P_3_lab[0] = kin3.p_3_lab*vx3;
  kin3.P_3_lab[1] = kin3.p_3_lab*vy3;
  kin3.P_3_lab[2] = kin3.p_3_lab*vz3;

  /*
  G4cout<<"kin3.p_3_lab:"<<kin3.p_3_lab<<G4endl;
  G4cout<<"kin3.E_3_lab:"<<kin3.E_3_lab<<G4endl;
  G4cout<<"kin3.E_3_lab-mass:"<<sqrt(pow(kin3.E_3_lab,2)-pow(m3,2))<<G4endl;
  */
  kin3.theta3 = Theta3;
  kin3.phi3 = Phi3;


  theta4 = -kin1.GetPhiLab();
  phi4= phi3;
  /*from three body original
  vx4 = cos(deg2rad(phi5))*cos(deg2rad(theta2))*sin(deg2rad(theta1)) +
     cos(deg2rad(theta1))*cos(deg2rad(phi5))*cos(deg2rad(phi3))*sin(deg2rad(theta2)) -
        sin(deg2rad(phi5))*sin(deg2rad(phi3))*sin(deg2rad(theta2));

  vy4 = -sin(deg2rad(phi5))*cos(deg2rad(theta2))*sin(deg2rad(theta1)) -
    cos(deg2rad(theta1))*sin(deg2rad(phi5))*cos(deg2rad(phi3))*sin(deg2rad(theta2)) -
    cos(deg2rad(phi5))*sin(deg2rad(phi3))*sin(deg2rad(theta2));
  vz4 = cos(deg2rad(theta2))*cos(deg2rad(theta1)) -
      sin(deg2rad(theta1))*cos(deg2rad(phi3))*sin(deg2rad(theta2));
  */

  ///////////from 3body, eular rotation

  vx4 = (cos(deg2rad(phi1))*cos(deg2rad(theta4))*sin(deg2rad(theta1)) +
    cos(deg2rad(theta1))*cos(deg2rad(phi1))*cos(deg2rad(phi4))*sin(deg2rad(theta4)) -
	 sin(deg2rad(phi1))*sin(deg2rad(phi4))*sin(deg2rad(theta4)));
  vy4 = (sin(deg2rad(phi1))*cos(deg2rad(theta4))*sin(deg2rad(theta1)) +
    cos(deg2rad(theta1))*sin(deg2rad(phi1))*cos(deg2rad(phi4))*sin(deg2rad(theta4)) +
	 cos(deg2rad(phi1))*sin(deg2rad(phi4))*sin(deg2rad(theta4)));
  vz4 = cos(deg2rad(theta4))*cos(deg2rad(theta1)) -
    sin(deg2rad(theta1))*cos(deg2rad(phi4))*sin(deg2rad(theta4));


  /*
    vx4 = sin(deg2rad(theta4))*cos(deg2rad(phi4));
    vy4 = sin(deg2rad(theta4))*sin(deg2rad(phi4));
    vz4 = cos(deg2rad(theta4));
  */

  //  G4cout<<"theta3,theta4:"<<theta3<<", "<<theta4<<G4endl;

  /*
//////////new cal by shhwang
  vx4=sin(deg2rad(theta1))*cos(deg2rad(phi1))*sin(deg2rad(theta4))*cos(deg2rad(phi4))
    +cos(deg2rad(phi1))*cos(deg2rad(theta1))*sin(deg2rad(theta4))*sin(deg2rad(phi4))
    -sin(deg2rad(phi1))*cos(deg2rad(theta4));
  vy4=sin(deg2rad(theta1))*sin(deg2rad(phi1))*sin(deg2rad(theta4))*cos(deg2rad(phi4))
    +sin(deg2rad(phi1))*cos(deg2rad(theta1))*sin(deg2rad(theta4))*sin(deg2rad(phi4))
    +cos(deg2rad(phi1))*cos(deg2rad(theta4));
  vz4=cos(deg2rad(theta1))*sin(deg2rad(theta4))*cos(deg2rad(phi4))
    -sin(deg2rad(theta1))*sin(deg2rad(theta4))*sin(deg2rad(phi4));
  */

  //  G4cout<<"|v4|:"<<sqrt(vx4*vx4+vy4*vy4+vz4*vz4)<<G4endl;

  //  G4cout<<"vx4, vy4, vz4:"<<vx4<<", "<<vy4<<", "<<vz4<<G4endl;
  CalcDistoribution(vx4, vy4, vz4, &Theta4, &Phi4);

  //  G4cout<<"vx4, vy4, vz4:"<<vx4<<", "<<vy4<<", "<<vz4<<G4endl;
  //  G4cout<<"atan2"<<atan2(vy5,vx5)/CLHEP::pi*180.<<G4endl;
  kin3.E_4_lab = kin1.GetEnergyLab(4);
  kin3.p_4_lab = kin1.GetMomentumLab(4);
  //  G4cout<<"kin3.p_4_lab:"<<kin3.p_4_lab<<G4endl;
  //  G4cout<<"kin3.E_4_lab:"<<kin3.E_4_lab<<G4endl;
  kin3.P_4_lab[0] = kin3.p_4_lab*vx4;
  kin3.P_4_lab[1] = kin3.p_4_lab*vy4;
  kin3.P_4_lab[2] = kin3.p_4_lab*vz4;
  kin3.theta4 = Theta4;
  kin3.phi4 = Phi4;

  //  G4cout<<"kin3.P_3_lab:kin3.P_4_lab"<<G4endl;
  //  G4cout<<kin3.P_3_lab[0]<<":"<<kin3.P_4_lab[0]<<G4endl;
  //  G4cout<<kin3.P_3_lab[1]<<":"<<kin3.P_4_lab[1]<<G4endl;
  //  G4cout<<kin3.P_3_lab[2]<<":"<<kin3.P_4_lab[2]<<G4endl;
  //  G4cout<<"-----------------------------------------------"<<G4endl;
  //Dump();
}

//_____________________________________________________________________________
KinemaFermi::~KinemaFermi( void )
{
}

//_____________________________________________________________________________
double
KinemaFermi::p2E( double p,double m )
{
  return sqrt(p*p + m*m);
}

//_____________________________________________________________________________
void
KinemaFermi::CalcDistoribution( double unitx, double unity, double /* unitz */,
				double *theta, double *phi )
{
  *theta = rag2deg(acos(unitx));
  *phi=rag2deg(atan2(unity,unitx));
  /*  if (unity>=0.0 && unitz>0.0)
    *phi = rag2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity<0.0 && unitz>=0.0)
    *phi = rag2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity<=0.0 && unitz<0.0)
    *phi = 360.0-rag2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity>0.0 && unitz<=0.0)
    *phi = 360.0-rag2deg(acos(unity/sin(deg2rad(*theta))));
  else {
    fprintf(stderr,
	  "KinemaFermi::CalcDistribution No such reagion unity=%f, unitz=%f\n",
	    unity, unitz);
    Dump();
    exit(1);
  }
  */
  return;
}

//_____________________________________________________________________________
double
KinemaFermi::deg2rad( double theta )
{
  return CLHEP::pi*theta/180.0;
}

//_____________________________________________________________________________
double
KinemaFermi::rag2deg( double rag )
{
  return 360.0 * rag/ (2.0 * CLHEP::pi);
}

//_____________________________________________________________________________
double
KinemaFermi::RandSin( void )
{
  int success=0;
  double x,fx;

  do {
    x = 180.0 * (double)CLHEP::RandFlat::shoot();
    //x = 180.0 * (double)rand()/(RAND_MAX+1.0);
    fx = sin(deg2rad(x));
    if (fx >= (double)CLHEP::RandFlat::shoot())
      success = 1;
  } while (success==0);

  return x;
}

//_____________________________________________________________________________
void
KinemaFermi::Dump( void )
{
  printf("======KinemaFermi Dump======\n");
  printf("--Particle1--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_1, kin3.p_1_lab[2], kin3.E_1_lab);
  printf("--Particle2--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_2, kin3.p_2_lab[2], kin3.E_2_lab);
  printf("--Particle3--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_3, kin3.p_3_lab, kin3.E_3_lab);
  printf("momentum=(%f, %f, %f), (theta, phi)=(%f, %f)\n",
	 kin3.P_3_lab[0], kin3.P_3_lab[1], kin3.P_3_lab[2],
	 kin3.theta3, kin3.phi3);
  printf("--Particle4--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_4, kin3.p_4_lab, kin3.E_4_lab);
  printf("momentum=(%f, %f, %f), (theta, phi)=(%f, %f)\n",
	 kin3.P_4_lab[0], kin3.P_4_lab[1], kin3.P_4_lab[2],
	 kin3.theta4, kin3.phi4);

  printf("Energy:E1+E2=%f, E3+E4=%f\n", kin3.E_1_lab+kin3.E_2_lab,
	 kin3.E_3_lab+kin3.E_4_lab);
  printf("Momentum: x-> p1+p2=%f, p3+p4=%f\n", kin3.p_1_lab[2]+kin3.p_2_lab[2],
	 kin3.P_3_lab[0]+kin3.P_4_lab[0]);
  printf("          y-> p3+p4=%f\n",
	 kin3.P_3_lab[1]+kin3.P_4_lab[1]);
  printf("          z-> p3+p4=%f\n",
	 kin3.P_3_lab[2]+kin3.P_4_lab[2]);

  return;
}

//_____________________________________________________________________________
double
KinemaFermi::GetEnergy( int i )
{
  switch (i) {
  case 1:
    return kin3.E_1_lab;
    break;
  case 2:
    return kin3.E_2_lab;
    break;
  case 3:
    return kin3.E_3_lab;
    break;
  case 4:
    return kin3.E_4_lab;
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetEnergy No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaFermi::GetMomentum( int i )
{
  switch (i) {
  case 1:
    return kin3.p_1_lab[2];
    break;
  case 2:
    return kin3.p_2_lab[2];
    break;
  case 3:
    return kin3.p_3_lab;
    break;
  case 4:
    return kin3.p_4_lab;
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
void
KinemaFermi::GetMomentum( int i, double *mom )
{
  switch (i) {
  case 1:
    mom[0] = kin3.p_1_lab[0];
    mom[1] = kin3.p_1_lab[1];
    mom[2] = kin3.p_1_lab[2];
    break;
  case 2:
    mom[0] = kin3.p_2_lab[0];
    mom[1] = kin3.p_2_lab[1];
    mom[2] = kin3.p_2_lab[2];
    break;
  case 3:
    mom[0] = kin3.P_3_lab[0];
    mom[1] = kin3.P_3_lab[1];
    mom[2] = kin3.P_3_lab[2];
    break;
  case 4:
    mom[0] = kin3.P_4_lab[0];
    mom[1] = kin3.P_4_lab[1];
    mom[2] = kin3.P_4_lab[2];
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaFermi::GetTheta( int i )
{
  switch (i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return kin3.theta3;
    break;
  case 4:
    return kin3.theta4;
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetTheta No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaFermi::GetPhi( int i )
{
  switch (i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return kin3.phi3;
    break;
  case 4:
    return kin3.phi4;
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetPhi No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaFermi::GetThetaCM( int i )
{
  switch (i) {
  case 1:
    return kin3.Theta1CM;
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetThetaCM index should be 1 or 2 ->%d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaFermi::GetPhiCM( int i )
{
  switch (i) {
  case 1:
    return kin3.Phi1;
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetPhiCM index should be 1 or 2 ->%d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
void
KinemaFermi::RotateMom( int i, double deg, double *mom )
{
  double Sin,Cos;

  Sin=sin(deg2rad(deg));
  Cos=cos(deg2rad(deg));
  switch (i) {
  case 3:
    mom[0] = Cos*kin3.P_3_lab[0] - Sin*kin3.P_3_lab[1];
    mom[1] = Sin*kin3.P_3_lab[0] + Cos*kin3.P_3_lab[1];
    mom[2] = kin3.P_3_lab[2];
    break;
  case 4:
    mom[0] = Cos*kin3.P_4_lab[0] - Sin*kin3.P_4_lab[1];
    mom[1] = Sin*kin3.P_4_lab[0] + Cos*kin3.P_4_lab[1];
    mom[2] = kin3.P_4_lab[2];
    break;
  default:
    fprintf(stderr, "KinemaFermi::RotateMom should be 3,4->%d\n",i);
    exit(1);
  }
}

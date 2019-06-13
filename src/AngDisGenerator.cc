//-----------------------------------------------------------
// AngDisGenerator.cc
// for the angular distribution for the E27 experiment
//-----------------------------------------------------------

#include "AngDisGenerator.hh"
#include "Randomize.hh"

#include <cmath>

G4double URand(){ return  G4UniformRand(); }


AngDisGenerator::AngDisGenerator( double cost1, double cost2 )
  : cost1_(cost1), cost2_(cost2)
{}

AGSWave::AGSWave( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGSWave::GenerateDirection( void ) const
{
  double cost=cost1_+URand()*(cost2_-cost1_);
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

AGPWaveFP::AGPWaveFP( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGPWaveFP::GenerateDirection( void ) const
{
  double a=2./(cost2_-cost1_)/(2.+cost2_+cost1_);
  double cost=-1.+sqrt((cost1_+1.)*(cost1_+1.)+2.*URand()/a);
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}


AGPWaveBP::AGPWaveBP( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGPWaveBP::GenerateDirection( void ) const
{
  double a=2./(cost2_-cost1_)/(2.-cost2_-cost1_);
  double cost=1.-sqrt((cost1_-1.)*(cost1_-1.)-2.*URand()/a);
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

AGDWave1::AGDWave1( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGDWave1::GenerateDirection( void ) const
{
  double cost=(cost2_-cost1_)*URand()+cost1_;
  double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

AGSigma1385Zero::AGSigma1385Zero( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGSigma1385Zero::GenerateDirection( void ) const
{

  //cost == cost1_ ~ cost2_
  double cost=(cost2_-cost1_)*URand()+cost1_;
  double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

AGSigma1385Plus::AGSigma1385Plus( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGSigma1385Plus::GenerateDirection( void ) const
{

  //cost == cost1_ ~ cost2_
  double cost=(cost2_-cost1_)*URand()+cost1_;
  double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

AGLambda1405::AGLambda1405( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGLambda1405::GenerateDirection( void ) const
{

  //cost == cost1_ ~ cost2_
  double cost=(cost2_-cost1_)*URand()+cost1_;
  double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}


AGLambda::AGLambda( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGLambda::GenerateDirection( void ) const
{

  //cost == cost1_ ~ cost2_
  double cost=(cost2_-cost1_)*URand()+cost1_;
  double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}



AGSigmaZ::AGSigmaZ( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}
G4ThreeVector AGSigmaZ::GenerateDirection( void ) const
{

  //cost == cost1_ ~ cost2_
  double cost=(cost2_-cost1_)*URand()+cost1_;
  double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

AGSigmaP::AGSigmaP( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGSigmaP::GenerateDirection( void ) const
{

  //cost == cost1_ ~ cost2_
  double cost=(cost2_-cost1_)*URand()+cost1_;
  double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

AGPol::AGPol( double cost1, double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGPol::GenerateDirection( void ) const
{
  double a=2./(cost2_-cost1_)/(2.+cost2_+cost1_);
  double cost=-1.+sqrt((cost1_+1.)*(cost1_+1.)+2.*URand()/a);
  double sint=sqrt(1.-cost*cost);
  double phi=URand()*acos(-1.)*2.;
  double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}




//-----------------------------------------------------------
// AngDisGenerator.hh
// for the angular distribution for the E27 experiment
//-----------------------------------------------------------

#ifndef AngDisGenerator_H
#define AngDisGenerator_H

#include "G4ThreeVector.hh"

class AngDisGenerator
{
public:
  AngDisGenerator( double cost1=1.0, double cost2=-1.0 );
  virtual ~AngDisGenerator() {};

  virtual G4ThreeVector GenerateDirection( void ) const = 0;
  virtual double GetDfuncVal( double x) const = 0;

protected:
  double cost1_, cost2_;

};

class AGSWave : public AngDisGenerator
{
  // D(x)=1/2 
public:
  AGSWave( double cost1=1.0, double cost2=-1.0 );
  ~AGSWave() {};

  G4ThreeVector GenerateDirection( void ) const;

  double GetDfuncVal( double x) const{ return Dfunc(x);}
private:
  inline double Dfunc( double x ) const { return 0.5; }

};

typedef AGSWave AGUniform;


class AGPWaveFP : public AngDisGenerator
{
  // D(x)=1/2*(1+x)
public:
  AGPWaveFP( double cost1=1.0, double cost2=-1.0 );
  ~AGPWaveFP() {};

  G4ThreeVector GenerateDirection( void ) const;

  double GetDfuncVal( double x) const{ return Dfunc(x);}
private:
  inline double Dfunc( double x ) const { return 0.5*(x+1.); }
};


class AGPWaveBP : public AngDisGenerator
{
  // D(x)=1/2*(1-x)
public:
  AGPWaveBP( double cost1=1.0, double cost2=-1.0 );
  ~AGPWaveBP() {};

  G4ThreeVector GenerateDirection( void ) const;
  
  double GetDfuncVal( double x) const{ return Dfunc(x);}
private:
  inline double Dfunc( double x ) const { return 0.5*(1.-x); }

};

class AGDWave1 : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGDWave1( double cost1=1.0, double cost2=-1.0 );
  ~AGDWave1() {};

  G4ThreeVector GenerateDirection( void ) const;
  double GetDfuncVal( double x) const{ return Dfunc(x);}

  
private:
  inline double Dfunc( double x ) const { return 0.5*(x*x+1.); }
};

class AGSigma1385Zero : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGSigma1385Zero( double cost1=1.0, double cost2=-1.0 );
  ~AGSigma1385Zero() {};

  G4ThreeVector GenerateDirection( void ) const;
  double GetDfuncVal( double x) const{ return Dfunc(x);}
  
private:
  inline double Dfunc( double x ) const { 
    return 
      (6.1+
       1.2*x+
       (-5.2)/2.*(3*x*x-1))/10.;
  }
};


class AGSigma1385Plus : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGSigma1385Plus( double cost1=1.0, double cost2=-1.0 );
  ~AGSigma1385Plus() {};

  G4ThreeVector GenerateDirection( void ) const;
  double GetDfuncVal( double x) const{ return Dfunc(x);}
  
private:
  inline double Dfunc( double x ) const { 
    return 
      (30.6394+
       (11.7159)*x+
       (-9.52849)/2.*(3*x*x-1)+
       (-13.9936)/2.*(5*x*x*x-3*x)+
       (-7.13241)/8.*(35*x*x*x*x-30*x*x+3)+
       (-3.61133)/8.*(63*x*x*x*x*x-70*x*x*x+15*x)+
       (-2.21394)/16.*(231*x*x*x*x*x*x-315*x*x*x*x+105*x*x-5))/51.;
    
 
  }
};



class AGLambda1405 : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGLambda1405( double cost1=1.0, double cost2=-1.0 );
  ~AGLambda1405() {};

  G4ThreeVector GenerateDirection( void ) const;
  double GetDfuncVal( double x) const{ return Dfunc(x);}
  
private:
  inline double Dfunc( double x ) const { 
    return
      
      (1.64+
      1.02*x+
      1.54/2.*(3*x*x-1)+
      0.96/2.*(5*x*x*x-3*x)+
      0.55/8.*(35*x*x*x*x-30*x*x+3)+
       (-0.56)/8.*(63*x*x*x*x*x-70*x*x*x+15*x))/6.;
    
  }
};


class AGLambda : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGLambda( double cost1=1.0, double cost2=-1.0 );
  ~AGLambda() {};

  G4ThreeVector GenerateDirection( void ) const;
  double GetDfuncVal( double x) const{ return Dfunc(x);}
  
private:
  inline double Dfunc( double x ) const { 
    return
      
      (13.9+
       9.8*x+
       20.1/2.*(3*x*x-1)-
       1.3/2.*(5*x*x*x-3*x)+
       12.8/8.*(35*x*x*x*x-30*x*x+3))/53.;
    
  }
};




class AGSigmaZ : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGSigmaZ( double cost1=1.0, double cost2=-1.0 );
  ~AGSigmaZ() {};

  G4ThreeVector GenerateDirection( void ) const;
  double GetDfuncVal( double x) const{ return Dfunc(x);}
  
private:
  inline double Dfunc( double x ) const { 
    return
      
      (9.6+
       5.3*x+
       8.0/2.*(3*x*x-1)+
       16.2/2.*(5*x*x*x-3*x)+
       6.7/8.*(35*x*x*x*x-30*x*x+3)+
       10.4/8.*(63*x*x*x*x*x-70*x*x*x+15*x)+
       8.0/16.*(231*x*x*x*x*x*x-315*x*x*x*x+105*x*x-5))/61.;
  }
};



class AGSigmaP : public AngDisGenerator
{
  // D(x)= 0.5*(x*x+1)
public:
  AGSigmaP( double cost1=1.0, double cost2=-1.0 );
  ~AGSigmaP() {};

  G4ThreeVector GenerateDirection( void ) const;
  double GetDfuncVal( double x) const{ return Dfunc(x);}
  
private:
  inline double Dfunc( double x ) const { 
    return
      
      (37.5+
       7.5*x+
       19.0/2.*(3*x*x-1)+
       23.1/2.*(5*x*x*x-3*x)+
       40.4/8.*(35*x*x*x*x-30*x*x+3)+
       14.6/8.*(63*x*x*x*x*x-70*x*x*x+15*x)+
       31.3/16.*(231*x*x*x*x*x*x-315*x*x*x*x+105*x*x-5))/190.;
  }
};


class AGPol : public AngDisGenerator
{
public:
  AGPol( double cost1=1.0, double cost2=-1.0 );
  ~AGPol() {};

  G4ThreeVector GenerateDirection( void ) const;

  double GetDfuncVal( double x) const{ return Dfunc(x);}
private:
  inline double Dfunc( double x ) const { return (1 - 0.4*x)/2.; }
};






#endif

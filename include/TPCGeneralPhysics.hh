// -*- C++ -*-

#ifndef TPC_GENERAL_PHYSICS_H
#define TPC_GENERAL_PHYSICS_H

#include <G4VPhysicsConstructor.hh>

//_____________________________________________________________________________
class TPCGeneralPhysics : public G4VPhysicsConstructor
{
public:
  TPCGeneralPhysics( const G4String& name="General physics" );
  virtual ~TPCGeneralPhysics( void );
  virtual void ConstructParticle( void );
  virtual void ConstructProcess( void );
};

#endif

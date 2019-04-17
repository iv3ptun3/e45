// ====================================================================
//   TPCGeneralPhysics.hh
//
// ====================================================================
#ifndef TPC_GENERAL_PHYSICS_H
#define TPC_GENERAL_PHYSICS_H

#include "G4VPhysicsConstructor.hh"

class TPCGeneralPhysics : public G4VPhysicsConstructor {
public: 
  TPCGeneralPhysics(const G4String& name = "General physics");
  virtual ~TPCGeneralPhysics();
  
  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
  
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();

};

#endif


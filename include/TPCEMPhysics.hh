// ====================================================================
//   TPCEMPhysics.hh
//
// ====================================================================
#ifndef TPC_EM_PHYSICS_H
#define TPC_EM_PHYSICS_H

#include "G4VPhysicsConstructor.hh"

class TPCEMPhysics : public G4VPhysicsConstructor {
public: 
  TPCEMPhysics(const G4String& name ="EM physics");
  virtual ~TPCEMPhysics();
  
  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
  
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();
};

#endif

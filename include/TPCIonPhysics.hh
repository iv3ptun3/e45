// ====================================================================
//   TPCIonPhysics.hh
//
// ====================================================================
#ifndef TPC_ION_PHYSICS_H
#define TPC_ION_PHYSICS_H

#include "G4VPhysicsConstructor.hh"

class TPCIonPhysics : public G4VPhysicsConstructor {
public: 
  TPCIonPhysics(const G4String& name ="Ion physics");
  virtual ~TPCIonPhysics();
  
  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
  
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();
};

#endif

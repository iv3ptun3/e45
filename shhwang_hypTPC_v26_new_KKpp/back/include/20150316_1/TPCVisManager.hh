// ====================================================================
//   TPCVisManager.hh
//
// ====================================================================
#ifndef TPC_VIS_MANAGER_H
#define TPC_VIS_MANAGER_H

#include "G4VisManager.hh"

class TPCVisManager : public G4VisManager {
public:
  TPCVisManager();
  virtual ~TPCVisManager();
  
private:
  virtual void RegisterGraphicsSystems();

};

#endif

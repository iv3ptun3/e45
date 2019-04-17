// ====================================================================
//   TPCDetectorConstruction.hh
//
// ====================================================================
#ifndef TPC_DETECTOR_CONSTRUCTION_H
#define TPC_DETECTOR_CONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"
class E42_Cham;
class TPCDetectorConstruction : public G4VUserDetectorConstruction {
public:
  TPCDetectorConstruction();
  ~TPCDetectorConstruction();
private:
  E42_Cham* par_cham;
  void ConstructForwardSpectrometer();

  // implement it.
  virtual G4VPhysicalVolume* Construct(); 
};

#endif

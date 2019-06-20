// -*- C++ -*-

#ifndef TPC_DETECTOR_CONSTRUCTION_HH
#define TPC_DETECTOR_CONSTRUCTION_HH

#include <G4VUserDetectorConstruction.hh>

//_____________________________________________________________________________
class TPCDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  TPCDetectorConstruction( void );
  ~TPCDetectorConstruction( void );

private:
  void ConstructElements( void );

private:
  void ConstructForwardSpectrometer( void );
  virtual G4VPhysicalVolume* Construct( void );
};

#endif

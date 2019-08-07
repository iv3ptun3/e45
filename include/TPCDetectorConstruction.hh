// -*- C++ -*-

#ifndef TPC_DETECTOR_CONSTRUCTION_HH
#define TPC_DETECTOR_CONSTRUCTION_HH

#include <map>

#include <G4VUserDetectorConstruction.hh>
#include <G4RotationMatrix.hh>
#include <G4String.hh>

class G4Element;
class G4Material;
class G4LogicalVolume;
class G4PVPlacement;
class TPCSDCSD;

//_____________________________________________________________________________
class TPCDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  static G4String ClassName( void );
  TPCDetectorConstruction( void );
  ~TPCDetectorConstruction( void );

private:
  G4int                           m_experiment;
  std::map<G4String, G4Element*>  m_element_map;
  std::map<G4String, G4Material*> m_material_map;
  G4LogicalVolume*                m_world_lv;
  G4LogicalVolume*                m_tpc_lv;
  G4double                        m_rotation_angle;
  G4RotationMatrix*               m_rotation_matrix;
  TPCSDCSD*                       m_sdc_sd;

private:
  // Materials
  void ConstructElements( void );
  void ConstructMaterials( void );
  // Detectors
  void ConstructFTOF( void );
  void ConstructHTOF( void );
  void ConstructHypTPC( void );
  void ConstructK18BeamlineSpectrometer( void );
  void ConstructKuramaMagnet( void );
  void ConstructNBAR( void );
  void ConstructPVAC2( void );
  void ConstructSCH( void );
  void ConstructSDC1( void );
  void ConstructSDC2( void );
  void ConstructSDC3( void );
  void ConstructShsMagnet( void );
  void ConstructTarget( void );
  void ConstructWC( void );

private:
  virtual G4VPhysicalVolume* Construct( void );
};

//_____________________________________________________________________________
inline G4String
TPCDetectorConstruction::ClassName( void )
{
  static G4String s_name("TPCDetectorConstruction");
  return s_name;
}

#endif

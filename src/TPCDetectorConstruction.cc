// -*- C++ -*-

#include "TPCDetectorConstruction.hh"

#include <G4Box.hh>
#include <G4ChordFinder.hh>
#include <G4Element.hh>
#include <G4FieldManager.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4Polyhedra.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4SDManager.hh>
#include <G4SubtractionSolid.hh>
#include <G4Transform3D.hh>
#include <G4TransportationManager.hh>
#include <G4Trd.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4UserLimits.hh>
#include <G4VisAttributes.hh>
#include <G4TwoVector.hh>
#include <G4GenericTrap.hh>

#include "BeamMan.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "TPCACSD.hh"
#include "TPCBH2SD.hh"
#include "TPCBVHSD.hh"
#include "TPCField.hh"
#include "TPCFTOFSD.hh"
#include "TPCHTOFSD.hh"
#include "TPCLACSD.hh"
#include "TPCNBARSD.hh"
#include "TPCPadSD.hh"
#include "TPCSCHSD.hh"
#include "TPCSDCSD.hh"
#include "TPCTargetSD.hh"
#include "TPCTargetVPSD.hh"
#include "TPCVPSD.hh"
#include "TPCWCSD.hh"
#include "TPCVC1SD.hh"
#include "TPCVC2SD.hh"
#include "padHelper.hh"

namespace
{
  using CLHEP::cm;
  using CLHEP::cm3;
  using CLHEP::deg;
  using CLHEP::g;
  using CLHEP::kelvin;
  using CLHEP::m;
  using CLHEP::mg;
  using CLHEP::mm;
  using CLHEP::mole;
  using CLHEP::STP_Temperature;
  using CLHEP::universe_mean_density;
  const auto& gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
  TPCField* myField = nullptr;
  // color
  const G4Colour AQUA( 0.247, 0.8, 1.0 );
  const G4Colour ORANGE( 1.0, 0.55, 0.0 );
  const G4Colour LAVENDER( 0.901, 0.901, 0.98 );
  const G4Colour MAROON( 0.5, 0.0, 0.0 );
  const G4Colour PINK( 1.0, 0.753, 0.796 );
 	const auto& DiscardData = gConf.Get<G4bool> ("DiscardData");
 	const auto& Generator = gConf.Get<G4int> ("Generator");
}

//_____________________________________________________________________________
TPCDetectorConstruction::TPCDetectorConstruction( void )
  : m_experiment( gConf.Get<Int_t>("Experiment") ),
    m_element_map(),
    m_material_map(),
    m_world_lv(),
    m_tpc_lv(),
    m_rotation_angle( gConf.Get<Double_t>("SpectrometerAngle")*deg ),
    m_rotation_matrix( new G4RotationMatrix ),
    m_sdc_sd()
{
  m_rotation_matrix->rotateY( - m_rotation_angle );
}

//_____________________________________________________________________________
TPCDetectorConstruction::~TPCDetectorConstruction( void )
{
}

//_____________________________________________________________________________
G4VPhysicalVolume*
TPCDetectorConstruction::Construct( void )
{
  ConstructElements();
  ConstructMaterials();
  auto world_solid = new G4Box( "WorldSolid", 10.*m/2, 6.*m/2, 16.*m/2 );
  m_world_lv = new G4LogicalVolume( world_solid, m_material_map["Air"],
				    "WorldLV" );
  m_world_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );
  auto world_pv = new G4PVPlacement( nullptr, G4ThreeVector(), m_world_lv,
				     "WorldPV", nullptr, false, 0 );

  myField = new TPCField;
  auto transMan = G4TransportationManager::GetTransportationManager();
  auto fieldMan = transMan->GetFieldManager();
  fieldMan->SetDetectorField( myField );
  fieldMan->CreateChordFinder( myField );
  // fieldMan->GetChordFinder()->SetDeltaChord( 1.e-3*mm );

  if( gConf.Get<G4int>("Generator") == 10 )
    ConstructK18BeamlineSpectrometer();

#if 1
  ConstructAreaTent();
#endif
#if 1 
  ConstructBC3();
  ConstructBC4();
#endif
#if 1
  ConstructBH2();
#endif

#if 1
  ConstructShsMagnet();
  ConstructHypTPC();
  ConstructHTOF();
#endif
#if 0
  ConstructPVAC2();
  ConstructNBAR();
#endif
#if 0 // at E42 1, modified by HeeJeong BYEON
  if( gConf.Get<G4int>("ConstructKurama") ){
    ConstructKuramaMagnet();
    ConstructSDC1();
    ConstructSDC2();
    ConstructSCH();
  	ConstructBVH();
    ConstructSDC3();
    ConstructSDC4();
    ConstructFTOF();
    ConstructLAC();
    ConstructWC();
  }
#endif
#if 1
  ConstructVC1();
  ConstructVC2();
#endif

  myField->Initialize();

  return world_pv;
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructElements( void )
{
  /* G4Element(name, symbol, Z, A) */
  G4String name, symbol;
  G4double Z, A;
  name = "Hydrogen";
  m_element_map[name] = new G4Element( name, symbol="H",  Z=1.,
				       A=1.00794 *g/mole );
  name = "Carbon";
  m_element_map[name] = new G4Element( name, symbol="C",  Z=6.,
				       A=12.011 *g/mole );
  name = "Nitrogen";
  m_element_map[name] = new G4Element( name, symbol="N",  Z=7.,
				       A=14.00674 *g/mole );
  name = "Oxygen";
  m_element_map[name] = new G4Element( name, symbol="O",  Z=8.,
				       A=15.9994 *g/mole );
  name = "Sodium";
  m_element_map[name] = new G4Element( name, symbol="Na", Z=11.,
				       A=22.989768 *g/mole );
  name = "Silicon";
  m_element_map[name] = new G4Element( name, symbol="Si", Z=14.,
				       A=28.0855 *g/mole );
  name = "Phoshorus";
  m_element_map[name] = new G4Element( name, symbol="P", Z=15.,
				       A=30.973762 *g/mole );
  name = "Sulfur";
  m_element_map[name] = new G4Element( name, symbol="S", Z=16.,
				       A=32.066 *g/mole );
  name = "Argon";
  m_element_map[name] = new G4Element( name, symbol="Ar", Z=18.,
				       A=39.948 *g/mole );
  name = "Chrominum";
  m_element_map[name] = new G4Element( name, symbol="Cr", Z=24.,
				       A=51.9961 *g/mole );
  name = "Manganese";
  m_element_map[name] = new G4Element( name, symbol="Mn", Z=25.,
				       A=54.93805 *g/mole );
  name = "Iron";
  m_element_map[name] = new G4Element( name, symbol="Fe", Z=26.,
				       A=55.847 *g/mole );
  name = "Nickel";
  m_element_map[name] = new G4Element( name, symbol="Ni", Z=28.,
				       A=58.69 *g/mole );
  name = "Iodine";
  m_element_map[name] = new G4Element( name, symbol="I",  Z=53.,
				       A=126.90447 *g/mole );
  name = "Cesium";
  m_element_map[name] = new G4Element( name, symbol="Cs", Z=55.,
				       A=132.90543 *g/mole );
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructMaterials( void )
{
  G4cout << FUNC_NAME << G4endl;
  /*
    G4Material( name, density, nelement, state, temperature, pressure );
    G4Material( name, z, a, density, state, temperature, pressure );
  */
  G4String name;
  G4double Z, A, density, massfraction;
  G4int natoms, nel, ncomponents;
  const G4double room_temp = STP_Temperature + 20.*kelvin;
  // Vacuum
  name = "Vacuum";
  m_material_map[name] = new G4Material( name, density=universe_mean_density,
					 nel=2 );
  m_material_map[name]->AddElement( m_element_map["Nitrogen"], 0.7 );
  m_material_map[name]->AddElement( m_element_map["Oxygen"], 0.3 );
  // Air
  name = "Air";
  m_material_map[name] = new G4Material( name, density=1.2929e-03*g/cm3,
					 nel=3, kStateGas, room_temp );
  G4double fracN  = 75.47;
  G4double fracO  = 23.20;
  G4double fracAr =  1.28;
  G4double denominator = fracN + fracO + fracAr;
  m_material_map[name]->AddElement( m_element_map["Nitrogen"],
				    massfraction=fracN/denominator );
  m_material_map[name]->AddElement( m_element_map["Oxygen"],
				    massfraction=fracO/denominator );
  m_material_map[name]->AddElement( m_element_map["Argon"],
				    massfraction=fracAr/denominator );
  // Water
  name = "Water";
  m_material_map[name] = new G4Material( name, density=1.*g/cm3, nel=2 );
  m_material_map[name]->AddElement( m_element_map["Hydrogen"], 2 );
  m_material_map[name]->AddElement( m_element_map["Oxygen"], 1 );
  // Aluminum
  name = "Aluminum";
  m_material_map[name] = new G4Material( name, Z=13., A=26.9815*g/mole,
					 density=2.699*g/cm3 );
  // Iron
  name = "Iron";
  m_material_map[name] = new G4Material( name, Z=26., A=55.85*g/mole,
					 density=7.87*g/cm3 );
  // Copper
  name = "Copper";
  m_material_map[name] = new G4Material( name, Z=29., A=63.546*g/mole,
					 density=8.96*g/cm3 );
  // Carbon
  name = "Carbon";
  m_material_map[name] = new G4Material( name, Z=6., A = 12.0107*g/mole,
					 density=2.265*g/cm3 );
  // Diamond
  name = "Diamond";
  m_material_map[name] = new G4Material( name, Z=6., A = 12.0107*g/mole,
					 density=3.34*g/cm3 );
  // LH2
  name = "LH2";
  m_material_map[name] = new G4Material( name, Z=1., A=1.008*g/mole,
					 density=70.99*mg/cm3 );
  // LD2
  name ="LD2";
  m_material_map[name] = new G4Material( name, Z=1., A=2.01410*g/mole,
					 density=166.0*mg/cm3 );
  // Ar gas
  name = "Argon";
  G4double densityAr = 1.782e-03 * g/cm3 * STP_Temperature / room_temp;
  density = densityAr;
  m_material_map[name] = new G4Material( name, Z=18., A=39.948*g/mole,
					 density, kStateGas, room_temp );
  // Ethane (C2H6)
  name = "Ethane";
  G4double densityEthane = 1.356e-3 *g/cm3 * STP_Temperature / room_temp;
  density = densityEthane;
  m_material_map[name] = new G4Material( name, density, nel=2,
					 kStateGas, room_temp );
  m_material_map[name]->AddElement( m_element_map["Carbon"], natoms=2 );
  m_material_map[name]->AddElement( m_element_map["Hydrogen"], natoms=6 );
  // Methane (CH4)
  name = "Methane";
  G4double densityMethane = 0.717e-3 *g/cm3 * STP_Temperature / room_temp;
  density = densityMethane;
  m_material_map[name] = new G4Material( name, density, nel=2,
					 kStateGas, room_temp );
  m_material_map[name]->AddElement( m_element_map["Carbon"], natoms=1 );
  m_material_map[name]->AddElement( m_element_map["Hydrogen"], natoms=4 );
  // Ar(50%) + Ethane(50%) mixture
  name = "ArEthane";
  density = 0.5*( densityAr + densityEthane );
  m_material_map[name] = new G4Material( name, density, ncomponents=2,
					 kStateGas, room_temp );
  m_material_map[name]->AddMaterial( m_material_map["Argon"],
				     massfraction=0.5*densityAr/density );
  m_material_map[name]->AddMaterial( m_material_map["Ethane"],
				     massfraction=0.5*densityEthane/density );
  // P10 gas Ar(90%) + Methane(10%) mixture
  name = "P10";
  density = 0.9*densityAr + 0.1*densityMethane;
  m_material_map[name] = new G4Material( name, density, nel=2,
					 kStateGas, room_temp );
  m_material_map[name]->AddMaterial( m_material_map["Argon"],
				     massfraction=0.9*densityAr/density );
  m_material_map[name]->AddMaterial( m_material_map["Methane"],
				     massfraction=0.1*densityMethane/density );
  // G10 epoxy glass
  name = "G10";
  m_material_map[name] = new G4Material( name, density=1.700*g/cm3,
					 ncomponents=4 );
  m_material_map[name]->AddElement( m_element_map["Silicon"], natoms=1 );
  m_material_map[name]->AddElement( m_element_map["Oxygen"] , natoms=2 );
  m_material_map[name]->AddElement( m_element_map["Carbon"] , natoms=3 );
  m_material_map[name]->AddElement( m_element_map["Hydrogen"] , natoms=3 );
  
  // Scintillator (Polystyene(C6H5CH=CH2))
  name = "Scintillator";
  m_material_map[name] = new G4Material( name, density=1.032*g/cm3, nel=2 );
  m_material_map[name]->AddElement( m_element_map["Carbon"], natoms=8 );
  m_material_map[name]->AddElement( m_element_map["Hydrogen"], natoms=8 );

 // EJ-320 added by Heejeong BYEON
  name = "EJ230";
  m_material_map[name] = new G4Material( name, density=1.032*g/cm3, nel=2 );
  m_material_map[name]->AddElement( m_element_map["Carbon"], natoms=9 );
  m_material_map[name]->AddElement( m_element_map["Hydrogen"], natoms=10 );

  
  // CH2 Polyethelene
  name = "CH2";
  m_material_map[name] = new G4Material( name, density=1.13*g/cm3, nel=2 );
  m_material_map[name]->AddElement( m_element_map["Carbon"], natoms=1 );
  m_material_map[name]->AddElement( m_element_map["Hydrogen"], natoms=2 );

  // Silica Aerogel for LAC
  name = "SilicaAerogelLAC";
  m_material_map[name] = new G4Material( name, density=0.18 *g/cm3, nel=2 );
  m_material_map[name]->AddElement( m_element_map["Silicon"], natoms=1 );
  m_material_map[name]->AddElement( m_element_map["Oxygen"],  natoms=2 );

  // Quartz (SiO2, crystalline)
  name = "Quartz";
  m_material_map[name] = new G4Material( name, density=2.64 *g/cm3, nel=2 );
  m_material_map[name]->AddElement( m_element_map["Silicon"], natoms=1 );
  m_material_map[name]->AddElement( m_element_map["Oxygen"],  natoms=2 );

  // Acrylic for WC
  name = "Acrylic";
  m_material_map[name] = new G4Material( name, density=1.18 *g/cm3, nel=3 );
  m_material_map[name]->AddElement( m_element_map["Carbon"], natoms=5 );
  m_material_map[name]->AddElement( m_element_map["Hydrogen"],  natoms=8 );
  m_material_map[name]->AddElement( m_element_map["Oxygen"],  natoms=2 );

  G4String target_material = gConf.Get<G4String>("TargetMaterial");
  G4cout << "   Target material : " << target_material << G4endl;
  if( target_material == "Carbon" ){
    m_material_map["Target"] = m_material_map["Carbon"];
  } else if( target_material == "Diamond" ){
    m_material_map["Target"] = m_material_map["Diamond"];
  } else if( target_material == "Copper" ){
    m_material_map["Target"] = m_material_map["Copper"];
  } else if( target_material == "LH2"  ){
    m_material_map["Target"] = m_material_map["LH2"];
  } else if( target_material == "LD2" ){
    m_material_map["Target"] = m_material_map["LD2"];
  } else if( target_material == "Vacuum" ){
    m_material_map["Target"] = m_material_map["Vacuum"];
  } else if( target_material == "Empty" ){
    m_material_map["Target"] = m_material_map["P10"];
  } else if( target_material == "CH2"  ){
    m_material_map["Target"] = m_material_map["CH2"];
  }
  else {
    std::string e(FUNC_NAME + " No target material : " + target_material );
    throw std::invalid_argument( e );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructAreaTent( void )
{
  //const G4ThreeVector size( 5.0*m/2, 5.0*m/2, 6.0*m/2 );
  const G4ThreeVector size( 4.8*m/2, 4.8*m/2, 6.3*m/2 );
  //const G4ThreeVector pos( 0., 0., -143.-1200.-170.*mm+size.z() );
  const G4ThreeVector pos( 0., 0., -1270.-119.*mm+size.z() );
  auto tent_out_solid = new G4Box( "AreaTentOutSolid",
				   size.x(), size.y(), size.z() );
  auto tent_in_solid = new G4Box( "AreaTentInSolid", size.x()-1.*mm,
				  size.y()-1.*mm, size.z()-1.*mm );
  auto tent_solid = new G4SubtractionSolid( "AreaTentSolid",
					    tent_out_solid, tent_in_solid );
  auto tent_lv = new G4LogicalVolume( tent_solid,
				      m_material_map["Air"],
				      "AreaTentLV" );
  tent_lv->SetVisAttributes( G4Color::Blue() );
  // tent_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );
  new G4PVPlacement( nullptr, pos, tent_lv,
		     "AreaTentPV", m_world_lv, false, 0 );
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructBC3( void )
{
  if( !m_sdc_sd ){
    m_sdc_sd = new TPCSDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector( m_sdc_sd );
  }

  const auto& bc3_pos = ( gGeom.GetGlobalPosition("KURAMA") +
			  ( gGeom.GetGlobalPosition("BC3-X1") +
			    gGeom.GetGlobalPosition("BC3-U2") ) * 0.5 );
  const auto& frame_size = gSize.GetSize( "Bc3Frame" ) * 0.5 * mm;
  const auto& drift_size = gSize.GetSize( "Bc3Drift" ) * 0.5 * mm;
  auto bc3_solid = new G4Box( "Bc3Solid", frame_size.x(),
			       frame_size.y(), frame_size.z() );
  auto bc3_lv = new G4LogicalVolume( bc3_solid, m_material_map["Argon"],
				      "Bc3LV", 0, 0, 0 );
  bc3_lv->SetVisAttributes( G4Colour::Green() );
  new G4PVPlacement( m_rotation_matrix, bc3_pos,
		     bc3_lv, "Bc3PV", m_world_lv, false, 0 );
  auto bc3pl_solid = new G4Box( "Bc3PlSolid", drift_size.x(),
				 drift_size.y(), drift_size.z() );
  G4String plane_name[] = { "Bc3X1", "Bc3X2", "Bc3V1",
			    "Bc3V2", "Bc3U1", "Bc3U2" };
  for( G4int i=0; i<NumOfLayersBcOut/2; ++i ){
    G4ThreeVector pos;
    switch (i) {
    case 0: pos.setZ( -22.*mm ); break;
    case 1: pos.setZ( -18.*mm ); break;
    case 2: pos.setZ(  -2.*mm ); break;
    case 3: pos.setZ(   2.*mm ); break;
    case 4: pos.setZ(  18.*mm ); break;
    case 5: pos.setZ(  22.*mm ); break;
    }
    auto bc3pl_lv = new G4LogicalVolume( bc3pl_solid,
					 m_material_map["Argon"],
					 plane_name[i] + "LV", 0, 0, 0 );
    bc3pl_lv->SetSensitiveDetector( m_sdc_sd );
    new G4PVPlacement( nullptr, pos, bc3pl_lv, plane_name[i] + "PV",
		       bc3_lv, false, 1001+i );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructBC4( void )
{
  if( !m_sdc_sd ){
    m_sdc_sd = new TPCSDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector( m_sdc_sd );
  }

  const auto& bc4_pos = ( gGeom.GetGlobalPosition("KURAMA") +
			  ( gGeom.GetGlobalPosition("BC4-U1") +
			    gGeom.GetGlobalPosition("BC4-X2") ) * 0.5 );
  const auto& frame_size = gSize.GetSize( "Bc4Frame" ) * 0.5 * mm;
  const auto& drift_size = gSize.GetSize( "Bc4Drift" ) * 0.5 * mm;
  auto bc4_solid = new G4Box( "Bc4Solid", frame_size.x(),
			       frame_size.y(), frame_size.z() );
  auto bc4_lv = new G4LogicalVolume( bc4_solid, m_material_map["Argon"],
				     "Bc4LV", 0, 0, 0 );
  bc4_lv->SetVisAttributes( G4Colour::Green() );
  new G4PVPlacement( m_rotation_matrix, bc4_pos,
		     bc4_lv, "Bc4PV", m_world_lv, false, 0 );
  auto bc4pl_solid = new G4Box( "Bc4PlSolid", drift_size.x(),
				 drift_size.y(), drift_size.z() );
  G4String plane_name[] = { "Bc4U1", "Bc4U2", "Bc4V1",
			    "Bc4V2", "Bc4X1", "Bc4X2" };
  for( G4int i=0; i<NumOfLayersBcOut/2; ++i ){
    G4ThreeVector pos;
    switch (i) {
    case 0: pos.setZ( -22.*mm ); break;
    case 1: pos.setZ( -18.*mm ); break;
    case 2: pos.setZ(  -2.*mm ); break;
    case 3: pos.setZ(   2.*mm ); break;
    case 4: pos.setZ(  18.*mm ); break;
    case 5: pos.setZ(  22.*mm ); break;
    }
    auto bc4pl_lv = new G4LogicalVolume( bc4pl_solid,
					 m_material_map["Argon"],
					 plane_name[i] + "LV", 0, 0, 0 );
    bc4pl_lv->SetSensitiveDetector( m_sdc_sd );
    new G4PVPlacement( nullptr, pos, bc4pl_lv, plane_name[i] + "PV",
		       bc4_lv, false, 2001+i );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructBH2( void )
{
  const auto& ra2 = gGeom.GetRotAngle2("BH2") * deg;
  const auto& half_size = gSize.GetSize("Bh2Seg") * 0.5 * mm;
  const G4double pitch = gGeom.GetWirePitch("BH2") *mm;
  auto bh2SD = new TPCBH2SD("/BH2");
  G4SDManager::GetSDMpointer()->AddNewDetector( bh2SD );
  // Mother
  auto mother_solid = new G4Box( "Bh2MotherSolid",
				 half_size.x()*NumOfSegBH2 + 50.*mm,
				 half_size.y() + 50.*mm,
				 half_size.z()*2 + 50.*mm );
  auto mother_lv = new G4LogicalVolume( mother_solid,
					m_material_map["Air"],
					"Bh2MotherLV" );
  auto rot = new G4RotationMatrix;
  rot->rotateY( - ra2 - m_rotation_angle );
  auto pos = ( gGeom.GetGlobalPosition("KURAMA") +
	       gGeom.GetGlobalPosition("BH2") );
  pos.rotateY( m_rotation_angle );
  new G4PVPlacement( rot, pos, mother_lv,
		     "Bh2MotherPV", m_world_lv, false, 0 );
  mother_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // Segment
  auto segment_solid = new G4Box( "Bh2SegmentSolid", half_size.x(),
				  half_size.y(), half_size.z() );
  auto segment_lv = new G4LogicalVolume( segment_solid,
					 m_material_map["Scintillator"],
					 "Bh2SegmentLV" );
    G4double edge=17.0 *mm;
  auto edge_segment_solid = new G4Box( "EdgeBh2SegmentSolid", edge/2.,
				       half_size.y(), half_size.z() );
  auto edge_segment_lv = new G4LogicalVolume( edge_segment_solid,
					      m_material_map["Scintillator"],
					      "EdgeBh2SegmentLV" );

  for( G4int i=0; i<NumOfSegBH2; ++i ){
    if(i==0){ //left edge
      pos = G4ThreeVector( ( -NumOfSegBH2/2 + i + 0.5 )*pitch - 1.5*mm, 0.*mm, 0.*mm );
      new G4PVPlacement( nullptr, pos, edge_segment_lv,
			 "Bh2SegmentPV", mother_lv, false, i );
    }
    else if(i==NumOfSegBH2-1){ //right edge
      pos = G4ThreeVector( ( -NumOfSegBH2/2 + i + 0.5 )*pitch + 1.5*mm, 0.*mm, 0.*mm );
      new G4PVPlacement( nullptr, pos, edge_segment_lv,
			 "Bh2SegmentPV", mother_lv, false, i );
    }
    else{ //center segments
      pos = G4ThreeVector( ( -NumOfSegBH2/2 + i + 0.5 )*pitch, 0.*mm, 0.*mm );
      new G4PVPlacement( nullptr, pos, segment_lv,
			 "Bh2SegmentPV", mother_lv, false, i );
    }
  }
  segment_lv->SetVisAttributes( G4Colour::Cyan() );
  segment_lv->SetSensitiveDetector( bh2SD );
  edge_segment_lv->SetVisAttributes( G4Colour::Cyan() );
  edge_segment_lv->SetSensitiveDetector( bh2SD );
}
//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructVC1( void )
{
  const auto& ra2 = gGeom.GetRotAngle2("VC1") * deg;
  const auto& half_size = gSize.GetSize("Vc1Seg") * 0.5 * mm;
  const G4double pitch = gGeom.GetWirePitch("VC1") * mm;

  auto vc1SD = new TPCVC1SD("/VC1");
  G4SDManager::GetSDMpointer()->AddNewDetector( vc1SD );

  // Mother Volume
  auto mother_solid = new G4Box( "Vc1MotherSolid",
                                 half_size.x()*10 + 50.*mm, // 10 segments
                                 half_size.y() + 50.*mm,
                                 half_size.z()*2 + 50.*mm );
  auto mother_lv = new G4LogicalVolume( mother_solid,
                                        m_material_map["Air"],
                                        "Vc1MotherLV" );

  auto rot = new G4RotationMatrix;
  rot->rotateY( - ra2 - m_rotation_angle );

  auto pos = ( gGeom.GetGlobalPosition("KURAMA") +
               gGeom.GetGlobalPosition("VC1") );
  pos.rotateY( m_rotation_angle );

  new G4PVPlacement( rot, pos, mother_lv,
                     "Vc1MotherPV", m_world_lv, false, 0 );

  mother_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );

  // Segment
  auto segment_solid = new G4Box( "Vc1SegmentSolid",
                                  half_size.x(),
                                  half_size.y(),
                                  half_size.z() );
  auto segment_lv = new G4LogicalVolume( segment_solid,
                                         m_material_map["EJ230"],
                                         "Vc1SegmentLV" );

  // set same segment 
  for( G4int i=0; i<10; ++i ){
    pos = G4ThreeVector( ( -10/2 + i + 0.5 ) * pitch, 0.*mm, 0.*mm );
    new G4PVPlacement( nullptr, pos, segment_lv,
                       "Vc1SegmentPV", mother_lv, false, i );
  }

  segment_lv->SetVisAttributes( G4Colour::Cyan() );
  segment_lv->SetSensitiveDetector( vc1SD );
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructVC2( void )
{
  const auto& ra2 = gGeom.GetRotAngle2("VC2") * deg;
  const auto& half_size = gSize.GetSize("Vc2Seg") * 0.5 * mm;
  const G4double pitch = gGeom.GetWirePitch("VC2") * mm;

  auto vc2SD = new TPCVC2SD("/VC2");
  G4SDManager::GetSDMpointer()->AddNewDetector( vc2SD );

  // Mother Volume
  auto mother_solid = new G4Box( "Vc2MotherSolid",
                                 half_size.x()*10 + 50.*mm, // 10 segments
                                 half_size.y() + 50.*mm,
                                 half_size.z()*2 + 50.*mm );
  auto mother_lv = new G4LogicalVolume( mother_solid,
                                        m_material_map["Air"],
                                        "Vc2MotherLV" );

  auto rot = new G4RotationMatrix;
  rot->rotateY( - ra2 - m_rotation_angle );

  auto pos = ( gGeom.GetGlobalPosition("KURAMA") +
               gGeom.GetGlobalPosition("VC2") );
  pos.rotateY( m_rotation_angle );

  new G4PVPlacement( rot, pos, mother_lv,
                     "Vc2MotherPV", m_world_lv, false, 0 );

  mother_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );

  // Segment
  auto segment_solid = new G4Box( "Vc2SegmentSolid",
                                  half_size.x(),
                                  half_size.y(),
                                  half_size.z() );
  auto segment_lv = new G4LogicalVolume( segment_solid,
                                         m_material_map["EJ230"],
                                         "Vc2SegmentLV" );

  // set same segment 
  for( G4int i=0; i<10; ++i ){
    pos = G4ThreeVector( ( -10/2 + i + 0.5 ) * pitch, 0.*mm, 0.*mm );
    new G4PVPlacement( nullptr, pos, segment_lv,
                       "Vc2SegmentPV", mother_lv, false, i );
  }

  segment_lv->SetVisAttributes( G4Colour::Cyan() );
  segment_lv->SetSensitiveDetector( vc2SD );
}


//_____________________________________________________________________________
	void
TPCDetectorConstruction::ConstructBVH( void )
{
  const auto& ra2 = gGeom.GetRotAngle2("BVH") * deg;
  const auto& half_size = gSize.GetSize("BVHSeg") * 0.5 * mm;
  const G4double pitch = gGeom.GetWirePitch("BVH") *mm;
  auto bvhSD = new TPCBVHSD("/BVH");
  G4SDManager::GetSDMpointer()->AddNewDetector( bvhSD );
  // Mother
  auto mother_solid = new G4Box( "BVHMotherSolid",
				 half_size.x()*NumOfSegBVH + 50.*mm,
				 half_size.y() + 50.*mm,
				 half_size.z()*2 + 50.*mm );
  auto mother_lv = new G4LogicalVolume( mother_solid,
					m_material_map["Air"],
					"BVHMotherLV" );
  auto rot = new G4RotationMatrix;
  rot->rotateY( - ra2 - m_rotation_angle );
  auto pos = ( gGeom.GetGlobalPosition("KURAMA") +
	       gGeom.GetGlobalPosition("BVH") );
  pos.rotateY( m_rotation_angle );
  new G4PVPlacement( rot, pos, mother_lv,
		     "BVHMotherPV", m_world_lv, false, 0 );
  mother_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // Segment
  auto segment_solid = new G4Box( "BVHSegmentSolid", half_size.x(),
				  half_size.y(), half_size.z() );
  auto segment_lv = new G4LogicalVolume( segment_solid,
					 m_material_map["Scintillator"],
					 "BVHSegmentLV" );

  for( G4int i=0; i<NumOfSegBVH; ++i ){
		pos = G4ThreeVector( ( -NumOfSegBVH/2 + i + 0.5 )*pitch, 0.*mm, 0.*mm );
		new G4PVPlacement( nullptr, pos, segment_lv,
				"BVHSegmentPV", mother_lv, false, i );
  }
  segment_lv->SetVisAttributes( G4Colour::Cyan() );
  segment_lv->SetSensitiveDetector( bvhSD );
}

//_____________________________________________________________________________

void
TPCDetectorConstruction::ConstructFTOF( void )
{
  const auto& ra2 = gGeom.GetRotAngle2("TOF") * deg;
  const auto& half_size = gSize.GetSize("FtofSeg") * 0.5 * mm;
  const G4double pitch = gGeom.GetWirePitch("TOF") * mm;
  auto ftofSD = new TPCFTOFSD("/FTOF");
  G4SDManager::GetSDMpointer()->AddNewDetector( ftofSD );
  // Mother
  auto mother_solid = new G4Box( "FtofMotherSolid",
				 half_size.x()*NumOfSegTOF + 50.*mm,
				 half_size.y() + 50.*mm,
				 half_size.z()*2 + 50.*mm );
  auto mother_lv = new G4LogicalVolume( mother_solid,
					m_material_map["Air"],
					"FtofMotherLV" );
  auto rot = new G4RotationMatrix;
  rot->rotateY( - ra2 - m_rotation_angle );
  auto pos = ( gGeom.GetGlobalPosition("KURAMA") +
	       gGeom.GetGlobalPosition("TOF") );
  pos.rotateY( m_rotation_angle );
  new G4PVPlacement( rot, pos, mother_lv,
		     "FtofMotherPV", m_world_lv, false, 0 );
  mother_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // Segment
  auto segment_solid = new G4Box( "FtofSegmentSolid", half_size.x(),
				  half_size.y(), half_size.z() );
  auto segment_lv = new G4LogicalVolume( segment_solid,
					 m_material_map["Scintillator"],
					 "FtofSegmentLV" );
  for( G4int i=0; i<NumOfSegTOF; ++i ){
    segment_lv->SetVisAttributes( G4Colour::Cyan() );
    segment_lv->SetSensitiveDetector( ftofSD );
    pos = G4ThreeVector( ( -NumOfSegTOF/2 + i )*pitch,
			 0.0,
			 2.*( - i%2 + 0.5 )*half_size.z() );
    new G4PVPlacement( nullptr, pos, segment_lv,
		       "FtofSegmentPV", mother_lv, false, i );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructHTOF( void )
{
  auto htof_sd = new TPCHTOFSD("/HTOF");
  G4SDManager::GetSDMpointer()->AddNewDetector( htof_sd );
  const auto& htof_pos = gGeom.GetGlobalPosition( "HTOF" );
  const auto& half_size = gSize.GetSize( "HtofSeg" ) * 0.5 * mm;
  const G4double L = gGeom.GetLocalZ( "HTOF" );
  const G4double dXdW = gGeom.GetWirePitch( "HTOF" );
  // Segment
  // Sintillator
  auto htof_scintilltor = new G4Box( "HtofScint", half_size.x(),
				     half_size.y(), half_size.z() );
  // Light-guides
  //const G4ThreeVector lg_size(35.0,100.0,4.0);
  const G4ThreeVector lg_size(17.5,50.0,2.0);
  // Upper one
  std::vector<G4TwoVector> upper_lg_vertices;
  upper_lg_vertices.push_back( G4TwoVector( -half_size.x() , -half_size.z() ) );
  upper_lg_vertices.push_back( G4TwoVector( -half_size.x() , half_size.z() ) );
  upper_lg_vertices.push_back( G4TwoVector( half_size.x() , half_size.z() ) );
  upper_lg_vertices.push_back( G4TwoVector( half_size.x() , -half_size.z() ) );
  upper_lg_vertices.push_back( G4TwoVector( -lg_size.x() , -half_size.z() ) );
  upper_lg_vertices.push_back( G4TwoVector( -lg_size.x() , -half_size.z() + lg_size.z() ) );
  upper_lg_vertices.push_back( G4TwoVector( lg_size.x() , -half_size.z() + lg_size.z() ) );
  upper_lg_vertices.push_back( G4TwoVector( lg_size.x() , -half_size.z() ) );
  auto htof_upper_lg = new G4GenericTrap("HtofLG_upper", lg_size.y(), upper_lg_vertices );
  auto rotM_upper_lg = new G4RotationMatrix;
  rotM_upper_lg->rotateX( 90.0 *deg );
  rotM_upper_lg->rotateZ( - 180.0 *deg );
  G4ThreeVector trans_upper_lg(0.*mm , half_size.y() + lg_size.y(),  0.*mm );

  // Lower one
  std::vector<G4TwoVector> lower_lg_vertices;
  lower_lg_vertices.push_back( G4TwoVector( -half_size.x() , -half_size.z() ) );
  lower_lg_vertices.push_back( G4TwoVector( -half_size.x() , half_size.z() ) );
  lower_lg_vertices.push_back( G4TwoVector( half_size.x() , half_size.z() ) );
  lower_lg_vertices.push_back( G4TwoVector( half_size.x() , -half_size.z() ) );
  lower_lg_vertices.push_back( G4TwoVector( -lg_size.x() , half_size.z() - lg_size.z() ) );
  lower_lg_vertices.push_back( G4TwoVector( -lg_size.x() , half_size.z() ) );
  lower_lg_vertices.push_back( G4TwoVector( lg_size.x() , half_size.z() ) );
  lower_lg_vertices.push_back( G4TwoVector( lg_size.x() , half_size.z() - lg_size.z() ) );
  auto htof_lower_lg = new G4GenericTrap("HtofLG_lower", lg_size.y(), lower_lg_vertices );
  auto rotM_lower_lg = new G4RotationMatrix;
  rotM_lower_lg->rotateX( - 90.0 *deg );
  rotM_lower_lg->rotateZ( - 180.0 *deg );
  G4ThreeVector trans_lower_lg(0.*mm , -half_size.y() - lg_size.y() , 0.*mm );
  auto solid_dummy = new G4UnionSolid( "dummysolid", htof_scintilltor, htof_upper_lg,
				       rotM_upper_lg, trans_upper_lg );

  auto htof_solid = new G4UnionSolid( "HtofSolid", solid_dummy,
				      htof_lower_lg, rotM_lower_lg,
				      trans_lower_lg );
  //Common slats
  auto htof_lv = new G4LogicalVolume( htof_solid, m_material_map["Scintillator"],
				      "HtofLV" );

  //HTOF beam-through part
  G4double HTOF_window=110.*mm;
  auto window_dummy = new G4Box( "window_dummy", half_size.x(), half_size.y()/2. - HTOF_window/4., half_size.z() );

  //Upper slats(Beam-through)
  G4ThreeVector trans_upper_window(0.*mm , half_size.y()/2. - HTOF_window/4. + lg_size.y(),  0.*mm );
  auto htof_solid_upper = new G4UnionSolid( "HtofSolid_upper", window_dummy, htof_upper_lg,
					    rotM_upper_lg, trans_upper_window );
  auto htof_upper_lv = new G4LogicalVolume( htof_solid_upper, m_material_map["Scintillator"],
					    "HtofUpperLV" );
  //Lower slats(Beam-through)
  G4ThreeVector trans_lower_window(0.*mm , -half_size.y()/2. + HTOF_window/4. - lg_size.y(),  0.*mm );
  auto htof_solid_lower = new G4UnionSolid( "HtofSolid_lower", window_dummy, htof_lower_lg,
					    rotM_lower_lg, trans_lower_window );
  auto htof_lower_lv = new G4LogicalVolume( htof_solid_lower, m_material_map["Scintillator"],
					    "HtofLowerLV" );

  for( G4int i=0; i<NumOfPlaneHTOF; ++i ){
    for( G4int j=0; j<NumOfSegHTOFOnePlane; ++j ){
      G4int seg = i*NumOfSegHTOFOnePlane + j;
      // ( lateral, height, radial )
      G4ThreeVector seg_pos( -dXdW * ( j - ( NumOfSegHTOFOnePlane - 1 )/2. ),
			     0.*mm,
			     -L );
      auto rotMOutP = new G4RotationMatrix;
      rotMOutP->rotateY( - i * 360./NumOfPlaneHTOF*deg );
      seg_pos.rotateY( i * 360./NumOfPlaneHTOF*deg );
      seg_pos += htof_pos;
      G4int copy_no = seg+2;

      G4ThreeVector window_pos( 0.*mm, half_size.y()/2. + HTOF_window/4., 0.*mm);
      //common slats
      if(i!=0 )	new G4PVPlacement( rotMOutP, seg_pos, htof_lv, Form("HtofPV%d", copy_no), m_world_lv, false, copy_no );
      else if(j==0) new G4PVPlacement( rotMOutP, seg_pos, htof_lv, Form("HtofPV%d", 0), m_world_lv, false, 0 );
      else if(j==3) new G4PVPlacement( rotMOutP, seg_pos, htof_lv, Form("HtofPV%d", 5), m_world_lv, false, 5 );
      //Beam-through slats
      else if(j==1){
	new G4PVPlacement( rotMOutP, seg_pos + window_pos, htof_upper_lv, Form("HtofPV%d", 1), m_world_lv, false, 1 );
	new G4PVPlacement( rotMOutP, seg_pos - window_pos, htof_lower_lv, Form("HtofPV%d", 2), m_world_lv, false, 2 );
      }
      else if(j==2){
	new G4PVPlacement( rotMOutP, seg_pos + window_pos, htof_upper_lv, Form("HtofPV%d", seg), m_world_lv, false, 3 );
	new G4PVPlacement( rotMOutP, seg_pos - window_pos, htof_lower_lv, Form("HtofPV%d", seg+31), m_world_lv, false, 4 );
      }
    }
  }

  htof_lv->SetSensitiveDetector( htof_sd );
  htof_upper_lv->SetSensitiveDetector( htof_sd );
  htof_lower_lv->SetSensitiveDetector( htof_sd );
  htof_lv->SetVisAttributes( G4Colour::Cyan() );
  htof_upper_lv->SetVisAttributes( G4Colour::Green() );
  htof_lower_lv->SetVisAttributes( G4Colour::Green() );

  // Supporting frame parts
  // Dummy for subtraction
  auto RingHoleOut = new G4Box( "RingHoleOut", 127.5*mm, 6.0*mm, 5.*mm );
  auto RingHoleIn = new G4Box( "RingHoleIn", 10.*mm, 6.*mm, 5.*mm );
  auto RingHole = new G4SubtractionSolid( "RingHole" , RingHoleOut, RingHoleIn );

  // Top Ring
  const G4double RingIn_zPlane[2] = { -5.*mm, 5.*mm };
  const G4double RingIn_rInner[2] = { 0.*mm, 0.*mm };
  const G4double RingIn_rOuter[2] =  { 334.13*mm, 334.13*mm } ;
  auto TopRingOutSolid = new G4Tubs( "TopRingOutSolid", 0.*mm, 391.97*mm,
				     5.*mm, 0.*deg, 360.*deg );
  auto TopRingInSolid = new G4Polyhedra( "TopRingInSolid", 22.5*deg, ( 360. + 22.5 )*deg,
					 8, 2,
					 RingIn_zPlane, RingIn_rInner, RingIn_rOuter );

  G4SubtractionSolid* TopRingSubSolid[8];
  for( G4int i=0; i<8; i++){
    auto rotMHole = new G4RotationMatrix;
    rotMHole->rotateZ( i * 45.0 * deg );
    G4ThreeVector trans_hole(0*mm, 348.13*mm, 0.*mm );
    trans_hole.rotateZ( - i * 45.0 * deg );
    if( i==0 ) TopRingSubSolid[0] = new G4SubtractionSolid( "TopRingSubSolid_0", TopRingOutSolid, RingHole, rotMHole, trans_hole );
    else TopRingSubSolid[i] = new G4SubtractionSolid( Form("TopRingSubSolid_%d", i ), TopRingSubSolid[i-1], RingHole, rotMHole, trans_hole );
  }

  auto TopRingSolid = new G4SubtractionSolid( "TopRingSoild" , TopRingSubSolid[7], TopRingInSolid );
  auto TopRing_lv = new G4LogicalVolume( TopRingSolid, m_material_map["Iron"], "TopRingLV" );
  TopRing_lv->SetVisAttributes( ORANGE );

  // Bottom Ring
  const G4double RingOut_zPlane[2] = { -5.*mm, 5.*mm };
  const G4double RingOut_rInner[2] = { 334.13*mm, 334.13*mm };
  const G4double RingOut_rOuter[2] = { 362.13*mm, 362.13*mm };
  auto BotRingOutSolid = new G4Polyhedra( "BotRingOutSolid", 22.5*deg, ( 360. + 22.5 )*deg,
					  8, 2,
					  RingOut_zPlane, RingOut_rInner, RingOut_rOuter );

  G4SubtractionSolid* BotRingSubSolid[8];
  for( G4int i=0; i<8; i++){
    auto rotMHole = new G4RotationMatrix;
    rotMHole->rotateZ( i * 45.0 * deg );
    G4ThreeVector trans_hole(0.*mm, 348.13*mm, 0.*mm );
    trans_hole.rotateZ( - i * 45.0 * deg );
    if( i==0 ) BotRingSubSolid[i] = new G4SubtractionSolid( Form("BotRingSubSolid_%d", i ) , BotRingOutSolid, RingHole, rotMHole, trans_hole );
    else BotRingSubSolid[i] = new G4SubtractionSolid( Form("BotRingSubSolid_%d", i ) , BotRingSubSolid[i-1], RingHole, rotMHole, trans_hole );
  }
  auto BotRing_lv = new G4LogicalVolume( BotRingSubSolid[7], m_material_map["Iron"],"BotRingLV");
  BotRing_lv->SetVisAttributes( G4Colour::Blue() );

  // Bracket
  auto BraInSolid = new G4Box( "BraInSolid", 137.5*mm, 75.*mm, 2.5*mm );
  auto BraInHole = new G4Box( "BraInHole", 127.5*mm, 30.*mm, 2.5*mm );
  G4ThreeVector brain_trans_hole(0.*mm, 25.*mm, 0.*mm );
  auto BraInSubSolid = new G4SubtractionSolid( "BraInSubSolid", BraInSolid, BraInHole,
					       0, brain_trans_hole );

  auto BraOutSolid = new G4Box( "BraOutSolid", 137.5*mm, 2.5*mm, 16.5*mm );
  G4ThreeVector braout_trans_hole(0.*mm, 0.*mm, 3.5*mm );
  auto BraOutSubSolid = new G4SubtractionSolid( "BraOutSubSolid", BraOutSolid, RingHole,
						0, braout_trans_hole );

  G4ThreeVector trans_braket(0.*mm, 72.5*mm, 14.*mm );
  auto BraSolid = new G4UnionSolid( "BraSolid", BraInSubSolid, BraOutSubSolid,
				    0, trans_braket );
  auto Bra_lv = new G4LogicalVolume( BraSolid, m_material_map["Iron"],
				     "BraLV" );
  Bra_lv->SetVisAttributes( G4Colour::Green() );

  // Preamp Support Frame
  const G4double PreFrameIn_zPlane[2] = { -37.*mm, 37.*mm };
  const G4double PreFrameIn_rInner[2] = { 334.13*mm, 334.13*mm };
  const G4double PreFrameIn_rOuter[2] = { 339.13*mm, 339.13*mm } ;
  auto PreFrameInSolid = new G4Polyhedra( "PreFrameInSolid", 22.5*deg, ( 360. + 22.5 )*deg,
					  8, 2,
					  PreFrameIn_zPlane,  PreFrameIn_rInner,  PreFrameIn_rOuter);

  const G4double PreFrameOut_zPlane[2] = { -32.5*mm, 32.5*mm };
  const G4double PreFrameOut_rInner[2] = { 354.13*mm, 354.13*mm };
  const G4double PreFrameOut_rOuter[2] = { 362.13*mm, 362.13*mm } ;
  auto PreFrameOutSolid = new G4Polyhedra( "PreFrameOutSolid", 22.5*deg, ( 360. + 22.5 )*deg,
					   8, 2,
					   PreFrameOut_zPlane,  PreFrameOut_rInner,  PreFrameOut_rOuter);
  auto PreFrameHole = new G4Box( "PreFrameHole", 3.5*mm, 4.0*mm, 32.5*mm );

  G4SubtractionSolid* PreFrameOutSubSolid[8];
  for( G4int i=0; i<8; i++){
    auto rotMHole = new G4RotationMatrix;
    rotMHole->rotateZ( i * 45.0 * deg );
    G4ThreeVector trans_hole(0.*mm, 348.13*mm, 0.*mm );
    trans_hole.rotateZ( - i * 45.0 * deg );
    if( i==0 ) PreFrameOutSubSolid[0] = new G4SubtractionSolid( "PreFrameOutSubSolid_0",
								PreFrameOutSolid, PreFrameHole );
    else PreFrameOutSubSolid[i] = new G4SubtractionSolid( Form( "PreFrameOutSubSolid_%d", i ),
							  PreFrameOutSubSolid[i-1], PreFrameHole,
							  rotMHole, trans_hole );
  }
  G4ThreeVector trans_preframe( 0.*mm, 0.*mm, -9.0*mm );
  auto PreFrameSolid = new G4UnionSolid( "PreFrameSolid", PreFrameInSolid, PreFrameOutSubSolid[7],
					 0, trans_preframe );
  auto PreFrame_lv = new G4LogicalVolume( PreFrameSolid, m_material_map["Aluminum"],
					 "PreFrameLV" );
  PreFrame_lv->SetVisAttributes( G4Colour::Red() );

  // Bar
  auto BarMainSolid = new G4Box( "BarMainSolid", 45.*mm, 618.25*mm, 5.*mm );
  auto BarSideSolid = new G4Box( "BarSideSolid", 5.*mm, 618.25*mm, 5.*mm );
  G4ThreeVector trans_bar_side( -40.*mm, 0.*mm, -10.*mm );
  auto BarUniSolid = new G4UnionSolid( "BarUniSolid", BarMainSolid, BarSideSolid,
				       0, trans_bar_side );

  auto BarBotSolid = new G4Box( "BarSideSolid", 45.*mm, 5.*mm, 5.*mm );
  G4ThreeVector trans_bar_bot( 0.*mm, -613.25*mm, -10.*mm );

  auto BarSolid = new G4UnionSolid( "BarSolid", BarUniSolid, BarBotSolid,
				    0, trans_bar_bot );
  auto Bar_lv = new G4LogicalVolume( BarSolid, m_material_map["Iron"],
				     "BarLV" );
  Bar_lv->SetVisAttributes( G4Colour::Blue() );

  //if(true){   // Htof Frame Placement
  if(false){   // Htof Frame Placement

    auto rotMOutRing = new G4RotationMatrix;
    rotMOutRing->rotateX( - 90. *deg );

    // Top Ring
    G4ThreeVector TopRing_pos( 0.*mm, 586.72*mm, 0.*mm );
    TopRing_pos += htof_pos;
    new G4PVPlacement( rotMOutRing, TopRing_pos, TopRing_lv, "TopRingPV", m_world_lv, false, 0 );

    // Bottom Ring
    G4ThreeVector BotRing_pos( 0.*mm, -586.72*mm, 0.*mm );
    BotRing_pos += htof_pos;
    new G4PVPlacement( rotMOutRing, BotRing_pos, BotRing_lv, "BotRingPV",
		       m_world_lv, false, 0 );

    // Preamp Support Frame
    for( G4int i=0; i<2; ++i ){
      auto rotMOutP_PreFrame = new G4RotationMatrix;
      rotMOutP_PreFrame->rotateX( ( 1 - 2 * i ) * 90.*deg );
      rotMOutP_PreFrame->rotateZ( i * 180.*deg );

      G4ThreeVector PreFrame_pos( 0.*mm, 453.22*mm, 0.*mm );
      PreFrame_pos.rotateZ( - i * 180.*deg );
      PreFrame_pos += htof_pos;
      new G4PVPlacement( rotMOutP_PreFrame, PreFrame_pos, PreFrame_lv,
			 Form("PreFramePV%d", i),
			 m_world_lv, false, i );
    }

    G4int Bar_seg=0;
    for( G4int i=0; i<NumOfPlaneHTOF; ++i ){
      //Bar
      auto rotMOutP_bar = new G4RotationMatrix;
      rotMOutP_bar->rotateY( - i * 45.*deg );
      G4ThreeVector Bar_pos( 0.*mm, - 36.53*mm, - 372.13*mm );
      Bar_pos.rotateY( i * 45.*deg );
      Bar_pos += htof_pos;
      if(i!=0 && i!=4){ //Beam through
	new G4PVPlacement( rotMOutP_bar, Bar_pos, Bar_lv,
			   Form("BarPV%d", Bar_seg),
			   m_world_lv, false, Bar_seg );
	Bar_seg++;
      }
      //Bracket
      for( G4int j=0; j<2; ++j ){
	auto rotMOutP_bra = new G4RotationMatrix;
	rotMOutP_bra->rotateY(  - ( 1 - 2 * j ) * i * 45.*deg );
	rotMOutP_bra->rotateZ(  - j * 180.*deg );
	G4ThreeVector Bra_pos( 0.*mm, 506.72*mm, - 364.63*mm );
	Bra_pos.rotateY( i * 45.*deg );
	Bra_pos.rotateZ( j * 180.*deg );
	Bra_pos += htof_pos;
	G4int Bra_seg = 2 * i + j;
	new G4PVPlacement( rotMOutP_bra, Bra_pos, Bra_lv,
			   Form("BraPV%d", Bra_seg),
			   m_world_lv, false, Bra_seg );
      }
    }
  }

}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructHypTPC( void )
{

  const auto& tpc_pos = gGeom.GetGlobalPosition("HypTPC");
  
  bool BeamData =  (DiscardData and ( abs(Generator) == 135 or abs(Generator) == 493 or abs(Generator) == 938 ));
  // Target
  auto target_sd = new TPCTargetSD( "/TGT" );
  G4SDManager::GetSDMpointer()->AddNewDetector( target_sd );
  auto tpc_sd = new TPCPadSD("/TPC");
	G4SDManager::GetSDMpointer()->AddNewDetector( tpc_sd );
	//
  auto vp_sd = new TPCTargetVPSD("/TGTVP");
  G4SDManager::GetSDMpointer()->AddNewDetector( vp_sd );
  const auto target_pos = gGeom.GetGlobalPosition( "SHSTarget" ) * mm;
  const auto target_size = gSize.GetSize( "Target" ) * 0.5 * mm;
  const auto holder_size = gSize.GetSize( "TargetHolder" ) * mm;
  G4ThreeVector holder_pos;
  G4VSolid* target_solid = nullptr;
  G4VSolid* holder_solid = nullptr;
	G4double maxStep = 0.3*mm;//~TPCres
	auto StepLimit = new G4UserLimits(maxStep);
  switch( m_experiment ){
  case 42: {
    target_solid = new G4Box( "Target", target_size.x(),
			      target_size.y(), target_size.z() );
    holder_solid = new G4Tubs( "TargetHolderSolid", holder_size.x() /* Rin */,
			       holder_size.y() /* Rout */,
			       holder_size.z()/2 /* DZ */, 0., 360.*deg );
    holder_pos = target_pos;
  }
    break;
  case 45: case 27: {
    //  G4Box* TargetSolid = new G4Box("target", 1.5*cm,0.25*cm,0.5*cm);
    G4double Target_r = gSize.Get( "Target", ThreeVector::X );
    // G4double Target_y = gSize.Get( "Target", ThreeVector::Y );
    G4double Target_z = gSize.Get( "Target", ThreeVector::Z )*0.5;
    target_solid = new G4Tubs( "TargetSolid", 0.*mm,
			       Target_r*mm,Target_z*mm, 0., 360*deg );
    //  G4Box* TargetSolid = new G4Box("target", 1.*cm,0.25*cm,1.*cm);
    holder_solid = new G4Tubs( "TargetHolderSolid", (Target_r + 15.)*mm,
			       (Target_r + 15.2)*mm, 200.*mm, 0., 360*deg );
    // holder_pos.set( 0*mm, 100.*mm, target_pos.z() );
  }
    break;
  default:
    return;
  }
  auto target_lv = new G4LogicalVolume( target_solid, m_material_map["Target"],
					"TargetLV" );
	target_lv->SetUserLimits(StepLimit);
  target_lv->SetSensitiveDetector( target_sd );
  target_lv->SetVisAttributes( G4Colour::Red() );
  new G4PVPlacement( nullptr, target_pos, target_lv, "TargetPV",
  		     m_world_lv, true, 0 );
 	double vp_th = 0.001*mm;
	auto vp_solid = new G4Box( "TGTVPSolid", 30.*cm/2, 20.*cm/2, vp_th/2 );
	auto vp_lv = new G4LogicalVolume( vp_solid, m_material_map["P10"], "VPLV" );
	vp_lv->SetSensitiveDetector( vp_sd );
	vp_lv->SetVisAttributes(G4Color::Blue());
	auto holder_lv = new G4LogicalVolume( holder_solid, m_material_map["P10"],
					"TargetHolderLV");
  holder_lv->SetVisAttributes( G4Colour::Blue() );
  // new G4PVPlacement( nullptr, holder_pos, holder_lv, "TargetHolderPV",
  // 		     m_world_lv, true, 0 );
  // HypTPC
  if(BeamData){
    auto vp_lv_dummy = new G4LogicalVolume( vp_solid, m_material_map["P10"], "TpcPadLV0" ); 
    vp_lv_dummy->SetSensitiveDetector( tpc_sd );
      new G4PVPlacement(nullptr,target_pos-G4ThreeVector(0,0,target_size.z()+3*vp_th/2),vp_lv_dummy,"TpcPadPV0",
          m_world_lv,0,0);
  }
		new G4PVPlacement(nullptr,target_pos-G4ThreeVector(0,0,target_size.z()+vp_th/2),vp_lv,"VPTargetBefore",
				m_world_lv,0,0);
		auto after_tgt = target_pos + G4ThreeVector(0,0,target_size.z()+2*mm) ; 
		new G4PVPlacement(nullptr,target_pos+G4ThreeVector(0,0,target_size.z()+vp_th/2),vp_lv,"VPTargetAfter",
				m_world_lv,0,1);
  if(BeamData)return;
  {
    const G4double Rin  = gSize.Get( "TpcRin" )*mm/2;
    const G4double Rout = gSize.Get( "TpcRout" )*mm/2;
    const G4double Dz   = gSize.Get( "TpcDz" )*mm;
    const G4int NumOfSide   = 8;
    const G4int NumOfZPlane = 2;
    const G4double zPlane[NumOfZPlane] = { -Dz, Dz };
    const G4double rInner[NumOfZPlane] = { 0.5*std::sqrt(3.)*Rin,
					   0.5*std::sqrt(3.)*Rin };
    const G4double rOuter[NumOfZPlane] = { Rout, Rout };
    auto tpc_out_solid = new G4Polyhedra( "TpcOutSolid",
					  22.5*deg, ( 360. + 22.5 )*deg,
					  NumOfSide, NumOfZPlane,
					  zPlane, rInner, rOuter );
    auto target_w_vp1 = new G4UnionSolid("Target_w_VP1",target_solid,vp_solid,nullptr,G4ThreeVector(0,0,-target_size.z()-vp_th/2));
    auto target_w_vp2 = new G4UnionSolid("Target_w_VP1",target_w_vp1,vp_solid,nullptr,G4ThreeVector(0,0,target_size.z()+vp_th/2));
    G4VSolid* Empty_volume = nullptr; 
    if(m_experiment == 42){
      G4ThreeVector Empty_Volume_Dim = G4ThreeVector( target_size.x()+2*mm, target_size.y()+2*mm, target_size.z()+2*mm );  
		  Empty_volume = new G4Box("hollow",Empty_Volume_Dim.x()/2,Empty_Volume_Dim.y()/2,Empty_Volume_Dim.z()/2);
    }
    else if (m_experiment == 27 or m_experiment == 45){
      G4double Target_r = gSize.Get( "Target", ThreeVector::X );
      Empty_volume = new G4Tubs("hollow",0,(Target_r+15.3)*mm,200.1*mm,0,360*deg);
    }
//    if(BeamData)Empty_Volume_Dim *= 2;
		
		auto pos = target_pos;
    pos.rotateX( 90.*deg );
		G4VSolid* tpc_solid = nullptr;
	/*
  	if(BeamData){
			tpc_solid = new G4SubtractionSolid( "TpcSolid",
					     tpc_out_solid, target_solid,
					     nullptr, pos );
				}
		else{
			tpc_solid = new G4SubtractionSolid( "TpcSolid",
					     tpc_out_solid, Empty_volume,
					     nullptr, pos );
		}
    */
			tpc_solid = new G4SubtractionSolid( "TpcSolid",
					     tpc_out_solid, Empty_volume,
					     nullptr, pos );
    auto rot = new G4RotationMatrix;
    rot->rotateX( 90.*deg );
    m_tpc_lv = new G4LogicalVolume( tpc_solid, m_material_map["P10"],
				    "TpcLV" );
		m_tpc_lv->SetUserLimits(StepLimit);
    new G4PVPlacement( rot, tpc_pos, m_tpc_lv, "TpcPV",
		       m_world_lv, false, 0 );
    m_tpc_lv->SetVisAttributes( G4Colour::White() );
		auto before_tgt = target_pos - G4ThreeVector(0,0,target_size.z()+2*mm); 
//		new G4PVPlacement(nullptr,before_tgt,vp_lv,"VPTargetBefore",
    //m_tpc_lv->SetSensitiveDetector( tpc_sd );
  }


  // Field Cage
  {
    const G4double Rin  = gSize.Get( "TpcRinFieldCage" )*mm*0.5;
    const G4double Rout = gSize.Get( "TpcRoutFieldCage" )*mm*0.5;
    const G4double Dz   = gSize.Get( "TpcDzFieldCage" )*mm;
    const G4int NumOfSide   = 8;
    const G4int NumOfZPlane = 2;
    const G4double zPlane[NumOfZPlane] = { -Dz, Dz };
    const G4double rInner[NumOfZPlane] = { Rin, Rin };
    const G4double rOuter[NumOfZPlane] = { Rout, Rout };
    auto fc_solid = new G4Polyhedra( "FcSolid", 22.5*deg, ( 360. + 22.5 )*deg,
				     NumOfSide, NumOfZPlane,
				     zPlane, rInner, rOuter );
    auto fc_lv = new G4LogicalVolume( fc_solid, m_material_map["P10"],
				      "FcLV" );
    new G4PVPlacement( nullptr, G4ThreeVector(), fc_lv,
		       "FieldCagePV", m_tpc_lv, false, 0 );
    fc_lv->SetVisAttributes( G4Colour::Red() );
  }
  // Virtual pads
  G4Tubs* pad_solid[NumOfPadTPC];
  G4LogicalVolume* pad_lv[NumOfPadTPC];
  G4double angle[NumOfPadTPC] = {};
  const G4double pad_center_z = gSize.Get("TpcPadCenterZ")*mm;
  G4double pad_in[NumOfPadTPC] = {};
  G4double pad_out[NumOfPadTPC] = {};
  G4double tpc_rad = 250;
  // out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  const G4double pad_length_in = gSize.Get("TpcPadLengthIn");
  const G4double pad_length_out = gSize.Get("TpcPadLengthOut");
  const G4double pad_gap = gSize.Get("TpcPadGap");
  const G4int pad_configure = gSize.Get("TpcPadConfigure");
  switch( pad_configure ){
  case 1:
    for( G4int i=0; i<NumOfPadTPC; ++i ){
      if( i<NumOfPadTPCIn ){
	pad_in[i]  = 10.+(pad_length_in+pad_gap)*i;
	pad_out[i] = 10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i]   = 360.;
      } else {
	pad_in[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn);
	pad_out[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn) + pad_length_out;
	angle[i] =
	  -std::acos( ( math::Sq( pad_out[i] ) +
			math::Sq( pad_center_z ) -
			math::Sq( tpc_rad ) ) /
		      ( 2*pad_out[i]*pad_center_z ) )*math::Rad2Deg() + 180.;
      }
    }
    break;
  case 2:
    for( G4int i=0; i<NumOfPadTPC; ++i ){
      if( i<NumOfPadTPCIn ){
	pad_in[i]  = 10.+(pad_length_in+pad_gap)*i;
	pad_out[i] = 10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i]   = 360.;
      }else {
	pad_in[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn);
	pad_out[i] = 10.+(pad_length_in+pad_gap)*NumOfPadTPCIn +
	  (pad_length_out+pad_gap)*(i-NumOfPadTPCIn) + pad_length_out;
      }
    }
    angle[10] = 180. - 155.35;
    angle[11] = 180. - 144.8;
    angle[12] = 180. - 138.;
    angle[13] = 180. - 116.73;
    angle[14] = 180. - 106.;
    angle[15] = 180. - 98.77;
    angle[16] = 180. - 94.29;
    angle[17] = 180. - 89.8;
    angle[18] = 180. - 87.18;
    angle[19] = 180. - 84.16;
    angle[20] = 180. - 81.48;
    angle[21] = 180. - 73.39;
    angle[22] = 180. - 65.51011;
    angle[23] = 180. - 60.19;
    angle[24] = 180. - 56.35239;
    angle[25] = 180. - 52.85;
    angle[26] = 180. - 50.14;
    angle[27] = 180. - 47.17;
    angle[28] = 180. - 41.24;
    angle[29] = 180. - 29.;
    angle[30] = 180. - 23.23;
    angle[31] = 180. - 18.69;
    break;


  case 3:
    //for tracking analysis
    //If you need the dE/dx information, it should be modified.
    //Thin sensitive detector is introduced.
    for( G4int i=0; i<NumOfPadTPC; ++i ){
      double pad_radius = padHelper::getRadius(i);
      pad_in[i] = pad_radius;
      pad_out[i] = pad_radius + 0.1*mm;
			/*
			double pad_halflength = padHelper::getLength(i)/2;
      pad_in[i] = pad_radius-pad_halflength;
      pad_out[i] = pad_radius + pad_halflength;
      */
			if( i<NumOfPadTPCIn ){
	angle[i]   = 360.;
      }
      else{
	angle[i] = padHelper::getsTheta(i);
      }
    }
    break;


  default:
    break;
  }

  G4double below_target = 0;
  if(m_experiment==42){
    below_target=0.2;
  }else if(m_experiment==45||m_experiment==27){
    below_target=60.4;
  }
  // Inner Pads
  for( G4int i=0; i<NumOfPadTPCIn; ++i ){
    if( pad_out[i]+pad_length_in/2<below_target ){
      if(m_experiment==42){
	pad_solid[i] = new G4Tubs( Form("TpcPadSolid%d", i), pad_in[i]*mm,
				   pad_out[i]*mm, 120.*mm, 0.,
				   angle[i]*deg );
      }else if(m_experiment==45||m_experiment==27){
	pad_solid[i] = new G4Tubs( Form("TpcPadSolid%d", i), pad_in[i]*mm,
				   pad_out[i]*mm, 200./2.*mm, 0.,
				   angle[i]*deg );
      }
      pad_lv[i]  = new G4LogicalVolume( pad_solid[i], m_material_map["P10"],
					Form("TpcPadLV%d", i) );
    } else {
      pad_solid[i] = new G4Tubs( "TpcPadSolid", pad_in[i]*mm, pad_out[i]*mm,
				 275.*mm, 0., angle[i]*deg );
      pad_lv[i]  = new G4LogicalVolume( pad_solid[i], m_material_map["P10"],
					Form("TpcPadLV%d", i) );
    }
    pad_lv[i]->SetVisAttributes( G4Colour::Blue() );
    if( pad_out[i] < below_target ){
      if( m_experiment == 42 ){
	new G4PVPlacement( nullptr, G4ThreeVector( 0., -pad_center_z, (-120.-25.)*mm ),
			   pad_lv[i], Form("TpcPadPV%d", i), m_tpc_lv, true, i );
      } else if( m_experiment == 45 || m_experiment == 27 ){
	new G4PVPlacement( nullptr, G4ThreeVector( 0., -pad_center_z, (-200.)*mm ),
			   pad_lv[i], Form("TpcPadPV%d", i), m_tpc_lv, true, i );
      }
    }else{
      new G4PVPlacement( nullptr, G4ThreeVector( 0., -pad_center_z, -25.*mm ),
			 pad_lv[i], Form("TpcPadPV%d", i), m_tpc_lv, true, i );
    }
  }
  G4ThreeVector padpos( 0., -pad_center_z, -25.*mm );
  for( G4int i=NumOfPadTPCIn; i<NumOfPadTPC; ++i ){
    pad_solid[i] = new G4Tubs( Form("TpcPadSolid%d", i), pad_in[i]*mm,
			       pad_out[i]*mm, 275.*mm, (90.+angle[i])*deg,
			       (360.-2.*angle[i])*deg );
    pad_lv[i]  = new G4LogicalVolume( pad_solid[i], m_material_map["P10"],
				     Form("TpcPadLV%d", i) );
    pad_lv[i]->SetVisAttributes( ORANGE );
      new G4PVPlacement( nullptr, padpos, pad_lv[i], Form("TpcPadPV%d", i),
		       m_tpc_lv, true, i );
  }
  // Dead area
  auto dead_solid = new G4Box( "DeadSolid", 5*mm, 250*mm, 0.001*mm );
  auto dead_lv = new G4LogicalVolume( dead_solid, m_material_map["Carbon"],
				      "DeadLV" );
  auto rotdead1 = new G4RotationMatrix;
  rotdead1->rotateZ( 45.*deg );
  new G4PVPlacement( rotdead1, G4ThreeVector( 0., 0.*mm, -300.1*mm ),
		     dead_lv, "DeadPV1", m_tpc_lv, true, 0 );
  auto rotdead2 = new G4RotationMatrix;
  rotdead2->rotateZ( -45.*deg );
  new G4PVPlacement( rotdead2, G4ThreeVector( 0., 0.*mm, -300.1*mm ),
		     dead_lv, "DeadPV2", m_tpc_lv, true, 1 );
  dead_lv->SetVisAttributes( G4Colour::Gray() );
  // Virtual pad
  G4Tubs* vpad_solid[NumOfPadTPC];
  G4LogicalVolume* vpad_lv[NumOfPadTPC];
  for( G4int i=0; i<NumOfPadTPCIn; ++i ){
    vpad_solid[i] = new G4Tubs( Form("TpcVPadSolid%d", i), pad_in[i]*mm,
				pad_out[i]*mm, 0.5*mm, 0., 360.*deg );
    vpad_lv[i] = new G4LogicalVolume( vpad_solid[i], m_material_map["P10"],
				      Form("TpcVPadLV%d", i) );
    vpad_lv[i]->SetVisAttributes( G4Colour::Blue() );
    new G4PVPlacement( nullptr,
		       G4ThreeVector( 0., -pad_center_z, -302.*mm ),
		       vpad_lv[i], Form("TpcVPadPV%d", i), m_tpc_lv, true, 0 );
  }
  for( G4int i=NumOfPadTPCIn; i<NumOfPadTPC; ++i ){
    vpad_solid[i] = new G4Tubs( Form("TpcVPadSolid%d", i),  pad_in[i]*mm,
				pad_out[i]*mm,
				0.5*mm, (90.+angle[i])*deg,
				(360.-2.*angle[i])*deg );
    vpad_lv[i] = new G4LogicalVolume( vpad_solid[i], m_material_map["P10"],
				      Form("TpcVPadLV%d", i) );
    vpad_lv[i]->SetVisAttributes( ORANGE );
    new G4PVPlacement( nullptr,
		       G4ThreeVector( 0., -pad_center_z, -302.*mm ),
		       vpad_lv[i], Form("TpcVPadPV%d", i), m_tpc_lv, true, 0 );
  }
  for( G4int i=0; i<NumOfPadTPC; ++i ){
    pad_lv[i]->SetSensitiveDetector( tpc_sd );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructK18BeamlineSpectrometer( void )
{
  G4cout << FUNC_NAME << G4endl;
  // D4 magnet
  const G4double D4Rho = gSize.Get( "D4Rho" )*mm;
  const G4double D4BendAngle = gSize.Get( "D4BendAngle" )*deg; // [deg]
  const auto& D4FieldSize = gSize.GetSize( "D4Field" )*mm;
  const G4double D4r1 = D4Rho - D4FieldSize.x()/2.;
  const G4double D4r2 = D4Rho + D4FieldSize.x()/2.;
  const G4double D4Width1 = 1700.*mm/2;
  // const auto& D4Coil1Size = gSize.GetSize( "D4Coil1" )*mm;
  // const G4double D4CoilWidth1 = 500.*mm/2;
  // const G4double D4CoilLength1 = 200.*mm/2;
  // const G4double D4CoilLength2 = 100.*mm/2;
  // const G4double D4Length = 4.468*m;
  const G4double D4B0 = 15.010222 * CLHEP::kilogauss;
  const G4double D4alpha = 23.5; // [deg]
  const G4double D4beta = 0; // [deg]
  // Q10 magnet
  const auto& Q10Size = gSize.GetSize( "Q10" )*mm;
  const G4double Q10B0 = -9.814096 * CLHEP::kilogauss;
  const G4double Q10a0 = 0.1*m;    // [m]
  // Q11 magnet
  const auto& Q11Size = gSize.GetSize( "Q11" )*mm;
  const G4double Q11B0 = 7.493605 * CLHEP::kilogauss;
  const G4double Q11a0 = 0.1*m;
  // Q12 magnet
  const auto& Q12Size = gSize.GetSize( "Q12" )*mm;
  const G4double Q12B0 = -5.87673 * CLHEP::kilogauss;
  const G4double Q12a0 = 0.1*m;
  // Q13 magnet
  const auto& Q13Size = gSize.GetSize( "Q13" )*mm;
  const G4double Q13B0 = 0. * CLHEP::kilogauss;
  const G4double Q13a0 = 0.1*m;      // [m]
  const G4double driftL0 = gSize.Get( "K18L0" )*mm; // VI-Q10
  const G4double driftL1 = gSize.Get( "K18L1" )*mm; // Q10-Q11
  const G4double driftL2 = gSize.Get( "K18L2" )*mm; // Q11-D4
  const G4double driftL3 = gSize.Get( "K18L3" )*mm; // D4-Q12
  const G4double driftL4 = gSize.Get( "K18L4" )*mm; // Q12-Q13
  const G4double driftL5 = gSize.Get( "K18L5" )*mm; // Q13-VO
  const G4double driftL6 = gSize.Get( "K18L6" )*mm; // VO-FF:1200, VO-HS
  // const G4double driftL7 = gSize.Get( "K18L7" )*mm; // HS-KURAMA
  const G4ThreeVector VC1Size( Q10Size.x(), Q10Size.y()/2, driftL1 );
  const G4ThreeVector VC2Size( Q10Size.x(), Q10Size.y()/2, driftL2 );
  const G4ThreeVector VC3Size( Q12Size.x(), Q12Size.y()/2, driftL3 );
  const G4ThreeVector VC4Size( Q13Size.x(), Q13Size.y()/2, driftL4 );

  G4double x,y,z;
  G4RotationMatrix rotZero;

  // VI
  x = ( D4Rho*( 1. - std::cos( D4BendAngle ) ) +
	( driftL0 + Q10Size.z() + VC1Size.z() + Q11Size.z() + VC2Size.z() ) *
	std::sin( D4BendAngle ) );
  y = 0.*m;
  z = ( - D4Rho*std::sin( D4BendAngle )
	- ( driftL0 + Q10Size.z() + VC1Size.z() + Q11Size.z() + VC2Size.z() ) *
	std::cos( D4BendAngle )
	- driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6 );
  BeamMan::GetInstance().SetVIPosition( G4ThreeVector( x, y, z ) );

  // Q10 magnet
  auto Q10Solid = new G4Box( "Q10Solid", Q10Size.x()/2,
			     Q10Size.y()/2, Q10Size.z()/2 );
  auto Q10LV = new G4LogicalVolume( Q10Solid,
				    m_material_map["Vacuum"],
				    "Q10LV" );
  Q10LV->SetVisAttributes( ORANGE );
  x = ( D4Rho*( 1. - std::cos( D4BendAngle ) ) +
	( Q10Size.z()/2 + VC1Size.z() + Q11Size.z() + VC2Size.z() ) *
	std::sin( D4BendAngle ) );
  y = 0.*m;
  z = ( - D4Rho*std::sin( D4BendAngle )
	- ( Q10Size.z()/2 + VC1Size.z() + Q11Size.z() + VC2Size.z() ) *
	std::cos( D4BendAngle )
	- driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6 );
  auto rotQ10 = new G4RotationMatrix;
  rotQ10->rotateY( D4BendAngle );
  new G4PVPlacement( rotQ10, G4ThreeVector( x, y, z ),
		     Q10LV, "Q10PV", m_world_lv, false, 0 );
  MagnetInfo Q10Info( "Q10" );
  Q10Info.type = MagnetInfo::kQuadrupole;
  Q10Info.b0 = Q10B0;
  Q10Info.a0 = Q10a0;
  Q10Info.pos = G4ThreeVector( x, y, z );
  Q10Info.size = Q10Size*0.5;
  Q10Info.ra1 = -D4BendAngle;
  myField->AddMagnetInfo( Q10Info );
  // Q11 magnet
  auto Q11Solid = new G4Box( "Q11Solid", Q11Size.x()/2,
			     Q11Size.y()/2, Q11Size.z()/2 );
  auto Q11LV = new G4LogicalVolume( Q11Solid,
				    m_material_map["Vacuum"],
				    "Q11LV" );
  Q11LV->SetVisAttributes( ORANGE );
  x = ( D4Rho*( 1. - std::cos( D4BendAngle ) ) +
	( VC2Size.z() + Q11Size.z()/2 ) * std::sin( D4BendAngle ) );
  y = 0.*m;
  z = ( - D4Rho*std::sin( D4BendAngle )
	- ( VC2Size.z() + Q11Size.z()/2 ) * std::cos( D4BendAngle )
	- driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6 );
  auto rotQ11 = new G4RotationMatrix;
  rotQ11->rotateY( D4BendAngle );
  new G4PVPlacement( rotQ11, G4ThreeVector( x, y, z ),
		     Q11LV, "Q11PV", m_world_lv, false, 0 );
  MagnetInfo Q11Info( "Q11" );
  Q11Info.type = MagnetInfo::kQuadrupole;
  Q11Info.b0 = Q11B0;
  Q11Info.a0 = Q11a0;
  Q11Info.pos = G4ThreeVector( x, y, z );
  Q11Info.size = Q11Size*0.5;
  Q11Info.ra1 = -D4BendAngle;
  myField->AddMagnetInfo( Q11Info );
  // D4 magnet
  G4cout << "   D4Rho = " << D4Rho << G4endl;
  G4cout << "   D4r1 = " << D4r1/m  << ", D4r2 = " << D4r2/m  << G4endl;
  auto D4Solid = new G4Tubs( "D4Solid", D4r1, D4r2, D4FieldSize.y()/2,
			     0.*deg, D4BendAngle );
  auto D4OutSolid = new G4Tubs( "D4OutSolid", D4Rho - D4Width1,
				D4Rho + D4Width1, Q13Size.y()/2,
				0.*deg, D4BendAngle );
  auto D4PoleSolid = new G4SubtractionSolid( "D4PoleSolid", D4OutSolid,
					     D4Solid );
  auto D4PoleLV = new G4LogicalVolume( D4PoleSolid,
				       m_material_map["Iron"],
				       "D4PoleLV" );
  D4PoleLV->SetVisAttributes( G4Colour::Green() );
  x = D4Rho;
  y = 0.*m;
  z = - driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6;
  auto D4Rot = new G4RotationMatrix;
  D4Rot->rotateX( 90.*deg );
  D4Rot->rotateZ( 180.*deg + D4BendAngle );
  new G4PVPlacement( D4Rot, G4ThreeVector( x, y, z ),
		     D4PoleLV, "D4PolePV", m_world_lv, false, 0 );
  auto D4LV = new G4LogicalVolume( D4Solid,
				   m_material_map["Vacuum"],
				   "D4LV" );
  D4LV->SetVisAttributes( G4Colour::White() );
  new G4PVPlacement( D4Rot, G4ThreeVector( x, y, z ),
		     D4LV, "D4PV", m_world_lv, false, 0 );
#if 0
  auto D4Coil1Solid = new G4Box( "D4Coil1Solid", D4Coil1Size.x()/2,
				 D4Coil1Size.y()/2, D4Coil1Size.z()/2 );
  auto D4Coil1LV = new G4LogicalVolume( D4Coil1Solid,
					m_material_map["Iron"],
					"D4Coil1LV" );
  D4Coil1LV->SetVisAttributes( G4Colour::Red() );
  x = 0.*m;
  y = 0.*m;
  z = ( D4Coil1Size.z()/2 - driftL3 - Q12Size.z() - driftL4 - Q13Size.z()
	- driftL5 - driftL6 );
  new G4PVPlacement( nullptr, G4ThreeVector( x, y, z ),
		     D4Coil1LV, "D4Coil1PV", m_world_lv, false, 0 );
#endif
  MagnetInfo D4Info( "D4" );
  D4Info.type = MagnetInfo::kDipole;
  D4Info.b0 = D4B0;
  D4Info.rho = D4Rho;
  D4Info.alpha = D4alpha;
  D4Info.beta = D4beta;
  D4Info.pos = G4ThreeVector( x, y, z );
  D4Info.size = D4FieldSize*0.5;
  D4Info.ra1 = -D4BendAngle;
  D4Info.bend = D4BendAngle/deg;
  myField->AddMagnetInfo( D4Info );
  // Q12 magnet
  auto Q12Solid = new G4Box( "Q12Solid", Q12Size.x()/2,
			     Q12Size.y()/2, Q12Size.z()/2 );
  auto Q12LV = new G4LogicalVolume( Q12Solid,
				    m_material_map["Vacuum"],
				    "Q12LV" );
  Q12LV->SetVisAttributes( ORANGE );
  x = 0.*m;
  y = 0.*m;
  z = - Q12Size.z()/2 - driftL4 - Q13Size.z() - driftL5 - driftL6;
  new G4PVPlacement( nullptr, G4ThreeVector( x, y, z ),
		     Q12LV, "Q12PV", m_world_lv, false, 0 );
  MagnetInfo Q12Info( "Q12" );
  Q12Info.type = MagnetInfo::kQuadrupole;
  Q12Info.b0 = Q12B0;
  Q12Info.a0 = Q12a0;
  Q12Info.pos = G4ThreeVector( x, y, z );
  Q12Info.size = Q12Size*0.5;
  Q12Info.ra1 = -D4BendAngle;
  myField->AddMagnetInfo( Q12Info );
  // Q13 magnet
  auto Q13Solid = new G4Box( "Q13Solid", Q13Size.x()/2,
			     Q13Size.y()/2, Q13Size.z()/2 );
  auto Q13LV = new G4LogicalVolume( Q13Solid,
				    m_material_map["Vacuum"],
				    "Q13LV" );
  Q13LV->SetVisAttributes( ORANGE );
  x = 0.*m;
  y = 0.*m;
  z =  - Q13Size.z()/2 - driftL5 - driftL6;
  new G4PVPlacement( nullptr, G4ThreeVector( x, y, z ),
		     Q13LV, "Q13PV", m_world_lv, false, 0 );
  MagnetInfo Q13Info( "Q13" );
  Q13Info.type = MagnetInfo::kQuadrupole;
  Q13Info.b0 = Q13B0;
  Q13Info.a0 = Q13a0;
  Q13Info.pos = G4ThreeVector( x, y, z );
  Q13Info.size = Q13Size*0.5;
  Q13Info.ra1 = -D4BendAngle;
  myField->AddMagnetInfo( Q13Info );
  // Vacuum Chamber 1
  auto VC1Solid = new G4Box( "VC1Solid", VC1Size.x()/2,
			     VC1Size.y()/2, VC1Size.z()/2 );
  auto VC1LV = new G4LogicalVolume( VC1Solid,
				    m_material_map["Vacuum"],
				    "VC1LV" );
  x = ( D4Rho*( 1. - std::cos( D4BendAngle ) ) +
	( VC1Size.z()/2 + Q11Size.z() + VC2Size.z() ) * std::sin( D4BendAngle ) );
  y = 0.*m;
  z = ( - D4Rho*std::sin( D4BendAngle )
	- ( VC1Size.z()/2 + Q11Size.z() + VC2Size.z() ) * std::cos( D4BendAngle )
	- driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6 );
  auto rotVC1 = new G4RotationMatrix;
  rotVC1->rotateY( D4BendAngle );
  new G4PVPlacement( rotVC1, G4ThreeVector( x, y, z ),
		     VC1LV, "VC1PV", m_world_lv, false, 0 );
  // Vacuum Chamber 2
  auto VC2Solid = new G4Box( "slidVC2", VC2Size.x()/2,
			     VC2Size.y()/2, VC2Size.z()/2 );
  auto VC2LV = new G4LogicalVolume( VC2Solid,
				    m_material_map["Vacuum"],
				    "VC2LV" );
  x = ( D4Rho*( 1. - std::cos( D4BendAngle ) ) +
	VC2Size.z()/2 * std::sin( D4BendAngle ) );
  y = 0.*m;
  z = ( - D4Rho*std::sin( D4BendAngle )
	- VC2Size.z()/2 * std::cos( D4BendAngle )
	- driftL3 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6 );
  auto rotVC2 = new G4RotationMatrix;
  rotVC2->rotateY( D4BendAngle );
  new G4PVPlacement( rotVC2, G4ThreeVector( x, y, z ),
		     VC2LV, "VC2PV", m_world_lv, false, 0 );
  // Vacuum Chamber 3
  auto VC3Solid = new G4Box( "VC3Solid", VC3Size.x()/2,
			     VC3Size.y()/2, VC3Size.z()/2 );
  auto VC3LV = new G4LogicalVolume( VC3Solid,
				    m_material_map["Vacuum"],
				    "VC3LV" );
  x = 0.*m;
  y = 0.*m;
  z = -driftL3/2 - Q12Size.z() - driftL4 - Q13Size.z() - driftL5 - driftL6;
  new G4PVPlacement( nullptr, G4ThreeVector( x, y, z ),
		     VC3LV, "VC3PV", m_world_lv, false, 0 );
  // Vacuum Chamber 4
  auto VC4Solid = new G4Box( "VC4Solid", VC4Size.x()/2, VC4Size.y()/2, VC4Size.z()/2 );
  auto VC4LV = new G4LogicalVolume( VC4Solid,
				    m_material_map["Vacuum"],
				    "VC4LV" );
  x = 0.*m;
  y = 0.*m;
  z = - driftL4/2 - Q13Size.z() - driftL5 - driftL6;
  new G4PVPlacement( nullptr, G4ThreeVector( x, y, z ),
		     VC4LV, "VC4PV", m_world_lv, false, 0 );
  myField->SetStatusK18Field( true );
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructKuramaMagnet( void )
{
  auto vp_sd = new TPCVPSD("/VP");
  G4SDManager::GetSDMpointer()->AddNewDetector( vp_sd );
  const auto field_size = gSize.GetSize("KuramaField") * 0.5 * mm;
  const auto coil_color = G4Colour::Red();
  // size
  const G4ThreeVector coil1_size( 900.*mm/2, 193.*mm/2, 280.*mm/2 );
  const G4ThreeVector coil2_size( 193.*mm/2, 117.*mm/2, 280.*mm/2 );
  const G4ThreeVector coil3_size( 193.*mm/2, 137.*mm/2, 280.*mm/2 );
  const G4ThreeVector coil4_size( 193.*mm/2, 280.*mm/2, 740.*mm/2 );
  const G4ThreeVector coil5_size( 900.*mm/2, 214.*mm/2, 280.*mm/2 );
  const G4ThreeVector yoke_inner_size( field_size.x() + 193.*mm,
				       field_size.y(), field_size.z() );
  const G4ThreeVector yoke_outer_size( 1800.*mm/2, 2200.*mm/2,
				       field_size.z() );
  const G4ThreeVector yoke_ud_size( 2200.*mm/2, 370.*mm/2, field_size.z() );
  const G4ThreeVector yoke_lr_size( 200.*mm/2, field_size.y(),
				    field_size.z() );
  const G4ThreeVector uguard_inner_size( 600.*mm/2, 300.*mm/2, 100.*mm/2 );
  const G4ThreeVector uguard_inner2_size( 50.*mm/2, 50.*mm/2, 100.*mm/2 );
  const G4ThreeVector uguard_outer_size( 1600.*mm/2, 2200.*mm/2, 100.*mm/2 );
  const G4ThreeVector dguard_inner_size( 1100.*mm/2, 1100.*mm/2, 100.*mm/2 );
  const G4ThreeVector dguard_outer_size( 1600.*mm/2, 2200.*mm/2, 100.*mm/2 );
  // position
  const G4ThreeVector field_pos = gGeom.GetGlobalPosition("KURAMA");
  const G4ThreeVector yoke_pos( field_pos );
  const G4ThreeVector uguard_pos( field_pos.x(), field_pos.y(),
				  field_pos.z() - ( 820. - 50. )*mm );
  const G4ThreeVector dguard_pos( field_pos.x(), field_pos.y(),
				  field_pos.z() + ( 820. - 50. )*mm );
  const G4ThreeVector coil1_pos =
    { field_pos.x(),
      field_size.y() + yoke_ud_size.y()*2. - ( coil1_size.y() + 20.*mm ),
      field_pos.z() - field_size.z() - ( coil1_size.z() + 20.*mm ) };
  const G4ThreeVector coil4l_pos =
    { field_pos.x() + field_size.x() + coil4_size.x(),
      field_size.y()/2,
      field_pos.z() };
  const G4ThreeVector coil4r_pos =
    { field_pos.x() - field_size.x() - coil4_size.x(),
      field_size.y()/2,
      field_pos.z() };
  const G4ThreeVector coil2l_pos =
    { field_pos.x() + field_size.x() + coil4_size.x(),
      ( coil4l_pos.y() + coil1_pos.y() + coil4_size.y() - coil1_size.y() )/2,
      field_pos.z() - field_size.z() - coil1_size.z() + 20.*mm };
  const G4ThreeVector coil2r_pos =
    { field_pos.x() - field_size.x() - coil4_size.x(),
      ( coil4r_pos.y() + coil1_pos.y() + coil4_size.y() - coil1_size.y() )/2,
      field_pos.z() - field_size.z() - coil1_size.z() + 20.*mm };
  const G4ThreeVector coil5_pos =
    { field_pos.x(),
      field_size.y() + yoke_ud_size.y()*2. - coil5_size.y(),
      field_pos.z() + field_size.z() + coil5_size.z() + 21.*mm };
  const G4ThreeVector coil3l_pos =
    { field_pos.x() + field_size.x() + coil4_size.x(),
      ( coil4l_pos.y() + coil5_pos.y() + coil4_size.y() - coil5_size.y() )/2,
      field_pos.z() + field_size.z() + coil5_size.z() + 21.*mm };
  const G4ThreeVector coil3r_pos =
    { field_pos.x() - field_size.x() - coil4_size.x(),
      ( coil4r_pos.y() + coil5_pos.y() + coil4_size.y() - coil5_size.y() )/2,
      field_pos.z() + field_size.z() + coil5_size.z() + 21.*mm };
  const G4ThreeVector yoke_u_pos( field_pos.x(),
				  field_size.y() + yoke_ud_size.y(),
				  field_pos.z() );
  const G4ThreeVector yoke_d_pos( field_pos.x(),
				  - field_size.y() - yoke_ud_size.y(),
				  field_pos.z() );
  const G4ThreeVector yoke_l_pos =
    { field_pos.x() + field_size.x() + yoke_lr_size.x() + 200.*mm,
      0.*mm,
      field_pos.z() };
  const G4ThreeVector yoke_r_pos =
    { field_pos.x() - field_size.x() + yoke_lr_size.x() - 200.*mm,
      0.*mm,
      field_pos.z() };
  // Construct KURAMA Magnet
  // auto kurama_solid = new G4Box( "KuramaSolid", 4.*m/2, 3.*m/2, 4.*m/2 );
  // auto kurama_lv = new G4LogicalVolume( kurama_solid, m_material_map["Air"],
  // 					"KuramaLV" );
  // kurama_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // auto kurama_pv = new G4PVPlacement( m_rotation_angle, G4ThreeVector(), "KuramaPV",
  // 				      kurama_lv, m_world_pv, false, 0 );

  // Coil1
  auto coil1_solid = new G4Box( "Coil1Solid", coil1_size.x(),
				coil1_size.y(), coil1_size.z() );
  auto coil1_lv = new G4LogicalVolume( coil1_solid, m_material_map["Copper"],
				       "Coil1LV" );
  coil1_lv->SetVisAttributes( coil_color );
  auto coil1u_pos = coil1_pos;
  coil1u_pos.rotateY( m_rotation_angle );
  new G4PVPlacement( m_rotation_matrix, coil1u_pos, coil1_lv,
		     "Coil1UPV", m_world_lv, false, 0 );
  auto coil1d_pos = G4ThreeVector( coil1_pos.x(),
				   -coil1_pos.y(),
				   coil1_pos.z() );
  coil1d_pos.rotateY( m_rotation_angle );
  new G4PVPlacement( m_rotation_matrix,	coil1d_pos, coil1_lv,
		     "Coil1DPV", m_world_lv, false, 0 );

  // Coil4
  auto coil4_solid = new G4Box( "Coil4Solid", coil4_size.x(),
				coil4_size.y(), coil4_size.z() );
  auto coil4_lv = new G4LogicalVolume( coil4_solid, m_material_map["Copper"],
				       "Coil4LV" );
  coil4_lv->SetVisAttributes( coil_color );
  auto coil4ur_pos = coil4l_pos;
  coil4ur_pos.rotateY( m_rotation_angle );
  new G4PVPlacement( m_rotation_matrix, coil4ur_pos, coil4_lv,
		     "Coil4URPV", m_world_lv, false, 0 );
  auto coil4ul_pos = coil4r_pos;
  coil4ul_pos.rotateY( m_rotation_angle );
  new G4PVPlacement( m_rotation_matrix, coil4ul_pos, coil4_lv,
		     "Coil4ULPV", m_world_lv, false, 0 );
  auto coil4dr_pos = G4ThreeVector( coil4l_pos.x(),
				    -coil4l_pos.y(),
				    coil4l_pos.z() ).rotateY( m_rotation_angle );
  new G4PVPlacement( m_rotation_matrix, coil4dr_pos, coil4_lv,
		     "Coil4DRPV", m_world_lv, false, 0 );
  auto coil4dl_pos = G4ThreeVector( coil4r_pos.x(),
				    -coil4r_pos.y(),
				    coil4r_pos.z() ).rotateY( m_rotation_angle );
  new G4PVPlacement( m_rotation_matrix, coil4ul_pos, coil4_lv,
		     "Coil4DLPV", m_world_lv, false, 0 );

  // Coil5
  auto coil5_solid = new G4Box( "Coil5Solid", coil5_size.x(),
				coil5_size.y(), coil5_size.z() );
  auto coil5_lv = new G4LogicalVolume( coil5_solid, m_material_map["Copper"],
				       "Coil5LV" );
  coil5_lv->SetVisAttributes( coil_color );
  auto coil5u_pos = G4ThreeVector( coil5_pos.x(),
				   coil5_pos.y(),
				   coil5_pos.z() ).rotateY( m_rotation_angle );
  new G4PVPlacement( m_rotation_matrix,	coil5u_pos, coil5_lv,
		     "Coil5UPV", m_world_lv, false, 0 );
  auto coil5d_pos = G4ThreeVector( coil5_pos.x(),
				   -coil5_pos.y(),
				   coil5_pos.z() ).rotateY( m_rotation_angle );
  new G4PVPlacement( m_rotation_matrix, coil5d_pos, coil5_lv,
		     "Coil5DPV", m_world_lv, false, 0 );

  // Coil6
  G4double size_COIL6[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL6[0] = 50.0*mm;
  size_COIL6[1] = coil1_size.y()*2+size_COIL6[0];
  //  size_COIL6[1] = 330.*mm;
  size_COIL6[2] = 280.*mm/2;
  size_COIL6[3] = 90.*deg;

  G4double pos_COIL6LU[3];
  G4double pos_COIL6RU[3];
  G4double pos_COIL6LD[3];
  G4double pos_COIL6RD[3];
  //LU
  pos_COIL6LU[ThreeVector::X] = field_pos.x() +field_size.x()-size_COIL6[0];
  pos_COIL6LU[ThreeVector::Y] = coil1_pos.y()  -(size_COIL6[0]+coil1_size.y());
  pos_COIL6LU[ThreeVector::Z] = field_pos.z() - field_size.z()-(size_COIL6[ThreeVector::Z]+21.*mm);
  //RU
  pos_COIL6RU[ThreeVector::X] = field_pos.x() -field_size.x()+size_COIL6[0];
  pos_COIL6RU[ThreeVector::Y] = coil1_pos.y()  -(size_COIL6[0]+coil1_size.y());
  pos_COIL6RU[ThreeVector::Z] = field_pos.z() - field_size.z()-(size_COIL6[ThreeVector::Z]+21.*mm);
  //LD
  pos_COIL6LD[ThreeVector::X] = field_pos.x()  +field_size.x()-size_COIL6[0];
  pos_COIL6LD[ThreeVector::Y] = -coil1_pos.y() +(size_COIL6[0]+coil1_size.y());
  pos_COIL6LD[ThreeVector::Z] = field_pos.z()  - field_size.z()-(size_COIL6[ThreeVector::Z]+21.*mm);
  //RD
  pos_COIL6RD[ThreeVector::X] = field_pos.x()  -field_size.x()+size_COIL6[0];
  pos_COIL6RD[ThreeVector::Y] = -coil1_pos.y() +(size_COIL6[0]+coil1_size.y());
  pos_COIL6RD[ThreeVector::Z] = field_pos.z()  - field_size.z()-(size_COIL6[ThreeVector::Z]+21.*mm);


  G4Tubs* Coil6_tub = new G4Tubs("Coil6_tubs",
				 size_COIL6[0],size_COIL6[1],size_COIL6[2],0.,size_COIL6[3]);
  G4LogicalVolume*  Coil6_log = new G4LogicalVolume(Coil6_tub, m_material_map["Copper"], "Coil6_log",0,0,0);
  Coil6_log->SetVisAttributes( coil_color );
  G4RotationMatrix *rotcoil6lu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6ru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6ld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6rd = new G4RotationMatrix();
  rotcoil6lu->rotateY( - m_rotation_angle );
  rotcoil6ru->rotateY( - m_rotation_angle );
  rotcoil6ld->rotateY( - m_rotation_angle );
  rotcoil6rd->rotateY( - m_rotation_angle );

  rotcoil6lu->rotateZ(0.*deg);
  new G4PVPlacement( rotcoil6lu,
		     G4ThreeVector( pos_COIL6LU[ThreeVector::X],
				    pos_COIL6LU[ThreeVector::Y],
				    pos_COIL6LU[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil6_log,
		     "Coil6LU_phys",
		     m_world_lv,
		     false,
		     0 );
  rotcoil6ru->rotateZ(-90.*deg);
  new G4PVPlacement( rotcoil6ru,
		     G4ThreeVector( pos_COIL6RU[ThreeVector::X],
				    pos_COIL6RU[ThreeVector::Y],
				    pos_COIL6RU[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil6_log,
		     "Coil6RU_phys",
		     m_world_lv,
		     false,
		     0 );
  rotcoil6ld->rotateZ(-180.*deg);
  new G4PVPlacement( rotcoil6ld,
		     G4ThreeVector( pos_COIL6RD[ThreeVector::X],
				    pos_COIL6RD[ThreeVector::Y],
				    pos_COIL6RD[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil6_log,
		     "Coil6RD_phys",
		     m_world_lv,
		     false,
		     0 );
  rotcoil6rd->rotateZ(-270.*deg);
  new G4PVPlacement( rotcoil6rd,
		     G4ThreeVector( pos_COIL6LD[ThreeVector::X],
				    pos_COIL6LD[ThreeVector::Y],
				    pos_COIL6LD[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil6_log,
		     "Coil6LD_phys",
		     m_world_lv,
		     false,
		     0 );
  //  m_rotation_matrix->rotateZ(90.*deg);
  //  m_rotation_matrix->rotateZ(90.*deg);

  //////////////coil8RLUD
  G4double size_COIL8[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL8[0] = 50.0*mm;
  size_COIL8[1] = coil5_size.y()*2.+size_COIL8[0];
  size_COIL8[2] = 280.*mm/2;
  size_COIL8[3] = 90.*deg;

  G4double pos_COIL8LU[3];
  G4double pos_COIL8RU[3];
  G4double pos_COIL8LD[3];
  G4double pos_COIL8RD[3];
  //LU
  pos_COIL8LU[ThreeVector::X] = field_pos.x() +field_size.x()-size_COIL8[0];
  pos_COIL8LU[ThreeVector::Y] = coil5_pos.y()  -(size_COIL8[0]+coil5_size.y());
  pos_COIL8LU[ThreeVector::Z] = field_pos.z() + field_size.z()+(size_COIL8[ThreeVector::Z]+21.*mm);
  //RU
  pos_COIL8RU[ThreeVector::X] = field_pos.x() -field_size.x()+size_COIL8[0];
  pos_COIL8RU[ThreeVector::Y] = coil5_pos.y()  -(size_COIL8[0]+coil5_size.y());
  pos_COIL8RU[ThreeVector::Z] = field_pos.z() + field_size.z()+(size_COIL8[ThreeVector::Z]+21.*mm);
  //LD
  pos_COIL8LD[ThreeVector::X] = field_pos.x()  +field_size.x()-size_COIL8[0];
  pos_COIL8LD[ThreeVector::Y] = -coil5_pos.y() +(size_COIL8[0]+coil5_size.y());
  pos_COIL8LD[ThreeVector::Z] = field_pos.z()  + field_size.z()+(size_COIL8[ThreeVector::Z]+21.*mm);
  //RD
  pos_COIL8RD[ThreeVector::X] = field_pos.x()  -field_size.x()+size_COIL8[0];
  pos_COIL8RD[ThreeVector::Y] = -coil5_pos.y() +(size_COIL8[0]+coil5_size.y());
  pos_COIL8RD[ThreeVector::Z] = field_pos.z()  + field_size.z()+(size_COIL8[ThreeVector::Z]+21.*mm);

  G4Tubs* Coil8_tub = new G4Tubs("Coil8_tubs",
				 size_COIL8[0],size_COIL8[1],size_COIL8[2],0.,size_COIL8[3]);
  G4LogicalVolume*  Coil8_log = new G4LogicalVolume(Coil8_tub, m_material_map["Copper"], "Coil8_log",0,0,0);
  Coil8_log->SetVisAttributes( coil_color );
  G4RotationMatrix *rotcoil8lu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8ru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8ld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8rd = new G4RotationMatrix();
  rotcoil8lu->rotateY( - m_rotation_angle );
  rotcoil8ru->rotateY( - m_rotation_angle );
  rotcoil8ld->rotateY( - m_rotation_angle );
  rotcoil8rd->rotateY( - m_rotation_angle );

  rotcoil8lu->rotateZ(0.*deg);
  new G4PVPlacement( rotcoil8lu,
		     G4ThreeVector( pos_COIL8LU[ThreeVector::X],
				    pos_COIL8LU[ThreeVector::Y],
				    pos_COIL8LU[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil8_log,
		     "Coil8LU_phys",
		     m_world_lv,
		     false,
		     0 );
  rotcoil8ru->rotateZ(-90.*deg);
  new G4PVPlacement( rotcoil8ru,
		     G4ThreeVector( pos_COIL8RU[ThreeVector::X],
				    pos_COIL8RU[ThreeVector::Y],
				    pos_COIL8RU[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil8_log,
		     "Coil8RU_phys",
		     m_world_lv,
		     false,
		     0 );
  rotcoil8ld->rotateZ(-180.*deg);
  new G4PVPlacement( rotcoil8ld,
		     G4ThreeVector( pos_COIL8RD[ThreeVector::X],
				    pos_COIL8RD[ThreeVector::Y],
				    pos_COIL8RD[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil8_log,
		     "Coil8RD_phys",
		     m_world_lv,
		     false,
		     0);
  rotcoil8rd->rotateZ(-270.*deg);
  new G4PVPlacement( rotcoil8rd,
		     G4ThreeVector( pos_COIL8LD[ThreeVector::X],
				    pos_COIL8LD[ThreeVector::Y],
				    pos_COIL8LD[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil8_log,
		     "Coil8LD_phys",
		     m_world_lv,
		     false,
		     0 );
  //  m_rotation_matrix->rotateZ(90.*deg);
  //  m_rotation_matrix->rotateZ(90.*deg);

  //////////////coil7RLUD
  G4double size_COIL7[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL7[0] = 50.0*mm;
  size_COIL7[1] = coil4_size.y()*2.+size_COIL7[0];
  size_COIL7[2] = coil4_size.x();
  size_COIL7[3] = 90.*deg;

  G4double pos_COIL7ULU[3];
  G4double pos_COIL7URU[3];
  G4double pos_COIL7ULD[3];
  G4double pos_COIL7URD[3];
  G4double pos_COIL7DLU[3];
  G4double pos_COIL7DRU[3];
  G4double pos_COIL7DLD[3];
  G4double pos_COIL7DRD[3];
  //ULU
  // pos_COIL7ULU[ThreeVector::X] = field_pos.x() +field_size.x()+size_COIL7[0];
  pos_COIL7ULU[ThreeVector::X] = coil4l_pos.x();
  pos_COIL7ULU[ThreeVector::Y] = coil4l_pos.y() +(size_COIL7[0]+coil4_size.y());
  pos_COIL7ULU[ThreeVector::Z] = field_pos.z() - coil4_size.z();
  //URU
  pos_COIL7URU[ThreeVector::X] = coil4r_pos.x();
  pos_COIL7URU[ThreeVector::Y] = coil4r_pos.y() +(size_COIL7[0]+coil4_size.y());
  pos_COIL7URU[ThreeVector::Z] = field_pos.z() - coil4_size.z();
  //ULD
  pos_COIL7ULD[ThreeVector::X] = coil4l_pos.x();
  pos_COIL7ULD[ThreeVector::Y] = -coil4l_pos.y() -(size_COIL7[0]+coil4_size.y());
  pos_COIL7ULD[ThreeVector::Z] = field_pos.z()  - coil4_size.z();
  //URD
  pos_COIL7URD[ThreeVector::X] = coil4r_pos.x();
  pos_COIL7URD[ThreeVector::Y] = -coil4r_pos.y() -(size_COIL7[0]+coil4_size.y());
  pos_COIL7URD[ThreeVector::Z] = field_pos.z() - coil4_size.z();


  //DLU
  //  pos_COIL7ULU[ThreeVector::X] = field_pos.x() +field_size.x()+size_COIL7[0];
  pos_COIL7DLU[ThreeVector::X] = coil4l_pos.x();
  pos_COIL7DLU[ThreeVector::Y] = coil4l_pos.y() +(size_COIL7[0]+coil4_size.y());
  pos_COIL7DLU[ThreeVector::Z] = field_pos.z() + coil4_size.z();
  //DRU
  pos_COIL7DRU[ThreeVector::X] = coil4r_pos.x();
  pos_COIL7DRU[ThreeVector::Y] = coil4r_pos.y() +(size_COIL7[0]+coil4_size.y());
  pos_COIL7DRU[ThreeVector::Z] = field_pos.z() + coil4_size.z();
  //DLD
  pos_COIL7DLD[ThreeVector::X] = coil4l_pos.x();
  pos_COIL7DLD[ThreeVector::Y] = -coil4l_pos.y() -(size_COIL7[0]+coil4_size.y());
  pos_COIL7DLD[ThreeVector::Z] = field_pos.z()  + coil4_size.z();
  //DRD
  pos_COIL7DRD[ThreeVector::X] = coil4r_pos.x();
  pos_COIL7DRD[ThreeVector::Y] = -coil4r_pos.y() -(size_COIL7[0]+coil4_size.y());
  pos_COIL7DRD[ThreeVector::Z] = field_pos.z() + coil4_size.z();


  G4Tubs* Coil7_tub = new G4Tubs("Coil7_tubs",
				 size_COIL7[0],size_COIL7[1],size_COIL7[2],0.,size_COIL7[3]);
  G4LogicalVolume*  Coil7_log = new G4LogicalVolume(Coil7_tub, m_material_map["Copper"], "Coil7_log",0,0,0);
  Coil7_log->SetVisAttributes( coil_color );
  G4RotationMatrix *rotcoil7ulu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7uru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7uld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7urd = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dlu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7drd = new G4RotationMatrix();
  rotcoil7ulu->rotateY( - m_rotation_angle );
  rotcoil7uru->rotateY( - m_rotation_angle );
  rotcoil7uld->rotateY( - m_rotation_angle );
  rotcoil7urd->rotateY( - m_rotation_angle );
  rotcoil7dlu->rotateY( - m_rotation_angle );
  rotcoil7dru->rotateY( - m_rotation_angle );
  rotcoil7dld->rotateY( - m_rotation_angle );
  rotcoil7drd->rotateY( - m_rotation_angle );

  rotcoil7ulu->rotateZ(0.*deg);
  rotcoil7ulu->rotateX(180.*deg);
  rotcoil7ulu->rotateY(90.*deg);
  new G4PVPlacement( rotcoil7ulu,
		     G4ThreeVector( pos_COIL7ULU[ThreeVector::X],
				    pos_COIL7ULU[ThreeVector::Y],
				    pos_COIL7ULU[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil7_log,
		     "Coil7ULU_phys",
		     m_world_lv,
		     false,
		     0 );
  rotcoil7uru->rotateZ(0.*deg);
  rotcoil7uru->rotateX(180.*deg);
  rotcoil7uru->rotateY(90.*deg);
  new G4PVPlacement( rotcoil7uru,
		     G4ThreeVector( pos_COIL7URU[ThreeVector::X],
				    pos_COIL7URU[ThreeVector::Y],
				    pos_COIL7URU[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil7_log,
		     "Coil7URU_phys",
		     m_world_lv,
		     false,
		     0);
  rotcoil7urd->rotateZ(0.*deg);
  rotcoil7urd->rotateX(90.*deg);
  rotcoil7urd->rotateY(90.*deg);
  new G4PVPlacement( rotcoil7urd,
		     G4ThreeVector( pos_COIL7URD[ThreeVector::X],
				    pos_COIL7URD[ThreeVector::Y],
				    pos_COIL7URD[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil7_log,
		     "Coil7URD_phys",
		     m_world_lv,
		     false,
		     0 );
  rotcoil7uld->rotateZ(0.*deg);
  rotcoil7uld->rotateX(90.*deg);
  rotcoil7uld->rotateY(90.*deg);
  new G4PVPlacement( rotcoil7uld,
		     G4ThreeVector( pos_COIL7ULD[ThreeVector::X],
				    pos_COIL7ULD[ThreeVector::Y],
				    pos_COIL7ULD[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil7_log,
		     "Coil7ULD_phys",
		     m_world_lv,
		     false,
		     0);
  rotcoil7dlu->rotateZ(0.*deg);
  rotcoil7dlu->rotateX(-90.*deg);
  rotcoil7dlu->rotateY(90.*deg);
  new G4PVPlacement( rotcoil7dlu,
		     G4ThreeVector( pos_COIL7DLU[ThreeVector::X],
				    pos_COIL7DLU[ThreeVector::Y],
				    pos_COIL7DLU[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil7_log,
		     "Coil7DLU_phys",
		     m_world_lv,
		     false,
		     0);
  rotcoil7dru->rotateZ(0.*deg);
  rotcoil7dru->rotateX(-90.*deg);
  rotcoil7dru->rotateY(90.*deg);
  new G4PVPlacement( rotcoil7dru,
		     G4ThreeVector( pos_COIL7DRU[ThreeVector::X],
				    pos_COIL7DRU[ThreeVector::Y],
				    pos_COIL7DRU[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil7_log,
		     "Coil7DRU_phys",
		     m_world_lv,
		     false,
		     0 );
  rotcoil7drd->rotateZ(0.*deg);
  rotcoil7drd->rotateX(0.*deg);
  rotcoil7drd->rotateY(90.*deg);
  new G4PVPlacement( rotcoil7drd,
		     G4ThreeVector( pos_COIL7DRD[ThreeVector::X],
				    pos_COIL7DRD[ThreeVector::Y],
				    pos_COIL7DRD[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil7_log,
		     "Coil7DRD_phys",
		     m_world_lv,
		     false,
		     0 );
  rotcoil7dld->rotateZ(0.*deg);
  rotcoil7dld->rotateX(0.*deg);
  rotcoil7dld->rotateY(90.*deg);
  new G4PVPlacement( rotcoil7dld,
		     G4ThreeVector( pos_COIL7DLD[ThreeVector::X],
				    pos_COIL7DLD[ThreeVector::Y],
				    pos_COIL7DLD[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil7_log,
		     "Coil7DLD_phys",
		     m_world_lv,
		     false,
		     0 );

  ///coil2
  //////////////coil2
  G4Box* Coil2_box = new G4Box( "Coil2_box",
				coil2_size.x(),
				coil2_size.y(),
				coil2_size.z() );
  G4LogicalVolume*  Coil2_log = new G4LogicalVolume(Coil2_box, m_material_map["Copper"], "Coil2_log",0,0,0);
  Coil2_log->SetVisAttributes( coil_color );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( coil2l_pos.x(),
				    coil2l_pos.y(),
				    coil2l_pos.z() ).rotateY( m_rotation_angle ),
		     Coil2_log,
		     "Coil2UL_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( coil2r_pos.x(),
				    coil2r_pos.y(),
				    coil2r_pos.z() ).rotateY( m_rotation_angle ),
		     Coil2_log,
		     "Coil2UR_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( coil2l_pos.x(),
				    -coil2l_pos.y(),
				    coil2l_pos.z() ).rotateY( m_rotation_angle ),
		     Coil2_log,
		     "Coil2DL_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( coil2r_pos.x(),
				    -coil2r_pos.y(),
				    coil2r_pos.z() ).rotateY( m_rotation_angle ),
		     Coil2_log,
		     "Coil2DR_phys",
		     m_world_lv,
		     false,
		     0 );

  ///coil3
  G4Box* Coil3_box = new G4Box("Coil3_box",
			       coil3_size.x(), coil3_size.y(), coil3_size.z() );
  G4LogicalVolume*  Coil3_log = new G4LogicalVolume(Coil3_box, m_material_map["Copper"], "Coil3_log",0,0,0);
  Coil3_log->SetVisAttributes( coil_color );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( coil3l_pos.x(),
				    coil3l_pos.y(),
				    coil3l_pos.z() ).rotateY( m_rotation_angle ),
		     Coil3_log,
		     "Coil3UL_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( coil3r_pos.x(),
				    coil3r_pos.y(),
				    coil3r_pos.z() ).rotateY( m_rotation_angle ),
		     Coil3_log,
		     "Coil3UR_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( coil3l_pos.x(),
				    -coil3l_pos.y(),
				    coil3l_pos.z() ).rotateY( m_rotation_angle ),
		     Coil3_log,
		     "Coil3DL_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( coil3r_pos.x(),
				    -coil3r_pos.y(),
				    coil3r_pos.z() ).rotateY( m_rotation_angle ),
		     Coil3_log,
		     "Coil3DR_phys",
		     m_world_lv,
		     false,
		     0 );
  // Upstream End Guard
  auto uguard_inner_solid = new G4Box( "UGuardInnerSolid",
				       uguard_inner_size.x(),
				       uguard_inner_size.y(),
				       uguard_inner_size.z() );
  auto uguard_inner2_solid = new G4Box( "UGuardInner2Solid",
					uguard_inner2_size.x(),
					uguard_inner2_size.y(),
					uguard_inner2_size.z() );
  auto uguard_outer_solid = new G4Box( "UGuardOuterSolid",
				       uguard_outer_size.x(),
				       uguard_outer_size.y(),
				       uguard_outer_size.z() );
  auto sub_solid = new G4SubtractionSolid( "SubSolid",
					   uguard_outer_solid,
					   uguard_inner_solid );
  G4ThreeVector uhole( -uguard_inner_size.x()-uguard_inner2_size.x(), 0., 0. );
  auto uguard_solid = new G4SubtractionSolid( "UGuardSolid",
					      sub_solid, uguard_inner2_solid,
					      nullptr, uhole );
  auto uguard_lv = new G4LogicalVolume( uguard_solid, m_material_map["Iron"],
					"UGuardLV" );
  new G4PVPlacement( m_rotation_matrix,
  		     G4ThreeVector( uguard_pos ).rotateY( m_rotation_angle ),
  		     uguard_lv, "UGuardPV", m_world_lv, false, 0 );
  // Yoke
  auto yoke_inner_solid = new G4Box( "YokeInnerSolid",
				     yoke_inner_size.x(),
				     yoke_inner_size.y(),
				     yoke_inner_size.z() );
  auto yoke_outer_solid = new G4Box( "YokeOuterSolid",
				     yoke_outer_size.x(),
				     yoke_outer_size.y(),
				     yoke_outer_size.z() );
  auto yoke_solid = new G4SubtractionSolid( "YokeSolid",
					    yoke_outer_solid,
					    yoke_inner_solid );
  auto yoke_lv = new G4LogicalVolume( yoke_solid, m_material_map["Iron"],
				      "YokeLV" );
  new G4PVPlacement( m_rotation_matrix,
  		     G4ThreeVector( yoke_pos ).rotateY( m_rotation_angle ),
  		     yoke_lv, "YokePV", m_world_lv, false, 0 );
  // Downstream End Guard
  auto dguard_inner_solid = new G4Box( "DGuardInnerSolid",
				       dguard_inner_size.x(),
				       dguard_inner_size.y(),
				       dguard_inner_size.z() );
  auto dguard_outer_solid = new G4Box( "DGuardOuterSolid",
				       dguard_outer_size.x(),
				       dguard_outer_size.y(),
				       dguard_outer_size.z() );
  auto dguard_solid = new G4SubtractionSolid( "DGuardSolid",
					      dguard_outer_solid,
					      dguard_inner_solid );
  auto dguard_lv = new G4LogicalVolume( dguard_solid, m_material_map["Iron"],
					"DGuardLV" );
  new G4PVPlacement( m_rotation_matrix,
  		     G4ThreeVector( dguard_pos ).rotateY( m_rotation_angle ),
  		     dguard_lv, "DGuardPV", m_world_lv, false, 0 );
  // Virtual Plane
  auto vp_solid = new G4Box( "VPSolid", 2.*m/2, 2.*m/2, 0.001*mm/2 );
  auto vp_lv = new G4LogicalVolume( vp_solid, m_material_map["Air"], "VPLV" );
  vp_lv->SetSensitiveDetector( vp_sd );
  vp_lv->SetVisAttributes(G4Color::Blue());;
			//G4VisAttributes::GetInvisible() );
  for( G4int i=0; i<NumOfSegVP; ++i ){
    auto vp_pos = ( gGeom.GetGlobalPosition( "KURAMA" ) +
		    gGeom.GetGlobalPosition( Form( "VP%d", i+1 ) ) );
    new G4PVPlacement( m_rotation_matrix, vp_pos, vp_lv,
		       Form( "VP%dPV", i+1 ), m_world_lv, false, i );
  }
  myField->SetStatusKuramaField( true );
  myField->SetKuramaFieldMap( gConf.Get<G4String>( "KURAMAFLDMAP" ) );
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructLAC( void )
{
  auto lac_sd = new TPCLACSD("/LAC");
  lac_sd->SetRefractiveIndex( 1.05 );
  G4SDManager::GetSDMpointer()->AddNewDetector( lac_sd );
  const auto& ra2 = gGeom.GetRotAngle2("LAC") * deg;
  const auto& frame_size = gSize.GetSize("LacFrame") * 0.5 * mm;
  const auto& radiator_size = gSize.GetSize("LacRadiator") * 0.5 * mm;
  // Mother
  auto mother_solid = new G4Box( "LacMotherSolid",
				 frame_size.x() + 50.*mm,
				 frame_size.y() + 50.*mm,
				 frame_size.z() + 50.*mm );
  auto mother_lv = new G4LogicalVolume( mother_solid,
					m_material_map["Air"],
					"LacMotherLV" );
  auto rot = new G4RotationMatrix;
  rot->rotateY( - ra2 - m_rotation_angle );
  auto pos = ( gGeom.GetGlobalPosition("KURAMA") +
	       gGeom.GetGlobalPosition("LAC") );
  pos.rotateY( m_rotation_angle );
  new G4PVPlacement( rot, pos, mother_lv,
		     "LacMotherPV", m_world_lv, false, 0 );
  mother_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // Frame
  auto frame_solid = new G4Box( "LacFrameSolid", frame_size.x(),
				frame_size.y(), frame_size.z() );
  auto frame_lv = new G4LogicalVolume( frame_solid,
				       m_material_map["Air"],
				       "LacFrameLV" );
  pos.setMag( 0. );
  new G4PVPlacement( nullptr, pos, frame_lv,
		     "LacFramePV", mother_lv, false, 0 );
  // Radiator
  auto radiator_solid = new G4Box( "LacRadiatorSolid", radiator_size.x(),
				   radiator_size.y(), radiator_size.z() );
  auto radiator_lv = new G4LogicalVolume( radiator_solid,
					  m_material_map["SilicaAerogelLAC"],
					  "LacRadiatorLV" );
  radiator_lv->SetSensitiveDetector( lac_sd );
  radiator_lv->SetVisAttributes( G4Color::Magenta() );
  pos.set( 0., 0., -frame_size.z() + radiator_size.z() + 10.*mm );
  new G4PVPlacement( nullptr, pos, radiator_lv,
		     "LacRadiatorPV", frame_lv, false, 0 );
  // Mirror
  const G4double mirror_thickness = 1.*mm/2.;
  const G4double mirror_space = 20.*mm;
  const G4ThreeVector triangle_size( 710.59*mm, frame_size.y(), 500.*mm );
  const G4double mirror_angle = std::atan2( triangle_size.z(),
					    triangle_size.x() );
  const G4ThreeVector mirror1_size( ( frame_size.x() - triangle_size.x() )/2.,
				    triangle_size.y(), mirror_thickness );
  const G4ThreeVector mirror2_size( std::hypot( triangle_size.x(),
						triangle_size.z() )/2.,
				    triangle_size.y(),
				    mirror_thickness );
  auto mirror1_solid = new G4Box( "LacMirror1Solid", mirror1_size.x(),
				  mirror1_size.y(), mirror1_size.z() );
  auto mirror1_lv = new G4LogicalVolume( mirror1_solid,
					 m_material_map["Aluminum"],
					"LacMirror1LV" );
  auto mirror2_solid = new G4Box( "LacMirror2Solid", mirror2_size.x(),
				  mirror2_size.y(), mirror2_size.z() );
  auto mirror2_lv = new G4LogicalVolume( mirror2_solid,
					 m_material_map["Aluminum"],
					"LacMirror2LV" );
  for( G4int i=0; i<2; ++i ){
    pos.set( ( triangle_size.x() + mirror1_size.x() ) * ( i*2 - 1 ),
	     0., frame_size.z() - mirror_space );
    new G4PVPlacement( nullptr, pos, mirror1_lv,
		       "LacMirrorPV", frame_lv, false, 0 );
    pos.set( triangle_size.x()/2 * ( i*2 - 1 ),
	     0., frame_size.z() - triangle_size.z()/2 - mirror_space );
    rot = new G4RotationMatrix;
    rot->rotateY( mirror_angle * ( i*2 - 1 ) );
    new G4PVPlacement( rot, pos, mirror2_lv,
		       "LacMirrorPV", frame_lv, false, 0 );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructNBAR( void )
{
  G4int n_bar_use = gConf.Get<G4int>( "UseNBar" );

  ////////////////////////n_bar
  const G4int NPHI_NBAR = 32;
  G4LogicalVolume* nbarLV[NPHI_NBAR];
  G4VisAttributes* nbarVisAtt[NPHI_NBAR];

  if(n_bar_use==1.){
    G4double n_bar_radius=379.5;
    const G4int num_one_plane_nbar=NPHI_NBAR/8.;
    const G4double DX_NBAR1 = 65.0/2.*mm;
    G4double plane_width_nbar = (n_bar_radius-DX_NBAR1)*tan(22.5*deg)*2;

    const G4double DZ_NBAR1 = (n_bar_radius-DX_NBAR1)*tan(22.5*deg)*2/num_one_plane_nbar*mm/2;
    const G4double DY_NBAR1 = 500.0*mm;
    G4Box* nbarSolid1= new G4Box("SIDE NBAR1", DX_NBAR1,
				 DY_NBAR1,  DZ_NBAR1);
    // 16 NBAR
    const G4double dangleOut_nbar = 22.5*2*deg;
    G4ThreeVector posMOut1_nbar(n_bar_radius*mm,0.*mm,(-DZ_NBAR1)); //x,z,y??
    G4ThreeVector posMOut2_nbar(n_bar_radius*mm,0.*mm,(+DZ_NBAR1)); //x,z,y??
    G4RotationMatrix* rotMOutP_nbar = new G4RotationMatrix;
    rotMOutP_nbar->rotateY(dangleOut_nbar*0.5-22.5*deg);
    posMOut1_nbar.rotateY(dangleOut_nbar*0.5-22.5*deg);
    posMOut2_nbar.rotateY(dangleOut_nbar*0.5-22.5*deg);
    for(G4int i=0;i<8; i++){
      for(G4int j=0;j<num_one_plane_nbar;j++){
	G4ThreeVector posMOut1( n_bar_radius*mm, 0.*mm,
				plane_width_nbar/2-DZ_NBAR1-DZ_NBAR1*2*j); //x,z,y??
	G4RotationMatrix* rotMOutP = new G4RotationMatrix;
	G4Transform3D transformMP1(rotMOutP->rotateY(dangleOut_nbar*(i)),
				   posMOut1.rotateY(dangleOut_nbar*(i)));
	nbarLV[i*num_one_plane_nbar+j] =
	  new G4LogicalVolume(nbarSolid1,
			      m_material_map["Scintillator"],
			      Form("NBARLV%d",i*num_one_plane_nbar+j));
	/////id 0 is start forward scintillator
	G4int copyno=i*num_one_plane_nbar+j+6;
	if(copyno>31) copyno=copyno-32;
	new G4PVPlacement( transformMP1, nbarLV[i*num_one_plane_nbar+j],
			   Form("NBARPV%d",i*num_one_plane_nbar+j), m_world_lv, false, copyno );
	nbarVisAtt[i*num_one_plane_nbar+j]= new G4VisAttributes(true, G4Colour(0.,0.8,1.));
	nbarLV[i*num_one_plane_nbar+j]->SetVisAttributes(nbarVisAtt[i*num_one_plane_nbar+j]);
      }
    }
  }

  TPCNBARSD* nbarSD = new TPCNBARSD("/NBAR");
  G4SDManager::GetSDMpointer()->AddNewDetector( nbarSD );
  for(G4int i = 0; i<32; i++){
    nbarLV[i]->SetSensitiveDetector(nbarSD);
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructPVAC2( void )
{
  G4int ac_use = gConf.Get<G4int>( "UseAC" );
  const G4int NPHI_PVC = 32;
  G4LogicalVolume* pvcLV[NPHI_PVC];
  G4VisAttributes* pvcVisAtt[NPHI_PVC];
  if(ac_use==1.){
    const G4double DX_PVC1 = 10.0/2.*mm;
    const G4double DZ_PVC1 = (340.-DX_PVC1)*tan(22.5*deg)*mm/2;
    const G4double DY_PVC1 = 440.0*mm;
    // const G4double DX_PVC2 = 10.0/2.*mm;
    const G4double DZ_PVC2 = (340.-DX_PVC1)*tan(22.5*deg)*mm/2;
    // const G4double DY_PVC2 = 440.0*mm;
    G4Box* pvcSolid1= new G4Box("SIDE PVC1", DX_PVC1, DY_PVC1,  DZ_PVC1);
    // G4Box* pvcSolid2= new G4Box("SIDE PVC2", DX_PVC2,
    //                                    DY_PVC2,  DZ_PVC2);

    // 16 PVC
    const G4double dangleOut_pvc = 22.5*2*deg;
    G4ThreeVector posMOut1_pvc(340.*mm,0.*mm,(-DZ_PVC2)); //x,z,y??
    G4ThreeVector posMOut2_pvc(340.*mm,0.*mm,(+DZ_PVC2)); //x,z,y??
    G4RotationMatrix* rotMOutP_pvc = new G4RotationMatrix;
    rotMOutP_pvc->rotateY(dangleOut_pvc*0.5-22.5*deg);
    posMOut1_pvc.rotateY(dangleOut_pvc*0.5-22.5*deg);
    posMOut2_pvc.rotateY(dangleOut_pvc*0.5-22.5*deg);
    for(G4int k=0;k<8; k++){
      G4Transform3D transformMP1_pvc(*rotMOutP_pvc, posMOut1_pvc);
      pvcLV[k*2] = new G4LogicalVolume(pvcSolid1, m_material_map["Scintillator"], Form("PVCLV%d",k*2));
      new G4PVPlacement( transformMP1_pvc, pvcLV[k*2], Form("PVCPV%d",k*2), m_world_lv, false, k*2 );
      G4Transform3D transformMP2_pvc(*rotMOutP_pvc, posMOut2_pvc);
      pvcLV[k*2+1] = new G4LogicalVolume(pvcSolid1, m_material_map["Scintillator"], Form("pvcLV%d",k*2+1));
      new G4PVPlacement( transformMP2_pvc, pvcLV[k*2], Form("pvcPV%d",k*2+1),
			 m_world_lv, false, k*2 );
      pvcVisAtt[k*2]= new G4VisAttributes(true, G4Colour(1.,0.0,1.));
      pvcLV[k*2]->SetVisAttributes(pvcVisAtt[k*2]);
      pvcVisAtt[k*2+1]= new G4VisAttributes(true, G4Colour(1.,0.0,1.));
      pvcLV[k*2+1]->SetVisAttributes(pvcVisAtt[k*2+1]);
      rotMOutP_pvc->rotateY(dangleOut_pvc);
      posMOut1_pvc.rotateY(dangleOut_pvc);
      posMOut2_pvc.rotateY(dangleOut_pvc);
    }
  }

  TPCACSD* acSD= new TPCACSD("/AC");
  G4SDManager::GetSDMpointer()->AddNewDetector(acSD);
  for(G4int i = 0; i<16; i++){
    pvcLV[i]->SetSensitiveDetector(acSD);
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructSCH( void )
{
  const auto& sch_pos = ( gGeom.GetGlobalPosition("KURAMA") +
			  gGeom.GetGlobalPosition("SCH") );
  const auto& sch_size = gSize.GetSize( "SchSeg" ) * 0.5 * mm;
  const G4double dXdW = gGeom.GetWirePitch( "SCH" ) * mm;
  G4LogicalVolume* sch_lv[NumOfSegSCH];
  auto sch_solid = new G4Box( "SchSolid", sch_size.x(),
			      sch_size.y(), sch_size.z() );
  for( G4int i=0; i<NumOfSegSCH; ++i ){
    sch_lv[i] = new G4LogicalVolume( sch_solid, m_material_map["Scintillator"],
				     Form( "SchSeg%dLV", i ), 0, 0, 0 );
    sch_lv[i]->SetVisAttributes( G4Colour::Green() );
    G4double ipos_x = dXdW * ( i - ( NumOfSegSCH - 1 )/2. );
    G4ThreeVector pos( sch_pos.x() + ipos_x,
		       sch_pos.y(),
		       sch_pos.z() + 1.*mm*( 1- 2*(i%2) ) );
    // pos.rotateY( m_rotation_angle );
    new G4PVPlacement( m_rotation_matrix, pos, sch_lv[i],
		       Form( "SchSeg%dPV", i ), m_world_lv, false, i );
  }
  auto schSD = new TPCSCHSD("/SCH");
  G4SDManager::GetSDMpointer()->AddNewDetector( schSD );
  for( G4int i = 0; i<NumOfSegSCH; ++i ){
    sch_lv[i]->SetSensitiveDetector( schSD );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructSDC1( void )
{
  if( !m_sdc_sd ){
    m_sdc_sd = new TPCSDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector( m_sdc_sd );
  }

  const auto& sdc1_pos = ( gGeom.GetGlobalPosition("KURAMA") +
			   ( gGeom.GetGlobalPosition("SDC1-V1") +
			     gGeom.GetGlobalPosition("SDC1-U2") ) * 0.5 );
  const auto& frame_size = gSize.GetSize( "Sdc1Frame" ) * 0.5 * mm;
  const auto& drift_size = gSize.GetSize( "Sdc1Drift" ) * 0.5 * mm;
	auto sdc1_solid = new G4Box( "Sdc1Solid", frame_size.x(),
			       frame_size.y(), frame_size.z() );
  auto sdc1_lv = new G4LogicalVolume( sdc1_solid, m_material_map["Argon"],
				      "Sdc1LV", 0, 0, 0 );
  sdc1_lv->SetVisAttributes( G4Colour::Green() );
  new G4PVPlacement( m_rotation_matrix, sdc1_pos,
		     sdc1_lv, "Sdc1PV", m_world_lv, false, 0 );
  auto sdc1pl_solid = new G4Box( "Sdc1PlSolid", drift_size.x(),
				 drift_size.y(), drift_size.z() );
  G4String plane_name[] = { "Sdc1V1", "Sdc1V2", "Sdc1X1",
			    "Sdc1X2", "Sdc1U1", "Sdc1U2" };
  for( G4int i=0; i<NumOfLayersSDC1; ++i ){
    G4ThreeVector pos;
    G4RotationMatrix* rot = new G4RotationMatrix;
		switch (i) {
    case 0:
      pos.setZ( -22.5985*mm );
//			rot->rotateZ(15.*deg);
      break;
    case 1:
      pos.setZ( -17.4015*mm );
//			rot->rotateZ(15.*deg);
      break;
    case 2:
      pos.setZ( -2.5985*mm );
      break;
    case 3:
      pos.setZ( 2.5985*mm );
      break;
    case 4:
      pos.setZ( 17.4015*mm );
//			rot->rotateZ(-15.*deg);
      break;
    case 5:
      pos.setZ( 22.5985*mm );
//			rot->rotateZ(-15.*deg);
      break;
    }
    auto sdc1pl_lv = new G4LogicalVolume( sdc1pl_solid,
					  m_material_map["Argon"],
					  plane_name[i] + "LV", 0, 0, 0 );
    sdc1pl_lv->SetSensitiveDetector( m_sdc_sd );

    new G4PVPlacement( rot, pos, sdc1pl_lv, plane_name[i] + "PV",
		       sdc1_lv, false, 101+i );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructSDC2( void )
{
  if( !m_sdc_sd ){
    m_sdc_sd = new TPCSDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector( m_sdc_sd );
  }

  const auto& sdc2_pos = ( gGeom.GetGlobalPosition("KURAMA") +
			   ( gGeom.GetGlobalPosition("SDC2-X1") +
			     gGeom.GetGlobalPosition("SDC2-Y2") ) * 0.5 );
  const auto& frame_size = gSize.GetSize( "Sdc2Frame" ) * 0.5 * mm;
  const auto& drift_size = gSize.GetSize( "Sdc2Drift" ) * 0.5 * mm;
  auto sdc2_solid = new G4Box( "Sdc2Solid", frame_size.x(),
			       frame_size.y(), frame_size.z() );
  auto sdc2_lv = new G4LogicalVolume( sdc2_solid, m_material_map["Argon"],
				      "Sdc2LV", 0, 0, 0 );
  sdc2_lv->SetVisAttributes( G4Colour::Green() );
  new G4PVPlacement( m_rotation_matrix, sdc2_pos,
		     sdc2_lv, "Sdc2PV", m_world_lv, false, 0 );
  auto sdc2pl_solid = new G4Box( "Sdc2PlSolid", drift_size.x(),
				 drift_size.y(), drift_size.z() );
  G4String plane_name[] = { "Sdc2X1", "Sdc2X2", "Sdc2Y1", "Sdc2Y2" };
  for( G4int i=0; i<NumOfLayersSDC2; ++i ){
    G4ThreeVector pos;
    switch (i) {
    case 0:
      pos.setZ( -23.76*mm );
      break;
    case 1:
      pos.setZ( -13.76*mm );
      break;
    case 2:
      pos.setZ( 13.76*mm );
      break;
    case 3:
      pos.setZ( 23.76*mm );
      break;
    }
    auto sdc2pl_lv = new G4LogicalVolume( sdc2pl_solid,
					  m_material_map["Argon"],
					  plane_name[i] + "LV", 0, 0, 0 );
    sdc2pl_lv->SetSensitiveDetector( m_sdc_sd );
    new G4PVPlacement( nullptr, pos, sdc2pl_lv, plane_name[i] + "PV",
		       sdc2_lv, false, 201+i );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructSDC3( void )
{
  if( !m_sdc_sd ){
    m_sdc_sd = new TPCSDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector( m_sdc_sd );
  }

  const auto& sdc3_pos = ( gGeom.GetGlobalPosition("KURAMA") +
			   ( gGeom.GetGlobalPosition("SDC3-X1") +
			     gGeom.GetGlobalPosition("SDC3-Y2") ) * 0.5 );
  const auto& frame_size = gSize.GetSize( "Sdc3Frame" ) * 0.5 * mm;
  const auto& drift_size = gSize.GetSize( "Sdc3Drift" ) * 0.5 * mm;
  auto sdc3_solid = new G4Box( "Sdc3Solid", frame_size.x(),
			       frame_size.y(), frame_size.z() );
  auto sdc3_lv = new G4LogicalVolume( sdc3_solid, m_material_map["Argon"],
				      "Sdc3LV", 0, 0, 0 );
  sdc3_lv->SetVisAttributes( G4Colour::Green() );
  new G4PVPlacement( m_rotation_matrix, sdc3_pos,
		     sdc3_lv, "Sdc3PV", m_world_lv, false, 0 );
  auto sdc3pl_solid = new G4Box( "Sdc3PlSolid", drift_size.x(),
				 drift_size.y(), drift_size.z() );
  G4String plane_name[] = { "Sdc3X1", "Sdc3X2", "Sdc3Y1", "Sdc3Y2" };
  for( G4int i=0; i<NumOfLayersSDC3; ++i ){
    G4ThreeVector pos;
    switch (i) {
    case 0:
      pos.setZ( -16.0*mm );
      break;
    case 1:
      pos.setZ( -8.206*mm );
      break;
    case 2:
      pos.setZ( 8.206*mm );
      break;
    case 3:
      pos.setZ( 16.0*mm );
      break;
    }
    auto sdc3pl_lv = new G4LogicalVolume( sdc3pl_solid,
					  m_material_map["Argon"],
					  plane_name[i] + "LV", 0, 0, 0 );
    sdc3pl_lv->SetSensitiveDetector( m_sdc_sd );
    new G4PVPlacement( nullptr, pos, sdc3pl_lv, plane_name[i] + "PV",
		       sdc3_lv, false, 301+i );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructSDC4( void )
{
  if( !m_sdc_sd ){
    m_sdc_sd = new TPCSDCSD("/SDC");
    G4SDManager::GetSDMpointer()->AddNewDetector( m_sdc_sd );
  }

  const auto& sdc4_pos = ( gGeom.GetGlobalPosition("KURAMA") +
			   ( gGeom.GetGlobalPosition("SDC4-Y1") +
			     gGeom.GetGlobalPosition("SDC4-X2") ) * 0.5 );
  const auto& frame_size = gSize.GetSize( "Sdc4Frame" ) * 0.5 * mm;
  const auto& drift_size = gSize.GetSize( "Sdc4Drift" ) * 0.5 * mm;
  auto sdc4_solid = new G4Box( "Sdc4Solid", frame_size.x(),
			       frame_size.y(), frame_size.z() );
  auto sdc4_lv = new G4LogicalVolume( sdc4_solid, m_material_map["Argon"],
				      "Sdc4LV", 0, 0, 0 );
  sdc4_lv->SetVisAttributes( G4Colour::Green() );
  new G4PVPlacement( m_rotation_matrix, sdc4_pos,
		     sdc4_lv, "Sdc4PV", m_world_lv, false, 0 );
  auto sdc4pl_solid = new G4Box( "Sdc4PlSolid", drift_size.x(),
				 drift_size.y(), drift_size.z() );
  G4String plane_name[] = { "Sdc4Y1", "Sdc4Y2", "Sdc4X1", "Sdc4X2" };
  for( G4int i=0; i<NumOfLayersSDC4; ++i ){
    G4ThreeVector pos;
    switch (i) {
    case 0:
      pos.setZ( -34.868*mm );
      break;
    case 1:
      pos.setZ( -17.547*mm );
      break;
    case 2:
      pos.setZ( 17.547*mm );
      break;
    case 3:
      pos.setZ( 34.868*mm );
      break;
    }
    auto sdc4pl_lv = new G4LogicalVolume( sdc4pl_solid,
					  m_material_map["Argon"],
					  plane_name[i] + "LV", 0, 0, 0 );
    sdc4pl_lv->SetSensitiveDetector( m_sdc_sd );
    new G4PVPlacement( nullptr, pos, sdc4pl_lv, plane_name[i] + "PV",
		       sdc4_lv, false, 401+i );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructShsMagnet( void )
{
#if 0 // Old version
  const auto& tpc_pos = gGeom.GetGlobalPosition("HypTPC");
  const G4double DPHI_TPC = 360.*deg;
  auto tube_solid = new G4Tubs( "TubeSolid", 45.*cm, 80.*cm, 135./2.*cm, 0.*deg, DPHI_TPC);
  auto hole_solid = new G4Box( "HoleSolid", 100./2.*cm, 100./2.*cm, 60./2.*cm );
  G4ThreeVector HTrans(0, -60.*cm, 0);
  G4RotationMatrix* yRot = new G4RotationMatrix;  // Rotates X and Z axes only
  G4Transform3D transform( *yRot, HTrans );
  auto coil_solid = new G4SubtractionSolid( "ShsMagnetCoilSolid",
					    tube_solid, hole_solid, transform );
  auto coil_lv = new G4LogicalVolume( coil_solid, m_material_map["Iron"],
				      "ShsMagnetCoilLV" );
  G4RotationMatrix *rotHelm = new G4RotationMatrix();
  rotHelm->rotateX( 90.*deg );
  rotHelm->rotateZ( - m_rotation_angle );
  new G4PVPlacement( rotHelm, tpc_pos, coil_lv, "ShsMagnetCoilPV",
  		     m_world_lv, false, 0 );
  coil_lv->SetVisAttributes( PINK );
#else // Current version
  const G4ThreeVector yoke_size( 1550./2.*mm, 950./2.*mm, 1200./2.*mm );
  // auto shs_solid = new G4Box( "ShsMagnetSolid", yoke_size.x(),
  // 			      yoke_size.y(), yoke_size.z() );
  // auto shs_lv = new G4LogicalVolume( shs_solid, m_material_map["Air"],
  // 				     "ShsMagnetLV" );
  // shs_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // auto shs_pv = new G4PVPlacement( nullptr, tpc_pos, shs_lv,
  // 				   "ShsMagnetPV", m_world_lv, false, 0 );
  const G4int NumOfParams = 4;
  G4double yoke_width[NumOfParams] = { 1550*mm, 1530*mm, 1480*mm, 1470*mm };
  G4double yoke_depth[NumOfParams] = { 1200*mm, 1180*mm, 1140*mm, 1130*mm };
  G4double yoke_height[NumOfParams] = { 950*mm, 890*mm, 732*mm, 712*mm };
  G4double coner_size[NumOfParams] = { 1449.569*mm, 1429.569*mm,
				       1357.645*mm, 1347.645*mm };
  G4double chamber_rad[NumOfParams] = { 800./2.*mm, 820./2.*mm,
					850./2.*mm, 860./2.*mm };
  G4double hole_gap[NumOfParams] = { 300*mm, 320*mm, 348*mm, 358*mm };
  G4double side_gap[NumOfParams];
  for( G4int i=0 ; i<NumOfParams; ++i ){
    side_gap[i] = yoke_depth[i]/2.0*std::sqrt(2.)*mm;
  }
  G4double dummy_size[3] = { 2000*mm, 2000*mm, 2000*mm };
  G4double hole_pos[NumOfParams];
  for( G4int i=0; i<NumOfParams; ++i ){
    hole_pos[i] = yoke_depth[i]/2.;
  }
  G4Box* box_solid[NumOfParams];
  G4Tubs* tube_solid[NumOfParams];
  G4Box* box_outer_solid[NumOfParams];
  G4Box* box_inner_solid[NumOfParams];
  G4Box* box_side_solid[NumOfParams];
  G4SubtractionSolid* octa_solid[NumOfParams];
  G4SubtractionSolid* corner_solid[NumOfParams];
  G4SubtractionSolid* side1_solid[NumOfParams];
  G4SubtractionSolid* side2_solid[NumOfParams];
  G4SubtractionSolid* chamber_solid[NumOfParams];
  G4RotationMatrix* rot_magnet  = new G4RotationMatrix;
  rot_magnet->rotateZ( 45.*deg );
  for( G4int i=0; i<NumOfParams; ++i ){
    box_solid[i] = new G4Box( Form( "Box%dSolid", i ), yoke_width[i]/2.,
			      yoke_depth[i]/2., yoke_height[i]/2. );
    box_outer_solid[i] = new G4Box( Form( "BoxOuter%dSolid", i ),
				    dummy_size[0], dummy_size[1],
				    dummy_size[2] );
    box_inner_solid[i] = new G4Box( Form( "BoxInner%dSolid", i ),
				    coner_size[i], coner_size[i],
				    dummy_size[2] );
    box_side_solid[i] = new G4Box( Form("BoxSide%dSolid", i ),
				   side_gap[i]/2., side_gap[i]/2.,
				   hole_gap[i]/2. );
    tube_solid[i] = new G4Tubs( Form("Tube%dSolid", i ), 0., chamber_rad[i],
				2000*mm, 0.*deg, 360.*deg );
    corner_solid[i] = new G4SubtractionSolid( Form( "Corner%dSolid", i ),
					      box_outer_solid[i],
					      box_inner_solid[i],
					      nullptr, G4ThreeVector() );
    octa_solid[i] = new G4SubtractionSolid( Form( "Octa%dSolid", i ),
					    box_solid[i], corner_solid[i],
					    rot_magnet, G4ThreeVector() );
    G4ThreeVector side1_size( 0., hole_pos[i], 0. );
    side1_solid[i] = new G4SubtractionSolid( Form( "Side1-%dSolid", i ),
					     octa_solid[i], box_side_solid[i],
					     rot_magnet, side1_size );
    G4ThreeVector side2_size( 0., -hole_pos[i], 0. );
    side2_solid[i] = new G4SubtractionSolid( Form("Side2-%dSolid", i ),
					     side1_solid[i], box_side_solid[i],
					     rot_magnet, side2_size );
    chamber_solid[i] = new G4SubtractionSolid( Form("Chamber%dSolid", i ),
					       side2_solid[i], tube_solid[i],
					       nullptr, G4ThreeVector() );
  }
  auto frame1_solid = new G4SubtractionSolid( "Frame1Solid", chamber_solid[0],
					      chamber_solid[1], nullptr,
					      G4ThreeVector() );
  auto frame2_solid = new G4SubtractionSolid( "Frame2Solid", chamber_solid[2],
					      chamber_solid[3], nullptr,
					      G4ThreeVector() );
  auto magnet_solid = new G4UnionSolid( "ShsMagnetSolid", frame1_solid,
					frame2_solid, nullptr,
					G4ThreeVector() );
  auto magnet_lv = new G4LogicalVolume( magnet_solid,
					m_material_map["Iron"],
					"ShsMagnetLV" );
  // magnet_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );
  magnet_lv->SetVisAttributes( PINK );
  G4RotationMatrix rot_frame;
  rot_frame.rotateX( 90.*deg );
  new G4PVPlacement( G4Transform3D( rot_frame, G4ThreeVector() ),
		      magnet_lv, "ShsMagnetPV", m_world_lv, false, 0 );
  // Coil Support
  G4double CoilSupPos_height = 250*mm;
  G4double RadIn = 445*mm;
  G4double RadOut= 545*mm;
  G4double GapRadIn = 465*mm;
  G4double GapRadOut = 545*mm;
  G4double SupHeight = 112*mm;
  G4double GapHeight = 62*mm;
  std::string fullNameCoilSup = "SCCoilSup";
  auto solidTube_Sup = new G4Tubs( "CoilSupMainSolid", RadIn, RadOut, SupHeight/2.,
				   0*deg, 360*deg );
  auto solidTube_Sub = new G4Tubs( "CoilSupSubSolid", GapRadIn, GapRadOut+20*mm,
				   GapHeight/2., 0*deg, 360*deg );
  auto solidDetectorCS = new G4SubtractionSolid( "CoilSupSolid", solidTube_Sup,
						 solidTube_Sub, nullptr,
						 G4ThreeVector() );
  auto logicDetectorCS = new G4LogicalVolume( solidDetectorCS,
					      m_material_map["Iron"],
					      "CoilSupLV" );
  logicDetectorCS->SetVisAttributes( G4Color::Green() );
  G4RotationMatrix rot_sup;
  rot_sup.rotateX( 90.*deg );
  new G4PVPlacement( G4Transform3D( rot_sup, G4ThreeVector( 0, CoilSupPos_height, 0 ) ),
		     logicDetectorCS, "CoilSupUpPV", m_world_lv, false, 0 );
  new G4PVPlacement( G4Transform3D( rot_sup, G4ThreeVector( 0, -CoilSupPos_height, 0 ) ),
		     logicDetectorCS, "CoilSupDwPV", m_world_lv, false, 0 );
  // Coil
  const G4double coilRad_in = 466*mm;
  const G4double coilRad_out = 535*mm;
  const G4double coilHeight = 60*mm;
  const G4ThreeVector coilu_pos( 0*mm, 250*mm, 0*mm );
  const G4ThreeVector coild_pos( 0*mm, -250*mm, 0*mm );
  auto solidDetectorCoil = new G4Tubs( "CoilSolid", coilRad_in, coilRad_out,
				       coilHeight/2., 0*deg, 360*deg );
  auto logicDetectorCoil = new G4LogicalVolume( solidDetectorCoil,
						m_material_map["Copper"],
						"CoilLV" );
  logicDetectorCoil->SetVisAttributes( G4Color::Yellow() );
  G4RotationMatrix rot_coil;
  rot_coil.rotateX( 90.*deg );
  new G4PVPlacement( G4Transform3D( rot_coil, coilu_pos ),
		     logicDetectorCoil, "CoilUpPV", m_world_lv, false, 0 );
  new G4PVPlacement( G4Transform3D( rot_coil, coild_pos ),
		     logicDetectorCoil, "CoilDwPV", m_world_lv, false, 0 );
#endif
  myField->SetStatusShsField( true );
  myField->SetShsFieldMap( gConf.Get<G4String>( "SHSFLDMAP" ) );
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructWC( void )
{
  const auto& ra2 = gGeom.GetRotAngle2("WC") * deg;
  const auto& half_size_In = gSize.GetSize("WcSegIn") * 0.5 * mm;
  const auto& half_size_Out = gSize.GetSize("WcSegOut") * 0.5 * mm;
  const G4double pitch = gGeom.GetWirePitch("WC");
  auto wcSD = new TPCWCSD("/WC");
  wcSD->SetRefractiveIndex( 1.33 );
  G4SDManager::GetSDMpointer()->AddNewDetector( wcSD );
  // Mother
  auto mother_solid = new G4Box( "WcMotherSolid",
				 half_size_Out.x()*NumOfSegWC + 200.*mm,
				 half_size_Out.y() + 200.*mm,
				 half_size_Out.z()*2 + 200.*mm );
				 // half_size_Out.x()*NumOfSegWC + 50.*mm,
				 // half_size_Out.y() + 50.*mm,
				 // half_size_Out.z()*2 + 50.*mm );

  auto mother_lv = new G4LogicalVolume( mother_solid,
					m_material_map["Air"],
					"WcMotherLV" );
  auto rot = new G4RotationMatrix;
  rot->rotateY( - ra2 - m_rotation_angle );
  auto pos = ( gGeom.GetGlobalPosition("KURAMA") +
	       gGeom.GetGlobalPosition("WC") );
  pos.rotateY( m_rotation_angle );
  new G4PVPlacement( rot, pos, mother_lv,
		     "WcMotherPV", m_world_lv, false, 0 );
  mother_lv->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // Segment
  auto segment_solid = new G4Box( "WcSegmentSolid", half_size_In.x(),
				  half_size_In.y(), half_size_In.z() );
  auto segment_lv = new G4LogicalVolume( segment_solid,
					 m_material_map["Water"],
					 "WcSegmentLV" );

  auto WCContainer     = new G4Box("WCContainer",
				   half_size_Out.x(),
				   half_size_Out.y(),
				   half_size_Out.z());
  auto WCContainer_gap = new G4Box("WCContainer_gap",
				   half_size_In.x(),
				   half_size_In.y(),
				   half_size_In.z());
  //G4RotationMatrix* rot_wccontainer_gap;
  auto rot_wccontainer_gap = new G4RotationMatrix;
  G4ThreeVector pos_wccontainer_gap(0.0, 0.0 ,0.0);
  auto solid_WCContainer
    = new G4SubtractionSolid("solid_WCContainer",
   			     WCContainer, WCContainer_gap,
   			     rot_wccontainer_gap,
			     pos_wccontainer_gap);

  auto logWCContainer = new G4LogicalVolume(solid_WCContainer,
					    m_material_map["Acrylic"],
					    "logWCContainer");

  for( G4int i=0; i<NumOfSegWC; ++i ){


    pos = G4ThreeVector( ( -NumOfSegWC/2 + i )*pitch,
			 0.0,
			 2.*( - i%2 + 0.5 )*half_size_Out.z() );
    //for Vessel
    //    logWCContainer->SetVisAttributes( G4Colour::White() );
    logWCContainer->SetVisAttributes( G4Colour::Cyan() );
    new G4PVPlacement( nullptr, pos, logWCContainer,
		       "WcSegmentContainerPV", mother_lv, false, i );
    //for Water
    segment_lv->SetVisAttributes( G4Colour::Cyan() );
    segment_lv->SetSensitiveDetector( wcSD );
    new G4PVPlacement( nullptr, pos, segment_lv,
		       "WcSegmentPV", mother_lv, false, i );

  }
 /*
  //Temporary!!!!!!
  //WC frame (test)
  //double Alcut_z =0.*mm;
  double Alcut_z =(1000.-(250.+403.))*mm;

  auto Alframe1     = new G4Box("Alframe1",
				80.*mm/2.,
				80.*mm/2.,
				2000.*mm/2. - Alcut_z/2.);

  auto Alframe2     = new G4Box("Alframe2",
				80.*mm/2.,
				80.*mm/2.,
				750.*mm/2.);



  auto Alframe1_lv = new G4LogicalVolume( Alframe1,
					  m_material_map["Aluminum"],
					  "Alframe1_LV" );
  auto Alframe2_lv = new G4LogicalVolume( Alframe2,
					  m_material_map["Aluminum"],
					  "Alframe2_LV" );

  auto pos_frame1_1 = G4ThreeVector( -pitch/2.-3000.*mm/2.,
				     -2000.*mm,
				     -1.*Alcut_z/2.);

  auto pos_frame1_2 = G4ThreeVector( -pitch/2.+3000.*mm/2.,
				     -2000.*mm,
				     -1.*Alcut_z/2.);

  Alframe1_lv->SetVisAttributes( G4Colour::White() );
  new G4PVPlacement( nullptr, pos_frame1_1, Alframe1_lv,
		     "WcFramePV", mother_lv, false, 0 );

  new G4PVPlacement( nullptr, pos_frame1_2, Alframe1_lv,
		     "WcFramePV", mother_lv, false, 1 );


  auto pos_frame2_1 = G4ThreeVector( -pitch/2.-3000.*mm/2.,
				     -1000.*mm,
				     750/2.);

  auto pos_frame2_2 = G4ThreeVector( -pitch/2.+3000.*mm/2.,
				     -1000.*mm,
				     750/2.);

  Alframe2_lv->SetVisAttributes( G4Colour::Red() );
  new G4PVPlacement( nullptr, pos_frame2_1, Alframe2_lv,
		     "WcFramePV", mother_lv, false, 3 );

  new G4PVPlacement( nullptr, pos_frame2_2, Alframe2_lv,
		     "WcFramePV", mother_lv, false, 4 );

  //Temporary!!!!!!
*/


}

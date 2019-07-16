// -*- C++ -*-

#include "TPCDetectorConstruction.hh"

#include <G4Box.hh>
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
#include <G4VisAttributes.hh>
#include <tools/mathd>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "TPCACSD.hh"
#include "TPCDCSD.hh"
#include "TPCField.hh"
#include "TPCFTOFSD.hh"
#include "TPCNBARSD.hh"
#include "TPCHTOFSD.hh"
#include "TPCPadSD.hh"
#include "TPCSCHSD.hh"
#include "TPCScintSD.hh"
#include "TPCTargetSD.hh"

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
  // color
  const G4Colour WHITE( 1., 1., 1. );
  const G4Colour BLACK( 0., 0., 0. );
  const G4Colour RED( 1.0, 0.0, 0.0 );
  const G4Colour GREEN( 0.0, 1.0, 0.0 );
  const G4Colour BLUE( 0.0, 0.0, 1.0 );
  const G4Colour CYAN( 0.0, 1.0, 1.0 );
  const G4Colour AQUA( 0.247, 0.8, 1.0 );
  const G4Colour MAGENTA( 1.0, 0.0, 1.0 );
  const G4Colour YELLOW( 1.0, 1.0, 0.0 );
  const G4Colour ORANGE( 1.0, 0.55, 0.0 );
  const G4Colour GRAY( 0.5, 0.5, 0.5 );
  const G4Colour LAVENDER( 0.901, 0.901, 0.98 );
  const G4Colour MAROON( 0.5, 0.0, 0.0 );
  const G4Colour PINK( 1.0, 0.753, 0.796 );
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
    m_dc_sd()
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
  auto world_solid = new G4Box( "WorldSolid", 8.*m/2, 6.*m/2, 16.*m/2 );
  m_world_lv = new G4LogicalVolume( world_solid, m_material_map["Air"],
				    "WorldLV" );
  auto world_vis_attr = new G4VisAttributes( false, BLACK );
  m_world_lv->SetVisAttributes( world_vis_attr );
  auto world_pv = new G4PVPlacement( nullptr, G4ThreeVector(), m_world_lv,
				     "WorldPV", nullptr, false, 0 );
  ConstructHypTPC();
  ConstructHTOF();
  ConstructTarget();
  // ConstructPVAC2();
  // ConstructNBAR();
  ConstructShsMagnet();
  if( gConf.Get<G4int>("ConstructKurama") ){
    ConstructKuramaMagnet();
    ConstructSDC1();
    ConstructSCH();
    ConstructSDC2();
    ConstructSDC3();
    ConstructFTOF();
  }
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
  name = "Argon";
  m_element_map[name] = new G4Element( name, symbol="Ar", Z=18.,
				       A=39.948 *g/mole );
  name = "Silicon";
  m_element_map[name] = new G4Element( name, symbol="Si", Z=14.,
				       A=28.0855 *g/mole );
  name = "Iodine";
  m_element_map[name] = new G4Element( name, symbol="I",  Z=53.,
				       A=126.90447 *g/mole );
  name = "Cesium";
  m_element_map[name] = new G4Element( name, symbol="Cs", Z=55.,
				       A=132.90543 *g/mole );
  name = "Sodium";
  m_element_map[name] = new G4Element( name, symbol="Na", Z=11.,
				       A=22.989768 *g/mole );
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructMaterials( void )
{
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
  // CH2 Polyethelene
  name = "CH2";
  m_material_map[name] = new G4Material( name, density=0.95*g/cm3, nel=2 );
  m_material_map[name]->AddElement( m_element_map["Carbon"], natoms=1 );
  m_material_map[name]->AddElement( m_element_map["Hydrogen"], natoms=2 );
  // Quartz (SiO2, crystalline)
  name = "Quartz";
  m_material_map[name] = new G4Material( name, density=2.64 *g/cm3, nel=2 );
  m_material_map[name]->AddElement( m_element_map["Silicon"], natoms=1 );
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
  } else {
    std::string e(FUNC_NAME + " No target material : " + target_material );
    throw std::invalid_argument( e );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructFTOF( void )
{
  const auto& ra2 = gGeom.GetRotAngle2("TOF") * deg;
  const auto& half_size = gSize.GetSize("FtofSeg") * 0.5 * mm;
  const G4double pitch = gGeom.GetWirePitch("TOF");
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
  mother_lv->SetVisAttributes( new G4VisAttributes( false, RED ) );
  // Segment
  auto segment_solid = new G4Box( "FtofSegmentSolid", half_size.x(),
				  half_size.y(), half_size.z() );
  auto segment_lv = new G4LogicalVolume( segment_solid,
					 m_material_map["Scintillator"],
					 "FtofSegmentLV" );
  for( G4int i=0; i<NumOfSegTOF; ++i ){
    segment_lv->SetVisAttributes( AQUA );
    // segment_lv->SetUserLimits( new G4UserLimits( maxStep ) );
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
  const auto& htof_pos = gGeom.GetGlobalPosition("HTOF");
  const auto& half_size = gSize.GetSize("HtofSeg") * 0.5 * mm;
  const G4double L = gGeom.GetLocalZ("HTOF");
  const G4double dXdW = gGeom.GetWirePitch("HTOF");
  auto htof_solid = new G4Box( "HtofSolid", half_size.x(),
			       half_size.y(), half_size.z() );
  auto htof_lv = new G4LogicalVolume( htof_solid, m_material_map["Scintillator"],
				      "HtofLV" );
  htof_lv->SetSensitiveDetector( htof_sd );
  htof_lv->SetVisAttributes( CYAN );
  for( G4int i=0; i<NumOfPlaneHTOF; ++i ){
    for( G4int j=0; j<NumOfSegHTOFOnePlane; ++j ){
      G4int seg = i*NumOfSegHTOFOnePlane + j;
      // ( lateral, height, radial )
      G4ThreeVector seg_pos( dXdW * ( j - ( NumOfSegHTOFOnePlane - 1 )/2. ),
			     0.*mm,
			     -L );
      auto rotMOutP = new G4RotationMatrix;
      rotMOutP->rotateY( - i * 360./NumOfPlaneHTOF*deg );
      seg_pos.rotateY( i * 360./NumOfPlaneHTOF*deg );
      seg_pos += htof_pos;
      G4int copy_no = seg + 6;
      if( copy_no > 31 )
	copy_no -= 32;
      new G4PVPlacement( rotMOutP, seg_pos, htof_lv, Form("HtofPV%d", seg),
			 m_world_lv, false, copy_no );
    }
  }

  // ==============================================================
  // light guide for Scintillators
  // ==============================================================
  /*
    const G4int NTOF_LG = 16;
    G4LogicalVolume* TOFLGLV;
    G4VPhysicalVolume* TOFLGPV[NTOF_LG];
    G4VisAttributes* TOFLGVisAtt;

    // ==============================================================
    // Side TOFLGillators (Outer)
    // ==============================================================
    //thickness 5mm
    //  const G4double R_TOFLG =  ROUT_TPC*0.5*sqrt(3.0)+2.5*1.0;
    const G4double Angle_TOFLG = 15.0*deg;
    const G4double DX_TOFLG1 = cos(15.*deg)*5.*mm;
    const G4double DZ_TOFLG1 = (320.-DX_TOFLG1)*tan(22.5*deg)/2*mm;
    const G4double DY_TOFLG1 = 100.0*mm;

    G4Trd* TOFLGSolid1= new G4Trd("SIDE TOFLG1", DX_TOFLG1, 30.*mm,
    DY_TOFLG1, DY_TOFLG1, DZ_TOFLG1);
    TOFLGLV = new G4LogicalVolume(htof_solid, Scinti, name1);
    // 16 side TOFLG.

    const G4double dangleTOFLG = 22.5*2*deg;
    G4ThreeVector posTOFLG1(520.*mm,0.*mm,(-DZ_TOFLG1)); //x,z,y??
    G4ThreeVector posTOFLG2(520.*mm,0.*mm,(+DZ_TOFLG1)); //x,z,y??
    G4RotationMatrix* rotTOFLGP = new G4RotationMatrix;
    rotTOFLGP->rotateY(dangleTOFLG*0.5-22.5*deg);
    posTOFLG1.rotateY(dangleTOFLG*0.5-22.5*deg);
    posTOFLG2.rotateY(dangleTOFLG*0.5-22.5*deg);

    //  for(G4int k=0;k<NPHI_TOFLG; k++){
    for(G4int k=0;k<8; k++){
    G4Transform3D transformMP1(*rotTOFLGP, posTOFLG1);
    TOFLGPV[k*2] = new G4PVPlacement(transformMP1,"TOFLGPV", TOFLGLV, m_world_pv, FALSE, 0);
    G4Transform3D transformMP2(*rotTOFLGP, posTOFLG2);
    TOFLGPV[k*2+1] = new G4PVPlacement(transformMP2,"TOFLGPV", TOFLGLV, m_world_pv, FALSE, 0);
    rotTOFLGP->rotateY(dangleTOFLG);
    posTOFLG1.rotateY(dangleTOFLG);
    posTOFLG2.rotateY(dangleTOFLG);
    }


    TOFLGVisAtt= new G4VisAttributes(true, G4Colour(0.,0.8,0.));
    TOFLGLV->SetVisAttributes(TOFLGVisAtt);
    //    TOFLGVisAtt[k*2+1]= new G4VisAttributes(true, AQUA );
    TOFLGVisAtt= new G4VisAttributes(true, G4Colour(0.,0.8,0.));
    TOFLGLV->SetVisAttributes(TOFLGVisAtt);

    // ==============================================================
    // end light guide
    // ==============================================================
    */

}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructHypTPC( void )
{
  auto tpc_sd = new TPCPadSD("/TPC");
  G4SDManager::GetSDMpointer()->AddNewDetector( tpc_sd );
  const auto& tpc_pos = gGeom.GetGlobalPosition("HypTPC");
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
    auto tpc_solid = new G4Polyhedra( "TpcSolid",
				      22.5*deg, ( 360. + 22.5 )*deg,
				      NumOfSide, NumOfZPlane,
				      zPlane, rInner, rOuter );
    m_tpc_lv = new G4LogicalVolume( tpc_solid, m_material_map["P10"],
				    "TpcLV" );
    auto rot = new G4RotationMatrix;
    rot->rotateX( 90.*deg );
    new G4PVPlacement( rot, tpc_pos, m_tpc_lv, "TpcPV",
		       m_world_lv, false, 0 );
    m_tpc_lv->SetVisAttributes( WHITE );
    m_tpc_lv->SetSensitiveDetector( tpc_sd );
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
		       "FcPV", m_tpc_lv, false, 0 );
    fc_lv->SetVisAttributes( RED );
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
    // numpads[10]=208.;
    // numpads[11]=218.;
    // numpads[12]=230.;
    // numpads[13]=214.;
    // numpads[14]=212.;
    // numpads[15]=214.;
    // numpads[16]=220.;
    // numpads[17]=224.;
    // numpads[18]=232.;
    // numpads[19]=238.;
    // numpads[20]=244.;
    // numpads[21]=232.;
    // numpads[22]=218.;
    // numpads[23]=210.;
    // numpads[24]=206.;
    // numpads[25]=202.;
    // numpads[26]=200.;
    // numpads[27]=196.;
    // numpads[28]=178.;
    // numpads[29]=130.;
    // numpads[30]=108.;
    // numpads[31]=90.;
    break;
  default:
    break;
  }

  G4double below_target = 0;
  if(m_experiment==42){
    below_target=32.4;
  }else if(m_experiment==45||m_experiment==27){
    below_target=60.4;
  }
  // Inner Pads
  for( G4int i=0; i<NumOfPadTPCIn; ++i ){
    if( pad_out[i]<below_target ){
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
    pad_lv[i]->SetVisAttributes( BLUE );
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
  dead_lv->SetVisAttributes( GRAY );
  // Virtual pad
  G4Tubs* vpad_solid[NumOfPadTPC];
  G4LogicalVolume* vpad_lv[NumOfPadTPC];
  for( G4int i=0; i<NumOfPadTPCIn; ++i ){
    vpad_solid[i] = new G4Tubs( Form("TpcVPadSolid%d", i), pad_in[i]*mm,
				pad_out[i]*mm, 0.5*mm, 0., 360.*deg );
    vpad_lv[i] = new G4LogicalVolume( vpad_solid[i], m_material_map["P10"],
				      Form("TpcVPadLV%d", i) );
    vpad_lv[i]->SetVisAttributes( BLUE );
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
TPCDetectorConstruction::ConstructKuramaMagnet( void )
{
  const auto field_size = gSize.GetSize("KuramaField") * 0.5 * mm;
  const auto coil_color = RED;
  const auto yoke_color = LAVENDER;

  // size
  const G4ThreeVector coil1_size( 900.*mm/2, 193.*mm/2, 280.*mm/2 );
  const G4ThreeVector coil2_size( 193.*mm/2, 117.*mm/2, 280.*mm/2 );
  const G4ThreeVector coil3_size( 193.*mm/2, 137.*mm/2, 280.*mm/2 );
  const G4ThreeVector coil4_size( 193.*mm/2, 280.*mm/2, 740.*mm/2 );
  const G4ThreeVector coil5_size( 900.*mm/2, 214.*mm/2, 280.*mm/2 );
  const G4ThreeVector yoke_ud_size( 2200.*mm/2, 370.*mm/2, field_size.z() );
  const G4ThreeVector yoke_lr_size( 200.*mm/2, field_size.y(), field_size.z() );
  const G4ThreeVector uguard_ud_size( 1900.*mm/2, 620.*mm/2, 100.*mm/2 );
  const G4ThreeVector uguard_lr_size( 800.*mm/2, 300.*mm/2, 100.*mm/2 );
  const G4ThreeVector dguard_ud_size( 1600.*mm/2, 420.*mm/2, 100.*mm/2 );
  const G4ThreeVector dguard_lr_size( 250.*mm/2,
				      field_size.y() + 800.*mm/2,
				      100.*mm/2 );
  // position
  const G4ThreeVector field_pos = gGeom.GetGlobalPosition("KURAMA");
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
  const G4ThreeVector uguard_u_pos =
    { 0.*mm,
      uguard_lr_size.y() + uguard_ud_size.y(),
      field_pos.z() - 820.*mm + uguard_ud_size.z() };
  const G4ThreeVector uguard_d_pos =
    { 0.*mm,
      - uguard_lr_size.y() - uguard_ud_size.y(),
      field_pos.z() - 820.*mm + uguard_ud_size.z() };
  const G4ThreeVector uguard_l_pos =
    { 300.0*mm + uguard_lr_size.x(), // 150 -->gap
      0.*mm,
      field_pos.z() - 820.*mm + uguard_lr_size.z() };
  const G4ThreeVector uguard_r_pos =
    { - 300.0*mm - uguard_lr_size.x(),
      0.*mm,
      field_pos.z() - 820.*mm + uguard_lr_size.z() };
  const G4ThreeVector dguard_u_pos =
    { field_pos.x(),
      dguard_lr_size.y() + dguard_ud_size.y() - 200.*mm,
      field_pos.z() + 820.*mm - dguard_ud_size.z() };
  const G4ThreeVector dguard_d_pos =
    { field_pos.x(),
      - dguard_lr_size.y() - dguard_ud_size.y() + 200.*mm,
      field_pos.z() + 820.*mm - dguard_ud_size.z() };
  const G4ThreeVector dguard_l_pos =
    { field_pos.x() + field_size.x() + 300.*mm + dguard_lr_size.x(),
      0.*mm,
      field_pos.z() + 820.*mm - dguard_lr_size.z() };
  const G4ThreeVector dguard_r_pos =
    { field_pos.x() - field_size.x() + 300.*mm + dguard_lr_size.x(),
      0.*mm,
      field_pos.z() + 820.*mm - dguard_lr_size.z() };

  // Construct KURAMA Magnet
  // auto kurama_solid = new G4Box( "KuramaSolid", 4.*m/2, 3.*m/2, 4.*m/2 );
  // auto kurama_lv = new G4LogicalVolume( kurama_solid, m_material_map["Air"],
  // 					"KuramaLV" );
  // kurama_lv->SetVisAttributes( new G4VisAttributes( false, BLACK ) );
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
  // maxStep=0.00001*mm;
  // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

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
  // maxStep=0.00001*mm;
  // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

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
  // maxStep=0.00001*mm;
  // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
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
  // maxStep=0.00001*mm;
  // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
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
  //////////////coil3
  G4Box* Coil3_box = new G4Box("Coil3_box",
			       coil3_size.x(), coil3_size.y(), coil3_size.z() );
  G4LogicalVolume*  Coil3_log = new G4LogicalVolume(Coil3_box, m_material_map["Copper"], "Coil3_log",0,0,0);
  Coil3_log->SetVisAttributes( coil_color );
  // maxStep=0.00001*mm;
  // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
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

  //-------------------- Upstraam End Guard
  if(m_experiment!=3 && m_experiment!=42 && m_experiment!=27){
    G4Box* upGuard_UD_box = new G4Box("upGuard_UD_box",
				      uguard_ud_size.x(),uguard_ud_size.y(),uguard_ud_size.z());
    G4LogicalVolume*  upGuard_U_log = new G4LogicalVolume(upGuard_UD_box, m_material_map["Iron"], "upGuard_U_log",0,0,0);
    upGuard_U_log->SetVisAttributes( yoke_color );
    // maxStep=0.00001*mm;
    // upGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));
    new G4PVPlacement( m_rotation_matrix,
		       G4ThreeVector( uguard_u_pos.x(),
				      uguard_u_pos.y(),
				      uguard_u_pos.z() ).rotateY( m_rotation_angle ),
		       upGuard_U_log,
		       "upGuard_U_phys",
		       m_world_lv,
		       false,
		       0 );

    G4LogicalVolume*  upGuard_D_log = new G4LogicalVolume(upGuard_UD_box, m_material_map["Iron"], "upGuard_D_log",0,0,0);
    upGuard_D_log->SetVisAttributes( yoke_color );
    // maxStep=0.00001*mm;
    // upGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
    new G4PVPlacement( m_rotation_matrix,
		       G4ThreeVector( uguard_d_pos.x(),
				      uguard_d_pos.y(),
				      uguard_d_pos.z() ).rotateY( m_rotation_angle ),
		       upGuard_D_log,
		       "upGuard_D_phys",
		       m_world_lv,
		       false,
		       0 );

    G4Box* upGuard_LR_box = new G4Box("upGuard_LR_box",
				      uguard_lr_size.x(),uguard_lr_size.y(),uguard_lr_size.z());
    G4LogicalVolume*  upGuard_L_log = new G4LogicalVolume(upGuard_LR_box, m_material_map["Iron"], "upGuard_L_log",0,0,0);
    upGuard_L_log->SetVisAttributes( yoke_color );
    // maxStep=0.00001*mm;
    // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
    new G4PVPlacement( m_rotation_matrix,
		       G4ThreeVector( uguard_l_pos.x(),
				      uguard_l_pos.y(),
				      uguard_l_pos.z() ).rotateY( m_rotation_angle ),
		       upGuard_L_log,
		       "upGuard_L_phys",
		       m_world_lv,
		       false,
		       0 );

    G4LogicalVolume*  upGuard_R_log = new G4LogicalVolume(upGuard_LR_box, m_material_map["Iron"], "upGuard_R_log",0,0,0);
    upGuard_R_log->SetVisAttributes( yoke_color );
    // maxStep=0.00001*mm;
    // upGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
    new G4PVPlacement( m_rotation_matrix,
		       G4ThreeVector( uguard_r_pos.x(),
				      uguard_r_pos.y(),
				      uguard_r_pos.z() ).rotateY( m_rotation_angle ),
		       upGuard_R_log,
		       "upGuard_R_phys",
		       m_world_lv,
		       false,
		       0 );
  }

  //-------------------- Yoke
  G4Box* Yoke_UD_box = new G4Box( "Yoke_UD_box",
				  yoke_ud_size.x(), yoke_ud_size.y(), yoke_ud_size.z() );
  G4LogicalVolume*  Yoke_U_log = new G4LogicalVolume(Yoke_UD_box, m_material_map["Iron"], "Yoke_U_log",0,0,0);

  Yoke_U_log->SetVisAttributes( yoke_color );
  // maxStep=0.00001*mm;
  // Yoke_U_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( yoke_u_pos.x(),
				    yoke_u_pos.y(),
				    yoke_u_pos.z() ).rotateY( m_rotation_angle ),
		     Yoke_U_log,
		     "Yoke_U_phys",
		     m_world_lv,
		     false,
		     0 );
  G4LogicalVolume*  Yoke_D_log = new G4LogicalVolume(Yoke_UD_box, m_material_map["Iron"], "Yoke_D_log",0,0,0);
  Yoke_D_log->SetVisAttributes( yoke_color );
  // maxStep=0.00001*mm;
  // Yoke_D_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( yoke_d_pos.x(),
				    yoke_d_pos.y(),
				    yoke_d_pos.z() ).rotateY( m_rotation_angle ),
		     Yoke_D_log,
		     "Yoke_D_phys",
		     m_world_lv,
		     false,
		     0 );
  G4Box* Yoke_LR_box = new G4Box( "Yoke_LR_box",
				  yoke_lr_size.x(), yoke_lr_size.y(), yoke_lr_size.z() );
  G4LogicalVolume*  Yoke_L_log = new G4LogicalVolume(Yoke_LR_box, m_material_map["Iron"], "Yoke_L_log",0,0,0);
  Yoke_L_log->SetVisAttributes( yoke_color );
  // maxStep=0.00001*mm;
  // Yoke_L_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( yoke_l_pos.x(),
				    yoke_l_pos.y(),
				    yoke_l_pos.z() ).rotateY( m_rotation_angle ),
		     Yoke_L_log,
		     "Yoke_L_phys",
		     m_world_lv,
		     false,
		     0 );
  G4LogicalVolume*  Yoke_R_log = new G4LogicalVolume(Yoke_LR_box, m_material_map["Iron"], "Yoke_R_log",0,0,0);
  Yoke_R_log->SetVisAttributes( yoke_color );
  // maxStep=0.00001*mm;
  // Yoke_R_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( yoke_r_pos.x(),
				    yoke_r_pos.y(),
				    yoke_r_pos.z() ).rotateY( m_rotation_angle ),
		     Yoke_R_log,
		     "Yoke_R_phys",
		     m_world_lv,
		     false,
		     0 );

  //-------------------- Downstream End Guard
  G4Box* downGuard_UD_box = new G4Box("downGuard_UD_box",
				      dguard_ud_size.x(),dguard_ud_size.y(),dguard_ud_size.z());
  G4LogicalVolume*  downGuard_U_log = new G4LogicalVolume(downGuard_UD_box,
							  m_material_map["Iron"],
							  "downGuard_U_log");
  downGuard_U_log->SetVisAttributes( yoke_color );
  // maxStep=0.00001*mm;
  // downGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( dguard_u_pos.x(),
				    dguard_u_pos.y(),
				    dguard_u_pos.z() ).rotateY( m_rotation_angle ),
		     downGuard_U_log,
		     "downGuard_U_phys",
		     m_world_lv,
		     false,
		     0 );
  G4LogicalVolume*  downGuard_D_log = new G4LogicalVolume(downGuard_UD_box,
							  m_material_map["Iron"],
							  "downGuard_D_log");
  downGuard_D_log->SetVisAttributes( yoke_color );
  // maxStep=0.00001*mm;
  // downGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( dguard_d_pos.x(),
				    dguard_d_pos.y(),
				    dguard_d_pos.z() ).rotateY( m_rotation_angle ),
		     downGuard_D_log,
		     "downGuard_D_phys",
		     m_world_lv,
		     false,
		     0 );
  G4Box* downGuard_LR_box = new G4Box("downGuard_LR_box",
				      dguard_lr_size.x(),
				      dguard_lr_size.y(),
				      dguard_lr_size.z());
  G4LogicalVolume*  downGuard_L_log = new G4LogicalVolume(downGuard_LR_box,
							  m_material_map["Iron"],
							  "downGuard_L_log");
  downGuard_L_log->SetVisAttributes( yoke_color );
  // maxStep=0.00001*mm;
  // downGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( dguard_l_pos.x(),
				    dguard_l_pos.y(),
				    dguard_l_pos.z() ).rotateY( m_rotation_angle ),
		     downGuard_L_log,
		     "downGuard_L_phys",
		     m_world_lv,
		     false,
		     0 );
  G4LogicalVolume*  downGuard_R_log = new G4LogicalVolume( downGuard_LR_box,
							   m_material_map["Iron"],
							   "downGuard_R_log" );
  downGuard_R_log->SetVisAttributes( yoke_color );
  // maxStep=0.00001*mm;
  // downGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( dguard_r_pos.x(),
				    dguard_r_pos.y(),
				    dguard_r_pos.z() ).rotateY( m_rotation_angle ),
		     downGuard_R_log,
		     "downGuard_R_phys",
		     m_world_lv,
		     false,
		     0 );
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
	G4ThreeVector posMOut1(n_bar_radius*mm,0.*mm,+plane_width_nbar/2-DZ_NBAR1-DZ_NBAR1*2*j); //x,z,y??
	G4RotationMatrix* rotMOutP = new G4RotationMatrix;
	G4Transform3D transformMP1(rotMOutP->rotateY(dangleOut_nbar*(i)), posMOut1.rotateY(dangleOut_nbar*(i)));
	nbarLV[i*num_one_plane_nbar+j] = new G4LogicalVolume(nbarSolid1, m_material_map["Scintillator"], Form("NBARLV%d",i*num_one_plane_nbar+j));
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
  G4LogicalVolume* SCH_log[NumOfSegSCH];
  G4double size_SCH[3];
  size_SCH[ThreeVector::X] = 11.5*mm;
  size_SCH[ThreeVector::Y] = 400.0*mm;
  size_SCH[ThreeVector::Z] = 2.0*mm;
  G4Box* SCH_box = new G4Box("SCH_box",size_SCH[ThreeVector::X]/2,size_SCH[ThreeVector::Y]/2,size_SCH[ThreeVector::Z]/2);
  G4double SCH_Overlap = 1.*mm;
  //  G4LogicalVolume* SCH_log[NumOfSegSCH];
  for (int i=0; i<NumOfSegSCH; i++) {
    SCH_log[i] = new G4LogicalVolume(SCH_box, m_material_map["Scintillator"], Form("SCH%d_log", i),
				    0, 0, 0 );
    SCH_log[i]->SetVisAttributes( AQUA );
    //    par_cham->calpc((double *)pos_SCH, SCHx, (double)(i+1));
    G4double ipos_x = -NumOfSegSCH/2.*(size_SCH[ThreeVector::X]-SCH_Overlap)+5.75*mm+(size_SCH[ThreeVector::X]-SCH_Overlap)*i;
    if(i%2==0){
      new G4PVPlacement( m_rotation_matrix,
			 G4ThreeVector( ipos_x,
					sch_pos[ThreeVector::Y],
					sch_pos[ThreeVector::Z] - 1.*mm ).rotateY( m_rotation_angle ),
			 SCH_log[i],
			 Form("SCH%d", i),
			 m_world_lv,
			 false,
			 i );
    }else if(i%2==1){
      new G4PVPlacement( m_rotation_matrix,
			 G4ThreeVector( ipos_x,
					sch_pos[ThreeVector::Y],
					sch_pos[ThreeVector::Z] + 1.*mm ).rotateY( m_rotation_angle ),
			 SCH_log[i],
			 Form("SCH%d", i),
			 m_world_lv,
			 false,
			 i);
    }
  }

  auto schSD = new TPCSCHSD("/SCH");
  G4SDManager::GetSDMpointer()->AddNewDetector( schSD );
  for( G4int i = 0; i<NumOfSegSCH; ++i ){
    SCH_log[i]->SetSensitiveDetector( schSD );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructSDC1( void )
{
  if( !m_dc_sd ){
    m_dc_sd = new TPCDCSD("/DC");
    G4SDManager::GetSDMpointer()->AddNewDetector( m_dc_sd );
  }

  const auto& sdc1_pos = ( gGeom.GetGlobalPosition("KURAMA") +
			   ( gGeom.GetGlobalPosition("SDC1-V1") +
			     gGeom.GetGlobalPosition("SDC1-U2") ) * 0.5 );
  G4LogicalVolume*    DC1Plane_log[6];
  G4double size_DC1[3];
  //  size_DC1[ThreeVector::X] = 289.0*2.*mm;
  //  size_DC1[ThreeVector::Y] = 215.0*2.*mm;
  //  size_DC1[ThreeVector::Z] = 146.0*mm;
  ///Tanida DC1 is to small
  size_DC1[ThreeVector::X] = 389.0*2.*mm;
  size_DC1[ThreeVector::Y] = 315.0*2.*mm;
  size_DC1[ThreeVector::Z] = 146.0*mm;

  G4double size_DC1Plane[3];
  size_DC1Plane[ThreeVector::X] = 6.0*127.0*0.5*mm;
  size_DC1Plane[ThreeVector::Y] = 6.0*97.0*0.5*mm;
  size_DC1Plane[ThreeVector::Z] = 0.0001*mm;

  G4Box* DC1_box = new G4Box("DC1_box",size_DC1[ThreeVector::X]/2.,size_DC1[ThreeVector::Y]/2,size_DC1[ThreeVector::Z]/2);
  G4LogicalVolume*  DC1_log = new G4LogicalVolume(DC1_box, m_material_map["Argon"], "DC1_log",0,0,0);
  DC1_log->SetVisAttributes( GREEN );
  // G4double sdc1_pos[3];
  // sdc1_pos[ThreeVector::X] = par_cham->get_DCPlaneCenter(DC1X, ThreeVector::X)*mm;
  // sdc1_pos[ThreeVector::Y] = par_cham->get_DCPlaneCenter(DC1X, ThreeVector::Y)*mm;
  // sdc1_pos[ThreeVector::Z] = (par_cham->get_DCPlaneCenter(DC1X, ThreeVector::Z)
  //		    +par_cham->get_DCPlaneCenter(DC1U, ThreeVector::Z))*0.5*mm;
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( sdc1_pos[ThreeVector::X],
				    sdc1_pos[ThreeVector::Y],
				    sdc1_pos[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     DC1_log,
		     "DC1_phys",
		     m_world_lv,
		     false,
		     0 );
  //---------DC1 Planes
  G4Box* DC1Plane_box = new G4Box("DC1Plane_box",
				  size_DC1Plane[ThreeVector::X],size_DC1Plane[ThreeVector::Y],size_DC1Plane[ThreeVector::Z]);
  // G4LogicalVolume*    DC1Plane_log[6];
  G4double pos_DC1Plane[3]={0.};
  G4String name1, name2;
  for (int i=0; i<6; i++) {
    switch (i) {
    case 0:
      name1 = "DC1U_log";
      name2 = "DC1U_phys";
      pos_DC1Plane[ThreeVector::X] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Z] = -50.0*mm;
      break;
    case 1:
      name1 = "DC1Up_log";
      name2 = "DC1Up_phys";
      pos_DC1Plane[ThreeVector::X] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Z] = -30.*mm;
      break;
    case 2:
      name1 = "DC1X_log";
      name2 = "DC1X_phys";
      pos_DC1Plane[ThreeVector::X] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Z] = -10.0*mm;
      break;
    case 3:
      name1 = "DC1Xp_log";
      name2 = "DC1Xp_phys";
      pos_DC1Plane[ThreeVector::X] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Z] = 10.0*mm;
      break;
    case 4:
      name1 = "DC1V_log";
      name2 = "DC1V_phys";
      pos_DC1Plane[ThreeVector::X] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Z] = 30.0*mm;
      break;
    case 5:
      name1 = "DC1Vp_log";
      name2 = "DC1Vp_phys";
      pos_DC1Plane[ThreeVector::X] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC1Plane[ThreeVector::Z] = 50.0*mm;
      break;
    }
    DC1Plane_log[i] = new G4LogicalVolume(DC1Plane_box, m_material_map["Argon"], name1,0,0,0);
    new G4PVPlacement( 0,
		       G4ThreeVector( pos_DC1Plane[ThreeVector::X],
				      pos_DC1Plane[ThreeVector::Y],
				      pos_DC1Plane[ThreeVector::Z] ),
		       DC1Plane_log[i],
		       name2,
		       DC1_log,
		       false,
		       101+i );
  }
  for(G4int i = 0; i<6; i++){
    DC1Plane_log[i]->SetSensitiveDetector( m_dc_sd );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructSDC2( void )
{
  if( !m_dc_sd ){
    m_dc_sd = new TPCDCSD("/DC");
    G4SDManager::GetSDMpointer()->AddNewDetector( m_dc_sd );
  }

  const auto& sdc2_pos = ( gGeom.GetGlobalPosition("KURAMA") +
			   ( gGeom.GetGlobalPosition("SDC2-X1") +
			     gGeom.GetGlobalPosition("SDC2-Y2") ) * 0.5 );
  G4LogicalVolume*    DC2Plane_log[6];
  G4double size_DC2[3];
  size_DC2[ThreeVector::X] = 1186.5*mm;//# of wires are 128
  size_DC2[ThreeVector::Y] = 1186.5*mm;//# of wires are 128
  //size_DC2[ThreeVector::Z] = 45.0*mm;
  size_DC2[ThreeVector::Z] = 100.0*mm;
  G4double size_DC2Plane[3];
  size_DC2Plane[ThreeVector::X] = 9.0*128.0*0.5*mm;
  size_DC2Plane[ThreeVector::Y] = 9.0*128.0*0.5*mm;
  size_DC2Plane[ThreeVector::Z] = 0.0001*mm;
  G4Box* DC2_box = new G4Box("DC2_box",size_DC2[ThreeVector::X]/2,size_DC2[ThreeVector::Y]/2,size_DC2[ThreeVector::Z]/2);
  G4LogicalVolume*  DC2_log = new G4LogicalVolume(DC2_box, m_material_map["Argon"], "DC2_log",0,0,0);
  DC2_log->SetVisAttributes( GREEN );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( sdc2_pos[ThreeVector::X],
				    sdc2_pos[ThreeVector::Y],
				    sdc2_pos[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     DC2_log,
		     "DC2_phys",
		     m_world_lv,
		     false,
		     0 );
  G4Box* DC2Plane_box = new G4Box("DC2Plane_box",
				  size_DC2Plane[ThreeVector::X],
				  size_DC2Plane[ThreeVector::Y],
				  size_DC2Plane[ThreeVector::Z]);
  //  G4LogicalVolume* DC2Plane_log[4];
  G4double pos_DC2Plane[3]={0.};
  G4String name1, name2;
  for (int i=0; i<4; i++) {
    switch (i) {
    case 0:
      name1 = "DC2X_log";
      name2 = "DC2X_phys";
      pos_DC2Plane[ThreeVector::X] = 0.0*mm;
      pos_DC2Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC2Plane[ThreeVector::Z] = -16.5*mm;
      break;
    case 1:
      name1 = "DC2Xp_log";
      name2 = "DC2Xp_phys";
      pos_DC2Plane[ThreeVector::X] = 0.0*mm;
      pos_DC2Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC2Plane[ThreeVector::Z] = -8.7*mm;
      break;
    case 2:
      name1 = "DC2Y_log";
      name2 = "DC2Y_phys";
      pos_DC2Plane[ThreeVector::X] = 0.0*mm;
      pos_DC2Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC2Plane[ThreeVector::Z] = 8.7*mm;
      break;
    case 3:
      name1 = "DC2Yp_log";
      name2 = "DC2Yp_phys";
      pos_DC2Plane[ThreeVector::X] = 0.0*mm;
      pos_DC2Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC2Plane[ThreeVector::Z] = 16.5*mm;
      break;
    }
    DC2Plane_log[i] = new G4LogicalVolume(DC2Plane_box, m_material_map["Argon"], name1,0,0,0);
    new G4PVPlacement( 0,
		       G4ThreeVector( pos_DC2Plane[ThreeVector::X],
				      pos_DC2Plane[ThreeVector::Y],
				      pos_DC2Plane[ThreeVector::Z] ),
		       DC2Plane_log[i],
		       name2,
		       DC2_log,
		       false,
		       121+i );
  }

  for( G4int i = 0; i<4; ++i ){
    DC2Plane_log[i]->SetSensitiveDetector( m_dc_sd );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructSDC3( void )
{
  if( !m_dc_sd ){
    m_dc_sd = new TPCDCSD("/DC");
    G4SDManager::GetSDMpointer()->AddNewDetector( m_dc_sd );
  }

  const auto& sdc3_pos = ( gGeom.GetGlobalPosition("KURAMA") +
			   ( gGeom.GetGlobalPosition("SDC3-Y1") +
			     gGeom.GetGlobalPosition("SDC3-X2") ) * 0.5 );
  G4LogicalVolume* DC3Plane_log[6];
  G4double size_DC3[3];
  size_DC3[ThreeVector::X] = 1900.0*mm;
  size_DC3[ThreeVector::Y] = 1280.0*mm;
  size_DC3[ThreeVector::Z] = 150.*mm;

  G4double size_DC3Plane[3];
  size_DC3Plane[ThreeVector::X] = 20.0*96.0*0.5*mm;
  size_DC3Plane[ThreeVector::Y] = 20.0*64.0*0.5*mm;
  size_DC3Plane[ThreeVector::Z] = 0.0001*mm;

  //--------------DC3
  G4Box* DC3_box = new G4Box( "DC3_box",
			      size_DC3[ThreeVector::X]/2,
			      size_DC3[ThreeVector::Y]/2,
			      size_DC3[ThreeVector::Z]/2 );
  G4LogicalVolume* DC3_log = new G4LogicalVolume( DC3_box, m_material_map["Argon"],
						  "DC3_log", 0, 0, 0 );
  DC3_log->SetVisAttributes( GREEN );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( sdc3_pos[ThreeVector::X],
				    sdc3_pos[ThreeVector::Y],
				    sdc3_pos[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     DC3_log,
		     "DC3_phys",
		     m_world_lv,
		     false,
		     0 );

  //---------DC3 Planes
  G4Box* DC3Plane_box = new G4Box("DC3Plane_box",
				  size_DC3Plane[ThreeVector::X],size_DC3Plane[ThreeVector::Y],size_DC3Plane[ThreeVector::Z]);
  //  G4LogicalVolume* DC3Plane_log[4];
  G4double pos_DC3Plane[3]={0.};
  G4String name1, name2;
  for (int i=0; i<4; i++) {
    switch (i) {
    case 0:
      name1 = "DC3X_log";
      name2 = "DC3X_phys";
      pos_DC3Plane[ThreeVector::X] = 0.0*mm;
      pos_DC3Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC3Plane[ThreeVector::Z] = -46.5*mm;
      break;
    case 1:
      name1 = "DC3Xp_log";
      name2 = "DC3Xp_phys";
      pos_DC3Plane[ThreeVector::X] = 0.0*mm;
      pos_DC3Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC3Plane[ThreeVector::Z] = -14.5*mm;
      break;
    case 2:
      name1 = "DC3Y_log";
      name2 = "DC3Y_phys";
      pos_DC3Plane[ThreeVector::X] = 0.0*mm;
      pos_DC3Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC3Plane[ThreeVector::Z] = 16.5*mm;
      break;
    case 3:
      name1 = "DC3Yp_log";
      name2 = "DC3Yp_phys";
      pos_DC3Plane[ThreeVector::X] = 0.0*mm;
      pos_DC3Plane[ThreeVector::Y] = 0.0*mm;
      pos_DC3Plane[ThreeVector::Z] = 46.5*mm;
      break;
    }
    DC3Plane_log[i] = new G4LogicalVolume(DC3Plane_box, m_material_map["Argon"], name1,0,0,0);
    new G4PVPlacement( 0,
		       G4ThreeVector( pos_DC3Plane[ThreeVector::X],
				      pos_DC3Plane[ThreeVector::Y],
				      pos_DC3Plane[ThreeVector::Z] ),
		       DC3Plane_log[i],
		       name2,
		       DC3_log,
		       false,
		       131+i );
  }
  for( G4int i = 0; i<4; ++i ){
    DC3Plane_log[i]->SetSensitiveDetector( m_dc_sd );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructShsMagnet( void )
{
  const auto& tpc_pos = gGeom.GetGlobalPosition("HypTPC");
  const G4double DPHI_TPC   = 360.*deg;
  G4Tubs* HelmSolid_t =
    new G4Tubs("HelmSolid_t",45.*cm,80.*cm,135./2.*cm,0.*deg,DPHI_TPC);

  G4Box* HelmDHole =
    new G4Box("HelmDHole",100./2.*cm,100./2.*cm,60./2.*cm);
  G4ThreeVector HTrans(0, -60.*cm, 0);
  G4RotationMatrix* yRot = new G4RotationMatrix;  // Rotates X and Z axes only

  G4Transform3D transform(*yRot, HTrans);


  G4SubtractionSolid* HelmSolid= new G4SubtractionSolid("HelmSolid",HelmSolid_t,HelmDHole,transform);

  //    new G4Tubs("HelmSolid",40.*cm,50.*cm,5*cm,0.*deg,DPHI_TPC);
  G4LogicalVolume* HelmLV=
    new G4LogicalVolume(HelmSolid, m_material_map["Iron"], "HelmLV");
  G4RotationMatrix *rotHelm = new G4RotationMatrix();
  rotHelm->rotateX(90.*deg);
  rotHelm->rotateZ( - m_rotation_angle );
  new G4PVPlacement( rotHelm, G4ThreeVector( 0.*mm, 0.*mm, tpc_pos.z()*mm ), HelmLV,
		     "HelmPV", m_world_lv, FALSE, 0 );

  //    new G4PVPlacement(rotHelm, G4ThreeVector(0,+20.*cm,0),
  //		      "HelmPV", HelmLV, m_world_pv, FALSE, 0);

  //  HelmPV[1] =
  //    new G4PVPlacement(rotHelm, G4ThreeVector(0,-20.*cm,0),
  //		      "HelmPV", HelmLV, m_world_pv, FALSE, 0);
  G4VisAttributes* HelmVisAtt= new G4VisAttributes( true, PINK );
  HelmLV->SetVisAttributes(HelmVisAtt);

  /*  /////// supporter (from bottom to helmholtz)
      G4Box* HelmSuSolid =
      new G4Box("HelmSuSolid",3.*cm,60.*cm,3*cm);
      G4LogicalVolume* HelmSuLV=
      new G4LogicalVolume(HelmSuSolid, m_material_map["Iron"], "HelmSuLV");
      G4PVPlacement* HelmSuPV[2];
      HelmSuPV[0] =
      new G4PVPlacement(0, G4ThreeVector(31*cm,-85.*cm,31*cm),
      "HelmSuPV", HelmSuLV, m_world_pv, FALSE, 0);
      HelmSuPV[1] =
      new G4PVPlacement(0, G4ThreeVector(-31*cm,-85.*cm,31*cm),
      "HelmSuPV", HelmSuLV, m_world_pv, FALSE, 0);
      HelmSuPV[2] =
      new G4PVPlacement(0, G4ThreeVector(-31*cm,-85.*cm,-31*cm),
      "HelmSuPV", HelmSuLV, m_world_pv, FALSE, 0);
      HelmSuPV[3] =
      new G4PVPlacement(0, G4ThreeVector(31*cm,-85.*cm,-31*cm),
      "HelmSuPV", HelmSuLV, m_world_pv, FALSE, 0);
      G4VisAttributes* HelmSuVisAtt= new G4VisAttributes(true, G4Colour(0.,0.0,1.));
      HelmSuLV->SetVisAttributes(HelmSuVisAtt);

      /////// supporter (Between helmholtz coils)
      G4Box* HelmSuBeSolid =
      new G4Box("HelmSu1Solid",1.*cm,15.*cm,1*cm);
      G4LogicalVolume* HelmSuBeLV=
      new G4LogicalVolume(HelmSuBeSolid, m_material_map["Iron"], "HelmSuBeLV");
      G4PVPlacement* HelmSuBePV[2];
      HelmSuBePV[0] =
      new G4PVPlacement(0, G4ThreeVector(31*cm,0.*cm,31*cm),
      "HelmSuBePV", HelmSuBeLV, m_world_pv, FALSE, 0);
      HelmSuBePV[1] =
      new G4PVPlacement(0, G4ThreeVector(-31*cm,0.*cm,31*cm),
      "HelmSuBePV", HelmSuBeLV, m_world_pv, FALSE, 0);
      HelmSuBePV[2] =
      new G4PVPlacement(0, G4ThreeVector(-31*cm,0.*cm,-31*cm),
      "HelmSuBePV", HelmSuBeLV, m_world_pv, FALSE, 0);
      HelmSuBePV[3] =
      new G4PVPlacement(0, G4ThreeVector(31*cm,0.*cm,-31*cm),
      "HelmSuBePV", HelmSuBeLV, m_world_pv, FALSE, 0);
      G4VisAttributes* HelmSuBeVisAtt= new G4VisAttributes(true, G4Colour(0.,0.0,1.));
      HelmSuBeLV->SetVisAttributes(HelmSuBeVisAtt);
  */

  auto myfield = new TPCField("helmholtz_field.dat", "KuramaMap80cm.dat");
  auto fieldMgr =
    G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetDetectorField( myfield );
  fieldMgr->CreateChordFinder( myfield );
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructTarget( void )
{
  auto target_pos = gGeom.GetGlobalPosition( "Target" );
  G4double target_pos_z = target_pos.z()*mm;
  G4LogicalVolume* TargetLV = nullptr;
  G4LogicalVolume* TargetHolderLV;
  switch( m_experiment ){
  case 42: {
    G4double Target_x = gSize.Get( "Target", ThreeVector::X );
    G4double Target_y = gSize.Get( "Target", ThreeVector::Y );
    G4double Target_z = gSize.Get( "Target", ThreeVector::Z );
    //  G4Box* TargetSolid = new G4Box("target", 1.5*cm,0.25*cm,0.5*cm); // 30 x 10 x 5

    G4Box* TargetSolid = new G4Box("target", Target_x*mm/2,Target_z*mm/2,Target_y*mm/2); // 30 x 10 x 15
    TargetLV=
      new G4LogicalVolume(TargetSolid, m_material_map["Target"], "TargetLV");
    new G4PVPlacement( 0, G4ThreeVector( 0*mm, target_pos_z, 0.*mm ),
		       TargetLV, "TargetPV", m_world_lv, true, 0 );
    G4VisAttributes* TargetVisAtt= new G4VisAttributes(true, G4Colour(1.,0.,0.));
    TargetLV->SetVisAttributes(TargetVisAtt);
    G4Tubs* TargetHolderSolid = new G4Tubs("TargetHolderSolid", 16.*mm, 16.2*mm,
					   155.*mm, 0., 360*deg);
    TargetHolderLV=
      new G4LogicalVolume(TargetHolderSolid, m_material_map["P10"], "TargetHolderLV");

    new G4PVPlacement( 0, G4ThreeVector( 0*mm, target_pos_z, 145.*mm ),
		       TargetHolderLV, "TargetHolderPV", m_world_lv, true, 0 );
    G4VisAttributes* TargetHolderVisAtt= new G4VisAttributes(true, G4Colour(0.1,0.1,0.));
    TargetHolderLV->SetVisAttributes(TargetHolderVisAtt);
  }
    break;
  case 45: case 27: {
    //  G4Box* TargetSolid = new G4Box("target", 1.5*cm,0.25*cm,0.5*cm); // 30 x 10 x 5
    G4double Target_r = gSize.Get( "Target", ThreeVector::X );
    // G4double Target_y = gSize.Get( "Target", ThreeVector::Y );
    G4double Target_z = gSize.Get( "Target", ThreeVector::Z );

    G4Tubs* TargetSolid = new G4Tubs("TargetSolid", 0.*mm,
				     Target_r*mm,Target_z*mm, 0., 360*deg);
    //  G4Box* TargetSolid = new G4Box("target", 1.*cm,0.25*cm,1.*cm);
    TargetLV=
      new G4LogicalVolume(TargetSolid, m_material_map["Target"], "TargetLV");
    new G4PVPlacement( 0, G4ThreeVector( 0*mm, target_pos_z, 0.*mm ),
		       TargetLV, "TargetPV", m_world_lv, true, 0 );
    G4VisAttributes* TargetVisAtt= new G4VisAttributes(true, G4Colour(1.,0.,0.));
    TargetLV->SetVisAttributes(TargetVisAtt);
    /////////////////////////////////////
    ///////// target holder  ////////////
    /////////////////////////////////////

    G4Tubs* TargetHolderSolid = new G4Tubs("TargetHolderSolid", (Target_r+15.)*mm, (Target_r+15.2)*mm,
					   200.*mm, 0., 360*deg);
    TargetHolderLV=
      new G4LogicalVolume(TargetHolderSolid, m_material_map["P10"], "TargetHolderLV");

    new G4PVPlacement( 0, G4ThreeVector( 0*mm, target_pos_z, 100.*mm ),
		       TargetHolderLV, "TargetHolderPV", m_world_lv, true, 0 );
    G4VisAttributes* TargetHolderVisAtt= new G4VisAttributes(true, G4Colour(0.1,0.1,0.));
    TargetHolderLV->SetVisAttributes(TargetHolderVisAtt);
  }
    break;
  default:
    return;
  }
  TPCTargetSD* tarSD= new TPCTargetSD("/TAR");
  G4SDManager::GetSDMpointer()->AddNewDetector(tarSD);
  TargetLV->SetSensitiveDetector(tarSD);
}

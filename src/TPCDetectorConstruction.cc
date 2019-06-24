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
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4VisAttributes.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "DetSizeMan.hh"
#include "ThreeVector.hh"
#include "TPCACSD.hh"
#include "TPCSCHSD.hh"
#include "TPCDCSD.hh"
#include "TPCField.hh"
#include "TPCFTOFSD.hh"
#include "TPCNBARSD.hh"
#include "TPCPadSD.hh"
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
  const G4Colour BLACK( 1., 1., 1. );
  const G4Colour RED( 1.0, 0.0, 0.0 );
  const G4Colour GREEN( 0.0, 1.0, 0.0 );
  const G4Colour BLUE( 0.0, 0.0, 1.0 );
  const G4Colour CYAN( 0.0, 1.0, 1.0 );
  const G4Colour AQUA( 0.247, 0.8, 1.0 );
  const G4Colour MAGENTA( 1.0, 0.0, 1.0 );
  const G4Colour YELLOW( 1.0, 1.0, 0.0 );
  const G4Colour GRAY( 0.5, 0.5, 0.5 );
  const G4Colour LAVENDER( 0.901, 0.901, 0.98 );
  const G4Colour MAROON( 0.5, 0.0, 0.0 );
}

//_____________________________________________________________________________
TPCDetectorConstruction::TPCDetectorConstruction( void )
  : m_experiment( gConf.Get<Int_t>("Experiment") ),
    m_element_map(),
    m_material_map(),
    m_world_lv(),
    m_world_pv(),
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
  auto world_solid = new G4Box("world_s", 4.0*m, 3.0*m, 8.0*m);
  m_world_lv = new G4LogicalVolume( world_solid, m_material_map["Air"],
				    "world_lv" );
  auto world_vis_attr = new G4VisAttributes( false, BLACK );
  m_world_lv->SetVisAttributes( world_vis_attr );
  m_world_pv = new G4PVPlacement( 0, G4ThreeVector(), "world_pv",
				  m_world_lv, 0, false, 0 );
  ConstructHypTPC();
  ConstructHTOF();
  ConstructTarget();
  ConstructPVAC2();
  ConstructNBAR();
  ConstructShsMagnet();
  if( gConf.Get<G4int>("ConstructKurama") ){
    ConstructKuramaMagnet();
    ConstructSDC1();
    ConstructSCH();
    ConstructSDC2();
    ConstructSDC3();
    ConstructFTOF();
  }
  return m_world_pv;
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
    G4cout << " * Wrong target material" << G4endl;
    exit(-1);
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructFTOF( void )
{
  const auto& tof_pos = gGeom.GetGlobalPosition("TOF");
  const auto& tof_ra2 = gGeom.GetRotAngle2("TOF") * deg;
  G4LogicalVolume* FTOF_log[NumOfSegTOF];
  G4double size_FTOF[3];
  size_FTOF[ThreeVector::X] = 40.0*mm;
  size_FTOF[ThreeVector::Y] = 900.0*mm;
  size_FTOF[ThreeVector::Z] = 15.0*mm;

  G4String name1, name2;
  G4Box* FTOF_mo_box = new G4Box("FTOF_mo_box",size_FTOF[ThreeVector::X]*24+50.*mm,size_FTOF[ThreeVector::Y]+50.*mm,size_FTOF[ThreeVector::Z]+0.1*mm);
  G4LogicalVolume* FTOF_mo_log = new G4LogicalVolume(FTOF_mo_box, m_material_map["Air"], "FTOF_mo_log",0,0,0);
  G4RotationMatrix* rot_FTOF = new G4RotationMatrix();
  rot_FTOF->rotateY(-tof_ra2 - m_rotation_angle );
  new G4PVPlacement( rot_FTOF,
		     G4ThreeVector( tof_pos[ThreeVector::X],
				    tof_pos[ThreeVector::Y],
				    tof_pos[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     FTOF_mo_log,
		     "FTOF_mo_phys",
		     m_world_lv,
		     false,
		     0 );
  G4VisAttributes* FTOF_mo_VisAtt= new G4VisAttributes(false, G4Colour(1.,0.,0.));
  FTOF_mo_log->SetVisAttributes(FTOF_mo_VisAtt);

  ////end mother volume of FTOF
  G4Box* FTOF_box = new G4Box("FTOF_box",size_FTOF[ThreeVector::X],size_FTOF[ThreeVector::Y],size_FTOF[ThreeVector::Z]);
  //  G4LogicalVolume* FTOF_log[NumOfSegTOF];

  //  G4RotationMatrix* rot_FTOF = new G4RotationMatrix();
  //  rot_FTOF->rotateY(-tof_ra2 - m_rotation_angle );
  G4double FTOF_Overlap=5.*mm;
  G4double pos_FTOF_bar[3];
  for (int i=0; i<NumOfSegTOF; i++) {
    FTOF_log[i] = new G4LogicalVolume(FTOF_box, m_material_map["Scintillator"],
				      Form("FTOF%d_log", i),0,0,0);
    FTOF_log[i]->SetVisAttributes( AQUA );
    // maxStep=0.1*mm;
    // FTOF_log[i]->SetUserLimits(new G4UserLimits(maxStep));
    pos_FTOF_bar[ThreeVector::X]=-NumOfSegTOF/2*(size_FTOF[ThreeVector::X]*2-FTOF_Overlap)+(size_FTOF[ThreeVector::X]*2-FTOF_Overlap)*i;
    //    pos_FTOF[ThreeVector::Z]=pos_FTOF[ThreeVector::Z]-sin((size_FTOF[ThreeVector::X]*2.-FTOF_Overlap)*i);
    if( i%2==0 ){
      new G4PVPlacement( 0,
			 G4ThreeVector( pos_FTOF_bar[ThreeVector::X],
					0.,
					size_FTOF[ThreeVector::Z] ),
			 FTOF_log[i],
			 Form("FTOF%d", i),
			 FTOF_mo_log,
			 false,
			 i);
    } else {
      new G4PVPlacement( 0,
			 G4ThreeVector( pos_FTOF_bar[ThreeVector::X],
					0.,
					-size_FTOF[ThreeVector::Z] ),
			 FTOF_log[i],
			 Form("FTOF%d", i),
			 FTOF_mo_log,
			 false,
			 i);
    }
  }
  TPCFTOFSD* ftofSD = new TPCFTOFSD("/FTOF");
  G4SDManager::GetSDMpointer()->AddNewDetector( ftofSD );
  for( G4int i = 0; i<24; ++i ){
    FTOF_log[i]->SetSensitiveDetector( ftofSD );
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructHTOF( void )
{
  // ==============================================================
  // Scintillators
  // ==============================================================
  G4LogicalVolume* scintLV[100];
  G4VisAttributes* scintVisAtt[100];

  // ==============================================================
  // Side Scintillators (Outer)
  // ==============================================================
  //thickness 5mm
  //  const G4double R_SCINT =  ROUT_TPC*0.5*sqrt(3.0)+2.5*1.0;
  const G4int num_one_plane=NumOfSegHTOF/8.;
  const G4double DX_SCINT1 = 10.0/2.*mm;
  //  const G4double DY_SCINT = 72.5*mm;
  //  const G4double DZ_SCINT1 = (350.-DX_SCINT1)*tan(22.5*deg)/2*mm;
  //  const G4double DZ_SCINT1 = (320.-DX_SCINT1)*tan(22.5*deg)*2/num_one_plane/2* mm;
  const G4double DZ_SCINT1 = (337.-DX_SCINT1)*tan(22.5*deg)*2/num_one_plane/2* mm;
  const G4double DY_SCINT1 = 400.0*mm;
  G4Box* scintSolid1= new G4Box("SIDE SCINT1", DX_SCINT1,
				DY_SCINT1,  DZ_SCINT1);
  //side scint.
  //  G4double plane_width = (320.-DX_SCINT1)*tan(22.5*deg)*2;
  G4double plane_width = (337.-DX_SCINT1)*tan(22.5*deg)*2;
  const G4double dangleOut = 22.5*2*deg;
  //    for(G4int m=0;m<num_one_place;m++){
  //    G4ThreeVector posMOut2(320.*mm,0.*mm,(+DZ_SCINT1)); //x,z,y??


  //  for(G4int k=0;k<8; k++){
  //    rotMOutP->rotateY(dangleOut*0.5-(360/8./2.)*deg);
  //    posMOut1.rotateY(dangleOut*0.5-(360/8./2.)*deg);
  //    G4cout<<DZ_SCINT1<<" "<<num_one_plane<<" "<<plane_width<<G4endl;

  for(G4int i=0;i<8; i++){
    for(G4int j=0;j<num_one_plane;j++){
      G4ThreeVector posMOut1(337.*mm,0.*mm,+plane_width/2-DZ_SCINT1-DZ_SCINT1*2*j); //x,z,y??
      G4RotationMatrix* rotMOutP = new G4RotationMatrix;
      G4Transform3D transformMP1(rotMOutP->rotateY(dangleOut*(i)), posMOut1.rotateY(dangleOut*(i)));

      scintLV[i*num_one_plane+j] = new G4LogicalVolume(scintSolid1, m_material_map["Air"],
						       Form("ScintLV%d",i*num_one_plane+j));//Scinti-->Air
      G4int copyno=i*num_one_plane+j+6;
      if(copyno>31) copyno=copyno-32;
      new G4PVPlacement( transformMP1, Form("ScintPV%d",i*num_one_plane+j), scintLV[i*num_one_plane+j], m_world_pv, FALSE, copyno );
      scintVisAtt[i*num_one_plane+j] = new G4VisAttributes(true, G4Colour(0.,0.8,1.));
      scintLV[i*num_one_plane+j]->SetVisAttributes(scintVisAtt[i*num_one_plane+j]);
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
    TOFLGLV = new G4LogicalVolume(scintSolid1, Scinti, name1);
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

  TPCScintSD* scintSD = new TPCScintSD("/SCINT");
  G4SDManager::GetSDMpointer()->AddNewDetector(scintSD);
  for(G4int i = 0; i<NumOfSegHTOF; i++){
    scintLV[i]->SetSensitiveDetector(scintSD);
  }
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructHypTPC( void )
{
  const G4double RIN_TPC    = gSize.Get( "TpcRin" )*mm*0.5;
  const G4double ROUT_TPC   = gSize.Get( "TpcRout" )*mm*0.5;
  const G4double DZ_TPC     = gSize.Get( "TpcDz" )*mm;
  const G4double DPHI_TPC   = 360.*deg;
  const G4double zz_o[2]    = { -DZ_TPC, DZ_TPC };
  const G4double r_in_o[2]  = { 0.5*sqrt(3.0)*RIN_TPC, 0.5*sqrt(3.0)*RIN_TPC };
  const G4double r_out_o[2] = { ROUT_TPC, ROUT_TPC};
  auto TPCSolid = new G4Polyhedra( "TPCHedra", 22.5*deg,DPHI_TPC + 22.5*deg,
				   8, 2, zz_o, r_in_o, r_out_o );
  m_tpc_lv = new G4LogicalVolume( TPCSolid, m_material_map["P10"], "TPC_LV" );
  auto rotTPC = new G4RotationMatrix;
  rotTPC->rotateX( 90.*deg );
  new G4PVPlacement( rotTPC, G4ThreeVector( 0, 0, 0 ),
		     "TPC_PV", m_tpc_lv, m_world_pv, false, -2 );
  auto TPCVisAtt = new G4VisAttributes( true, G4Colour( 1., 1., 1. ) );
  m_tpc_lv->SetVisAttributes( TPCVisAtt );

  // ==============================================================
  // Field cage
  // ==============================================================
  const G4double RIN_FC        = gSize.Get( "TpcRinFieldCage" )*mm*0.5;
  const G4double ROUT_FC       = gSize.Get( "TpcRoutFieldCage" )*mm*0.5;
  const G4double DZ_FC         = gSize.Get( "TpcDzFieldCage" )*mm;
  const G4double DPHI_FC       = 360.*deg;
  const G4double zz_o_fc[2]    = { -DZ_FC, DZ_FC };
  const G4double r_in_o_fc[2]  = { RIN_FC, RIN_FC };
  const G4double r_out_o_fc[2] = { ROUT_FC, ROUT_FC };
  auto FCSolid = new G4Polyhedra( "FCHedra", 22.5*deg, DPHI_FC + 22.5*deg,
				  8, 2, zz_o_fc, r_in_o_fc, r_out_o_fc );
  auto FCLV = new G4LogicalVolume( FCSolid, m_material_map["P10"], "FC_LV" );
  auto rotFC = new G4RotationMatrix;
  rotFC->rotateX( 90.*deg );
  new G4PVPlacement( rotFC, G4ThreeVector( 0, 0, 0 ),
		     "FC_PV", FCLV, m_world_pv, false, -2 );
  auto FCVisAtt = new G4VisAttributes( true, G4Colour( 1., 0., 0. ) );
  FCLV->SetVisAttributes( FCVisAtt );

  // ==============================================================
  // Virtual pads
  // ==============================================================
  char name5[30];
  char name6[30];
  char name7[30];
  G4Tubs* padSolid[50];
  G4LogicalVolume* padLV[50];
  //  G4VPhysicalVolume* padPV[40];

  G4double angle[40]={0};
  // G4int numpads[40] = {};//number of pads in the layers
  ////we can change very easy the pad structure.
  G4double pad_center=-143.;
  G4double cen_diff=fabs(pad_center)*mm;
  //  G4double cen_diff=abs(-143.)*mm;// pad position is fixed
  G4double pad_in[40]={0};
  G4double pad_out[40]={0};
  G4double tpc_rad=250;

  // out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  G4double pad_length_in = gConf.Get<G4double>("PadLengthIn");
  G4double pad_length_out = gConf.Get<G4double>("PadLengthOut"); //gap 1mm
  G4double pad_gap = gConf.Get<G4double>("PadGap");
  G4int pad_in_num = gConf.Get<G4int>("PadNumIn");
  G4int pad_out_num = gConf.Get<G4int>("PadNumOut");
  G4int pad_configure = gConf.Get<G4int>("PadConfigure");
  const int NUM_PAD = gConf.Get<G4int>("PadNumIn") + gConf.Get<G4int>("PadNumOut");

  switch( pad_configure ){
  case 1:
    for( G4int i=0; i<pad_in_num+pad_out_num; ++i ){
      if( i<pad_in_num ){
	pad_in[i] = 10.+(pad_length_in+pad_gap)*i;
	pad_out[i] = 10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i] = 360.;
      } else {
	pad_in[i] = 10.+(pad_length_in+pad_gap)*pad_in_num +
	  (pad_length_out+pad_gap)*(i-pad_in_num);
	pad_out[i] = 10.+(pad_length_in+pad_gap)*pad_in_num +
	  (pad_length_out+pad_gap)*(i-pad_in_num) + pad_length_out;
	angle[i] = 180.-acos((pow(pad_out[i],2)+pow(cen_diff,2)-pow(tpc_rad,2))/(2*pad_out[i]*cen_diff))*180/3.141592654;
      }
    }
    break;
  case 2:
    for( G4int i=0; i<pad_in_num+pad_out_num; ++i ){
      if( i<pad_in_num ){
	pad_in[i] = 10.+(pad_length_in+pad_gap)*i;
	pad_out[i] = 10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i] = 360.;
	// if( i==0 ){
	//   numpads[i]=48.;
	// } else if( i<pad_in_num ){
	//   numpads[i] = 24.*2.*(i+1.)/2.;
	// }
      }else {
	pad_in[i] = 10.+(pad_length_in+pad_gap)*pad_in_num +
	  (pad_length_out+pad_gap)*(i-pad_in_num);
	pad_out[i] = 10.+(pad_length_in+pad_gap)*pad_in_num +
	  (pad_length_out+pad_gap)*(i-pad_in_num) + pad_length_out;
      }
    }
    angle[10]=180.-155.35;
    angle[11]=180.-144.8;
    angle[12]=180.-138.;
    angle[13]=180.-116.73;
    angle[14]=180.-106.;
    angle[15]=180.-98.77;
    angle[16]=180.-94.29;
    angle[17]=180.-89.8;
    angle[18]=180.-87.18;
    angle[19]=180.-84.16;
    angle[20]=180.-81.48;
    angle[21]=180.-73.39;
    angle[22]=180.-65.51011;
    angle[23]=180.-60.19;
    angle[24]=180.-56.35239;
    angle[25]=180.-52.85;
    angle[26]=180.-50.14;
    angle[27]=180.-47.17;
    angle[28]=180.-41.24;
    angle[29]=180.-29.;
    angle[30]=180.-23.23;
    angle[31]=180.-18.69;
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

  G4VisAttributes* padVisAtt= new G4VisAttributes(false, G4Colour(0.,0.,1.));
  //  G4cout<<"pad install"<<G4endl;
  G4RotationMatrix *rotPad = new G4RotationMatrix();
  // G4RotationMatrix *rotPad_in = new G4RotationMatrix();
  G4double below_target = 0;
  if(m_experiment==42){
    below_target=32.4;
  }else if(m_experiment==45||m_experiment==27){
    below_target=60.4;
  }


  for(G4int i=0;i<pad_in_num;i++ ){
    sprintf(name5, "padSolid[%d]", i);
    sprintf(name6, "PadLV%d", i);
    sprintf(name7, "PadPV%d", i);

    if(pad_out[i]<below_target){
      if(m_experiment==42){
	padSolid[i] = new G4Tubs("TPC pad", pad_in[i]*mm, pad_out[i]*mm,
				 120.*mm, 0., angle[i]*deg);
      }else if(m_experiment==45||m_experiment==27){
	padSolid[i] = new G4Tubs("TPC pad", pad_in[i]*mm, pad_out[i]*mm,
				 200./2.*mm, 0., angle[i]*deg);
      }
      padLV[i]  = new G4LogicalVolume(padSolid[i],m_material_map["P10"],name6);
      //    G4cout<<"1-------------"<<G4endl;

    }else{

      padSolid[i] = new G4Tubs("TPC pad", pad_in[i]*mm, pad_out[i]*mm,
			       275.*mm, 0., angle[i]*deg);
      padLV[i]  = new G4LogicalVolume(padSolid[i],m_material_map["P10"],name6);
    }
    //    padLV[i]  = new G4LogicalVolume(padSolid[i],m_material_map["P10"],name6);
    padLV[i]->SetVisAttributes(padVisAtt);

    if(pad_out[i]<below_target){
      if( m_experiment == 42 ){
	new G4PVPlacement( rotPad, G4ThreeVector( 0., fabs(pad_center)*mm, (-120.-25.)*mm ),
			   padLV[i], name7, m_tpc_lv, true, i );
      } else if( m_experiment == 45 || m_experiment == 27 ){
	new G4PVPlacement( rotPad, G4ThreeVector( 0., fabs(pad_center)*mm, (-200.)*mm ),
			   padLV[i], name7, m_tpc_lv, true, i );
      }
    }else{
      new G4PVPlacement( rotPad, G4ThreeVector( 0., fabs(pad_center)*mm, -25.*mm ),
			 padLV[i], name7, m_tpc_lv, true, i );
    }
  }

  // G4RotationMatrix* rotpad = new G4RotationMatrix;
  G4ThreeVector padpos(0.,fabs(pad_center)*mm,-25.*mm);
  //  rotpad.rotateZ(-90.*deg);
  //  G4VisAttributes* padVisAtt1= new G4VisAttributes(false, G4Colour(1.,0.5,0.));
  G4VisAttributes* padVisAtt1= new G4VisAttributes(false, G4Colour(1.,0.5,0.));

  for(G4int i=pad_in_num;i<(pad_in_num+pad_out_num);i++ ){
    sprintf(name5, "padSolid[%d]", i);
    sprintf(name6, "PadLV%d", i);
    sprintf(name7, "PadPV%d", i);
    padSolid[i] = new G4Tubs("TPC pad",  pad_in[i]*mm, pad_out[i]*mm,
			     275.*mm,(90.+angle[i])*deg, (360.-2.*angle[i])*deg);
    padLV[i]  = new G4LogicalVolume(padSolid[i],m_material_map["P10"],name6);
    padLV[i]->SetVisAttributes(padVisAtt1);
    new G4PVPlacement( rotPad, padpos, padLV[i], name7, m_tpc_lv, true, i );
  }

  //////////////// dead layers
  G4Box* Deadlayer_solid = new G4Box("Deadlayer_solid", 5*mm,250*mm,0.001*mm); // 30 x 10 x 10
  G4LogicalVolume* DeadlayerLV=
    new G4LogicalVolume(Deadlayer_solid, m_material_map["Carbon"], "DeadlayerLV1");
  G4RotationMatrix *rotdead = new G4RotationMatrix();
  rotdead->rotateZ(45.*deg);
  new G4PVPlacement(rotdead, G4ThreeVector(0.,0.*mm,-300.1*mm),
		    DeadlayerLV,"DeadLayer1",m_tpc_lv, true, 0);

  G4RotationMatrix *rotdead1 = new G4RotationMatrix();
  rotdead1->rotateZ(-45.*deg);
  new G4PVPlacement(rotdead1, G4ThreeVector(0.,0.*mm,-300.1*mm),
		    DeadlayerLV,"DeadLayer2",m_tpc_lv, true, 1);

  G4VisAttributes* DeadVisAtt= new G4VisAttributes(true, G4Colour(0.1,0.1,0.));
  DeadlayerLV->SetVisAttributes(DeadVisAtt);

  //////////////// virtual pad
  G4VisAttributes* padVVisAtt= new G4VisAttributes(true, G4Colour(0.,0.,1.));
  G4Tubs* padVSolid[40];
  G4LogicalVolume* padVLV[40];
  for(G4int i=0;i<pad_in_num;i++ ){
    sprintf(name5, "padVSolid[%d]", i);
    sprintf(name6, "PadVLV%d", i);
    sprintf(name7, "PadVPV%d", i);
    padVSolid[i] = new G4Tubs("TPC pad",  pad_in[i]*mm, pad_out[i]*mm,
			      0.5*mm, 0., 360.*deg);
    padVLV[i]  = new G4LogicalVolume(padVSolid[i],m_material_map["P10"],name6);
    padVLV[i]->SetVisAttributes(padVVisAtt);
    new G4PVPlacement( rotPad,
		       G4ThreeVector( 0., fabs(pad_center)*mm, -302.*mm ),
		       padVLV[i], name7, m_tpc_lv, true, 0 );
  }

  G4VisAttributes* padVVisAtt1= new G4VisAttributes(true, G4Colour(1.,0.5,0.));
  for(G4int i=pad_in_num;i<pad_in_num+pad_out_num;i++ ){
    sprintf(name5, "padVSolid[%d]", i);
    sprintf(name6, "PadVLV%d", i);
    sprintf(name7, "PadVPV%d", i);
    padVSolid[i] = new G4Tubs("TPC pad",  pad_in[i]*mm, pad_out[i]*mm,
			      0.5*mm,(90.+angle[i])*deg, (360.-2.*angle[i])*deg);
    padVLV[i]  = new G4LogicalVolume(padVSolid[i],m_material_map["P10"],name6);
    padVLV[i]->SetVisAttributes(padVVisAtt1);
    new G4PVPlacement( rotPad,
		       G4ThreeVector( 0., fabs(pad_center)*mm, -302.*mm ),
		       padVLV[i], name7, m_tpc_lv, true, 0 );
  }

  auto padSD = new TPCPadSD("/TPC");
  G4SDManager::GetSDMpointer()->AddNewDetector(padSD);
  for(G4int i = 0;i<NUM_PAD;i++){
    padLV[i]->SetSensitiveDetector(padSD);
  }
  m_tpc_lv->SetSensitiveDetector(padSD);
}

//_____________________________________________________________________________
void
TPCDetectorConstruction::ConstructKuramaMagnet( void )
{
  //out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  G4double size_MFIELD[ThreeVector::SIZE];
  size_MFIELD[ThreeVector::X] = 500.0*mm;
  size_MFIELD[ThreeVector::Y] = 400.0*mm;
  size_MFIELD[ThreeVector::Z] = 400.0*mm;

  G4double size_COIL1[3];
  size_COIL1[ThreeVector::X] = 900./2*mm;
  size_COIL1[ThreeVector::Y] = 193./2*mm;
  size_COIL1[ThreeVector::Z] = 280./2*mm;

  G4double size_COIL2[3];
  size_COIL2[ThreeVector::X] = 193./2.*mm;
  size_COIL2[ThreeVector::Y] = 58.5*mm;
  size_COIL2[ThreeVector::Z] = 280.0/2.*mm;

  G4double size_COIL3[3];
  size_COIL3[ThreeVector::X] = 193.0/2.*mm;
  size_COIL3[ThreeVector::Y] = 68.5*mm;
  size_COIL3[ThreeVector::Z] = 280.0/2.*mm;

  G4double size_COIL4[3];
  size_COIL4[ThreeVector::X] = 193./2*mm;
  size_COIL4[ThreeVector::Y] = 280./2*mm;
  size_COIL4[ThreeVector::Z] = 740.0/2*mm;

  G4double size_COIL5[3];
  size_COIL5[ThreeVector::X] = 900.0/2*mm;
  size_COIL5[ThreeVector::Y] = 214./2*mm;
  size_COIL5[ThreeVector::Z] = 280./2*mm;

  G4double size_YOKE_UD[3];
  size_YOKE_UD[ThreeVector::X] = 2200.0/2*mm;
  size_YOKE_UD[ThreeVector::Y] = 370.0/2.*mm;
  size_YOKE_UD[ThreeVector::Z] = size_MFIELD[ThreeVector::Z]*mm;

  G4double size_YOKE_LR[3];
  size_YOKE_LR[ThreeVector::X] = 200.0*mm;
  size_YOKE_LR[ThreeVector::Y] = size_MFIELD[ThreeVector::Y]*mm;
  size_YOKE_LR[ThreeVector::Z] = size_MFIELD[ThreeVector::Z]*mm;

  G4double size_UGUARD_UD[3];
  size_UGUARD_UD[ThreeVector::X] = 950.0*mm;
  size_UGUARD_UD[ThreeVector::Y] = 310.0*mm;
  size_UGUARD_UD[ThreeVector::Z] = 50.0*mm;

  G4double size_UGUARD_LR[3];
  size_UGUARD_LR[ThreeVector::X] = 400.0*mm;
  size_UGUARD_LR[ThreeVector::Y] = 150.0*mm;
  size_UGUARD_LR[ThreeVector::Z] = 50.0*mm;

  G4double size_DGUARD_UD[3];
  size_DGUARD_UD[ThreeVector::X] = 800.0*mm;
  size_DGUARD_UD[ThreeVector::Y] = 210.0*mm;
  size_DGUARD_UD[ThreeVector::Z] = 50.0*mm;

  G4double size_DGUARD_LR[3];
  size_DGUARD_LR[ThreeVector::X] = 125.0*mm;
  size_DGUARD_LR[ThreeVector::Y] = size_MFIELD[ThreeVector::Y]+400.0*mm;
  size_DGUARD_LR[ThreeVector::Z] = 50.0*mm;

  G4double pos_MFIELD[3];
  pos_MFIELD[ThreeVector::X] = 0.0*mm; // no shift
  pos_MFIELD[ThreeVector::Y] = 0.0*mm;
  G4double env_kurama_pos_z   = 1500.*mm;
  pos_MFIELD[ThreeVector::Z] = env_kurama_pos_z*mm;
  //  pos_MFIELD[ThreeVector::Z] = (1105.+400.)*mm;
  //  pos_MFIELD[ThreeVector::Z] = 0.*mm;


  G4double pos_COIL1[3];
  pos_COIL1[ThreeVector::X] = pos_MFIELD[ThreeVector::X];
  pos_COIL1[ThreeVector::Y] = size_MFIELD[ThreeVector::Y]+size_YOKE_UD[ThreeVector::Y]*2. -(size_COIL1[ThreeVector::Y]+20.*mm);
  pos_COIL1[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]-size_MFIELD[ThreeVector::Z]-(size_COIL1[ThreeVector::Z]+20.*mm);

  G4double pos_COIL4L[3];
  pos_COIL4L[ThreeVector::X] = pos_MFIELD[ThreeVector::X]+size_MFIELD[ThreeVector::X]+size_COIL4[ThreeVector::X];
  pos_COIL4L[ThreeVector::Y] = size_MFIELD[ThreeVector::Y]/2;
  pos_COIL4L[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z];


  G4double pos_COIL4R[3];
  pos_COIL4R[ThreeVector::X] = pos_MFIELD[ThreeVector::X]-size_MFIELD[ThreeVector::X]-size_COIL4[ThreeVector::X];
  pos_COIL4R[ThreeVector::Y] = size_MFIELD[ThreeVector::Y]/2;
  pos_COIL4R[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z];


  G4double pos_COIL2L[3];
  pos_COIL2L[ThreeVector::X] = pos_MFIELD[ThreeVector::X]+size_MFIELD[ThreeVector::X]+size_COIL4[ThreeVector::X];
  pos_COIL2L[ThreeVector::Y] = (pos_COIL4L[ThreeVector::Y]+pos_COIL1[ThreeVector::Y])/2+(size_COIL4[ThreeVector::Y]-size_COIL1[ThreeVector::Y])/2;
  pos_COIL2L[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]-size_MFIELD[ThreeVector::Z]-(size_COIL1[ThreeVector::Z]+20.*mm);

  G4double pos_COIL2R[3];
  pos_COIL2R[ThreeVector::X] = pos_MFIELD[ThreeVector::X]-size_MFIELD[ThreeVector::X]-size_COIL4[ThreeVector::X];
  pos_COIL2R[ThreeVector::Y] = (pos_COIL4R[ThreeVector::Y]+pos_COIL1[ThreeVector::Y])/2+(size_COIL4[ThreeVector::Y]-size_COIL1[ThreeVector::Y])/2;
  pos_COIL2R[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]-size_MFIELD[ThreeVector::Z]-(size_COIL1[ThreeVector::Z]+20.*mm);


  G4double pos_COIL5[3];
  pos_COIL5[ThreeVector::X] = pos_MFIELD[ThreeVector::X];
  pos_COIL5[ThreeVector::Y] = size_MFIELD[ThreeVector::Y]+size_YOKE_UD[ThreeVector::Y]*2. -(size_COIL5[ThreeVector::Y]);
  pos_COIL5[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] + size_MFIELD[ThreeVector::Z]+(size_COIL5[ThreeVector::Z]+21.*mm);


  G4double pos_COIL3L[3];
  pos_COIL3L[ThreeVector::X] = pos_MFIELD[ThreeVector::X]+size_MFIELD[ThreeVector::X]+size_COIL4[ThreeVector::X];
  pos_COIL3L[ThreeVector::Y] = (pos_COIL4L[ThreeVector::Y]+pos_COIL5[ThreeVector::Y])/2+(size_COIL4[ThreeVector::Y]-size_COIL5[ThreeVector::Y])/2;
  pos_COIL3L[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]+size_MFIELD[ThreeVector::Z]+(size_COIL5[ThreeVector::Z]+21.*mm);

  G4double pos_COIL3R[3];
  pos_COIL3R[ThreeVector::X] = pos_MFIELD[ThreeVector::X]-size_MFIELD[ThreeVector::X]-size_COIL4[ThreeVector::X];
  pos_COIL3R[ThreeVector::Y] = (pos_COIL4R[ThreeVector::Y]+pos_COIL5[ThreeVector::Y])/2+(size_COIL4[ThreeVector::Y]-size_COIL5[ThreeVector::Y])/2;
  pos_COIL3R[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]+size_MFIELD[ThreeVector::Z]+(size_COIL5[ThreeVector::Z]+21.*mm);



  G4double pos_YOKE_U[3];
  pos_YOKE_U[ThreeVector::X] = pos_MFIELD[ThreeVector::X];
  pos_YOKE_U[ThreeVector::Y] = size_MFIELD[ThreeVector::Y]+size_YOKE_UD[ThreeVector::Y];
  pos_YOKE_U[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z];

  G4double pos_YOKE_D[3];
  pos_YOKE_D[ThreeVector::X] = pos_MFIELD[ThreeVector::X];
  pos_YOKE_D[ThreeVector::Y] = -(size_MFIELD[ThreeVector::Y]+size_YOKE_UD[ThreeVector::Y]);
  pos_YOKE_D[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z];

  G4double pos_YOKE_L[3];
  pos_YOKE_L[ThreeVector::X] = pos_MFIELD[ThreeVector::X] + size_MFIELD[ThreeVector::X]+size_YOKE_LR[ThreeVector::X]+200.*mm;
  pos_YOKE_L[ThreeVector::Y] = 0.0*mm;
  pos_YOKE_L[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z];

  G4double pos_YOKE_R[3];
  pos_YOKE_R[ThreeVector::X] = pos_MFIELD[ThreeVector::X] - (size_MFIELD[ThreeVector::X]+size_YOKE_LR[ThreeVector::X])-200.*mm;
  pos_YOKE_R[ThreeVector::Y] = 0.0*mm;
  pos_YOKE_R[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z];

  ///up guard
  G4double pos_UGUARD_U[3];
  pos_UGUARD_U[ThreeVector::X] = 0.0*mm;
  pos_UGUARD_U[ThreeVector::Y] = size_UGUARD_LR[ThreeVector::Y]+size_UGUARD_UD[ThreeVector::Y];
  pos_UGUARD_U[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] - 820.0*mm + size_UGUARD_UD[ThreeVector::Z];

  G4double pos_UGUARD_D[3];
  pos_UGUARD_D[ThreeVector::X] = 0.0*mm;
  pos_UGUARD_D[ThreeVector::Y] = -(size_UGUARD_LR[ThreeVector::Y]+size_UGUARD_UD[ThreeVector::Y]);
  pos_UGUARD_D[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]- 820.0*mm + size_UGUARD_UD[ThreeVector::Z];

  G4double pos_UGUARD_L[3];
  pos_UGUARD_L[ThreeVector::X] = 0.0*mm + (300.0*mm + size_UGUARD_LR[ThreeVector::X]); // 150 -->gap
  pos_UGUARD_L[ThreeVector::Y] = 0.0*mm;
  pos_UGUARD_L[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]- 820.0*mm + size_UGUARD_LR[ThreeVector::Z];

  G4double pos_UGUARD_R[3];
  pos_UGUARD_R[ThreeVector::X] = 0.0*mm - (300.0*mm + size_UGUARD_LR[ThreeVector::X]);
  pos_UGUARD_R[ThreeVector::Y] = 0.0*mm;
  pos_UGUARD_R[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]- 820.0*mm + size_UGUARD_LR[ThreeVector::Z];

  G4double pos_DGUARD_U[3];
  pos_DGUARD_U[ThreeVector::X] = pos_MFIELD[ThreeVector::X];
  pos_DGUARD_U[ThreeVector::Y] = size_DGUARD_LR[ThreeVector::Y]+size_DGUARD_UD[ThreeVector::Y]-200.;
  pos_DGUARD_U[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] + 820.0*mm - size_DGUARD_UD[ThreeVector::Z];

  G4double pos_DGUARD_D[3];
  pos_DGUARD_D[ThreeVector::X] = pos_MFIELD[ThreeVector::X];
  pos_DGUARD_D[ThreeVector::Y] = -(size_DGUARD_LR[ThreeVector::Y]+size_DGUARD_UD[ThreeVector::Y])+200.;
  pos_DGUARD_D[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]+ 820.0*mm - size_DGUARD_UD[ThreeVector::Z];

  G4double pos_DGUARD_L[3];
  pos_DGUARD_L[ThreeVector::X] = pos_MFIELD[ThreeVector::X] + (size_MFIELD[ThreeVector::X]+300.0*mm + size_DGUARD_LR[ThreeVector::X]);
  pos_DGUARD_L[ThreeVector::Y] = 0.0*mm;
  pos_DGUARD_L[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]+ 820.0*mm - size_DGUARD_LR[ThreeVector::Z];

  G4double pos_DGUARD_R[3];
  pos_DGUARD_R[ThreeVector::X] = pos_MFIELD[ThreeVector::X] - (size_MFIELD[ThreeVector::X]+300.0*mm + size_DGUARD_LR[ThreeVector::X]);
  pos_DGUARD_R[ThreeVector::Y] = 0.0*mm;
  pos_DGUARD_R[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]+ 820.0*mm - size_DGUARD_LR[ThreeVector::Z];

  /////////////
  // Construct KURAMA Magnet
  //////////////coil1U
  G4Box* Coil1_box = new G4Box("Coil1_box",
			       size_COIL1[ThreeVector::X],size_COIL1[ThreeVector::Y],size_COIL1[ThreeVector::Z]);
  G4LogicalVolume*  Coil1_log = new G4LogicalVolume(Coil1_box, m_material_map["Copper"], "Coil1_log",0,0,0);
  Coil1_log->SetVisAttributes( RED );
  // G4double maxStep = 0.*mm;
  // maxStep = 0.00001*mm;
  // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL1[ThreeVector::X],
				    pos_COIL1[ThreeVector::Y],
				    pos_COIL1[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil1_log,
		     "Coil1U_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL1[ThreeVector::X],
				    -pos_COIL1[ThreeVector::Y],
				    pos_COIL1[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil1_log,
		     "Coil1D_phys",
		     m_world_lv,
		     false,
		     0 );

  //////////////coil4RLUD
  G4Box* Coil4_box = new G4Box("Coil4_box",
			       size_COIL4[ThreeVector::X],size_COIL4[ThreeVector::Y],size_COIL4[ThreeVector::Z]);
  G4LogicalVolume*  Coil4_log = new G4LogicalVolume(Coil4_box, m_material_map["Copper"], "Coil4_log",0,0,0);
  Coil4_log->SetVisAttributes( RED );
  // maxStep=0.00001*mm;
  // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL4L[ThreeVector::X],
				    pos_COIL4L[ThreeVector::Y],
				    pos_COIL4L[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil4_log,
		     "Coil4UR_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL4R[ThreeVector::X],
				    pos_COIL4R[ThreeVector::Y],
				    pos_COIL4R[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil4_log,
		     "Coil4UL_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL4L[ThreeVector::X],
				    -pos_COIL4L[ThreeVector::Y],
				    pos_COIL4L[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil4_log,
		     "Coil4DR_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL4R[ThreeVector::X],
				    -pos_COIL4R[ThreeVector::Y],
				    pos_COIL4R[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil4_log,
		     "Coil4UL_phys",
		     m_world_lv,
		     false,
		     0 );

  //////////////coil5UD
  G4Box* Coil5_box = new G4Box("Coil5_box",
			       size_COIL5[ThreeVector::X],size_COIL5[ThreeVector::Y],size_COIL5[ThreeVector::Z]);
  G4LogicalVolume*  Coil5_log = new G4LogicalVolume(Coil5_box, m_material_map["Copper"], "Coil5_log",0,0,0);
  Coil5_log->SetVisAttributes( RED );
  // maxStep=0.00001*mm;
  // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL5[ThreeVector::X],
				    pos_COIL5[ThreeVector::Y],
				    pos_COIL5[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil5_log,
		     "Coil5U_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL5[ThreeVector::X],
				    -pos_COIL5[ThreeVector::Y],
				    pos_COIL5[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil5_log,
		     "Coil5D_phys",
		     m_world_lv,
		     false,
		     0 );

  //    padSolid[i] = new G4Tubs("TPC pad", pad_in[i]*mm, pad_out[i]*mm,
  //			     250.*mm, 0., angle[i]*deg);

  //////////////coil6RLUD
  G4double size_COIL6[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL6[0] = 50.0*mm;
  size_COIL6[1] = size_COIL1[1]*2+size_COIL6[0];
  //  size_COIL6[1] = 330.*mm;
  size_COIL6[2] = 280./2*mm;
  size_COIL6[3] = 90.*deg;

  G4double pos_COIL6LU[3];
  G4double pos_COIL6RU[3];
  G4double pos_COIL6LD[3];
  G4double pos_COIL6RD[3];
  //LU
  pos_COIL6LU[ThreeVector::X] = pos_MFIELD[ThreeVector::X] +size_MFIELD[ThreeVector::X]-size_COIL6[0];
  pos_COIL6LU[ThreeVector::Y] = pos_COIL1[ThreeVector::Y]  -(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6LU[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] - size_MFIELD[ThreeVector::Z]-(size_COIL6[ThreeVector::Z]+21.*mm);
  //RU
  pos_COIL6RU[ThreeVector::X] = pos_MFIELD[ThreeVector::X] -size_MFIELD[ThreeVector::X]+size_COIL6[0];
  pos_COIL6RU[ThreeVector::Y] = pos_COIL1[ThreeVector::Y]  -(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6RU[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] - size_MFIELD[ThreeVector::Z]-(size_COIL6[ThreeVector::Z]+21.*mm);
  //LD
  pos_COIL6LD[ThreeVector::X] = pos_MFIELD[ThreeVector::X]  +size_MFIELD[ThreeVector::X]-size_COIL6[0];
  pos_COIL6LD[ThreeVector::Y] = -pos_COIL1[ThreeVector::Y] +(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6LD[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]  - size_MFIELD[ThreeVector::Z]-(size_COIL6[ThreeVector::Z]+21.*mm);
  //RD
  pos_COIL6RD[ThreeVector::X] = pos_MFIELD[ThreeVector::X]  -size_MFIELD[ThreeVector::X]+size_COIL6[0];
  pos_COIL6RD[ThreeVector::Y] = -pos_COIL1[ThreeVector::Y] +(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6RD[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]  - size_MFIELD[ThreeVector::Z]-(size_COIL6[ThreeVector::Z]+21.*mm);


  G4Tubs* Coil6_tub = new G4Tubs("Coil6_tubs",
				 size_COIL6[0],size_COIL6[1],size_COIL6[2],0.,size_COIL6[3]);
  G4LogicalVolume*  Coil6_log = new G4LogicalVolume(Coil6_tub, m_material_map["Copper"], "Coil6_log",0,0,0);
  Coil6_log->SetVisAttributes( RED );
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
  size_COIL8[1] = size_COIL5[1]*2.+size_COIL8[0];
  size_COIL8[2] = 280./2*mm;
  size_COIL8[3] = 90.*deg;

  G4double pos_COIL8LU[3];
  G4double pos_COIL8RU[3];
  G4double pos_COIL8LD[3];
  G4double pos_COIL8RD[3];
  //LU
  pos_COIL8LU[ThreeVector::X] = pos_MFIELD[ThreeVector::X] +size_MFIELD[ThreeVector::X]-size_COIL8[0];
  pos_COIL8LU[ThreeVector::Y] = pos_COIL5[ThreeVector::Y]  -(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8LU[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] + size_MFIELD[ThreeVector::Z]+(size_COIL8[ThreeVector::Z]+21.*mm);
  //RU
  pos_COIL8RU[ThreeVector::X] = pos_MFIELD[ThreeVector::X] -size_MFIELD[ThreeVector::X]+size_COIL8[0];
  pos_COIL8RU[ThreeVector::Y] = pos_COIL5[ThreeVector::Y]  -(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8RU[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] + size_MFIELD[ThreeVector::Z]+(size_COIL8[ThreeVector::Z]+21.*mm);
  //LD
  pos_COIL8LD[ThreeVector::X] = pos_MFIELD[ThreeVector::X]  +size_MFIELD[ThreeVector::X]-size_COIL8[0];
  pos_COIL8LD[ThreeVector::Y] = -pos_COIL5[ThreeVector::Y] +(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8LD[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]  + size_MFIELD[ThreeVector::Z]+(size_COIL8[ThreeVector::Z]+21.*mm);
  //RD
  pos_COIL8RD[ThreeVector::X] = pos_MFIELD[ThreeVector::X]  -size_MFIELD[ThreeVector::X]+size_COIL8[0];
  pos_COIL8RD[ThreeVector::Y] = -pos_COIL5[ThreeVector::Y] +(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8RD[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]  + size_MFIELD[ThreeVector::Z]+(size_COIL8[ThreeVector::Z]+21.*mm);

  G4Tubs* Coil8_tub = new G4Tubs("Coil8_tubs",
				 size_COIL8[0],size_COIL8[1],size_COIL8[2],0.,size_COIL8[3]);
  G4LogicalVolume*  Coil8_log = new G4LogicalVolume(Coil8_tub, m_material_map["Copper"], "Coil8_log",0,0,0);
  Coil8_log->SetVisAttributes( RED );
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
  size_COIL7[1] = size_COIL4[1]*2.+size_COIL7[0];
  size_COIL7[2] = size_COIL4[0];
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
  // pos_COIL7ULU[ThreeVector::X] = pos_MFIELD[ThreeVector::X] +size_MFIELD[ThreeVector::X]+size_COIL7[0];
  pos_COIL7ULU[ThreeVector::X] = pos_COIL4L[ThreeVector::X];
  pos_COIL7ULU[ThreeVector::Y] = pos_COIL4L[ThreeVector::Y] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7ULU[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] - size_COIL4[ThreeVector::Z];
  //URU
  pos_COIL7URU[ThreeVector::X] = pos_COIL4R[ThreeVector::X];
  pos_COIL7URU[ThreeVector::Y] = pos_COIL4R[ThreeVector::Y] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7URU[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] - size_COIL4[ThreeVector::Z];
  //ULD
  pos_COIL7ULD[ThreeVector::X] = pos_COIL4L[ThreeVector::X];
  pos_COIL7ULD[ThreeVector::Y] = -pos_COIL4L[ThreeVector::Y] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7ULD[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]  - size_COIL4[ThreeVector::Z];
  //URD
  pos_COIL7URD[ThreeVector::X] = pos_COIL4R[ThreeVector::X];
  pos_COIL7URD[ThreeVector::Y] = -pos_COIL4R[ThreeVector::Y] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7URD[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] - size_COIL4[ThreeVector::Z];


  //DLU
  //  pos_COIL7ULU[ThreeVector::X] = pos_MFIELD[ThreeVector::X] +size_MFIELD[ThreeVector::X]+size_COIL7[0];
  pos_COIL7DLU[ThreeVector::X] = pos_COIL4L[ThreeVector::X];
  pos_COIL7DLU[ThreeVector::Y] = pos_COIL4L[ThreeVector::Y] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DLU[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] + size_COIL4[ThreeVector::Z];
  //DRU
  pos_COIL7DRU[ThreeVector::X] = pos_COIL4R[ThreeVector::X];
  pos_COIL7DRU[ThreeVector::Y] = pos_COIL4R[ThreeVector::Y] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DRU[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] + size_COIL4[ThreeVector::Z];
  //DLD
  pos_COIL7DLD[ThreeVector::X] = pos_COIL4L[ThreeVector::X];
  pos_COIL7DLD[ThreeVector::Y] = -pos_COIL4L[ThreeVector::Y] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DLD[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z]  + size_COIL4[ThreeVector::Z];
  //DRD
  pos_COIL7DRD[ThreeVector::X] = pos_COIL4R[ThreeVector::X];
  pos_COIL7DRD[ThreeVector::Y] = -pos_COIL4R[ThreeVector::Y] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DRD[ThreeVector::Z] = pos_MFIELD[ThreeVector::Z] + size_COIL4[ThreeVector::Z];


  G4Tubs* Coil7_tub = new G4Tubs("Coil7_tubs",
				 size_COIL7[0],size_COIL7[1],size_COIL7[2],0.,size_COIL7[3]);
  G4LogicalVolume*  Coil7_log = new G4LogicalVolume(Coil7_tub, m_material_map["Copper"], "Coil7_log",0,0,0);
  Coil7_log->SetVisAttributes( RED );
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
				size_COIL2[ThreeVector::X],
				size_COIL2[ThreeVector::Y],
				size_COIL2[ThreeVector::Z] );
  G4LogicalVolume*  Coil2_log = new G4LogicalVolume(Coil2_box, m_material_map["Copper"], "Coil2_log",0,0,0);
  Coil2_log->SetVisAttributes( RED );
  // maxStep=0.00001*mm;
  // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL2L[ThreeVector::X],
				    pos_COIL2L[ThreeVector::Y],
				    pos_COIL2L[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil2_log,
		     "Coil2UL_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL2R[ThreeVector::X],
				    pos_COIL2R[ThreeVector::Y],
				    pos_COIL2R[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil2_log,
		     "Coil2UR_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL2L[ThreeVector::X],
				    -pos_COIL2L[ThreeVector::Y],
				    pos_COIL2L[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil2_log,
		     "Coil2DL_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL2R[ThreeVector::X],
				    -pos_COIL2R[ThreeVector::Y],
				    pos_COIL2R[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil2_log,
		     "Coil2DR_phys",
		     m_world_lv,
		     false,
		     0 );

  ///coil3
  //////////////coil3
  G4Box* Coil3_box = new G4Box("Coil3_box",
			       size_COIL3[ThreeVector::X],size_COIL3[ThreeVector::Y],size_COIL3[ThreeVector::Z]);
  G4LogicalVolume*  Coil3_log = new G4LogicalVolume(Coil3_box, m_material_map["Copper"], "Coil3_log",0,0,0);
  Coil3_log->SetVisAttributes( RED );
  // maxStep=0.00001*mm;
  // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL3L[ThreeVector::X],
				    pos_COIL3L[ThreeVector::Y],
				    pos_COIL3L[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil3_log,
		     "Coil3UL_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL3R[ThreeVector::X],
				    pos_COIL3R[ThreeVector::Y],
				    pos_COIL3R[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil3_log,
		     "Coil3UR_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL3L[ThreeVector::X],
				    -pos_COIL3L[ThreeVector::Y],
				    pos_COIL3L[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil3_log,
		     "Coil3DL_phys",
		     m_world_lv,
		     false,
		     0 );
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_COIL3R[ThreeVector::X],
				    -pos_COIL3R[ThreeVector::Y],
				    pos_COIL3R[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Coil3_log,
		     "Coil3DR_phys",
		     m_world_lv,
		     false,
		     0 );

  //-------------------- Upstraam End Guard
  if(m_experiment!=3 && m_experiment!=42 && m_experiment!=27){
    G4Box* upGuard_UD_box = new G4Box("upGuard_UD_box",
				      size_UGUARD_UD[ThreeVector::X],size_UGUARD_UD[ThreeVector::Y],size_UGUARD_UD[ThreeVector::Z]);
    G4LogicalVolume*  upGuard_U_log = new G4LogicalVolume(upGuard_UD_box, m_material_map["Iron"], "upGuard_U_log",0,0,0);
    upGuard_U_log->SetVisAttributes( BLUE );
    // maxStep=0.00001*mm;
    // upGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));
    new G4PVPlacement( m_rotation_matrix,
		       G4ThreeVector( pos_UGUARD_U[ThreeVector::X],
				      pos_UGUARD_U[ThreeVector::Y],
				      pos_UGUARD_U[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		       upGuard_U_log,
		       "upGuard_U_phys",
		       m_world_lv,
		       false,
		       0 );

    G4LogicalVolume*  upGuard_D_log = new G4LogicalVolume(upGuard_UD_box, m_material_map["Iron"], "upGuard_D_log",0,0,0);
    upGuard_D_log->SetVisAttributes( BLUE );
    // maxStep=0.00001*mm;
    // upGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
    new G4PVPlacement( m_rotation_matrix,
		       G4ThreeVector( pos_UGUARD_D[ThreeVector::X],
				      pos_UGUARD_D[ThreeVector::Y],
				      pos_UGUARD_D[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		       upGuard_D_log,
		       "upGuard_D_phys",
		       m_world_lv,
		       false,
		       0 );

    G4Box* upGuard_LR_box = new G4Box("upGuard_LR_box",
				      size_UGUARD_LR[ThreeVector::X],size_UGUARD_LR[ThreeVector::Y],size_UGUARD_LR[ThreeVector::Z]);
    G4LogicalVolume*  upGuard_L_log = new G4LogicalVolume(upGuard_LR_box, m_material_map["Iron"], "upGuard_L_log",0,0,0);
    upGuard_L_log->SetVisAttributes( BLUE );
    // maxStep=0.00001*mm;
    // upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
    new G4PVPlacement( m_rotation_matrix,
		       G4ThreeVector( pos_UGUARD_L[ThreeVector::X],
				      pos_UGUARD_L[ThreeVector::Y],
				      pos_UGUARD_L[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		       upGuard_L_log,
		       "upGuard_L_phys",
		       m_world_lv,
		       false,
		       0 );

    G4LogicalVolume*  upGuard_R_log = new G4LogicalVolume(upGuard_LR_box, m_material_map["Iron"], "upGuard_R_log",0,0,0);
    upGuard_R_log->SetVisAttributes( BLUE );
    // maxStep=0.00001*mm;
    // upGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
    new G4PVPlacement( m_rotation_matrix,
		       G4ThreeVector( pos_UGUARD_R[ThreeVector::X],
				      pos_UGUARD_R[ThreeVector::Y],
				      pos_UGUARD_R[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		       upGuard_R_log,
		       "upGuard_R_phys",
		       m_world_lv,
		       false,
		       0 );
  }

  //-------------------- Yoke
  G4Box* Yoke_UD_box = new G4Box("Yoke_UD_box",
				 size_YOKE_UD[ThreeVector::X],size_YOKE_UD[ThreeVector::Y],size_YOKE_UD[ThreeVector::Z]);
  G4LogicalVolume*  Yoke_U_log = new G4LogicalVolume(Yoke_UD_box, m_material_map["Iron"], "Yoke_U_log",0,0,0);

  Yoke_U_log->SetVisAttributes( BLUE );
  // maxStep=0.00001*mm;
  // Yoke_U_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_YOKE_U[ThreeVector::X],
				    pos_YOKE_U[ThreeVector::Y],
				    pos_YOKE_U[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Yoke_U_log,
		     "Yoke_U_phys",
		     m_world_lv,
		     false,
		     0 );
  G4LogicalVolume*  Yoke_D_log = new G4LogicalVolume(Yoke_UD_box, m_material_map["Iron"], "Yoke_D_log",0,0,0);
  Yoke_D_log->SetVisAttributes( BLUE );
  // maxStep=0.00001*mm;
  // Yoke_D_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_YOKE_D[ThreeVector::X],
				    pos_YOKE_D[ThreeVector::Y],
				    pos_YOKE_D[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Yoke_D_log,
		     "Yoke_D_phys",
		     m_world_lv,
		     false,
		     0 );
  G4Box* Yoke_LR_box = new G4Box("Yoke_LR_box",
				 size_YOKE_LR[ThreeVector::X],size_YOKE_LR[ThreeVector::Y],size_YOKE_LR[ThreeVector::Z]);
  G4LogicalVolume*  Yoke_L_log = new G4LogicalVolume(Yoke_LR_box, m_material_map["Iron"], "Yoke_L_log",0,0,0);
  Yoke_L_log->SetVisAttributes( BLUE );
  // maxStep=0.00001*mm;
  // Yoke_L_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_YOKE_L[ThreeVector::X],
				    pos_YOKE_L[ThreeVector::Y],
				    pos_YOKE_L[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Yoke_L_log,
		     "Yoke_L_phys",
		     m_world_lv,
		     false,
		     0 );
  G4LogicalVolume*  Yoke_R_log = new G4LogicalVolume(Yoke_LR_box, m_material_map["Iron"], "Yoke_R_log",0,0,0);
  Yoke_R_log->SetVisAttributes( BLUE );
  // maxStep=0.00001*mm;
  // Yoke_R_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_YOKE_R[ThreeVector::X],
				    pos_YOKE_R[ThreeVector::Y],
				    pos_YOKE_R[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     Yoke_R_log,
		     "Yoke_R_phys",
		     m_world_lv,
		     false,
		     0 );

  //-------------------- Downstream End Guard
  G4Box* downGuard_UD_box = new G4Box("downGuard_UD_box",
				      size_DGUARD_UD[ThreeVector::X],size_DGUARD_UD[ThreeVector::Y],size_DGUARD_UD[ThreeVector::Z]);
  G4LogicalVolume*  downGuard_U_log = new G4LogicalVolume(downGuard_UD_box, m_material_map["Iron"], "downGuard_U_log",0,0,0);
  downGuard_U_log->SetVisAttributes( BLUE );
  // maxStep=0.00001*mm;
  // downGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_DGUARD_U[ThreeVector::X],
				    pos_DGUARD_U[ThreeVector::Y],
				    pos_DGUARD_U[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     downGuard_U_log,
		     "downGuard_U_phys",
		     m_world_lv,
		     false,
		     0 );
  G4LogicalVolume*  downGuard_D_log = new G4LogicalVolume(downGuard_UD_box, m_material_map["Iron"], "downGuard_D_log",0,0,0);
  downGuard_D_log->SetVisAttributes( BLUE );
  // maxStep=0.00001*mm;
  // downGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_DGUARD_D[ThreeVector::X],
				    pos_DGUARD_D[ThreeVector::Y],
				    pos_DGUARD_D[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     downGuard_D_log,
		     "downGuard_D_phys",
		     m_world_lv,
		     false,
		     0 );
  G4Box* downGuard_LR_box = new G4Box("downGuard_LR_box",
				      size_DGUARD_LR[ThreeVector::X],size_DGUARD_LR[ThreeVector::Y],size_DGUARD_LR[ThreeVector::Z]);
  G4LogicalVolume*  downGuard_L_log = new G4LogicalVolume(downGuard_LR_box, m_material_map["Iron"], "downGuard_L_log",0,0,0);
  downGuard_L_log->SetVisAttributes( BLUE );
  // maxStep=0.00001*mm;
  // downGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_DGUARD_L[ThreeVector::X],
				    pos_DGUARD_L[ThreeVector::Y],
				    pos_DGUARD_L[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     downGuard_L_log,
		     "downGuard_L_phys",
		     m_world_lv,
		     false,
		     0 );
  G4LogicalVolume*  downGuard_R_log = new G4LogicalVolume(downGuard_LR_box, m_material_map["Iron"], "downGuard_R_log",0,0,0);
  downGuard_R_log->SetVisAttributes( BLUE );
  // maxStep=0.00001*mm;
  // downGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
  new G4PVPlacement( m_rotation_matrix,
		     G4ThreeVector( pos_DGUARD_R[ThreeVector::X],
				    pos_DGUARD_R[ThreeVector::Y],
				    pos_DGUARD_R[ThreeVector::Z] ).rotateY( m_rotation_angle ),
		     downGuard_R_log,
		     "downGuard_R_phys",
		     m_world_lv,
		     false,
		     0 );

  // std::vector<G4LogicalVolume*> fsLV( 10 );
  // G4VisAttributes* fsVisAtt= new G4VisAttributes(true, G4Colour(1.,0.,0.));

  // //////////////////////////////shhwang checking position. hole
  // G4Box* fs_box = new G4Box("fs_box",800.* mm,450.*mm,0.0001*mm);
  // fsLV[0] =
  //   new G4LogicalVolume(fs_box, m_material_map["Air"], "fsLV0");
  // new G4PVPlacement( m_rotation_matrix, G4ThreeVector(0,0.*mm,450.1*mm).rotateY(m_rotation_angle),
  // 		     fsLV[0], "ScintPV32",m_world_lv, true, 32 );
  // fsLV[0]->SetVisAttributes(fsVisAtt);
  // /////////////////////////////////

  // //////////////////////////////shhwang checking position.  ch
  // G4Box* fs_box1 = new G4Box("fs_box1",800.* mm,450.*mm,0.0001*mm);
  // fsLV[1] = new G4LogicalVolume(fs_box1, m_material_map["Air"], "fsLV1");
  // new G4PVPlacement( m_rotation_matrix, G4ThreeVector( sch_pos[ThreeVector::X],
  // 						  0.*mm,
  // 						  sch_pos[ThreeVector::Z] + 0.1*mm ).rotateY( m_rotation_angle ),
  // 		     fsLV[1], "ScintPV33",m_world_lv, true, 33 );
  // fsLV[1]->SetVisAttributes(fsVisAtt);
  // /////////////////////////////////
  // //////////////////////////////shhwang checking position. dc1
  // G4Box* fs_box2 = new G4Box("fs_box2",800.* mm,450.*mm,0.0001*mm);
  // fsLV[2] = new G4LogicalVolume(fs_box2, m_material_map["Air"], "fsLV2");
  // new G4PVPlacement( m_rotation_matrix, G4ThreeVector( 70., 0.*mm, sdc1_pos[2]+0.1+146./2. ).rotateY( m_rotation_angle ),
  // 		     fsLV[2], "ScintPV34",m_world_lv, true, 34 );
  // fsLV[2]->SetVisAttributes(fsVisAtt);
  // /////////////////////////////////
  // //////////////////////////////shhwang checking position. enterance of Kurama spectrometer
  // G4Box* fs_box3 = new G4Box("fs_box3",800.* mm,450.*mm,0.0001*mm);
  // fsLV[3] = new G4LogicalVolume(fs_box3, m_material_map["Air"], "fsLV3");
  // new G4PVPlacement( m_rotation_matrix, G4ThreeVector(pos_MFIELD[0],0.*mm,pos_MFIELD[2]-400.-0.1).rotateY(m_rotation_angle),
  // 		     fsLV[3], "ScintPV35",m_world_lv, true, 35);
  // fsLV[3]->SetVisAttributes(fsVisAtt);
  // /////////////////////////////////
  // //////////////////////////////shhwang checking position. DG
  // G4Box* fs_box4 = new G4Box("fs_box4",1200.* mm,800.*mm,0.0001*mm);
  // fsLV[4] =	new G4LogicalVolume(fs_box4, m_material_map["Air"], "fsLV4");
  // new G4PVPlacement(m_rotation_matrix, G4ThreeVector(pos_MFIELD[0],0.*mm,pos_DGUARD_D[2]+50.1).rotateY(m_rotation_angle),
  // 		    fsLV[4], "ScintPV36",m_world_lv, true, 36);
  // fsLV[4]->SetVisAttributes(fsVisAtt);
  // /////////////////////////////////
  // //////////////////////////////shhwang checking position. DC2
  // G4Box* fs_box5 = new G4Box("fs_box5",1500.* mm,800.*mm,0.0001*mm);
  // fsLV[5] =	new G4LogicalVolume(fs_box5, m_material_map["Air"], "fsLV5");
  // new G4PVPlacement(m_rotation_matrix, G4ThreeVector(sdc2_pos[0],0.*mm,(sdc2_pos[2]+50.1)*mm).rotateY(m_rotation_angle),
  // 		    fsLV[5], "ScintPV37",m_world_lv, true, 37);
  // fsLV[5]->SetVisAttributes(fsVisAtt);
  // /////////////////////////////////
  // //////////////////////////////shhwang checking position. DC3
  // G4Box* fs_box6 = new G4Box("fs_box6",1500.* mm,800.*mm,0.0001*mm);
  // fsLV[6] =	new G4LogicalVolume(fs_box6, m_material_map["Air"], "fsLV6");
  // new G4PVPlacement(m_rotation_matrix, G4ThreeVector(sdc3_pos[0],0.*mm,sdc3_pos[2]+150./2.+0.1).rotateY(m_rotation_angle),
  // 		    fsLV[6], "ScintPV38",m_world_lv, true, 38);
  // fsLV[6]->SetVisAttributes(fsVisAtt);
  // /////////////////////////////////
  // //////////////////////////////shhwang checking position. ftof
  // G4RotationMatrix* rot_fs = new G4RotationMatrix();
  // rot_fs->rotateY(-tof_ra2 - m_rotation_angle );

  // G4Box* fs_box7 = new G4Box("fs_box7",1500.* mm,1000.*mm,0.0001*mm);
  // fsLV[7] =	new G4LogicalVolume(fs_box7, m_material_map["Air"], "fsLV7");
  // new G4PVPlacement( rot_fs, G4ThreeVector(tof_pos[0],0.*mm,tof_pos[2]+35.0).rotateY(m_rotation_angle),
  // 		     fsLV[7], "ScintPV39",m_world_lv, true, 39);
  // fsLV[7]->SetVisAttributes(fsVisAtt);
  // /////////////////////////////////
  // for( auto&& lv : fsLV ){
  //   if( lv )
  //     lv->SetSensitiveDetector(scintSD);
  // }
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

    const G4double DZ_NBAR1 = (n_bar_radius-DX_NBAR1)*tan(22.5*deg)*2/num_one_plane_nbar/2*mm;
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
	new G4PVPlacement( transformMP1, Form("NBARPV%d",i*num_one_plane_nbar+j), nbarLV[i*num_one_plane_nbar+j], m_world_pv, FALSE, copyno );
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
    const G4double DZ_PVC1 = (340.-DX_PVC1)*tan(22.5*deg)/2*mm;
    const G4double DY_PVC1 = 440.0*mm;
    // const G4double DX_PVC2 = 10.0/2.*mm;
    const G4double DZ_PVC2 = (340.-DX_PVC1)*tan(22.5*deg)/2*mm;
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
      new G4PVPlacement( transformMP1_pvc, Form("PVCPV%d",k*2), pvcLV[k*2], m_world_pv, FALSE, k*2 );
      G4Transform3D transformMP2_pvc(*rotMOutP_pvc, posMOut2_pvc);
      pvcLV[k*2+1] = new G4LogicalVolume(pvcSolid1, m_material_map["Scintillator"], Form("pvcLV%d",k*2+1));
      new G4PVPlacement( transformMP2_pvc, Form("pvcPV%d",k*2+1), pvcLV[k*2], m_world_pv, FALSE, k*2 );
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
  const auto& sch_pos = gGeom.GetGlobalPosition("SCH");
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

  const auto& sdc1_pos = ( gGeom.GetGlobalPosition("SDC1-V1") +
			   gGeom.GetGlobalPosition("SDC1-U2") ) * 0.5;
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

  const auto& sdc2_pos = ( gGeom.GetGlobalPosition("SDC2-X1") +
			   gGeom.GetGlobalPosition("SDC2-Y2") ) * 0.5;
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

  const auto& sdc3_pos = ( gGeom.GetGlobalPosition("SDC3-Y1") +
			   gGeom.GetGlobalPosition("SDC3-X2") ) * 0.5;
  G4LogicalVolume*    DC3Plane_log[6];
  G4double size_DC3[3];
  size_DC3[ThreeVector::X] = 1900.0*mm;
  size_DC3[ThreeVector::Y] = 1280.0*mm;
  size_DC3[ThreeVector::Z] = 150.*mm;

  G4double size_DC3Plane[3];
  size_DC3Plane[ThreeVector::X] = 20.0*96.0*0.5*mm;
  size_DC3Plane[ThreeVector::Y] = 20.0*64.0*0.5*mm;
  size_DC3Plane[ThreeVector::Z] = 0.0001*mm;

  //--------------DC3
  G4Box* DC3_box = new G4Box("DC3_box",size_DC3[ThreeVector::X]/2,size_DC3[ThreeVector::Y]/2,size_DC3[ThreeVector::Z]/2);
  G4LogicalVolume*  DC3_log = new G4LogicalVolume(DC3_box, m_material_map["Argon"], "DC3_log",0,0,0);
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
  new G4PVPlacement( rotHelm, G4ThreeVector(0,0.*cm,0),
		     "HelmPV", HelmLV, m_world_pv, FALSE, 0 );

  //    new G4PVPlacement(rotHelm, G4ThreeVector(0,+20.*cm,0),
  //		      "HelmPV", HelmLV, m_world_pv, FALSE, 0);

  //  HelmPV[1] =
  //    new G4PVPlacement(rotHelm, G4ThreeVector(0,-20.*cm,0),
  //		      "HelmPV", HelmLV, m_world_pv, FALSE, 0);
  G4VisAttributes* HelmVisAtt= new G4VisAttributes(true, G4Colour(1.,0.5,0));
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

    G4Box* TargetSolid = new G4Box("target", Target_x/2*mm,Target_z/2*mm,Target_y/2*mm); // 30 x 10 x 15
    TargetLV=
      new G4LogicalVolume(TargetSolid, m_material_map["Target"], "TargetLV");
    new G4PVPlacement( 0, G4ThreeVector( 0*mm, target_pos_z, 0.*mm ),
		       TargetLV, "TargetPV",m_tpc_lv, true, 0 );
    G4VisAttributes* TargetVisAtt= new G4VisAttributes(true, G4Colour(1.,0.,0.));
    TargetLV->SetVisAttributes(TargetVisAtt);
    G4Tubs* TargetHolderSolid = new G4Tubs("TargetHolderSolid", 16.*mm, 16.2*mm,
					   155.*mm, 0., 360*deg);
    TargetHolderLV=
      new G4LogicalVolume(TargetHolderSolid, m_material_map["P10"], "TargetHolderLV");

    new G4PVPlacement( 0, G4ThreeVector( 0*mm, target_pos_z, 145.*mm ),
		       TargetHolderLV, "TargetHolderPV",m_tpc_lv, true, 0 );
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
		       TargetLV, "TargetPV",m_tpc_lv, true, 0 );
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
		       TargetHolderLV, "TargetHolderPV",m_tpc_lv, true, 0 );
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

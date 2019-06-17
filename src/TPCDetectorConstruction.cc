// -*- C++ -*-

#include "TPCDetectorConstruction.hh"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

#include "G4Transform3D.hh"

#include "Spectrometer_par.hh"

// For magnetic field
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "TPCField.hh"

// For Sensitive Detector
#include "TPCACSD.hh"
#include "TPCNBARSD.hh"
#include "TPCPadSD.hh"
#include "TPCScintSD.hh"
#include "TPCTargetSD.hh"
#include "TPCFTOFSD.hh"
#include "TPCDCSD.hh"
#include "TPCCHSD.hh"

#include "common.hh"

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
}

//////////////////////////////////////////////////
TPCDetectorConstruction::TPCDetectorConstruction()
//////////////////////////////////////////////////
{
}

///////////////////////////////////////////////////
TPCDetectorConstruction::~TPCDetectorConstruction()
///////////////////////////////////////////////////
{
}

/////////////////////////////////////////////////////////
G4VPhysicalVolume* TPCDetectorConstruction::Construct()
/////////////////////////////////////////////////////////
{
  // ==============================================================
  // elements
  // ==============================================================
  G4double A, Z;
  G4String name, symbol;
  const double inch = 2.54*cm;
  char name1[30],name2[30];
  char name11[30],name22[30];


  A= 1.00794 *g/mole;
  G4Element* elH= new G4Element(name="Hydrogen", symbol="H", Z=1., A);

  A= 12.011 *g/mole;
  G4Element* elC= new G4Element(name="Carbon", symbol="C", Z=6., A);

  A= 14.00674 *g/mole;
  G4Element* elN= new G4Element(name="Nitrogen", symbol="N", Z=7., A);

  A= 15.9994 *g/mole;
  G4Element* elO= new G4Element(name="Oxygen", symbol="O", Z=8., A);

  A= 39.948 *g/mole;
  G4Element* elAr= new G4Element(name="Argon", symbol="Ar", Z=18., A);

  A= 28.0855 *g/mole;
  G4Element* elSi= new G4Element(name="Silicon", symbol="Si", Z=14., A);

  A= 126.90447 *g/mole;
  G4Element* elI= new G4Element(name="Iodine", symbol="I", Z=53., A);

  A= 132.90543 *g/mole;
  G4Element* elCs= new G4Element(name="Cesium", symbol="Cs", Z=55., A);

  A= 22.989768 *g/mole;
  G4Element* elNa= new G4Element(name="Sodium", symbol="Na", Z=11., A);

  // ==============================================================
  // materials
  // ==============================================================
  G4double density, massfraction;
  G4int natoms, nel;

  // temperature of experimental hall is controlled at 20 degree.
  const G4double expTemp= STP_Temperature+20.*kelvin;

  // vacuum
  density= universe_mean_density;
  G4Material* Vacuum= new G4Material(name="Vacuum", density, nel=2);
  Vacuum-> AddElement(elN, .7);
  Vacuum-> AddElement(elO, .3);


  // air
  density= 1.2929e-03 *g/cm3;  // at 20 degree
  G4Material* Air= new G4Material(name="Air", density, nel=3,
                                  kStateGas, expTemp);
  G4double ttt= 75.47+23.20+1.28;
  Air-> AddElement(elN,  massfraction= 75.47/ttt);
  Air-> AddElement(elO,  massfraction= 23.20/ttt);
  Air-> AddElement(elAr, massfraction=  1.28/ttt);


  //----------------Iron
  A = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Iron", Z=26., A, density);

  //----------------Copper
  A = 63.546*g/mole;
  density = 8.96*g/cm3;
  G4Material* Cu = new G4Material(name="Copper", Z=29., A, density);

  //----------------Carbon
  A = 12.0107*g/mole;
  //  density = 2.265*g/cm3; //--> not diamond
  density = 3.53*g/cm3; //--> diamond in wiki
  G4Material* C = new G4Material(name="Carbon", Z=6., A, density);

  //----------------LH2
  A = 1.008*g/mole;
  density = 70.99*mg/cm3; //--> shogun geant4
  G4Material* LH2 = new G4Material(name="LH2", Z=1., A, density);


  //----------------LD2
  A = 2.01410*g/mole;
  density = 166.0*mg/cm3; //--> shogun geant4
  G4Material* LD2 = new G4Material(name="LD2", Z=1., A, density);




  ///////////lets select target material from setenv
  char* env_target_material = getenv("Target_Material");
  //  G4cout<<vars()['env_target_material']<<G4endl;
  G4cout<<"<---------- target material ------------->"<<G4endl;
  G4cout<<"Target material : "<<env_target_material<<G4endl;
  G4cout<<"<---------------------------------------->"<<G4endl;
  G4Material* target_mat;
  if(strcmp(env_target_material, "C")==0){
    //----------------Carbon
    A = 12.0107*g/mole;
    density = 3.34*g/cm3; //--> from  H. Takahashi thesis
    target_mat = new G4Material(name="Carbon", Z=6., A, density);
    G4cout<<"Target material : "<<env_target_material<<G4endl;
  }else if(strcmp(env_target_material, "Cu")==0){
    //----------------Copper
    A = 63.546*g/mole;
    density = 8.96*g/cm3;
    target_mat = new G4Material(name="Copper", Z=29., A, density);
    G4cout<<"Target material : "<<env_target_material<<G4endl;
  }else if(strcmp(env_target_material, "LH2")==0){
    //----------------LH2
    A = 1.008*g/mole;
    density = 70.99*mg/cm3; //--> shogun geant4
    target_mat = new G4Material(name="LH2", Z=1., A, density);
    G4cout<<"Target material : "<<env_target_material<<G4endl;
  }else if(strcmp(env_target_material, "LD2")==0){
    //----------------LD2
    A = 2.0141*g/mole;
    density = 166.0*mg/cm3; //--> shogun geant4
    target_mat = new G4Material(name="LD2", Z=1., A, density);
    G4cout<<"Target material : "<<env_target_material<<G4endl;
  }else {
    G4cout<<"Wrong target material"<<G4endl;
    exit(-1);
  }



  // Ar gas
  A= 39.948 *g/mole;
  const G4double denAr= 1.782e-03 *g/cm3 * STP_Temperature/expTemp;
  G4Material* Ar= new G4Material(name="ArgonGas", Z=18., A, denAr,
                                 kStateGas, expTemp);



  // ethane (C2H6)
  const G4double denEthane= 1.356e-3 *g/cm3 * STP_Temperature/expTemp;
  G4Material* Ethane= new G4Material(name="Ethane", denEthane, nel=2,
                                     kStateGas, expTemp);
  Ethane-> AddElement(elC, natoms=2);
  Ethane-> AddElement(elH, natoms=6);

  // methane (CH4)
  const G4double denMethane= 0.717e-3 *g/cm3 * STP_Temperature/expTemp;
  G4Material* Methane= new G4Material(name="Methane", denMethane, nel=2,
                                     kStateGas, expTemp);
  Methane-> AddElement(elC, natoms=1);
  Methane-> AddElement(elH, natoms=4);

  // Ar(50%) + ethane(50%) mixture
  density=  (denAr+denEthane)/2.;
  G4Material* ArEthane= new G4Material(name="ArEthane", density, nel=2,
                                       kStateGas, expTemp);
  ArEthane-> AddMaterial(Ar, massfraction= denAr/2./density);
  ArEthane-> AddMaterial(Ethane, massfraction= denEthane/2./density);


  // P10 gas Ar(90%) + methane(10%) mixture
  density=  0.9*denAr+0.1*denMethane;
  G4Material* P10= new G4Material(name="P10", density, nel=2,
                                       kStateGas, expTemp);
  P10-> AddMaterial(Ar, massfraction= 0.9*denAr/density);
  P10-> AddMaterial(Methane, massfraction= 0.1*denMethane/density);

  // G10 epoxy glass
  G4int ncomponents;
  density = 1.700*g/cm3;
  G4Material* G10 = new G4Material(name="NemaG10", density, ncomponents=4);
  G10->AddElement(elSi, natoms=1);
  G10->AddElement(elO , natoms=2);
  G10->AddElement(elC , natoms=3);
  G10->AddElement(elH , natoms=3);

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // scintillator (Polystyene(C6H5CH=CH2))
  // implement here
  density = 1.032 *g/cm3;
  G4Material* Scinti = new G4Material(name="Scinti", density, nel=2);
  Scinti->AddElement(elC,natoms=8);
  Scinti->AddElement(elH,natoms=8);

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // CH2 Polyethelene
  //
  density = 0.95 *g/cm3;
  G4Material* CH2 = new G4Material(name="CH2", density, nel=2);
  CH2->AddElement(elC,natoms=1);
  CH2->AddElement(elH,natoms=2);

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // quartz (SiO2, crystalline)
  density= 2.64 *g/cm3;
  G4Material* Quartz= new G4Material(name="Quartz", density, nel= 2);
  Quartz-> AddElement(elSi, natoms=1);
  Quartz-> AddElement(elO,  natoms=2);

  // ==============================================================
  // color
  // ==============================================================
  G4Colour red(1.0, 0.0, 0.0);
  G4Colour green(0.0, 1.0, 0.0);
  G4Colour blue(0.0, 0.0, 1.0);
  G4Colour cyan(0.0, 1.0, 1.0);
  G4Colour aqua(0.247, 0.8, 1.0);
  G4Colour magenta(1.0, 0.0, 1.0);
  G4Colour yellow(1.0, 1.0, 0.0);
  G4Colour gray(0.5, 0.5, 0.5);
  G4Colour lavender(0.901, 0.901, 0.98);
  G4Colour maroon(0.5, 0.0, 0.0);

  // ==============================================================
  // Read position parameters
  // ==============================================================
  /////////////////
  //check experimental number
  G4String env_experiment_num = getenv("Experiment_NUM");
  G4int experiment_num=atoi(env_experiment_num.c_str() );

  G4String env_with_kurama = getenv("With_KURAMA");
  G4int with_kurama=atoi(env_with_kurama.c_str() );
  G4cout<<"--------------------------------------"<<G4endl;
  G4cout<<"With kurama:"<<with_kurama<<G4endl;
  G4cout<<"--------------------------------------"<<G4endl;
  G4RotationMatrix *rotForwardSp = new G4RotationMatrix();
  G4String spectrometer_ang = getenv("Spectrometer_ang");
  G4double fSpectrometerAngle=atof( spectrometer_ang.c_str())*deg;//out side less 100 mm. 10+5*x < 100 mm is pad_in_num

  G4String kurama_pos_z = getenv("Kurama_pos_z");
  G4double env_kurama_pos_z=atof( kurama_pos_z.c_str());//out side less 100 mm. 10+5*x < 100 mm is pad_in_num

  //  G4double fSpectrometerAngle=10*deg;
  rotForwardSp->rotateY(-fSpectrometerAngle);


  ////// first I would like to use MIWA-san's code
  ///////////////////////// read parameter from data file
  G4double pos_DC1[3]={0.}, shift_DC1[3]={0}, rotangle_DC1[3]={0};
  G4double pos_DC2[3]={0.}, shift_DC2[3]={0.}, rotangle_DC2[3]={0.};
  G4double pos_DC3[3]={0.}, shift_DC3[3]={0.}, rotangle_DC3[3]={0.};
  G4double pos_CH[3]={0.}, shift_CH[3]={0.}, rotangle_CH[3]={0.};
  G4double pos_FTOF[3]={0.}, shift_FTOF[3]={0.}, rotangle_FTOF[3]={0.};
  G4double pos_PVAC[3]={0.}, shift_PVAC[3]={0.}, rotangle_PVAC[3]={0.};
  G4double pos_FAC[3]={0.}, shift_FAC[3]={0.}, rotangle_FAC[3]={0.};
  G4double pos_BVAC[3]={0.}, shift_BVAC[3]={0.}, rotangle_BVAC[3]={0.};


  G4int i, id, tokennum;
  G4double orgn[3]={0.},shift[3]={0.},rotangle[3]={0.};
  char dummy[256],cham[5],namedet[10];
  char fnam[30];
  FILE *fchampos;

  if(experiment_num==3){
    sprintf(fnam, "champosE03.d");
  }else  if(experiment_num==7){
    sprintf(fnam, "champosE07.d");
  }else  if(experiment_num==42 || experiment_num==45 || experiment_num ==27){
    if(with_kurama==1){
      sprintf(fnam, "champosE42.d");
    }else  if(with_kurama==0){
      sprintf(fnam, "champosE42_wo_kurama.d");
    }
  }
  /*
  if ((fchampos = fopen(fnam, "r")) == NULL) {
    fprintf(stderr, "Cannot open parameter(posE42) file: %s\n", fnam);
    exit(-1);
  }

  printf("Let's read the parameters of each position ... \n");
  for (i=0;;) {
    //    printf(".");
    if (fgets(dummy,sizeof(dummy),fchampos)==NULL)
      break;
    if (dummy[0] == '#' || dummy[0] == '\n')
      continue;
    // orginal position
    tokennum = sscanf(dummy, "%s %d %s %lf %lf %lf",
		      cham, &id, namedet, &orgn[0], &orgn[1], &orgn[2]);
    printf("<------- %s %d ------->\n",cham,id);
    printf("orgn=(%f,%f,%f)\n",orgn[0], orgn[1], orgn[2]);
    if (tokennum != 6) {
      fprintf(stderr, "ChamPOSinit::invarid data structure\n");
      fprintf(stderr, "%s\n", dummy);
      exit(-1);
    }
    // position shift
    //example "%s %d %s %lf %lf %lf",  DC 1 origin:    11.0000000      1.05000000      782.9000000
    fgets(dummy, sizeof(dummy), fchampos);
    sscanf(dummy, "%s %d %s %lf %lf %lf",
	   cham, &id, namedet, &shift[0], &shift[1], &shift[2]);
    printf("shift=(%f,%f,%f)\n", shift[0], shift[1], shift[2]);
    if (tokennum != 6) {
      fprintf(stderr, "ChamPOSinit::invarid data structure\n");
      fprintf(stderr, "%s\n", dummy);
      exit(-1);
    }

    // rotation
    fgets(dummy, sizeof(dummy), fchampos);
    sscanf(dummy, "%s %d %s %lf %lf %lf",
	   cham, &id, namedet, &rotangle[0], &rotangle[1], &rotangle[2]);
    printf("rotangle=(%f,%f,%f)\n", rotangle[0], rotangle[1], rotangle[2]);

    if (tokennum != 6) {
      fprintf(stderr, "ChamPOSinit::invalid data structure\n");
      fprintf(stderr, "%s\n", dummy);
      exit(-1);
    }

    //      G4cout<<"shhwang:"<<dummy<<G4endl;
    if ( strcmp(cham, "DC") == 0 ) {
      if(id==1){
	for(G4int kk=0;kk<3;kk++){
	  pos_DC1[kk]=orgn[kk];
	  shift_DC1[kk]=shift[kk];
	  rotangle_DC1[kk]=rotangle[kk];
	}
      }
      if(id==2){
	for(G4int kk=0;kk<3;kk++){
	  pos_DC2[kk]=orgn[kk];
	  shift_DC2[kk]=shift[kk];
	  rotangle_DC2[kk]=rotangle[kk];
	}
      }
      if(id==3){
	for(G4int kk=0;kk<3;kk++){
	  pos_DC3[kk]=orgn[kk];
	  shift_DC3[kk]=shift[kk];
	  rotangle_DC3[kk]=rotangle[kk];
	}
      }
      //      G4cout<<"shhwang:"<<pos_DC1<<":"<<shift_DC1<<":"<<G4endl;
      //      G4cout<<"shhwang:"<<id<<G4endl;

      if (id < 1 || id > 4) {
	fprintf(stderr, "ChamPOSinit::no such DC:id=%d\n", id);
	exit(-1);
      }
    } else if ( strcmp(cham, "FTOF") == 0 ) {
      if(id==1){
	for(G4int kk=0;kk<3;kk++){
	  pos_FTOF[kk]=orgn[kk];
	  shift_FTOF[kk]=shift[kk];
	  rotangle_FTOF[kk]=rotangle[kk];
	}
      }
      if (id =! 1) {
	fprintf(stderr, "ChamPOSinit::no such FTOF:id=%d\n", id);
	exit(-1);
      }
    } else if ( strcmp(cham, "CH") == 0 ) {
      if(id==1.){
	for(G4int kk=0;kk<3;kk++){
	  pos_CH[kk]=orgn[kk];
	  shift_CH[kk]=shift[kk];
	  rotangle_CH[kk]=rotangle[kk];
	}
      }
      if (id =! 1 ) {
	fprintf(stderr, "ChamPOSinit::no such CH:id=%d\n", id);
	exit(-1);
      }
    } else if ( strcmp(cham, "PVAC") == 0 ) {
      if(id==1.){
	for(G4int kk=0;kk<3;kk++){
	  pos_PVAC[kk]=orgn[kk];
	  shift_PVAC[kk]=shift[kk];
	  rotangle_PVAC[kk]=rotangle[kk];
	}
      }
      if (id =! 1 ) {
	fprintf(stderr, "ChamPOSinit::no such PVAC:id=%d\n", id);
	exit(-1);
      }
    } else if ( strcmp(cham, "FAC") == 0 ) {
      if(id==1.){
	for(G4int kk=0;kk<3;kk++){
	  pos_FAC[kk]=orgn[kk];
	  shift_FAC[kk]=shift[kk];
	  rotangle_FAC[kk]=rotangle[kk];
	}
      }
      if (id =! 1 ) {
	fprintf(stderr, "ChamPOSinit::no such FAC:id=%d\n", id);
	exit(-1);
      }
    } else if ( strcmp(cham, "BVAC") == 0 ) {
      if(id==1.){
	for(G4int kk=0;kk<3;kk++){
	  pos_BVAC[kk]=orgn[kk];
	  shift_BVAC[kk]=shift[kk];
	  rotangle_BVAC[kk]=rotangle[kk];
	}
      }
      if (id =! 1 ) {
	fprintf(stderr, "ChamPOSinit::no such BVAC:id=%d\n", id);
	exit(-1);
      }
    }
  }

  //  pos_DC1[ZCOORD]=  pos_DC1[ZCOORD];
  //  pos_DC2[ZCOORD]=  pos_DC2[ZCOORD];
  //  pos_DC3[ZCOORD]=  pos_DC3[ZCOORD];
  //  pos_CH[ZCOORD]=  pos_CH[ZCOORD];
  //  pos_FTOF[ZCOORD]=  pos_FTOF[ZCOORD];
  */

  printf("done !\n");
  // fclose(fchampos);
  // ==============================================================
  // geometry
  // ==============================================================

  // ==============================================================
  // Experimental Hall (world)
  // ==============================================================
  //  const G4double R_EXPHALL=   3.0*m;
  //  const G4double DZ_EXPHALL=  10.0*m;
  G4double X_EXPHALL=4.0*m;
  G4double Y_EXPHALL=3.0*m;
  G4double Z_EXPHALL=8.0*m;
  G4Box* expHallSolid=
    new G4Box("EXP_HALL", X_EXPHALL, Y_EXPHALL, Z_EXPHALL);
    //    new G4Box("EXP_HALL", R_EXPHALL, DZ_EXPHALL, 0., 360.*deg);

  G4LogicalVolume* expHallLV=
    new G4LogicalVolume(expHallSolid, Air, "expHallLV");//EXP_HALL_LV");

  // visualization attributes
  G4VisAttributes* expHallVisAtt=
    new G4VisAttributes(false, G4Colour(1., 1., 1.));
  expHallLV-> SetVisAttributes(expHallVisAtt);

  G4PVPlacement* expHall= new G4PVPlacement(0, G4ThreeVector(), "EXP_HALL_PV",
                                            expHallLV, 0, false, 0);

  // ==============================================================
  //  Time projection chamber
  // ==============================================================
  ///check experiment

  // ==============================================================
  // Frame
  // ==============================================================


  //  const G4double ROUT_TPC=  260.*mm;// previous design

  const G4double ROUT_TPC=  574./2.*mm;//real design
  const G4double RIN_TPC=  0.0*mm;
  const G4double DZ_TPC=    302.*mm;
  const G4double DZ_TPC_OFFSET = 0.*mm;
  const G4double DPHI_TPC=  360.*deg;

  const G4double zz_o[2]={-DZ_TPC, DZ_TPC};
  const G4double r_in_o[2]={0.5*sqrt(3.0)*RIN_TPC, 0.5*sqrt(3.0)*RIN_TPC};
  //  const G4double r_out_o[2]={2/sqrt(3.0)*ROUT_TPC, 2/sqrt(3.0)*ROUT_TPC};
  const G4double r_out_o[2]={ROUT_TPC, ROUT_TPC};

  G4Polyhedra* TPCSolid =
    //    new G4Polyhedra("TPCHedra",22.5,DPHI_TPC,8,2,zz_o,r_in_o,r_out_o);
    new G4Polyhedra("TPCHedra",22.5*deg,DPHI_TPC+22.5*deg,8,2,zz_o,r_in_o,r_out_o);
  //  G4Tubs* TPCSolid =
  //    new G4Tubs("TPCHedra",0.,ROUT_TPC,DZ_TPC,0.*deg,DPHI_TPC);
  G4LogicalVolume* TPCLV=
    new G4LogicalVolume(TPCSolid, P10, "TPC_LV");
  G4RotationMatrix *rotTPC = new G4RotationMatrix();
  rotTPC->rotateX(90.*deg);
  //  rotTPC->rotateZ(90.*deg);

  G4PVPlacement* TPCPV =
    new G4PVPlacement(rotTPC, G4ThreeVector(0,0,DZ_TPC_OFFSET),
		      "TPC_PV", TPCLV, expHall, FALSE, -2);
  G4VisAttributes* TPCVisAtt= new G4VisAttributes(true, G4Colour(1.,1.,1.));
  TPCLV-> SetVisAttributes(TPCVisAtt);

  // ==============================================================
  // Field cage
  // ==============================================================


  //  const G4double ROUT_FC=  310.*mm;
  //  const G4double RIN_FC=  260.0*mm;
  //  const G4double ROUT_FC=  (574./2.+60.)*mm;
  const G4double ROUT_FC=  (654./2.)*mm;
  const G4double RIN_FC=  574.0/2.*mm;
  const G4double DZ_FC=    400.*mm;
  const G4double DZ_FC_OFFSET = 0.*mm;
  const G4double DPHI_FC=  360.*deg;

  const G4double zz_o_fc[2]={-DZ_FC, DZ_FC};
  //  const G4double r_in_o_fc[2]={0.5*sqrt(3.0)*RIN_FC, 0.5*sqrt(3.0)*RIN_FC};
  const G4double r_in_o_fc[2]={RIN_FC,RIN_FC};
  const G4double r_out_o_fc[2]={ROUT_FC, ROUT_FC};

  G4Polyhedra* FCSolid =
    //    new G4Polyhedra("FCHedra",22.5,DPHI_FC,8,2,zz_o_fc,r_in_o_fc,r_out_o_fc);
    new G4Polyhedra("FCHedra",22.5*deg,DPHI_FC+22.5*deg,8,2,zz_o_fc,r_in_o_fc,r_out_o_fc);
  //  G4Tubs* FCSolid =
  //    new G4Tubs("FCHedra",0.,ROUT_FC,DZ_FC,0.*deg,DPHI_FC);
  G4LogicalVolume* FCLV=
    new G4LogicalVolume(FCSolid, P10, "FC_LV");
  G4RotationMatrix *rotFC = new G4RotationMatrix();
  rotFC->rotateX(90.*deg);

  G4PVPlacement* FCPV =
    new G4PVPlacement(rotFC, G4ThreeVector(0,0,DZ_FC_OFFSET),
		      "FC_PV", FCLV, expHall, FALSE, -2);
  G4VisAttributes* FCVisAtt= new G4VisAttributes(true, G4Colour(1.,0.0,0.0));
  FCLV-> SetVisAttributes(FCVisAtt);


  // ==============================================================
  // Virtual pads
  // ==============================================================

  ///////////////////

  char name5[30];
  char name6[30];
  char name7[30];
  G4Tubs* padSolid[50];
  G4LogicalVolume* padLV[50];
  //  G4VPhysicalVolume* padPV[40];

  G4double angle[40]={0};
  G4int numpads[40]={0};//number of pads in the layers
  ////we can change very easy the pad structure.
  G4String env_target_pos_z = getenv("Target_Pos_z");
  G4double target_pos_z=atof( env_target_pos_z.c_str() );
  G4double pad_center=-143.;
  G4double cen_diff=fabs(pad_center)*mm;
  //  G4double cen_diff=abs(-143.)*mm;// pad position is fixed
  G4double pad_in[40]={0};
  G4double pad_out[40]={0};
  G4double tpc_rad=250;





  G4String env_pad_length_in = getenv("pad_length_in");
  G4double pad_length_in=atof( env_pad_length_in.c_str() );//out side less 100 mm. 10+5*x < 100 mm is pad_in_num

  G4String env_pad_length_out = getenv("pad_length_out");
  G4double pad_length_out=atof( env_pad_length_out.c_str() );//gap 1mm

  G4String env_pad_gap = getenv("pad_gap");
  G4double pad_gap=atof( env_pad_gap.c_str() );

  G4String env_pad_num_in = getenv("pad_num_in");
  G4int pad_in_num=atoi( env_pad_num_in.c_str() );
  G4String env_pad_num_out = getenv("pad_num_out");
  G4int pad_out_num=atoi( env_pad_num_out.c_str() );

  //  G4cout<<"test :"<<pad_in_num<<":"<<pad_out_num<<G4endl;
  //  G4cout<<"test1 :"<<pad_length_in<<":"<<pad_length_out<<G4endl;

  /*
  G4double pad_length_in=9.;//out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  G4double pad_length_out=13.;//gap 1mm
  G4double pad_gap=1.;

  G4int pad_in_num=10;
  G4int pad_out_num=20;
  */
  G4String env_pad_configure = getenv("pad_configure");
  if( atoi(env_pad_configure.c_str()) ==1 ){
    for(G4int i=0;i<pad_in_num+pad_out_num;i++){
      if(i<pad_in_num){
	pad_in[i]=10.+(pad_length_in+pad_gap)*i;
	pad_out[i]=10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i]=360.;

      }else {
	pad_in[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num);
	pad_out[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num) + pad_length_out;
	angle[i]=180.-acos((pow(pad_out[i],2)+pow(cen_diff,2)-pow(tpc_rad,2))/(2*pad_out[i]*cen_diff))*180/3.141592654;
      }
    }

  }else if( atoi(env_pad_configure.c_str()) ==2 ){
    for(G4int i=0;i<pad_in_num+pad_out_num;i++){
      if(i<pad_in_num){
	pad_in[i]=10.+(pad_length_in+pad_gap)*i;
	pad_out[i]=10.+(pad_length_in+pad_gap)*i+pad_length_in;
	angle[i]=360.;
	if(i==0){
	  numpads[i]=48.;
	}else if(i<pad_in_num){
	  numpads[i]=24.*2.*(i+1.)/2.;
	}
      }else {
	pad_in[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num);
	pad_out[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num) + pad_length_out;
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

	  numpads[10]=208.;
	  numpads[11]=218.;
	  numpads[12]=230.;
	  numpads[13]=214.;
	  numpads[14]=212.;
	  numpads[15]=214.;
	  numpads[16]=220.;
	  numpads[17]=224.;
	  numpads[18]=232.;
	  numpads[19]=238.;
	  numpads[20]=244.;
	  numpads[21]=232.;
	  numpads[22]=218.;
	  numpads[23]=210.;
	  numpads[24]=206.;
	  numpads[25]=202.;
	  numpads[26]=200.;
	  numpads[27]=196.;
	  numpads[28]=178.;
	  numpads[29]=130.;
	  numpads[30]=108.;
	  numpads[31]=90.;
	  //	  G4cout<<pad_in[i]<<G4endl;
	  //	  G4cout<<pad_out[i]<<G4endl;

  }
    //    G4cout<<i<<":"<<angle[i]<<G4endl;
      //  }else if( atoi(env_pad_configure.c_str()) ==2 ){
      //    for(G4int k=0.;k<pad_in_num;k++){

      //  }


  G4VisAttributes* padVisAtt= new G4VisAttributes(false, G4Colour(0.,0.,1.));
  //  G4cout<<"pad install"<<G4endl;
  G4PVPlacement* padPV[50];
  G4RotationMatrix *rotPad = new G4RotationMatrix();
  G4RotationMatrix *rotPad_in = new G4RotationMatrix();
  G4double below_target;
  if(experiment_num==42){
    below_target=32.4;
  }else if(experiment_num==45||experiment_num==27){
    below_target=60.4;
  }


  for(G4int i=0;i<pad_in_num;i++ ){
    sprintf(name5, "padSolid[%d]", i);
    sprintf(name6, "PadLV%d", i);
    sprintf(name7, "PadPV%d", i);

    if(pad_out[i]<below_target){
      if(experiment_num==42){
	padSolid[i] = new G4Tubs("TPC pad", pad_in[i]*mm, pad_out[i]*mm,
				 120.*mm, 0., angle[i]*deg);
      }else if(experiment_num==45||experiment_num==27){
	padSolid[i] = new G4Tubs("TPC pad", pad_in[i]*mm, pad_out[i]*mm,
				 200./2.*mm, 0., angle[i]*deg);
      }
      padLV[i]  = new G4LogicalVolume(padSolid[i],P10,name6);
      //    G4cout<<"1-------------"<<G4endl;

    }else{

      padSolid[i] = new G4Tubs("TPC pad", pad_in[i]*mm, pad_out[i]*mm,
			       275.*mm, 0., angle[i]*deg);
      padLV[i]  = new G4LogicalVolume(padSolid[i],P10,name6);
    }
    //    padLV[i]  = new G4LogicalVolume(padSolid[i],P10,name6);
    padLV[i]-> SetVisAttributes(padVisAtt);

    if(pad_out[i]<below_target){
      if(experiment_num==42){
	padPV[i] =
	  new G4PVPlacement(rotPad, G4ThreeVector(0.,fabs(pad_center)*mm,(-120.-25.)*mm),padLV[i],name7,TPCLV, true, i);


      }else if(experiment_num==45||experiment_num==27){
	padPV[i] =
	  new G4PVPlacement(rotPad, G4ThreeVector(0.,fabs(pad_center)*mm,(-200.)*mm),padLV[i],name7,TPCLV, true, i);
      }
    }else{
      padPV[i] =
	new G4PVPlacement(rotPad, G4ThreeVector(0.,fabs(pad_center)*mm,-25.*mm),padLV[i],name7,TPCLV, true, i);
      //	G4cout<<"333-------------:"<<i<<G4endl;
      //        G4cout<<angle[i]<<G4endl;
      //        G4cout<<numpads[i]<<G4endl;
    }

  }



  G4RotationMatrix* rotpad = new G4RotationMatrix();
  G4ThreeVector padpos(0.,fabs(pad_center)*mm,-25.*mm);
  //  rotpad.rotateZ(-90.*deg);
  //  G4VisAttributes* padVisAtt1= new G4VisAttributes(false, G4Colour(1.,0.5,0.));
  G4VisAttributes* padVisAtt1= new G4VisAttributes(false, G4Colour(1.,0.5,0.));
  //	G4cout<<"1234-------------:"<<G4endl;
	//        G4cout<<angle[i]<<G4endl;
	//        G4cout<<numpads[i]<<G4endl;

  for(G4int i=pad_in_num;i<(pad_in_num+pad_out_num);i++ ){
    sprintf(name5, "padSolid[%d]", i);
    sprintf(name6, "PadLV%d", i);
    sprintf(name7, "PadPV%d", i);
    padSolid[i] = new G4Tubs("TPC pad",  pad_in[i]*mm, pad_out[i]*mm,
			     275.*mm,(90.+angle[i])*deg, (360.-2.*angle[i])*deg);
    padLV[i]  = new G4LogicalVolume(padSolid[i],P10,name6);
    padLV[i]-> SetVisAttributes(padVisAtt1);
    padPV[i] =
      new G4PVPlacement(rotPad, padpos, padLV[i],name7,TPCLV, true, i);


    //	G4cout<<"123-------------:"<<i<<G4endl;
    //        G4cout<<angle[i]<<G4endl;
    //        G4cout<<numpads[i]<<G4endl;


  }

  //////////////// dead layers
  G4Box* Deadlayer_solid = new G4Box("Deadlayer_solid", 5*mm,250*mm,0.001*mm); // 30 x 10 x 10
  G4LogicalVolume* DeadlayerLV=
    new G4LogicalVolume(Deadlayer_solid, C, "DeadlayerLV1");
  G4RotationMatrix *rotdead = new G4RotationMatrix();
  rotdead->rotateZ(45.*deg);
  G4PVPlacement* DeadlayerPV1 =
    new G4PVPlacement(rotdead, G4ThreeVector(0.,0.*mm,-300.1*mm),DeadlayerLV,"DeadLayer1",TPCLV, true, 0);

  G4RotationMatrix *rotdead1 = new G4RotationMatrix();
  rotdead1->rotateZ(-45.*deg);
  G4PVPlacement* DeadlayerPV2 =
    new G4PVPlacement(rotdead1, G4ThreeVector(0.,0.*mm,-300.1*mm),DeadlayerLV,"DeadLayer2",TPCLV, true, 1);

  G4VisAttributes* DeadVisAtt= new G4VisAttributes(true, G4Colour(0.1,0.1,0.));
  DeadlayerLV-> SetVisAttributes(DeadVisAtt);




  //////////////// virtual pad

  G4VisAttributes* padVVisAtt= new G4VisAttributes(true, G4Colour(0.,0.,1.));
  G4Tubs* padVSolid[40];
  G4LogicalVolume* padVLV[40];
  G4PVPlacement* padVPV[40];
  for(G4int i=0;i<pad_in_num;i++ ){
    sprintf(name5, "padVSolid[%d]", i);
    sprintf(name6, "PadVLV%d", i);
    sprintf(name7, "PadVPV%d", i);
    padVSolid[i] = new G4Tubs("TPC pad",  pad_in[i]*mm, pad_out[i]*mm,
			     0.5*mm, 0., 360.*deg);
    padVLV[i]  = new G4LogicalVolume(padVSolid[i],P10,name6);
    padVLV[i]-> SetVisAttributes(padVVisAtt);
    padVPV[i] =
      new G4PVPlacement(rotPad, G4ThreeVector(0.,fabs(pad_center)*mm,-302.*mm),padVLV[i],name7,TPCLV, true, 0);
  }

  G4VisAttributes* padVVisAtt1= new G4VisAttributes(true, G4Colour(1.,0.5,0.));
  for(G4int i=pad_in_num;i<pad_in_num+pad_out_num;i++ ){
    sprintf(name5, "padVSolid[%d]", i);
    sprintf(name6, "PadVLV%d", i);
    sprintf(name7, "PadVPV%d", i);
    padVSolid[i] = new G4Tubs("TPC pad",  pad_in[i]*mm, pad_out[i]*mm,
			      0.5*mm,(90.+angle[i])*deg, (360.-2.*angle[i])*deg);
    padVLV[i]  = new G4LogicalVolume(padVSolid[i],P10,name6);
    padVLV[i]-> SetVisAttributes(padVVisAtt1);
    padVPV[i] =
      new G4PVPlacement(rotPad, G4ThreeVector(0.,fabs(pad_center)*mm,-302.*mm), padVLV[i],name7,TPCLV, true, 0);
  }

  const int NUM_PAD = pad_in_num+pad_out_num;
  ///////////////////

  //  G4int ipad = 0;
  //  printf("%d pads are installed\n",ipad);

  // ==============================================================
  //  Target
  // ==============================================================
  G4LogicalVolume* TargetLV;
  G4PVPlacement* TargetPV ;
  G4LogicalVolume* TargetHolderLV;
  G4PVPlacement* TargetHolderPV;
  if(experiment_num==42){
  //  G4Box* TargetSolid = new G4Box("target", 1.5*cm,0.25*cm,0.5*cm); // 30 x 10 x 5

  G4String env_target_size_x = getenv("Target_Size_x");
  G4String env_target_size_y = getenv("Target_Size_y");
  G4String env_target_size_z = getenv("Target_Size_z");

  G4double Target_x=atof( env_target_size_x.c_str() );
  G4double Target_y=atof( env_target_size_y.c_str() );
  G4double Target_z=atof( env_target_size_z.c_str() );

  char* env_target_material = getenv("Target_Material");
  //  G4cout<<env_target_material<<G4endl;

  G4Box* TargetSolid = new G4Box("target", Target_x/2*mm,Target_z/2*mm,Target_y/2*mm); // 30 x 10 x 15
  TargetLV=
    new G4LogicalVolume(TargetSolid, target_mat, "TargetLV");
    //    new G4LogicalVolume(TargetSolid, *env_target_material, "TargetLV");
  TargetPV =
    new G4PVPlacement(0, G4ThreeVector(0,fabs(target_pos_z)*mm,0.*mm),
                      TargetLV, "TargetPV",TPCLV, true, 0);
  G4VisAttributes* TargetVisAtt= new G4VisAttributes(true, G4Colour(1.,0.,0.));
  TargetLV-> SetVisAttributes(TargetVisAtt);
  /////////////////////////////////////
  ///////// target holder  ////////////
  /////////////////////////////////////
  G4Tubs* TargetHolderSolid = new G4Tubs("TargetHolderSolid", 16.*mm, 16.2*mm,
			     155.*mm, 0., 360*deg);
  TargetHolderLV=
    new G4LogicalVolume(TargetHolderSolid, P10, "TargetHolderLV");

 TargetHolderPV =
    new G4PVPlacement(0, G4ThreeVector(0,fabs(target_pos_z)*mm,145.*mm),
                      TargetHolderLV, "TargetHolderPV",TPCLV, true, 0);
  G4VisAttributes* TargetHolderVisAtt= new G4VisAttributes(true, G4Colour(0.1,0.1,0.));
  TargetHolderLV-> SetVisAttributes(TargetHolderVisAtt);
  }else if(experiment_num==45||experiment_num==27){
  //  G4Box* TargetSolid = new G4Box("target", 1.5*cm,0.25*cm,0.5*cm); // 30 x 10 x 5

    G4String env_target_size_r = getenv("Target_Size_x");//r
    G4String env_target_size_z = getenv("Target_Size_z");//z

  G4double Target_r=atof( env_target_size_r.c_str() );
  G4double Target_z=atof( env_target_size_z.c_str() );

  //  G4cout<<env_target_material<<G4endl;
  G4Tubs* TargetSolid = new G4Tubs("TargetSolid", 0.*mm,
					 Target_r*mm,Target_z*mm, 0., 360*deg);
  //  G4Box* TargetSolid = new G4Box("target", 1.*cm,0.25*cm,1.*cm);
  TargetLV=
    new G4LogicalVolume(TargetSolid, target_mat, "TargetLV");
    //    new G4LogicalVolume(TargetSolid, *env_target_material, "TargetLV");

  TargetPV =
    new G4PVPlacement(0, G4ThreeVector(0,fabs(target_pos_z)*mm,0.*mm),
                      TargetLV, "TargetPV",TPCLV, true, 0);
  G4VisAttributes* TargetVisAtt= new G4VisAttributes(true, G4Colour(1.,0.,0.));
  TargetLV-> SetVisAttributes(TargetVisAtt);
  /////////////////////////////////////
  ///////// target holder  ////////////
  /////////////////////////////////////

  G4Tubs* TargetHolderSolid = new G4Tubs("TargetHolderSolid", (Target_r+15.)*mm, (Target_r+15.2)*mm,
					 200.*mm, 0., 360*deg);
  TargetHolderLV=
    new G4LogicalVolume(TargetHolderSolid, P10, "TargetHolderLV");

 TargetHolderPV =
    new G4PVPlacement(0, G4ThreeVector(0,fabs(target_pos_z)*mm,100.*mm),
                      TargetHolderLV, "TargetHolderPV",TPCLV, true, 0);
  G4VisAttributes* TargetHolderVisAtt= new G4VisAttributes(true, G4Colour(0.1,0.1,0.));
  TargetHolderLV-> SetVisAttributes(TargetHolderVisAtt);
  }


  // ==============================================================
  //  Solenoid Magnet
  // ==============================================================


  // ==============================================================
  // Scintillators
  // ==============================================================
  G4String env_tof_segment = getenv("TOF_Segment");
  G4int tof_segment = atoi( env_tof_segment.c_str() );
  const G4int NPHI_SCINT = tof_segment;
  //  const G4int NPHI_SCINT = 32.;

  //  G4LogicalVolume* scintLV[NPHI_SCINT];
  //  G4VPhysicalVolume* scintPV[NPHI_SCINT];
  //  G4VisAttributes* scintVisAtt[NPHI_SCINT];


  G4LogicalVolume* scintLV[100];
  G4VPhysicalVolume* scintPV[100];
  G4VisAttributes* scintVisAtt[100];

  // ==============================================================
  // Side Scintillators (Outer)
  // ==============================================================
  //thickness 5mm
  //  const G4double R_SCINT =  ROUT_TPC*0.5*sqrt(3.0)+2.5*1.0;
  const G4int num_one_plane=NPHI_SCINT/8.;
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
    G4int m=0;
  for(G4int k=0;k<8; k++){
    for(G4int m=0;m<num_one_plane;m++){
      G4ThreeVector posMOut1(337.*mm,0.*mm,+plane_width/2-DZ_SCINT1-DZ_SCINT1*2*m); //x,z,y??
    G4RotationMatrix* rotMOutP = new G4RotationMatrix;
    G4Transform3D transformMP1(rotMOutP->rotateY(dangleOut*(k)), posMOut1.rotateY(dangleOut*(k)));

      sprintf(name1, "ScintLV%d",k*num_one_plane+m);
      sprintf(name2, "ScintPV%d",k*num_one_plane+m);
      scintLV[k*num_one_plane+m] = new G4LogicalVolume(scintSolid1, Air, name1);//Scinti-->Air
      G4int copyno=k*num_one_plane+m+6;
      if(copyno>31) copyno=copyno-32;
      scintPV[k*num_one_plane+m] = new G4PVPlacement(transformMP1, name2, scintLV[k*num_one_plane+m], expHall, FALSE, copyno);
      //      G4cout<<k*num_one_plane+m<<G4endl;
    //    G4Transform3D transformMP2(*rotMOutP, posMOut1);
    //    sprintf(name11, "ScintLV%d",k*num_one_plane+1);
    //    sprintf(name22, "ScintPV%d",k*num_one_plane+1);
    //    scintLV[k*2+1] = new G4LogicalVolume(scintSolid1, Scinti, name11);
    //    scintPV[k*2+1] = new G4PVPlacement(transformMP2, name22, scintLV[k*2], expHall, FALSE, k*2);
    //    scintVisAtt[k*2]= new G4VisAttributes(true, aqua);
      scintVisAtt[k*num_one_plane+m]= new G4VisAttributes(true, G4Colour(0.,0.8,1.));
      scintLV[k*num_one_plane+m]-> SetVisAttributes(scintVisAtt[k*num_one_plane+m]);
    //    scintVisAtt[k*2+1]= new G4VisAttributes(true, aqua);
      //    scintVisAtt[k*2+1]= new G4VisAttributes(true, G4Colour(0.,0.8,1.));
      //    scintLV[k*2+1]-> SetVisAttributes(scintVisAtt[k*2+1]);
      //      rotMOutP->rotateY(dangleOut);
      //      posMOut1.rotateY(dangleOut*(k));
      //      posMOut1.rotateY(dangleOut);
      //      posMOut1.rotateY(dangleOut);
          }
  }



  // ==============================================================
  // light guide for Scintillators
  // ==============================================================
  const G4int NTOF_LG = 16;
  G4LogicalVolume* TOFLGLV;
  G4VPhysicalVolume* TOFLGPV[NTOF_LG];
  G4VisAttributes* TOFLGVisAtt;

  /*
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
    TOFLGPV[k*2] = new G4PVPlacement(transformMP1,"TOFLGPV", TOFLGLV, expHall, FALSE, 0);
    G4Transform3D transformMP2(*rotTOFLGP, posTOFLG2);
    TOFLGPV[k*2+1] = new G4PVPlacement(transformMP2,"TOFLGPV", TOFLGLV, expHall, FALSE, 0);
    rotTOFLGP->rotateY(dangleTOFLG);
    posTOFLG1.rotateY(dangleTOFLG);
    posTOFLG2.rotateY(dangleTOFLG);
  }


  TOFLGVisAtt= new G4VisAttributes(true, G4Colour(0.,0.8,0.));
  TOFLGLV-> SetVisAttributes(TOFLGVisAtt);
    //    TOFLGVisAtt[k*2+1]= new G4VisAttributes(true, aqua);
  TOFLGVisAtt= new G4VisAttributes(true, G4Colour(0.,0.8,0.));
  TOFLGLV-> SetVisAttributes(TOFLGVisAtt);

  // ==============================================================
  // end light guide
  // ==============================================================
  */

  // ==============================================================
  // Chrenkov detector for detecting proton
  // ==============================================================

  G4String env_ac_use = getenv("AC_Use");
  G4int ac_use = atoi( env_ac_use.c_str() );

  G4String env_n_bar_use = getenv("N_Bar_Use");
  G4int n_bar_use = atoi( env_n_bar_use.c_str() );

  const G4int NPHI_PVC = 32;

  G4LogicalVolume* pvcLV[NPHI_PVC];
  G4VPhysicalVolume* pvcPV[NPHI_PVC];
  G4VisAttributes* pvcVisAtt[NPHI_PVC];
  if(ac_use==1.){
  const G4double DX_PVC1 = 10.0/2.*mm;
  const G4double DZ_PVC1 = (340.-DX_PVC1)*tan(22.5*deg)/2*mm;
  const G4double DY_PVC1 = 440.0*mm;

  const G4double DX_PVC2 = 10.0/2.*mm;
  const G4double DZ_PVC2 = (340.-DX_PVC1)*tan(22.5*deg)/2*mm;
  const G4double DY_PVC2 = 440.0*mm;


  G4Box* pvcSolid1= new G4Box("SIDE PVC1", DX_PVC1,
                                     DY_PVC1,  DZ_PVC1);
  G4Box* pvcSolid2= new G4Box("SIDE PVC2", DX_PVC2,
                                     DY_PVC2,  DZ_PVC2);

  // 16 PVC

  const G4double dangleOut_pvc = 22.5*2*deg;
  G4ThreeVector posMOut1_pvc(340.*mm,0.*mm,(-DZ_PVC2)); //x,z,y??
  G4ThreeVector posMOut2_pvc(340.*mm,0.*mm,(+DZ_PVC2)); //x,z,y??
  G4RotationMatrix* rotMOutP_pvc = new G4RotationMatrix;
  rotMOutP_pvc->rotateY(dangleOut_pvc*0.5-22.5*deg);
  posMOut1_pvc.rotateY(dangleOut_pvc*0.5-22.5*deg);
  posMOut2_pvc.rotateY(dangleOut_pvc*0.5-22.5*deg);
  const G4double dangleOut = 22.5*2*deg;
  for(G4int k=0;k<8; k++){
    G4Transform3D transformMP1_pvc(*rotMOutP_pvc, posMOut1_pvc);
    sprintf(name1, "PVCLV%d",k*2);
    sprintf(name2, "PVCPV%d",k*2);
    pvcLV[k*2] = new G4LogicalVolume(pvcSolid1, Scinti, name1);
    pvcPV[k*2] = new G4PVPlacement(transformMP1_pvc, name2, pvcLV[k*2], expHall, FALSE, k*2);

    G4Transform3D transformMP2_pvc(*rotMOutP_pvc, posMOut2_pvc);
    sprintf(name11, "pvcLV%d",k*2+1);
    sprintf(name22, "pvcPV%d",k*2+1);
    pvcLV[k*2+1] = new G4LogicalVolume(pvcSolid1, Scinti, name11);
    pvcPV[k*2+1] = new G4PVPlacement(transformMP2_pvc, name22, pvcLV[k*2], expHall, FALSE, k*2);


    pvcVisAtt[k*2]= new G4VisAttributes(true, G4Colour(1.,0.0,1.));
    pvcLV[k*2]-> SetVisAttributes(pvcVisAtt[k*2]);
    pvcVisAtt[k*2+1]= new G4VisAttributes(true, G4Colour(1.,0.0,1.));
    pvcLV[k*2+1]-> SetVisAttributes(pvcVisAtt[k*2+1]);

    rotMOutP_pvc->rotateY(dangleOut);
    posMOut1_pvc.rotateY(dangleOut);
    posMOut2_pvc.rotateY(dangleOut);
  }
  }



  ////////////////////////n_bar
  const G4int NPHI_NBAR = 32;
  G4LogicalVolume* nbarLV[NPHI_NBAR];
  G4VPhysicalVolume* nbarPV[NPHI_NBAR];
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

    for(G4int k=0;k<8; k++){
      for(G4int m=0;m<num_one_plane;m++){
	G4ThreeVector posMOut1(n_bar_radius*mm,0.*mm,+plane_width_nbar/2-DZ_NBAR1-DZ_NBAR1*2*m); //x,z,y??
	G4RotationMatrix* rotMOutP = new G4RotationMatrix;
	G4Transform3D transformMP1(rotMOutP->rotateY(dangleOut_nbar*(k)), posMOut1.rotateY(dangleOut*(k)));

      sprintf(name1, "NBARLV%d",k*num_one_plane+m);
      sprintf(name2, "NBARPV%d",k*num_one_plane+m);
      nbarLV[k*num_one_plane+m] = new G4LogicalVolume(nbarSolid1, Scinti, name1);

      /////id 0 is start forward scintillator
      G4int copyno=k*num_one_plane+m+6;
      if(copyno>31) copyno=copyno-32;

      nbarPV[k*num_one_plane+m] = new G4PVPlacement(transformMP1, name2, nbarLV[k*num_one_plane+m], expHall, FALSE, copyno);
      nbarVisAtt[k*num_one_plane+m]= new G4VisAttributes(true, G4Colour(0.,0.8,1.));
      nbarLV[k*num_one_plane+m]-> SetVisAttributes(nbarVisAtt[k*num_one_plane+m]);
      }
    }
  }///n bar end



  // ==============================================================
  // helmholtz coil
  // ==============================================================

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
    new G4LogicalVolume(HelmSolid, Fe, "HelmLV");
  G4RotationMatrix *rotHelm = new G4RotationMatrix();
  rotHelm->rotateX(90.*deg);
  rotHelm->rotateZ(-fSpectrometerAngle);
  G4PVPlacement* HelmPV[2];
  HelmPV[0] =
    new G4PVPlacement(rotHelm, G4ThreeVector(0,0.*cm,0),
    		      "HelmPV", HelmLV, expHall, FALSE, 0);



    //    new G4PVPlacement(rotHelm, G4ThreeVector(0,+20.*cm,0),
    //		      "HelmPV", HelmLV, expHall, FALSE, 0);

  //  HelmPV[1] =
  //    new G4PVPlacement(rotHelm, G4ThreeVector(0,-20.*cm,0),
  //		      "HelmPV", HelmLV, expHall, FALSE, 0);
  G4VisAttributes* HelmVisAtt= new G4VisAttributes(true, G4Colour(1.,0.5,0));
  HelmLV-> SetVisAttributes(HelmVisAtt);




  /*  /////// supporter (from bottom to helmholtz)
  G4Box* HelmSuSolid =
    new G4Box("HelmSuSolid",3.*cm,60.*cm,3*cm);
  G4LogicalVolume* HelmSuLV=
    new G4LogicalVolume(HelmSuSolid, Fe, "HelmSuLV");
  G4PVPlacement* HelmSuPV[2];
  HelmSuPV[0] =
    new G4PVPlacement(0, G4ThreeVector(31*cm,-85.*cm,31*cm),
		      "HelmSuPV", HelmSuLV, expHall, FALSE, 0);
  HelmSuPV[1] =
    new G4PVPlacement(0, G4ThreeVector(-31*cm,-85.*cm,31*cm),
		      "HelmSuPV", HelmSuLV, expHall, FALSE, 0);
  HelmSuPV[2] =
    new G4PVPlacement(0, G4ThreeVector(-31*cm,-85.*cm,-31*cm),
		      "HelmSuPV", HelmSuLV, expHall, FALSE, 0);
  HelmSuPV[3] =
    new G4PVPlacement(0, G4ThreeVector(31*cm,-85.*cm,-31*cm),
		      "HelmSuPV", HelmSuLV, expHall, FALSE, 0);
  G4VisAttributes* HelmSuVisAtt= new G4VisAttributes(true, G4Colour(0.,0.0,1.));
  HelmSuLV-> SetVisAttributes(HelmSuVisAtt);

  /////// supporter (Between helmholtz coils)
  G4Box* HelmSuBeSolid =
    new G4Box("HelmSu1Solid",1.*cm,15.*cm,1*cm);
  G4LogicalVolume* HelmSuBeLV=
    new G4LogicalVolume(HelmSuBeSolid, Fe, "HelmSuBeLV");
  G4PVPlacement* HelmSuBePV[2];
  HelmSuBePV[0] =
    new G4PVPlacement(0, G4ThreeVector(31*cm,0.*cm,31*cm),
		      "HelmSuBePV", HelmSuBeLV, expHall, FALSE, 0);
  HelmSuBePV[1] =
    new G4PVPlacement(0, G4ThreeVector(-31*cm,0.*cm,31*cm),
		      "HelmSuBePV", HelmSuBeLV, expHall, FALSE, 0);
  HelmSuBePV[2] =
    new G4PVPlacement(0, G4ThreeVector(-31*cm,0.*cm,-31*cm),
		      "HelmSuBePV", HelmSuBeLV, expHall, FALSE, 0);
  HelmSuBePV[3] =
    new G4PVPlacement(0, G4ThreeVector(31*cm,0.*cm,-31*cm),
		      "HelmSuBePV", HelmSuBeLV, expHall, FALSE, 0);
  G4VisAttributes* HelmSuBeVisAtt= new G4VisAttributes(true, G4Colour(0.,0.0,1.));
  HelmSuBeLV-> SetVisAttributes(HelmSuBeVisAtt);
  */

  ///for SD manager
 G4LogicalVolume*    DC1Plane_log[6];
 G4LogicalVolume*    DC2Plane_log[6];
 G4LogicalVolume*    DC3Plane_log[6];
 G4LogicalVolume* FTOF_log[FTOFMAX];
 G4LogicalVolume* CH_log[CHMAX];
 G4LogicalVolume* fsLV[10];
 G4PVPlacement* fsPV[10];
 G4VisAttributes* fsVisAtt= new G4VisAttributes(true, G4Colour(1.,0.,0.));



 if(experiment_num==42||experiment_num==27){

  // ==============================================================
  // ==============================================================
  // spectrometrr --> consist of kurama, SC, DCs
  // ==============================================================
  // ==============================================================


    if(with_kurama==1){
  G4String env_kurama_gap=getenv("Kurama_gap");
  G4double kurama_gap=atof( env_kurama_gap.c_str());

  G4double size_MFIELD[3];
  size_MFIELD[XCOORD] = 500.0*mm;
  size_MFIELD[YCOORD] = kurama_gap*mm;
  size_MFIELD[ZCOORD] = 400.0*mm;

  G4double size_COIL1[3];
  size_COIL1[XCOORD] = 900./2*mm;
  size_COIL1[YCOORD] = 193./2*mm;
  size_COIL1[ZCOORD] = 280./2*mm;

  G4double size_COIL2[3];
  size_COIL2[XCOORD] = 193./2.*mm;
  size_COIL2[YCOORD] = 58.5*mm;
  size_COIL2[ZCOORD] = 280.0/2.*mm;

  G4double size_COIL3[3];
  size_COIL3[XCOORD] = 193.0/2.*mm;
  size_COIL3[YCOORD] = 68.5*mm;
  size_COIL3[ZCOORD] = 280.0/2.*mm;

  G4double size_COIL4[3];
  size_COIL4[XCOORD] = 193./2*mm;
  size_COIL4[YCOORD] = 280./2*mm;
  size_COIL4[ZCOORD] = 740.0/2*mm;

  G4double size_COIL5[3];
  size_COIL5[XCOORD] = 900.0/2*mm;
  size_COIL5[YCOORD] = 214./2*mm;
  size_COIL5[ZCOORD] = 280./2*mm;



  G4double size_YOKE_UD[3];
  size_YOKE_UD[XCOORD] = 2200.0/2*mm;
  size_YOKE_UD[YCOORD] = 370.0/2.*mm;
  size_YOKE_UD[ZCOORD] = size_MFIELD[ZCOORD]*mm;

  G4double size_YOKE_LR[3];
  size_YOKE_LR[XCOORD] = 200.0*mm;
  size_YOKE_LR[YCOORD] = size_MFIELD[YCOORD]*mm;
  size_YOKE_LR[ZCOORD] = size_MFIELD[ZCOORD]*mm;

  G4double size_UGUARD_UD[3];
  size_UGUARD_UD[XCOORD] = 950.0*mm;
  size_UGUARD_UD[YCOORD] = 310.0*mm;
  size_UGUARD_UD[ZCOORD] = 50.0*mm;

  G4double size_UGUARD_LR[3];
  size_UGUARD_LR[XCOORD] = 400.0*mm;
  size_UGUARD_LR[YCOORD] = 150.0*mm;
  size_UGUARD_LR[ZCOORD] = 50.0*mm;

  G4double size_DGUARD_UD[3];
  size_DGUARD_UD[XCOORD] = 800.0*mm;
  size_DGUARD_UD[YCOORD] = 210.0*mm;
  size_DGUARD_UD[ZCOORD] = 50.0*mm;

  G4double size_DGUARD_LR[3];
  size_DGUARD_LR[XCOORD] = 125.0*mm;
  size_DGUARD_LR[YCOORD] = size_MFIELD[YCOORD]+400.0*mm;
  size_DGUARD_LR[ZCOORD] = 50.0*mm;


  G4double size_CH[3];
  size_CH[XCOORD] = 11.5*mm;
  size_CH[YCOORD] = 400.0*mm;
  size_CH[ZCOORD] = 2.0*mm;

  G4double size_FTOF[3];
  size_FTOF[XCOORD] = 40.0*mm;
  size_FTOF[YCOORD] = 900.0*mm;
  size_FTOF[ZCOORD] = 15.0*mm;

  G4double size_DC1[3];
  //  size_DC1[XCOORD] = 289.0*2.*mm;
  //  size_DC1[YCOORD] = 215.0*2.*mm;
  //  size_DC1[ZCOORD] = 146.0*mm;
  ///Tanida DC1 is to small
  size_DC1[XCOORD] = 389.0*2.*mm;
  size_DC1[YCOORD] = 315.0*2.*mm;
  size_DC1[ZCOORD] = 146.0*mm;


  G4double size_DC1Plane[3];
  size_DC1Plane[XCOORD] = 6.0*127.0*0.5*mm;
  size_DC1Plane[YCOORD] = 6.0*97.0*0.5*mm;
  size_DC1Plane[ZCOORD] = 0.0001*mm;

  G4double size_DC2[3];
  size_DC2[XCOORD] = 1186.5*mm;//# of wires are 128
  size_DC2[YCOORD] = 1186.5*mm;//# of wires are 128
  //size_DC2[ZCOORD] = 45.0*mm;
  size_DC2[ZCOORD] = 100.0*mm;

  G4double size_DC2Plane[3];
  size_DC2Plane[XCOORD] = 9.0*128.0*0.5*mm;
  size_DC2Plane[YCOORD] = 9.0*128.0*0.5*mm;
  size_DC2Plane[ZCOORD] = 0.0001*mm;

  G4double size_DC3[3];
  size_DC3[XCOORD] = 1900.0*mm;
  size_DC3[YCOORD] = 1280.0*mm;
  size_DC3[ZCOORD] = 150.*mm;

  G4double size_DC3Plane[3];
  size_DC3Plane[XCOORD] = 20.0*96.0*0.5*mm;
  size_DC3Plane[YCOORD] = 20.0*64.0*0.5*mm;
  size_DC3Plane[ZCOORD] = 0.0001*mm;

  G4String kurama_move_x=getenv("Kurama_move_x");
  G4double k_move_x=atof( kurama_move_x.c_str());

  G4double pos_MFIELD[3];
  pos_MFIELD[XCOORD] = k_move_x*mm;
  pos_MFIELD[YCOORD] = 0.0*mm;
  pos_MFIELD[ZCOORD] = env_kurama_pos_z*mm;
  //  pos_MFIELD[ZCOORD] = (1105.+400.)*mm;
  //  pos_MFIELD[ZCOORD] = 0.*mm;


  G4double pos_COIL1[3];
  pos_COIL1[XCOORD] = pos_MFIELD[XCOORD];
  pos_COIL1[YCOORD] = size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD]*2. -(size_COIL1[YCOORD]+20.*mm);
  pos_COIL1[ZCOORD] = pos_MFIELD[ZCOORD]-size_MFIELD[ZCOORD]-(size_COIL1[ZCOORD]+20.*mm);

  G4double pos_COIL4L[3];
  pos_COIL4L[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+size_COIL4[XCOORD];
  pos_COIL4L[YCOORD] = size_MFIELD[YCOORD]/2;
  pos_COIL4L[ZCOORD] = pos_MFIELD[ZCOORD];


  G4double pos_COIL4R[3];
  pos_COIL4R[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-size_COIL4[XCOORD];
  pos_COIL4R[YCOORD] = size_MFIELD[YCOORD]/2;
  pos_COIL4R[ZCOORD] = pos_MFIELD[ZCOORD];


  G4double pos_COIL2L[3];
  pos_COIL2L[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+size_COIL4[XCOORD];
  pos_COIL2L[YCOORD] = (pos_COIL4L[YCOORD]+pos_COIL1[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL1[YCOORD])/2;
  pos_COIL2L[ZCOORD] = pos_MFIELD[ZCOORD]-size_MFIELD[ZCOORD]-(size_COIL1[ZCOORD]+20.*mm);

  G4double pos_COIL2R[3];
  pos_COIL2R[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-size_COIL4[XCOORD];
  pos_COIL2R[YCOORD] = (pos_COIL4R[YCOORD]+pos_COIL1[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL1[YCOORD])/2;
  pos_COIL2R[ZCOORD] = pos_MFIELD[ZCOORD]-size_MFIELD[ZCOORD]-(size_COIL1[ZCOORD]+20.*mm);


  G4double pos_COIL5[3];
  pos_COIL5[XCOORD] = pos_MFIELD[XCOORD];
  pos_COIL5[YCOORD] = size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD]*2. -(size_COIL5[YCOORD]);
  pos_COIL5[ZCOORD] = pos_MFIELD[ZCOORD] + size_MFIELD[ZCOORD]+(size_COIL5[ZCOORD]+21.*mm);


  G4double pos_COIL3L[3];
  pos_COIL3L[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+size_COIL4[XCOORD];
  pos_COIL3L[YCOORD] = (pos_COIL4L[YCOORD]+pos_COIL5[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL5[YCOORD])/2;
  pos_COIL3L[ZCOORD] = pos_MFIELD[ZCOORD]+size_MFIELD[ZCOORD]+(size_COIL5[ZCOORD]+21.*mm);

  G4double pos_COIL3R[3];
  pos_COIL3R[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-size_COIL4[XCOORD];
  pos_COIL3R[YCOORD] = (pos_COIL4R[YCOORD]+pos_COIL5[YCOORD])/2+(size_COIL4[YCOORD]-size_COIL5[YCOORD])/2;
  pos_COIL3R[ZCOORD] = pos_MFIELD[ZCOORD]+size_MFIELD[ZCOORD]+(size_COIL5[ZCOORD]+21.*mm);



  G4double pos_YOKE_U[3];
  pos_YOKE_U[XCOORD] = pos_MFIELD[XCOORD];
  pos_YOKE_U[YCOORD] = size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD];
  pos_YOKE_U[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_D[3];
  pos_YOKE_D[XCOORD] = pos_MFIELD[XCOORD];
  pos_YOKE_D[YCOORD] = -(size_MFIELD[YCOORD]+size_YOKE_UD[YCOORD]);
  pos_YOKE_D[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_L[3];
  pos_YOKE_L[XCOORD] = pos_MFIELD[XCOORD] + size_MFIELD[XCOORD]+size_YOKE_LR[XCOORD]+200.*mm;
  pos_YOKE_L[YCOORD] = 0.0*mm;
  pos_YOKE_L[ZCOORD] = pos_MFIELD[ZCOORD];

  G4double pos_YOKE_R[3];
  pos_YOKE_R[XCOORD] = pos_MFIELD[XCOORD] - (size_MFIELD[XCOORD]+size_YOKE_LR[XCOORD])-200.*mm;
  pos_YOKE_R[YCOORD] = 0.0*mm;
  pos_YOKE_R[ZCOORD] = pos_MFIELD[ZCOORD];

  ///up guard
  G4double pos_UGUARD_U[3];
  pos_UGUARD_U[XCOORD] = 0.0*mm;
  pos_UGUARD_U[YCOORD] = size_UGUARD_LR[YCOORD]+size_UGUARD_UD[YCOORD];
  pos_UGUARD_U[ZCOORD] = pos_MFIELD[ZCOORD] - 820.0*mm + size_UGUARD_UD[ZCOORD];

  G4double pos_UGUARD_D[3];
  pos_UGUARD_D[XCOORD] = 0.0*mm;
  pos_UGUARD_D[YCOORD] = -(size_UGUARD_LR[YCOORD]+size_UGUARD_UD[YCOORD]);
  pos_UGUARD_D[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_UD[ZCOORD];

  G4double pos_UGUARD_L[3];
  pos_UGUARD_L[XCOORD] = 0.0*mm + (300.0*mm + size_UGUARD_LR[XCOORD]); // 150 --> gap
  pos_UGUARD_L[YCOORD] = 0.0*mm;
  pos_UGUARD_L[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];

  G4double pos_UGUARD_R[3];
  pos_UGUARD_R[XCOORD] = 0.0*mm - (300.0*mm + size_UGUARD_LR[XCOORD]);
  pos_UGUARD_R[YCOORD] = 0.0*mm;
  pos_UGUARD_R[ZCOORD] = pos_MFIELD[ZCOORD]- 820.0*mm + size_UGUARD_LR[ZCOORD];

  G4double pos_DGUARD_U[3];
  pos_DGUARD_U[XCOORD] = pos_MFIELD[XCOORD];
  pos_DGUARD_U[YCOORD] = size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD]-200.;
  pos_DGUARD_U[ZCOORD] = pos_MFIELD[ZCOORD] + 820.0*mm - size_DGUARD_UD[ZCOORD];

  G4double pos_DGUARD_D[3];
  pos_DGUARD_D[XCOORD] = pos_MFIELD[XCOORD];
  pos_DGUARD_D[YCOORD] = -(size_DGUARD_LR[YCOORD]+size_DGUARD_UD[YCOORD])+200.;
  pos_DGUARD_D[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_UD[ZCOORD];

  G4double pos_DGUARD_L[3];
  pos_DGUARD_L[XCOORD] = pos_MFIELD[XCOORD] + (size_MFIELD[XCOORD]+300.0*mm + size_DGUARD_LR[XCOORD]);
  pos_DGUARD_L[YCOORD] = 0.0*mm;
  pos_DGUARD_L[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_LR[ZCOORD];

  G4double pos_DGUARD_R[3];
  pos_DGUARD_R[XCOORD] = pos_MFIELD[XCOORD] - (size_MFIELD[XCOORD]+300.0*mm + size_DGUARD_LR[XCOORD]);
  pos_DGUARD_R[YCOORD] = 0.0*mm;
  pos_DGUARD_R[ZCOORD] = pos_MFIELD[ZCOORD]+ 820.0*mm - size_DGUARD_LR[ZCOORD];

  //  char name1[30], name2[30];
  G4double maxStep;


  /////////////
  // Construct KURAMA Magnet
  //////////////coil1U
  G4Box* Coil1_box = new G4Box("Coil1_box",
				    size_COIL1[XCOORD],size_COIL1[YCOORD],size_COIL1[ZCOORD]);
  G4LogicalVolume*  Coil1_log = new G4LogicalVolume(Coil1_box, Cu, "Coil1_log",0,0,0);
  Coil1_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4PVPlacement* Coil1U_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL1[XCOORD],pos_COIL1[YCOORD],pos_COIL1[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil1_log,
				     "Coil1U_phys",
				     expHallLV,
				     false,
				     0);
  G4PVPlacement* Coil1D_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL1[XCOORD],-pos_COIL1[YCOORD],pos_COIL1[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil1_log,
				     "Coil1D_phys",
				     expHallLV,
				     false,
				     0);

  //////////////coil4RLUD
  G4Box* Coil4_box = new G4Box("Coil4_box",
				    size_COIL4[XCOORD],size_COIL4[YCOORD],size_COIL4[ZCOORD]);
  G4LogicalVolume*  Coil4_log = new G4LogicalVolume(Coil4_box, Cu, "Coil4_log",0,0,0);
  Coil4_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4PVPlacement* Coil4UR_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL4L[XCOORD],pos_COIL4L[YCOORD],pos_COIL4L[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil4_log,
				     "Coil4UR_phys",
				     expHallLV,
				     false,
				     0);
  G4PVPlacement* Coil4UL_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL4R[XCOORD],pos_COIL4R[YCOORD],pos_COIL4R[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil4_log,
				     "Coil4UL_phys",
				     expHallLV,
				     false,
				     0);
  G4PVPlacement* Coil4DR_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL4L[XCOORD],-pos_COIL4L[YCOORD],pos_COIL4L[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil4_log,
				     "Coil4DR_phys",
				     expHallLV,
				     false,
				     0);
  G4PVPlacement* Coil4DL_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL4R[XCOORD],-pos_COIL4R[YCOORD],pos_COIL4R[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil4_log,
				     "Coil4UL_phys",
				     expHallLV,
				     false,
				     0);



  //////////////coil5UD
  G4Box* Coil5_box = new G4Box("Coil5_box",
				    size_COIL5[XCOORD],size_COIL5[YCOORD],size_COIL5[ZCOORD]);
  G4LogicalVolume*  Coil5_log = new G4LogicalVolume(Coil5_box, Cu, "Coil5_log",0,0,0);
  Coil5_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4PVPlacement* Coil5U_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL5[XCOORD],pos_COIL5[YCOORD],pos_COIL5[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil5_log,
				     "Coil5U_phys",
				     expHallLV,
				     false,
				     0);
  G4PVPlacement* Coil5D_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL5[XCOORD],-pos_COIL5[YCOORD],pos_COIL5[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil5_log,
				     "Coil5D_phys",
				     expHallLV,
				     false,
				     0);

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
  pos_COIL6LU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]-size_COIL6[0];
  pos_COIL6LU[YCOORD] = pos_COIL1[YCOORD]  -(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6LU[ZCOORD] = pos_MFIELD[ZCOORD] - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);
  //RU
  pos_COIL6RU[XCOORD] = pos_MFIELD[XCOORD] -size_MFIELD[XCOORD]+size_COIL6[0];
  pos_COIL6RU[YCOORD] = pos_COIL1[YCOORD]  -(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6RU[ZCOORD] = pos_MFIELD[ZCOORD] - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);
  //LD
  pos_COIL6LD[XCOORD] = pos_MFIELD[XCOORD]  +size_MFIELD[XCOORD]-size_COIL6[0];
  pos_COIL6LD[YCOORD] = -pos_COIL1[YCOORD] +(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6LD[ZCOORD] = pos_MFIELD[ZCOORD]  - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);
  //RD
  pos_COIL6RD[XCOORD] = pos_MFIELD[XCOORD]  -size_MFIELD[XCOORD]+size_COIL6[0];
  pos_COIL6RD[YCOORD] = -pos_COIL1[YCOORD] +(size_COIL6[0]+size_COIL1[1]);
  pos_COIL6RD[ZCOORD] = pos_MFIELD[ZCOORD]  - size_MFIELD[ZCOORD]-(size_COIL6[ZCOORD]+21.*mm);


  G4Tubs* Coil6_tub = new G4Tubs("Coil6_tubs",
				size_COIL6[0],size_COIL6[1],size_COIL6[2],0.,size_COIL6[3]);
  G4LogicalVolume*  Coil6_log = new G4LogicalVolume(Coil6_tub, Cu, "Coil6_log",0,0,0);
  Coil6_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4RotationMatrix *rotcoil6lu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6ru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6ld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6rd = new G4RotationMatrix();
  rotcoil6lu->rotateY(-fSpectrometerAngle);
  rotcoil6ru->rotateY(-fSpectrometerAngle);
  rotcoil6ld->rotateY(-fSpectrometerAngle);
  rotcoil6rd->rotateY(-fSpectrometerAngle);

  rotcoil6lu->rotateZ(0.*deg);
  G4PVPlacement* Coil6LU_phys
    = new G4PVPlacement(rotcoil6lu,
			G4ThreeVector(pos_COIL6LU[XCOORD],
				      pos_COIL6LU[YCOORD],
				      pos_COIL6LU[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil6_log,
			"Coil6LU_phys",
			expHallLV,
			false,
			0);
    rotcoil6ru->rotateZ(-90.*deg);
  G4PVPlacement* Coil6RU_phys
    = new G4PVPlacement(rotcoil6ru,
			G4ThreeVector(pos_COIL6RU[XCOORD],
				      pos_COIL6RU[YCOORD],
				      pos_COIL6RU[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil6_log,
			"Coil6RU_phys",
			expHallLV,
			false,
			0);

  rotcoil6ld->rotateZ(-180.*deg);
  G4PVPlacement* Coil6RD_phys
    = new G4PVPlacement(rotcoil6ld,
			G4ThreeVector(pos_COIL6RD[XCOORD],
				      pos_COIL6RD[YCOORD],
				      pos_COIL6RD[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil6_log,
			"Coil6RD_phys",
			expHallLV,
			false,
			0);
  rotcoil6rd->rotateZ(-270.*deg);
  G4PVPlacement* Coil6LD_phys
    = new G4PVPlacement(rotcoil6rd,
			G4ThreeVector(pos_COIL6LD[XCOORD],
				      pos_COIL6LD[YCOORD],
				      pos_COIL6LD[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil6_log,
			"Coil6LD_phys",
			expHallLV,
			false,
			0);
  //  rotForwardSp->rotateZ(90.*deg);
  //  rotForwardSp->rotateZ(90.*deg);



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
  pos_COIL8LU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]-size_COIL8[0];
  pos_COIL8LU[YCOORD] = pos_COIL5[YCOORD]  -(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8LU[ZCOORD] = pos_MFIELD[ZCOORD] + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);
  //RU
  pos_COIL8RU[XCOORD] = pos_MFIELD[XCOORD] -size_MFIELD[XCOORD]+size_COIL8[0];
  pos_COIL8RU[YCOORD] = pos_COIL5[YCOORD]  -(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8RU[ZCOORD] = pos_MFIELD[ZCOORD] + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);
  //LD
  pos_COIL8LD[XCOORD] = pos_MFIELD[XCOORD]  +size_MFIELD[XCOORD]-size_COIL8[0];
  pos_COIL8LD[YCOORD] = -pos_COIL5[YCOORD] +(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8LD[ZCOORD] = pos_MFIELD[ZCOORD]  + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);
  //RD
  pos_COIL8RD[XCOORD] = pos_MFIELD[XCOORD]  -size_MFIELD[XCOORD]+size_COIL8[0];
  pos_COIL8RD[YCOORD] = -pos_COIL5[YCOORD] +(size_COIL8[0]+size_COIL5[1]);
  pos_COIL8RD[ZCOORD] = pos_MFIELD[ZCOORD]  + size_MFIELD[ZCOORD]+(size_COIL8[ZCOORD]+21.*mm);


  G4Tubs* Coil8_tub = new G4Tubs("Coil8_tubs",
				size_COIL8[0],size_COIL8[1],size_COIL8[2],0.,size_COIL8[3]);
  G4LogicalVolume*  Coil8_log = new G4LogicalVolume(Coil8_tub, Cu, "Coil8_log",0,0,0);
  Coil8_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4RotationMatrix *rotcoil8lu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8ru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8ld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8rd = new G4RotationMatrix();
  rotcoil8lu->rotateY(-fSpectrometerAngle);
  rotcoil8ru->rotateY(-fSpectrometerAngle);
  rotcoil8ld->rotateY(-fSpectrometerAngle);
  rotcoil8rd->rotateY(-fSpectrometerAngle);

  rotcoil8lu->rotateZ(0.*deg);
  G4PVPlacement* Coil8LU_phys
    = new G4PVPlacement(rotcoil8lu,
			G4ThreeVector(pos_COIL8LU[XCOORD],
				      pos_COIL8LU[YCOORD],
				      pos_COIL8LU[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil8_log,
			"Coil8LU_phys",
			expHallLV,
			false,
			0);
    rotcoil8ru->rotateZ(-90.*deg);
  G4PVPlacement* Coil8RU_phys
    = new G4PVPlacement(rotcoil8ru,
			G4ThreeVector(pos_COIL8RU[XCOORD],
				      pos_COIL8RU[YCOORD],
				      pos_COIL8RU[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil8_log,
			"Coil8RU_phys",
			expHallLV,
			false,
			0);

  rotcoil8ld->rotateZ(-180.*deg);
  G4PVPlacement* Coil8RD_phys
    = new G4PVPlacement(rotcoil8ld,
			G4ThreeVector(pos_COIL8RD[XCOORD],
				      pos_COIL8RD[YCOORD],
				      pos_COIL8RD[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil8_log,
			"Coil8RD_phys",
			expHallLV,
			false,
			0);
  rotcoil8rd->rotateZ(-270.*deg);
  G4PVPlacement* Coil8LD_phys
    = new G4PVPlacement(rotcoil8rd,
			G4ThreeVector(pos_COIL8LD[XCOORD],
				      pos_COIL8LD[YCOORD],
				      pos_COIL8LD[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil8_log,
			"Coil8LD_phys",
			expHallLV,
			false,
			0);
  //  rotForwardSp->rotateZ(90.*deg);
  //  rotForwardSp->rotateZ(90.*deg);



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
  //  pos_COIL7ULU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]+size_COIL7[0];
  pos_COIL7ULU[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7ULU[YCOORD] = pos_COIL4L[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7ULU[ZCOORD] = pos_MFIELD[ZCOORD] - size_COIL4[ZCOORD];
  //URU
  pos_COIL7URU[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7URU[YCOORD] = pos_COIL4R[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7URU[ZCOORD] = pos_MFIELD[ZCOORD] - size_COIL4[ZCOORD];
  //ULD
  pos_COIL7ULD[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7ULD[YCOORD] = -pos_COIL4L[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7ULD[ZCOORD] = pos_MFIELD[ZCOORD]  - size_COIL4[ZCOORD];
  //URD
  pos_COIL7URD[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7URD[YCOORD] = -pos_COIL4R[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7URD[ZCOORD] = pos_MFIELD[ZCOORD] - size_COIL4[ZCOORD];


  //DLU
  //  pos_COIL7ULU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]+size_COIL7[0];
  pos_COIL7DLU[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7DLU[YCOORD] = pos_COIL4L[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DLU[ZCOORD] = pos_MFIELD[ZCOORD] + size_COIL4[ZCOORD];
  //DRU
  pos_COIL7DRU[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7DRU[YCOORD] = pos_COIL4R[YCOORD] +(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DRU[ZCOORD] = pos_MFIELD[ZCOORD] + size_COIL4[ZCOORD];
  //DLD
  pos_COIL7DLD[XCOORD] = pos_COIL4L[XCOORD];
  pos_COIL7DLD[YCOORD] = -pos_COIL4L[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DLD[ZCOORD] = pos_MFIELD[ZCOORD]  + size_COIL4[ZCOORD];
  //DRD
  pos_COIL7DRD[XCOORD] = pos_COIL4R[XCOORD];
  pos_COIL7DRD[YCOORD] = -pos_COIL4R[YCOORD] -(size_COIL7[0]+size_COIL4[1]);
  pos_COIL7DRD[ZCOORD] = pos_MFIELD[ZCOORD] + size_COIL4[ZCOORD];


  G4Tubs* Coil7_tub = new G4Tubs("Coil7_tubs",
				size_COIL7[0],size_COIL7[1],size_COIL7[2],0.,size_COIL7[3]);
  G4LogicalVolume*  Coil7_log = new G4LogicalVolume(Coil7_tub, Cu, "Coil7_log",0,0,0);
  Coil7_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4RotationMatrix *rotcoil7ulu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7uru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7uld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7urd = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dlu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7drd = new G4RotationMatrix();
  rotcoil7ulu->rotateY(-fSpectrometerAngle);
  rotcoil7uru->rotateY(-fSpectrometerAngle);
  rotcoil7uld->rotateY(-fSpectrometerAngle);
  rotcoil7urd->rotateY(-fSpectrometerAngle);
  rotcoil7dlu->rotateY(-fSpectrometerAngle);
  rotcoil7dru->rotateY(-fSpectrometerAngle);
  rotcoil7dld->rotateY(-fSpectrometerAngle);
  rotcoil7drd->rotateY(-fSpectrometerAngle);


  rotcoil7ulu->rotateZ(0.*deg);
  rotcoil7ulu->rotateX(180.*deg);
  rotcoil7ulu->rotateY(90.*deg);
  G4PVPlacement* Coil7ULU_phys
    = new G4PVPlacement(rotcoil7ulu,
			G4ThreeVector(pos_COIL7ULU[XCOORD],
				      pos_COIL7ULU[YCOORD],
				      pos_COIL7ULU[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil7_log,
			"Coil7ULU_phys",
			expHallLV,
			false,
			0);


  rotcoil7uru->rotateZ(0.*deg);
  rotcoil7uru->rotateX(180.*deg);
  rotcoil7uru->rotateY(90.*deg);
  G4PVPlacement* Coil7URU_phys
    = new G4PVPlacement(rotcoil7uru,
			G4ThreeVector(pos_COIL7URU[XCOORD],
				      pos_COIL7URU[YCOORD],
				      pos_COIL7URU[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil7_log,
			"Coil7URU_phys",
			expHallLV,
			false,
			0);

  rotcoil7urd->rotateZ(0.*deg);
  rotcoil7urd->rotateX(90.*deg);
  rotcoil7urd->rotateY(90.*deg);
  G4PVPlacement* Coil7URD_phys
    = new G4PVPlacement(rotcoil7urd,
			G4ThreeVector(pos_COIL7URD[XCOORD],
				      pos_COIL7URD[YCOORD],
				      pos_COIL7URD[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil7_log,
			"Coil7URD_phys",
			expHallLV,
			false,
			0);

  rotcoil7uld->rotateZ(0.*deg);
  rotcoil7uld->rotateX(90.*deg);
  rotcoil7uld->rotateY(90.*deg);
  G4PVPlacement* Coil7ULD_phys
    = new G4PVPlacement(rotcoil7uld,
			G4ThreeVector(pos_COIL7ULD[XCOORD],
				      pos_COIL7ULD[YCOORD],
				      pos_COIL7ULD[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil7_log,
			"Coil7ULD_phys",
			expHallLV,
			false,
			0);


  ///down


  rotcoil7dlu->rotateZ(0.*deg);
  rotcoil7dlu->rotateX(-90.*deg);
  rotcoil7dlu->rotateY(90.*deg);
  G4PVPlacement* Coil7DLU_phys
    = new G4PVPlacement(rotcoil7dlu,
			G4ThreeVector(pos_COIL7DLU[XCOORD],
				      pos_COIL7DLU[YCOORD],
				      pos_COIL7DLU[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil7_log,
			"Coil7DLU_phys",
			expHallLV,
			false,
			0);


  rotcoil7dru->rotateZ(0.*deg);
  rotcoil7dru->rotateX(-90.*deg);
  rotcoil7dru->rotateY(90.*deg);
  G4PVPlacement* Coil7DRU_phys
    = new G4PVPlacement(rotcoil7dru,
			G4ThreeVector(pos_COIL7DRU[XCOORD],
				      pos_COIL7DRU[YCOORD],
				      pos_COIL7DRU[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil7_log,
			"Coil7DRU_phys",
			expHallLV,
			false,
			0);

  rotcoil7drd->rotateZ(0.*deg);
  rotcoil7drd->rotateX(0.*deg);
  rotcoil7drd->rotateY(90.*deg);
  G4PVPlacement* Coil7DRD_phys
    = new G4PVPlacement(rotcoil7drd,
			G4ThreeVector(pos_COIL7DRD[XCOORD],
				      pos_COIL7DRD[YCOORD],
				      pos_COIL7DRD[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil7_log,
			"Coil7DRD_phys",
			expHallLV,
			false,
			0);

  rotcoil7dld->rotateZ(0.*deg);
  rotcoil7dld->rotateX(0.*deg);
  rotcoil7dld->rotateY(90.*deg);
  G4PVPlacement* Coil7DLD_phys
    = new G4PVPlacement(rotcoil7dld,
			G4ThreeVector(pos_COIL7DLD[XCOORD],
				      pos_COIL7DLD[YCOORD],
				      pos_COIL7DLD[ZCOORD]).rotateY(fSpectrometerAngle),
			Coil7_log,
			"Coil7DLD_phys",
			expHallLV,
			false,
			0);

  ///coil2
  //////////////coil2
  G4Box* Coil2_box = new G4Box("Coil2_box",
				    size_COIL2[XCOORD],size_COIL2[YCOORD],size_COIL2[ZCOORD]);
  G4LogicalVolume*  Coil2_log = new G4LogicalVolume(Coil2_box, Cu, "Coil2_log",0,0,0);
  Coil2_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4PVPlacement* Coil2UL_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL2L[XCOORD],pos_COIL2L[YCOORD],pos_COIL2L[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil2_log,
				     "Coil2UL_phys",
				     expHallLV,
				     false,
				     0);
  G4PVPlacement* Coil2UR_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL2R[XCOORD],pos_COIL2R[YCOORD],pos_COIL2R[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil2_log,
				     "Coil2UR_phys",
				     expHallLV,
				     false,
				     0);


  G4PVPlacement* Coil2DL_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL2L[XCOORD],-pos_COIL2L[YCOORD],pos_COIL2L[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil2_log,
				     "Coil2DL_phys",
				     expHallLV,
				     false,
				     0);
  G4PVPlacement* Coil2DR_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL2R[XCOORD],-pos_COIL2R[YCOORD],pos_COIL2R[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil2_log,
				     "Coil2DR_phys",
				     expHallLV,
				     false,
				     0);


  ///coil3
  //////////////coil3
  G4Box* Coil3_box = new G4Box("Coil3_box",
				    size_COIL3[XCOORD],size_COIL3[YCOORD],size_COIL3[ZCOORD]);
  G4LogicalVolume*  Coil3_log = new G4LogicalVolume(Coil3_box, Cu, "Coil3_log",0,0,0);
  Coil3_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4PVPlacement* Coil3UL_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL3L[XCOORD],pos_COIL3L[YCOORD],pos_COIL3L[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil3_log,
				     "Coil3UL_phys",
				     expHallLV,
				     false,
				     0);
  G4PVPlacement* Coil3UR_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL3R[XCOORD],pos_COIL3R[YCOORD],pos_COIL3R[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil3_log,
				     "Coil3UR_phys",
				     expHallLV,
				     false,
				     0);


  G4PVPlacement* Coil3DL_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL3L[XCOORD],-pos_COIL3L[YCOORD],pos_COIL3L[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil3_log,
				     "Coil3DL_phys",
				     expHallLV,
				     false,
				     0);
  G4PVPlacement* Coil3DR_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_COIL3R[XCOORD],-pos_COIL3R[YCOORD],pos_COIL3R[ZCOORD]).rotateY(fSpectrometerAngle),
				     Coil3_log,
				     "Coil3DR_phys",
				     expHallLV,
				     false,
				     0);




  //-------------------- Upstraam End Guard
  if(experiment_num!=3 && experiment_num!=42 && experiment_num!=27){
  G4Box* upGuard_UD_box = new G4Box("upGuard_UD_box",
				    size_UGUARD_UD[XCOORD],size_UGUARD_UD[YCOORD],size_UGUARD_UD[ZCOORD]);
  G4LogicalVolume*  upGuard_U_log = new G4LogicalVolume(upGuard_UD_box, Fe, "upGuard_U_log",0,0,0);
   upGuard_U_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //upGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));

  G4PVPlacement*  upGuard_U_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_UGUARD_U[XCOORD],pos_UGUARD_U[YCOORD],pos_UGUARD_U[ZCOORD]).rotateY(fSpectrometerAngle),
				     upGuard_U_log,
				     "upGuard_U_phys",
				     expHallLV,
				     false,
				     0);

G4LogicalVolume*  upGuard_D_log = new G4LogicalVolume(upGuard_UD_box, Fe, "upGuard_D_log",0,0,0);
  upGuard_D_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //upGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));

G4PVPlacement*  upGuard_D_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_UGUARD_D[XCOORD],pos_UGUARD_D[YCOORD],pos_UGUARD_D[ZCOORD]).rotateY(fSpectrometerAngle),
				     upGuard_D_log,
				     "upGuard_D_phys",
				     expHallLV,
				     false,
				     0);

  G4Box* upGuard_LR_box = new G4Box("upGuard_LR_box",
				    size_UGUARD_LR[XCOORD],size_UGUARD_LR[YCOORD],size_UGUARD_LR[ZCOORD]);
G4LogicalVolume*  upGuard_L_log = new G4LogicalVolume(upGuard_LR_box, Fe, "upGuard_L_log",0,0,0);
  upGuard_L_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
G4PVPlacement* upGuard_L_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_UGUARD_L[XCOORD],pos_UGUARD_L[YCOORD],pos_UGUARD_L[ZCOORD]).rotateY(fSpectrometerAngle),
				     upGuard_L_log,
				     "upGuard_L_phys",
				     expHallLV,
				     false,
				     0);

G4LogicalVolume*  upGuard_R_log = new G4LogicalVolume(upGuard_LR_box, Fe, "upGuard_R_log",0,0,0);
  upGuard_R_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //upGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
G4PVPlacement* upGuard_R_phys = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_UGUARD_R[XCOORD],pos_UGUARD_R[YCOORD],pos_UGUARD_R[ZCOORD]).rotateY(fSpectrometerAngle),
				     upGuard_R_log,
				     "upGuard_R_phys",
				     expHallLV,
				     false,
				     0);
  }



  //-------------------- Yoke
  G4Box* Yoke_UD_box = new G4Box("Yoke_UD_box",
				 size_YOKE_UD[XCOORD],size_YOKE_UD[YCOORD],size_YOKE_UD[ZCOORD]);
G4LogicalVolume*  Yoke_U_log = new G4LogicalVolume(Yoke_UD_box, Fe, "Yoke_U_log",0,0,0);

  Yoke_U_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //Yoke_U_log->SetUserLimits(new G4UserLimits(maxStep));
G4PVPlacement* Yoke_U_phys = new G4PVPlacement(rotForwardSp,
				  G4ThreeVector(pos_YOKE_U[XCOORD],pos_YOKE_U[YCOORD],pos_YOKE_U[ZCOORD]).rotateY(fSpectrometerAngle),
				  Yoke_U_log,
				  "Yoke_U_phys",
				  expHallLV,
				  false,
				  0);

G4LogicalVolume*  Yoke_D_log = new G4LogicalVolume(Yoke_UD_box, Fe, "Yoke_D_log",0,0,0);
  Yoke_D_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //Yoke_D_log->SetUserLimits(new G4UserLimits(maxStep));
G4PVPlacement* Yoke_D_phys = new G4PVPlacement(rotForwardSp,
				  G4ThreeVector(pos_YOKE_D[XCOORD],pos_YOKE_D[YCOORD],pos_YOKE_D[ZCOORD]).rotateY(fSpectrometerAngle),
				  Yoke_D_log,
				  "Yoke_D_phys",
				  expHallLV,
				  false,
				  0);

  G4Box* Yoke_LR_box = new G4Box("Yoke_LR_box",
				 size_YOKE_LR[XCOORD],size_YOKE_LR[YCOORD],size_YOKE_LR[ZCOORD]);
G4LogicalVolume*  Yoke_L_log = new G4LogicalVolume(Yoke_LR_box, Fe, "Yoke_L_log",0,0,0);
  Yoke_L_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //Yoke_L_log->SetUserLimits(new G4UserLimits(maxStep));
G4PVPlacement* Yoke_L_phys = new G4PVPlacement(rotForwardSp,
				  G4ThreeVector(pos_YOKE_L[XCOORD],pos_YOKE_L[YCOORD],pos_YOKE_L[ZCOORD]).rotateY(fSpectrometerAngle),
				  Yoke_L_log,
				  "Yoke_L_phys",
				  expHallLV,
				  false,
				  0);

G4LogicalVolume*  Yoke_R_log = new G4LogicalVolume(Yoke_LR_box, Fe, "Yoke_R_log",0,0,0);
  Yoke_R_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //Yoke_R_log->SetUserLimits(new G4UserLimits(maxStep));
G4PVPlacement* Yoke_R_phys = new G4PVPlacement(rotForwardSp,
				  G4ThreeVector(pos_YOKE_R[XCOORD],pos_YOKE_R[YCOORD],pos_YOKE_R[ZCOORD]).rotateY(fSpectrometerAngle),
				  Yoke_R_log,
				  "Yoke_R_phys",
				  expHallLV,
				  false,
				  0);


  //-------------------- Downstream End Guard
  G4Box* downGuard_UD_box = new G4Box("downGuard_UD_box",
				      size_DGUARD_UD[XCOORD],size_DGUARD_UD[YCOORD],size_DGUARD_UD[ZCOORD]);
G4LogicalVolume*  downGuard_U_log = new G4LogicalVolume(downGuard_UD_box, Fe, "downGuard_U_log",0,0,0);
  downGuard_U_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //downGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));
G4PVPlacement* downGuard_U_phys = new G4PVPlacement(rotForwardSp,
				       G4ThreeVector(pos_DGUARD_U[XCOORD],pos_DGUARD_U[YCOORD],pos_DGUARD_U[ZCOORD]).rotateY(fSpectrometerAngle),
				       downGuard_U_log,
				       "downGuard_U_phys",
				       expHallLV,
				       false,
				       0);

G4LogicalVolume*  downGuard_D_log = new G4LogicalVolume(downGuard_UD_box, Fe, "downGuard_D_log",0,0,0);
  downGuard_D_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //downGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
G4PVPlacement* downGuard_D_phys = new G4PVPlacement(rotForwardSp,
				       G4ThreeVector(pos_DGUARD_D[XCOORD],pos_DGUARD_D[YCOORD],pos_DGUARD_D[ZCOORD]).rotateY(fSpectrometerAngle),
				       downGuard_D_log,
				       "downGuard_D_phys",
				       expHallLV,
				       false,
				       0);

  G4Box* downGuard_LR_box = new G4Box("downGuard_LR_box",
				      size_DGUARD_LR[XCOORD],size_DGUARD_LR[YCOORD],size_DGUARD_LR[ZCOORD]);
G4LogicalVolume*  downGuard_L_log = new G4LogicalVolume(downGuard_LR_box, Fe, "downGuard_L_log",0,0,0);
  downGuard_L_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //downGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
G4PVPlacement* downGuard_L_phys = new G4PVPlacement(rotForwardSp,
				       G4ThreeVector(pos_DGUARD_L[XCOORD],pos_DGUARD_L[YCOORD],pos_DGUARD_L[ZCOORD]).rotateY(fSpectrometerAngle),
				       downGuard_L_log,
				       "downGuard_L_phys",
				       expHallLV,
				       false,
				       0);

G4LogicalVolume*  downGuard_R_log = new G4LogicalVolume(downGuard_LR_box, Fe, "downGuard_R_log",0,0,0);
  downGuard_R_log->SetVisAttributes(blue);
  maxStep=0.00001*mm;
  //downGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
G4PVPlacement*  downGuard_R_phys = new G4PVPlacement(rotForwardSp,
				       G4ThreeVector(pos_DGUARD_R[XCOORD],pos_DGUARD_R[YCOORD],pos_DGUARD_R[ZCOORD]).rotateY(fSpectrometerAngle),
				       downGuard_R_log,
				       "downGuard_R_phys",
				       expHallLV,
				       false,
				       0);



  ///////////////////end kurama magnet

  //--------------DC1
 G4Box* DC1_box = new G4Box("DC1_box",size_DC1[XCOORD]/2.,size_DC1[YCOORD]/2,size_DC1[ZCOORD]/2);
 G4LogicalVolume*  DC1_log = new G4LogicalVolume(DC1_box, Ar, "DC1_log",0,0,0);
 DC1_log->SetVisAttributes(green);
 // G4double pos_DC1[3];
 // pos_DC1[XCOORD] = par_cham->get_DCPlaneCenter(DC1X, XCOORD)*mm;
 // pos_DC1[YCOORD] = par_cham->get_DCPlaneCenter(DC1X, YCOORD)*mm;
 // pos_DC1[ZCOORD] = (par_cham->get_DCPlaneCenter(DC1X, ZCOORD)
 //		    +par_cham->get_DCPlaneCenter(DC1U, ZCOORD))*0.5*mm;
 pos_DC1[ZCOORD]=pos_DC1[ZCOORD];
 G4PVPlacement* DC1_phys = new G4PVPlacement(rotForwardSp,
					     G4ThreeVector(pos_DC1[XCOORD],pos_DC1[YCOORD],pos_DC1[ZCOORD]).rotateY(fSpectrometerAngle),
					     DC1_log,
					     "DC1_phys",
					     expHallLV,
					     false,
					     0);

  //---------DC1 Planes
 G4Box* DC1Plane_box = new G4Box("DC1Plane_box",
				 size_DC1Plane[XCOORD],size_DC1Plane[YCOORD],size_DC1Plane[ZCOORD]);
 // G4LogicalVolume*    DC1Plane_log[6];
 G4PVPlacement*   DC1Plane_phys[6];
 G4double pos_DC1Plane[3]={0.};
 for (int i=0; i<6; i++) {
   switch (i) {
   case 0:
     sprintf(name1, "DC1U_log");
     sprintf(name2, "DC1U_phys");
     pos_DC1Plane[XCOORD] = 0.0*mm;
     pos_DC1Plane[YCOORD] = 0.0*mm;
     pos_DC1Plane[ZCOORD] = -50.0*mm;
     break;
   case 1:
     sprintf(name1, "DC1Up_log");
     sprintf(name2, "DC1Up_phys");
     pos_DC1Plane[XCOORD] = 0.0*mm;
     pos_DC1Plane[YCOORD] = 0.0*mm;
     pos_DC1Plane[ZCOORD] = -30.*mm;
     break;
   case 2:
     sprintf(name1, "DC1X_log");
     sprintf(name2, "DC1X_phys");
     pos_DC1Plane[XCOORD] = 0.0*mm;
     pos_DC1Plane[YCOORD] = 0.0*mm;
     pos_DC1Plane[ZCOORD] = -10.0*mm;
     break;
   case 3:
     sprintf(name1, "DC1Xp_log");
     sprintf(name2, "DC1Xp_phys");
     pos_DC1Plane[XCOORD] = 0.0*mm;
     pos_DC1Plane[YCOORD] = 0.0*mm;
     pos_DC1Plane[ZCOORD] = 10.0*mm;
     break;
   case 4:
     sprintf(name1, "DC1V_log");
     sprintf(name2, "DC1V_phys");
     pos_DC1Plane[XCOORD] = 0.0*mm;
     pos_DC1Plane[YCOORD] = 0.0*mm;
     pos_DC1Plane[ZCOORD] = 30.0*mm;
     break;
   case 5:
     sprintf(name1, "DC1Vp_log");
     sprintf(name2, "DC1Vp_phys");
     pos_DC1Plane[XCOORD] = 0.0*mm;
     pos_DC1Plane[YCOORD] = 0.0*mm;
     pos_DC1Plane[ZCOORD] = 50.0*mm;
     break;
   }
   DC1Plane_log[i] = new G4LogicalVolume(DC1Plane_box, Ar, name1,0,0,0);
   DC1Plane_phys[i] = new G4PVPlacement(0,
					G4ThreeVector(pos_DC1Plane[XCOORD],pos_DC1Plane[YCOORD],pos_DC1Plane[ZCOORD]),
					DC1Plane_log[i],
					name2,
					DC1_log,
					false,
					101+i);
 }


  //////////////////////////////shhwang checking position. hole
  G4Box* fs_box = new G4Box("fs_box",800.* mm,450.*mm,0.0001*mm);
  fsLV[0]=
    new G4LogicalVolume(fs_box, Air, "fsLV0");
  fsPV[0] =
    new G4PVPlacement(rotForwardSp, G4ThreeVector(0,0.*mm,450.1*mm).rotateY(fSpectrometerAngle),
                      fsLV[0], "ScintPV32",expHallLV, true, 32);
  fsLV[0]-> SetVisAttributes(fsVisAtt);
  /////////////////////////////////

  //////////////////////////////shhwang checking position.  ch
  G4Box* fs_box1 = new G4Box("fs_box1",800.* mm,450.*mm,0.0001*mm);
  fsLV[1]=
    new G4LogicalVolume(fs_box1, Air, "fsLV1");
  fsPV[1] =
    new G4PVPlacement(rotForwardSp, G4ThreeVector(pos_CH[0],0.*mm,pos_CH[2]+0.1*mm).rotateY(fSpectrometerAngle),
                      fsLV[1], "ScintPV33",expHallLV, true, 33);
  fsLV[1]-> SetVisAttributes(fsVisAtt);
  /////////////////////////////////
  //////////////////////////////shhwang checking position. dc1
  G4Box* fs_box2 = new G4Box("fs_box2",800.* mm,450.*mm,0.0001*mm);
  fsLV[2]=
    new G4LogicalVolume(fs_box2, Air, "fsLV2");
  fsPV[2] =
    //    new G4PVPlacement(rotForwardSp, G4ThreeVector(pos_DC1[0],0.*mm,pos_DC1[2]+0.1+146./2.).rotateY(fSpectrometerAngle),

    new G4PVPlacement(rotForwardSp, G4ThreeVector(70.,0.*mm,pos_DC1[2]+0.1+146./2.).rotateY(fSpectrometerAngle),
                      fsLV[2], "ScintPV34",expHallLV, true, 34);
  fsLV[2]-> SetVisAttributes(fsVisAtt);
  /////////////////////////////////
  //////////////////////////////shhwang checking position. enterance of Kurama spectrometer
  G4Box* fs_box3 = new G4Box("fs_box3",800.* mm,450.*mm,0.0001*mm);
  fsLV[3]=
    new G4LogicalVolume(fs_box3, Air, "fsLV3");
  fsPV[3] =
    new G4PVPlacement(rotForwardSp, G4ThreeVector(pos_MFIELD[0],0.*mm,pos_MFIELD[2]-400.-0.1).rotateY(fSpectrometerAngle),
                      fsLV[3], "ScintPV35",expHallLV, true, 35);
  fsLV[3]-> SetVisAttributes(fsVisAtt);
  /////////////////////////////////
  //////////////////////////////shhwang checking position. DG
  G4Box* fs_box4 = new G4Box("fs_box4",1200.* mm,800.*mm,0.0001*mm);
  fsLV[4]=
    new G4LogicalVolume(fs_box4, Air, "fsLV4");
  fsPV[4] =
    new G4PVPlacement(rotForwardSp, G4ThreeVector(pos_MFIELD[0],0.*mm,pos_DGUARD_D[2]+50.1).rotateY(fSpectrometerAngle),
                      fsLV[4], "ScintPV36",expHallLV, true, 36);
  fsLV[4]-> SetVisAttributes(fsVisAtt);
  /////////////////////////////////
  //////////////////////////////shhwang checking position. DC2
  G4Box* fs_box5 = new G4Box("fs_box5",1500.* mm,800.*mm,0.0001*mm);
  fsLV[5]=
    new G4LogicalVolume(fs_box5, Air, "fsLV5");
  fsPV[5] =
    new G4PVPlacement(rotForwardSp, G4ThreeVector(pos_DC2[0],0.*mm,(pos_DC2[2]+50.1)*mm).rotateY(fSpectrometerAngle),
                      fsLV[5], "ScintPV37",expHallLV, true, 37);
  fsLV[5]-> SetVisAttributes(fsVisAtt);
  /////////////////////////////////
  //////////////////////////////shhwang checking position. DC3
  G4Box* fs_box6 = new G4Box("fs_box6",1500.* mm,800.*mm,0.0001*mm);
  fsLV[6]=
    new G4LogicalVolume(fs_box6, Air, "fsLV6");
  fsPV[6] =
    new G4PVPlacement(rotForwardSp, G4ThreeVector(pos_DC3[0],0.*mm,pos_DC3[2]+150./2.+0.1).rotateY(fSpectrometerAngle),
                      fsLV[6], "ScintPV38",expHallLV, true, 38);
  fsLV[6]-> SetVisAttributes(fsVisAtt);
  /////////////////////////////////
  //////////////////////////////shhwang checking position. ftof
  G4RotationMatrix* rot_fs = new G4RotationMatrix();
  rot_fs->rotateY(-rotangle_FTOF[1]*deg-fSpectrometerAngle);

  G4Box* fs_box7 = new G4Box("fs_box7",1500.* mm,1000.*mm,0.0001*mm);
  fsLV[7]=
    new G4LogicalVolume(fs_box7, Air, "fsLV7");
  fsPV[7] =
    new G4PVPlacement(rot_fs, G4ThreeVector(pos_FTOF[0],0.*mm,pos_FTOF[2]+35.0).rotateY(fSpectrometerAngle),
                      fsLV[7], "ScintPV39",expHallLV, true, 39);
  fsLV[7]-> SetVisAttributes(fsVisAtt);
  /////////////////////////////////





  //--------------CH
  G4Box* CH_box = new G4Box("CH_box",size_CH[XCOORD]/2,size_CH[YCOORD]/2,size_CH[ZCOORD]/2);
  G4double CH_Overlap = 1.*mm;
  //  G4LogicalVolume* CH_log[CHMAX];
  G4PVPlacement* CH_phys[CHMAX];
  for (int i=0; i<CHMAX; i++) {
    sprintf(name1, "CH%d_log", i);
    CH_log[i] = new G4LogicalVolume(CH_box, Scinti, name1,0,0,0);
    CH_log[i]->SetVisAttributes(aqua);
    //    par_cham->calpc((double *)pos_CH, CHx, (double)(i+1));
    pos_CH[XCOORD]=-CHMAX/2.*(size_CH[XCOORD]-CH_Overlap)+5.75*mm+(size_CH[XCOORD]-CH_Overlap)*i;
    sprintf(name1, "CH%d", i);
    if(i%2==0){
      CH_phys[i] = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_CH[XCOORD],pos_CH[YCOORD],pos_CH[ZCOORD]-1.*mm).rotateY(fSpectrometerAngle),
				     CH_log[i],
				     name1,
				     expHallLV,
				     false,
				     i);
    }else if(i%2==1){
      CH_phys[i] = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_CH[XCOORD],pos_CH[YCOORD],pos_CH[ZCOORD]+1.*mm).rotateY(fSpectrometerAngle),
				     CH_log[i],
				     name1,
				     expHallLV,
				     false,
				     i);
    }
  }



  //--------------DC2
  G4Box* DC2_box = new G4Box("DC2_box",size_DC2[XCOORD]/2,size_DC2[YCOORD]/2,size_DC2[ZCOORD]/2);
G4LogicalVolume*  DC2_log = new G4LogicalVolume(DC2_box, Ar, "DC2_log",0,0,0);
  DC2_log->SetVisAttributes(green);
  //  pos_DC2[XCOORD] = par_cham->get_DCPlaneCenter(DC2X, XCOORD)*mm;
  //  pos_DC2[YCOORD] = par_cham->get_DCPlaneCenter(DC2X, YCOORD)*mm;
  //  pos_DC2[ZCOORD] = (par_cham->get_DCPlaneCenter(DC2X, ZCOORD)
  //		     +par_cham->get_DCPlaneCenter(DC2Yp, ZCOORD))*0.5*mm;
 pos_DC2[ZCOORD]=pos_DC2[ZCOORD];
G4PVPlacement* DC2_phys = new G4PVPlacement(rotForwardSp,
			       G4ThreeVector(pos_DC2[XCOORD],pos_DC2[YCOORD],pos_DC2[ZCOORD]).rotateY(fSpectrometerAngle),
			       DC2_log,
			       "DC2_phys",
			       expHallLV,
			       false,
			       0);

  //---------DC2 Planes
  G4Box* DC2Plane_box = new G4Box("DC2Plane_box",
				  size_DC2Plane[XCOORD],size_DC2Plane[YCOORD],size_DC2Plane[ZCOORD]);

  //  G4LogicalVolume* DC2Plane_log[4];
  G4PVPlacement* DC2Plane_phys[4];
 G4double pos_DC2Plane[3]={0.};
  for (int i=0; i<4; i++) {
    switch (i) {
    case 0:
      sprintf(name1, "DC2X_log");
      sprintf(name2, "DC2X_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = -16.5*mm;
      break;
    case 1:
      sprintf(name1, "DC2Xp_log");
      sprintf(name2, "DC2Xp_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = -8.7*mm;
      break;
    case 2:
      sprintf(name1, "DC2Y_log");
      sprintf(name2, "DC2Y_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = 8.7*mm;
      break;
    case 3:
      sprintf(name1, "DC2Yp_log");
      sprintf(name2, "DC2Yp_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = 16.5*mm;
      break;
    }
    DC2Plane_log[i] = new G4LogicalVolume(DC2Plane_box, Ar, name1,0,0,0);
    DC2Plane_phys[i] = new G4PVPlacement(0,
					 G4ThreeVector(pos_DC2Plane[XCOORD],pos_DC2Plane[YCOORD],pos_DC2Plane[ZCOORD]),
					 DC2Plane_log[i],
					 name2,
					 DC2_log,
					 false,
					 121+i);
  }



  //--------------DC3
  G4Box* DC3_box = new G4Box("DC3_box",size_DC3[XCOORD]/2,size_DC3[YCOORD]/2,size_DC3[ZCOORD]/2);
  G4LogicalVolume*  DC3_log = new G4LogicalVolume(DC3_box, Ar, "DC3_log",0,0,0);
  DC3_log->SetVisAttributes(green);
  //  pos_DC3[XCOORD] = par_cham->get_DCPlaneCenter(DC3X, XCOORD)*mm;
  //  pos_DC3[YCOORD] = par_cham->get_DCPlaneCenter(DC3X, YCOORD)*mm;
  //  pos_DC3[ZCOORD] = (par_cham->get_DCPlaneCenter(DC3X, ZCOORD)
  //		     +par_cham->get_DCPlaneCenter(DC3Yp, ZCOORD))*0.5*mm;
 pos_DC3[ZCOORD]=pos_DC3[ZCOORD];
G4PVPlacement* DC3_phys = new G4PVPlacement(rotForwardSp,
			       G4ThreeVector(pos_DC3[XCOORD],pos_DC3[YCOORD],pos_DC3[ZCOORD]).rotateY(fSpectrometerAngle),
			       DC3_log,
			       "DC3_phys",
			       expHallLV,
			       false,
			       0);

  //---------DC3 Planes
  G4Box* DC3Plane_box = new G4Box("DC3Plane_box",
				  size_DC3Plane[XCOORD],size_DC3Plane[YCOORD],size_DC3Plane[ZCOORD]);
  //  G4LogicalVolume* DC3Plane_log[4];
  G4PVPlacement* DC3Plane_phys[4];
 G4double pos_DC3Plane[3]={0.};
  for (int i=0; i<4; i++) {
    switch (i) {
    case 0:
      sprintf(name1, "DC3X_log");
      sprintf(name2, "DC3X_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = -46.5*mm;
      break;
    case 1:
      sprintf(name1, "DC3Xp_log");
      sprintf(name2, "DC3Xp_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = -14.5*mm;
      break;
    case 2:
      sprintf(name1, "DC3Y_log");
      sprintf(name2, "DC3Y_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = 16.5*mm;
      break;
    case 3:
      sprintf(name1, "DC3Yp_log");
      sprintf(name2, "DC3Yp_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = 46.5*mm;
      break;
    }
    DC3Plane_log[i] = new G4LogicalVolume(DC3Plane_box, Ar, name1,0,0,0);
    DC3Plane_phys[i] = new G4PVPlacement(0,
					 G4ThreeVector(pos_DC3Plane[XCOORD],pos_DC3Plane[YCOORD],pos_DC3Plane[ZCOORD]),
					 DC3Plane_log[i],
					 name2,
					 DC3_log,
					 false,
					 131+i);
  }


  //--------------FTOF
  G4Box* FTOF_mo_box = new G4Box("FTOF_mo_box",size_FTOF[XCOORD]*24+50.*mm,size_FTOF[YCOORD]+50.*mm,size_FTOF[ZCOORD]+0.1*mm);
  G4LogicalVolume* FTOF_mo_log = new G4LogicalVolume(FTOF_mo_box, Air, "FTOF_mo_log",0,0,0);
 G4RotationMatrix* rot_FTOF = new G4RotationMatrix();
 rot_FTOF->rotateY(-rotangle_FTOF[1]*deg-fSpectrometerAngle);
 G4PVPlacement* FTOF_mo_phys = new G4PVPlacement(rot_FTOF,
						 G4ThreeVector(pos_FTOF[XCOORD],pos_FTOF[YCOORD],pos_FTOF[ZCOORD]).rotateY(fSpectrometerAngle),
					     FTOF_mo_log,
					     "FTOF_mo_phys",
					     expHallLV,
					     false,
					     0);
 //.rotateY(-fSpectrometerAngle*deg),
  G4VisAttributes* FTOF_mo_VisAtt= new G4VisAttributes(false, G4Colour(1.,0.,0.));
  FTOF_mo_log-> SetVisAttributes(FTOF_mo_VisAtt);

 ////end mother volume of FTOF
  G4Box* FTOF_box = new G4Box("FTOF_box",size_FTOF[XCOORD],size_FTOF[YCOORD],size_FTOF[ZCOORD]);
  G4PVPlacement* FTOF_phys[FTOFMAX];
  //  G4LogicalVolume* FTOF_log[FTOFMAX];

  //  G4RotationMatrix* rot_FTOF = new G4RotationMatrix();
  //  rot_FTOF->rotateY(-rotangle_FTOF[1]*deg-fSpectrometerAngle);
  G4double FTOF_Overlap=5.*mm;
  G4double pos_FTOF_bar[3];
  for (int i=0; i<FTOFMAX; i++) {
    sprintf(name1, "FTOF%d_log", i);
    FTOF_log[i] = new G4LogicalVolume(FTOF_box, Scinti, name1,0,0,0);
    FTOF_log[i]->SetVisAttributes(aqua);
    maxStep=0.1*mm;
    //    FTOF_log[i]->SetUserLimits(new G4UserLimits(maxStep));
    //    par_cham->calpc((double *)pos_FTOF, FTOFx, (double)(i+1));
    pos_FTOF_bar[XCOORD]=-FTOFMAX/2*(size_FTOF[XCOORD]*2-FTOF_Overlap)+(size_FTOF[XCOORD]*2-FTOF_Overlap)*i;
    //    pos_FTOF[ZCOORD]=pos_FTOF[ZCOORD]-sin((size_FTOF[XCOORD]*2.-FTOF_Overlap)*i);
    sprintf(name1, "FTOF%d", i);
    if(i%2==0){
      FTOF_phys[i] = new G4PVPlacement(0,
				       G4ThreeVector(pos_FTOF_bar[XCOORD],0.,size_FTOF[ZCOORD]),
				       FTOF_log[i],
				       name1,
				       FTOF_mo_log,
				       false,
				       i);
    }else if(i%2==1){
      FTOF_phys[i] = new G4PVPlacement(0,
				       G4ThreeVector(pos_FTOF_bar[XCOORD],0.,-size_FTOF[ZCOORD]),
				       FTOF_log[i],
				       name1,
				       FTOF_mo_log,
				       false,
				       i);
    }
  }

    }//with or w/o Kurama spectrometer
    else if(with_kurama==0){

  G4double size_CH[3];
  size_CH[XCOORD] = 11.5*mm;
  size_CH[YCOORD] = 400.0*mm;
  size_CH[ZCOORD] = 2.0*mm;

  G4double size_FTOF[3];
  size_FTOF[XCOORD] = 40.0*mm;
  size_FTOF[YCOORD] = 900.0*mm;
  size_FTOF[ZCOORD] = 15.0*mm;

  G4double size_DC1[3];
  //  size_DC1[XCOORD] = 289.0*2.*mm;
  //  size_DC1[YCOORD] = 215.0*2.*mm;
  //  size_DC1[ZCOORD] = 146.0*mm;
  ///Tanida DC1 is to small
  size_DC1[XCOORD] = 389.0*2.*mm;
  size_DC1[YCOORD] = 315.0*2.*mm;
  size_DC1[ZCOORD] = 146.0*mm;


  G4double size_DC1Plane[3];
  size_DC1Plane[XCOORD] = 6.0*127.0*0.5*mm;
  size_DC1Plane[YCOORD] = 6.0*97.0*0.5*mm;
  size_DC1Plane[ZCOORD] = 0.0001*mm;

  G4double size_DC2[3];
  size_DC2[XCOORD] = 1186.5*mm;//# of wires are 128
  size_DC2[YCOORD] = 1186.5*mm;//# of wires are 128
  //size_DC2[ZCOORD] = 45.0*mm;
  size_DC2[ZCOORD] = 100.0*mm;

  G4double size_DC2Plane[3];
  size_DC2Plane[XCOORD] = 9.0*128.0*0.5*mm;
  size_DC2Plane[YCOORD] = 9.0*128.0*0.5*mm;
  size_DC2Plane[ZCOORD] = 0.0001*mm;

  G4double size_DC3[3];
  size_DC3[XCOORD] = 1900.0*mm;
  size_DC3[YCOORD] = 1280.0*mm;
  size_DC3[ZCOORD] = 150.*mm;

  G4double size_DC3Plane[3];
  size_DC3Plane[XCOORD] = 20.0*96.0*0.5*mm;
  size_DC3Plane[YCOORD] = 20.0*64.0*0.5*mm;
  size_DC3Plane[ZCOORD] = 0.0001*mm;


  //--------------CH
  G4Box* CH_box = new G4Box("CH_box",size_CH[XCOORD]/2,size_CH[YCOORD]/2,size_CH[ZCOORD]/2);
  G4double CH_Overlap = 1.*mm;
  //  G4LogicalVolume* CH_log[CHMAX];
  G4PVPlacement* CH_phys[CHMAX];
  for (int i=0; i<CHMAX; i++) {
    sprintf(name1, "CH%d_log", i);
    CH_log[i] = new G4LogicalVolume(CH_box, Scinti, name1,0,0,0);
    CH_log[i]->SetVisAttributes(aqua);
    //    par_cham->calpc((double *)pos_CH, CHx, (double)(i+1));
    pos_CH[XCOORD]=-CHMAX/2.*(size_CH[XCOORD]-CH_Overlap)+5.75*mm+(size_CH[XCOORD]-CH_Overlap)*i;
    sprintf(name1, "CH%d", i);
    if(i%2==0){
      CH_phys[i] = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_CH[XCOORD],pos_CH[YCOORD],pos_CH[ZCOORD]-1.*mm).rotateY(fSpectrometerAngle),
				     CH_log[i],
				     name1,
				     expHallLV,
				     false,
				     i);
    }else if(i%2==1){
      CH_phys[i] = new G4PVPlacement(rotForwardSp,
				     G4ThreeVector(pos_CH[XCOORD],pos_CH[YCOORD],pos_CH[ZCOORD]+1.*mm).rotateY(fSpectrometerAngle),
				     CH_log[i],
				     name1,
				     expHallLV,
				     false,
				     i);
    }
  }



  //--------------DC2
  G4Box* DC2_box = new G4Box("DC2_box",size_DC2[XCOORD]/2,size_DC2[YCOORD]/2,size_DC2[ZCOORD]/2);
G4LogicalVolume*  DC2_log = new G4LogicalVolume(DC2_box, Ar, "DC2_log",0,0,0);
  DC2_log->SetVisAttributes(green);
  //  pos_DC2[XCOORD] = par_cham->get_DCPlaneCenter(DC2X, XCOORD)*mm;
  //  pos_DC2[YCOORD] = par_cham->get_DCPlaneCenter(DC2X, YCOORD)*mm;
  //  pos_DC2[ZCOORD] = (par_cham->get_DCPlaneCenter(DC2X, ZCOORD)
  //		     +par_cham->get_DCPlaneCenter(DC2Yp, ZCOORD))*0.5*mm;
 pos_DC2[ZCOORD]=pos_DC2[ZCOORD];
G4PVPlacement* DC2_phys = new G4PVPlacement(rotForwardSp,
			       G4ThreeVector(pos_DC2[XCOORD],pos_DC2[YCOORD],pos_DC2[ZCOORD]).rotateY(fSpectrometerAngle),
			       DC2_log,
			       "DC2_phys",
			       expHallLV,
			       false,
			       0);

  //---------DC2 Planes
  G4Box* DC2Plane_box = new G4Box("DC2Plane_box",
				  size_DC2Plane[XCOORD],size_DC2Plane[YCOORD],size_DC2Plane[ZCOORD]);

  //  G4LogicalVolume* DC2Plane_log[4];
  G4PVPlacement* DC2Plane_phys[4];
 G4double pos_DC2Plane[3]={0.};
  for (int i=0; i<4; i++) {
    switch (i) {
    case 0:
      sprintf(name1, "DC2X_log");
      sprintf(name2, "DC2X_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = -16.5*mm;
      break;
    case 1:
      sprintf(name1, "DC2Xp_log");
      sprintf(name2, "DC2Xp_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = -8.7*mm;
      break;
    case 2:
      sprintf(name1, "DC2Y_log");
      sprintf(name2, "DC2Y_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = 8.7*mm;
      break;
    case 3:
      sprintf(name1, "DC2Yp_log");
      sprintf(name2, "DC2Yp_phys");
      pos_DC2Plane[XCOORD] = 0.0*mm;
      pos_DC2Plane[YCOORD] = 0.0*mm;
      pos_DC2Plane[ZCOORD] = 16.5*mm;
      break;
    }
    DC2Plane_log[i] = new G4LogicalVolume(DC2Plane_box, Ar, name1,0,0,0);
    DC2Plane_phys[i] = new G4PVPlacement(0,
					 G4ThreeVector(pos_DC2Plane[XCOORD],pos_DC2Plane[YCOORD],pos_DC2Plane[ZCOORD]),
					 DC2Plane_log[i],
					 name2,
					 DC2_log,
					 false,
					 121+i);
  }



  //--------------DC3
  G4Box* DC3_box = new G4Box("DC3_box",size_DC3[XCOORD]/2,size_DC3[YCOORD]/2,size_DC3[ZCOORD]/2);
  G4LogicalVolume*  DC3_log = new G4LogicalVolume(DC3_box, Ar, "DC3_log",0,0,0);
  DC3_log->SetVisAttributes(green);
  //  pos_DC3[XCOORD] = par_cham->get_DCPlaneCenter(DC3X, XCOORD)*mm;
  //  pos_DC3[YCOORD] = par_cham->get_DCPlaneCenter(DC3X, YCOORD)*mm;
  //  pos_DC3[ZCOORD] = (par_cham->get_DCPlaneCenter(DC3X, ZCOORD)
  //		     +par_cham->get_DCPlaneCenter(DC3Yp, ZCOORD))*0.5*mm;
 pos_DC3[ZCOORD]=pos_DC3[ZCOORD];
G4PVPlacement* DC3_phys = new G4PVPlacement(rotForwardSp,
			       G4ThreeVector(pos_DC3[XCOORD],pos_DC3[YCOORD],pos_DC3[ZCOORD]).rotateY(fSpectrometerAngle),
			       DC3_log,
			       "DC3_phys",
			       expHallLV,
			       false,
			       0);

  //---------DC3 Planes
  G4Box* DC3Plane_box = new G4Box("DC3Plane_box",
				  size_DC3Plane[XCOORD],size_DC3Plane[YCOORD],size_DC3Plane[ZCOORD]);
  //  G4LogicalVolume* DC3Plane_log[4];
  G4PVPlacement* DC3Plane_phys[4];
 G4double pos_DC3Plane[3]={0.};
  for (int i=0; i<4; i++) {
    switch (i) {
    case 0:
      sprintf(name1, "DC3X_log");
      sprintf(name2, "DC3X_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = -46.5*mm;
      break;
    case 1:
      sprintf(name1, "DC3Xp_log");
      sprintf(name2, "DC3Xp_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = -14.5*mm;
      break;
    case 2:
      sprintf(name1, "DC3Y_log");
      sprintf(name2, "DC3Y_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = 16.5*mm;
      break;
    case 3:
      sprintf(name1, "DC3Yp_log");
      sprintf(name2, "DC3Yp_phys");
      pos_DC3Plane[XCOORD] = 0.0*mm;
      pos_DC3Plane[YCOORD] = 0.0*mm;
      pos_DC3Plane[ZCOORD] = 46.5*mm;
      break;
    }
    DC3Plane_log[i] = new G4LogicalVolume(DC3Plane_box, Ar, name1,0,0,0);
    DC3Plane_phys[i] = new G4PVPlacement(0,
					 G4ThreeVector(pos_DC3Plane[XCOORD],pos_DC3Plane[YCOORD],pos_DC3Plane[ZCOORD]),
					 DC3Plane_log[i],
					 name2,
					 DC3_log,
					 false,
					 131+i);
  }


  //--------------FTOF
  G4Box* FTOF_mo_box = new G4Box("FTOF_mo_box",size_FTOF[XCOORD]*24+50.*mm,size_FTOF[YCOORD]+50.*mm,size_FTOF[ZCOORD]+0.1*mm);
  G4LogicalVolume* FTOF_mo_log = new G4LogicalVolume(FTOF_mo_box, Air, "FTOF_mo_log",0,0,0);
 G4RotationMatrix* rot_FTOF = new G4RotationMatrix();
 rot_FTOF->rotateY(-rotangle_FTOF[1]*deg-fSpectrometerAngle);
 G4PVPlacement* FTOF_mo_phys = new G4PVPlacement(rot_FTOF,
						 G4ThreeVector(pos_FTOF[XCOORD],pos_FTOF[YCOORD],pos_FTOF[ZCOORD]).rotateY(fSpectrometerAngle),
					     FTOF_mo_log,
					     "FTOF_mo_phys",
					     expHallLV,
					     false,
					     0);
 //.rotateY(-fSpectrometerAngle*deg),
  G4VisAttributes* FTOF_mo_VisAtt= new G4VisAttributes(false, G4Colour(1.,0.,0.));
  FTOF_mo_log-> SetVisAttributes(FTOF_mo_VisAtt);

 ////end mother volume of FTOF
  G4Box* FTOF_box = new G4Box("FTOF_box",size_FTOF[XCOORD],size_FTOF[YCOORD],size_FTOF[ZCOORD]);
  G4PVPlacement* FTOF_phys[FTOFMAX];
  //  G4LogicalVolume* FTOF_log[FTOFMAX];

  //  G4RotationMatrix* rot_FTOF = new G4RotationMatrix();
  //  rot_FTOF->rotateY(-rotangle_FTOF[1]*deg-fSpectrometerAngle);
  G4double FTOF_Overlap=5.*mm;
  G4double pos_FTOF_bar[3];
  for (int i=0; i<FTOFMAX; i++) {
    sprintf(name1, "FTOF%d_log", i);
    FTOF_log[i] = new G4LogicalVolume(FTOF_box, Scinti, name1,0,0,0);
    FTOF_log[i]->SetVisAttributes(aqua);
    G4double maxStep=0.1*mm;
    //    FTOF_log[i]->SetUserLimits(new G4UserLimits(maxStep));
    //    par_cham->calpc((double *)pos_FTOF, FTOFx, (double)(i+1));
    pos_FTOF_bar[XCOORD]=-FTOFMAX/2*(size_FTOF[XCOORD]*2-FTOF_Overlap)+(size_FTOF[XCOORD]*2-FTOF_Overlap)*i;
    //    pos_FTOF[ZCOORD]=pos_FTOF[ZCOORD]-sin((size_FTOF[XCOORD]*2.-FTOF_Overlap)*i);
    sprintf(name1, "FTOF%d", i);
    if(i%2==0){
      FTOF_phys[i] = new G4PVPlacement(0,
				       G4ThreeVector(pos_FTOF_bar[XCOORD],0.,size_FTOF[ZCOORD]),
				       FTOF_log[i],
				       name1,
				       FTOF_mo_log,
				       false,
				       i);
    }else if(i%2==1){
      FTOF_phys[i] = new G4PVPlacement(0,
				       G4ThreeVector(pos_FTOF_bar[XCOORD],0.,-size_FTOF[ZCOORD]),
				       FTOF_log[i],
				       name1,
				       FTOF_mo_log,
				       false,
				       i);
    }
  }


    }

  }///E42 kurama setup


  TPCField* myfield = new TPCField("helmholtz_field.dat","KuramaMap80cm.dat");
  myfield->TPCField_Set();
 // myfield->TPCField_Set(on_off, on_off2, spectrometer_angle,k_gap,k_move_x, h_field, k_field);

  G4FieldManager* fieldMgr =
    G4TransportationManager::GetTransportationManager()->GetFieldManager();


  fieldMgr->SetDetectorField(myfield);
  fieldMgr->CreateChordFinder(myfield);


  //  G4cout << "test shhwang" << G4endl;
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // ==============================================================
  // define detector sensitivity
  // ==============================================================
  // sensitive Detectors

  G4SDManager* SDman= G4SDManager::GetSDMpointer();
  G4String SDname;
  //  G4cout << "test shhwang1" << G4endl;
  // ==============================================================
  // TPC Pad
  // ==============================================================

  TPCPadSD* padSD= new TPCPadSD(SDname="/TPC");
  padSD->TPCPadSD_Set();

  SDman-> AddNewDetector(padSD);
  //  for(G4int i = 0;i<NUM_PAD;i++){
  for(G4int i = 0;i<NUM_PAD;i++){
    padLV[i]->SetSensitiveDetector(padSD);
  }
  TPCLV->SetSensitiveDetector(padSD);
  //  TargetLV->SetSensitiveDetector(padSD);
  TPCTargetSD* tarSD= new TPCTargetSD(SDname="/TAR");
  SDman-> AddNewDetector(tarSD);
  TargetLV->SetSensitiveDetector(tarSD);

  //  G4cout << "test shhwang2" << G4endl;
  // ==============================================================
  // TPC scintillator (16 side + 2 downstream)
  // ==============================================================

  TPCScintSD* scintSD= new TPCScintSD(SDname="/SCINT");
  SDman-> AddNewDetector(scintSD);
  for(G4int i = 0; i<NPHI_SCINT; i++){
  //  for(G4int i = 0; i<1; i++){
    scintLV[i]->SetSensitiveDetector(scintSD);
  }

  if(with_kurama==1 && (experiment_num==42||experiment_num==27)){
  fsLV[0]->SetSensitiveDetector(scintSD);
  fsLV[1]->SetSensitiveDetector(scintSD);
  fsLV[2]->SetSensitiveDetector(scintSD);
  fsLV[3]->SetSensitiveDetector(scintSD);
  fsLV[4]->SetSensitiveDetector(scintSD);
  fsLV[5]->SetSensitiveDetector(scintSD);
  fsLV[6]->SetSensitiveDetector(scintSD);
  fsLV[7]->SetSensitiveDetector(scintSD);
  }

  //  fsLV[0]->SetSensitiveDetector(scintSD);
  //  fsLV[0]->SetSensitiveDetector(scintSD);
  //  fsLV[0]->SetSensitiveDetector(scintSD);
  //  fsLV[0]->SetSensitiveDetector(scintSD);
  //  fsLV[0]->SetSensitiveDetector(scintSD);

  if(ac_use==1){
    TPCACSD* acSD= new TPCACSD(SDname="/AC");
    SDman-> AddNewDetector(acSD);
    for(G4int i = 0; i<16; i++){
      pvcLV[i]->SetSensitiveDetector(acSD);
    }
  }


  if(n_bar_use==1){
    TPCNBARSD* nbarSD= new TPCNBARSD(SDname="/NBAR");
    SDman-> AddNewDetector(nbarSD);
    for(G4int i = 0; i<32; i++){
      nbarLV[i]->SetSensitiveDetector(nbarSD);
    }
  }


  if(experiment_num==42 || experiment_num==27){

  ////////dc
  TPCDCSD* dcSD= new TPCDCSD(SDname="/DC");
  SDman-> AddNewDetector(dcSD);
  if(with_kurama==1){
    for(G4int i = 0; i<6; i++){
      DC1Plane_log[i]->SetSensitiveDetector(dcSD);
    }
  }
  for(G4int i = 0; i<4; i++){
    DC2Plane_log[i]->SetSensitiveDetector(dcSD);
  }
  for(G4int i = 0; i<4; i++){
    DC3Plane_log[i]->SetSensitiveDetector(dcSD);
  }

  ////////////////FTOF
  TPCFTOFSD* ftofSD= new TPCFTOFSD(SDname="/FTOF");
  SDman-> AddNewDetector(ftofSD);
  for(G4int i = 0; i<24; i++){
    FTOF_log[i]->SetSensitiveDetector(ftofSD);
  }

  ////////////////CH
  TPCCHSD* chSD= new TPCCHSD(SDname="/CH");
  SDman-> AddNewDetector(chSD);
  for(G4int i = 0; i<CHMAX; i++){
    CH_log[i]->SetSensitiveDetector(chSD);
  }
  }
  //    G4cout << "test shhwang3" << G4endl;
  // ==============================================================
  // FDC
  // ==============================================================


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // ==============================================================
  // finally return the world volume
  // ==============================================================

  return expHall;

}

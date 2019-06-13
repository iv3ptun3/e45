#include "TPCField.hh"
#include "globals.hh"
#include <stdlib.h>
#include <fstream>
//#include "TPCParameters_getenv.hh"
const G4double Deg2Rad = acos(-1.)/180.;
const G4double Rad2Deg = 180./acos(-1.);

TPCField::TPCField(char* file1, char* file2)
{
  //, G4int on_off, G4int on_off2,G4double spec_ang, G4double k_gap, G4double k_move, G4double h_field, G4double k_field)
  //  G4cout<<TPCPameters_getenv::TPCParmeters_getenv::sh_test<<G4endl;
  //  G4String helmholtz_on_off=getenv("Helmholtz_fieldmap_on_off");
  //  G4int on_off=atoi( helmholtz_on_off.c_str()  );
  if(TPCField::env_helmholtz_on_off==1){
    std::ifstream ifs(file1);
    printf("Reading OPERA3D file by J.Y. Park \n");
    for(int ix= 0; ix<MAX_DIM_X_OPERA3D; ix++){
      for(int iy= 0; iy<MAX_DIM_Y_OPERA3D; iy++){
	for(int iz= 0; iz<MAX_DIM_Z_OPERA3D; iz++){
	  ifs >> xOPERA3D[ix] >> yOPERA3D[iy] >> zOPERA3D[iz] >>
	    bOPERA3D[0][ix][iy][iz] >> bOPERA3D[1][ix][iy][iz] >> bOPERA3D[2][ix][iy][iz];
	  bOPERA3D[0][ix][iy][iz] *=tesla;
	  bOPERA3D[1][ix][iy][iz] *=tesla;
	  bOPERA3D[2][ix][iy][iz] *=tesla;
	}
      }	
    }
     printf("Finish reading OPERA3D file\n");
  }
  //  G4cout<<TPCParameters_getenv::env_Generator<<G4endl; 
  //  G4String kurama_on_off=getenv("Kurama_fieldmap_on_off");
  //  G4int on_off2=atoi( kurama_on_off.c_str()  );
  if(TPCField::env_kurama_on_off==1 && env_with_kurama==1){
    std::ifstream ifs(file2);
    printf("Reading OPERA3D file by J.Y. Park \n");
    for(int ix= 0; ix<MAX_KURAMA_X_OPERA3D; ix++){
      for(int iy= 0; iy<MAX_KURAMA_Y_OPERA3D; iy++){
	for(int iz= 0; iz<MAX_KURAMA_Z_OPERA3D; iz++){
	  ifs >> xOPERA3D[ix] >> yOPERA3D[iy] >> zOPERA3D[iz] >>
	    bOPERA3D[0][ix][iy][iz] >> bOPERA3D[1][ix][iy][iz] >> bOPERA3D[2][ix][iy][iz];
	  bOPERA3D[0][ix][iy][iz] *=tesla;
	  bOPERA3D[1][ix][iy][iz] *=tesla;
	  bOPERA3D[2][ix][iy][iz] *=tesla;
	}
      }	
    }
     printf("Finish reading OPERA3D file\n");
  }

}  

//void TPCField::TPCField_Set(G4int on_off, G4int on_off2,G4double spec_ang, G4double k_gap, G4double k_move, G4double h_field, G4double k_field){
void TPCField::TPCField_Set(){
 G4String helmholtz_on_off=getenv("Helmholtz_fieldmap_on_off");
 G4int on_off=atoi( helmholtz_on_off.c_str()  );
 G4String kurama_on_off=getenv("Kurama_fieldmap_on_off");
 G4int on_off2=atoi( kurama_on_off.c_str()  );
 G4String spectrometer_ang = getenv("Spectrometer_ang");
 G4double spectrometer_angle=atof( spectrometer_ang.c_str()  )*deg;
 G4String env_kurama_gap=getenv("Kurama_gap");
 G4double k_gap=atof( env_kurama_gap.c_str())*mm;
 G4String kurama=getenv("Kurama_field");
 G4double k_field=atof( kurama.c_str()  );
 G4String helmholtz=getenv("Helmholtz_field");
 G4double h_field=atof( helmholtz.c_str()  );
  G4String env_kurama_pos_z = getenv("Kurama_pos_z");
  G4double kurama_pos_z = atof( env_kurama_pos_z.c_str())*mm;
  G4String kurama_move_x=getenv("Kurama_move_x");
  G4double k_move_x=atof( kurama_move_x.c_str());

  G4String with_kurama=getenv("With_KURAMA");
  G4int with_k=atoi( with_kurama.c_str());


  G4cout<<"-------------^____^--------------"<<G4endl;
  G4cout<<"Reading of B field parameters!!"<<G4endl;
  G4cout<<"---------------------------------"<<G4endl;
  env_with_kurama=with_k;
  env_k_gap=k_gap;
  env_k_move=k_move_x;
  env_helmholtz_on_off=on_off;
  env_kurama_on_off=on_off2;
  env_spec_angle=spectrometer_angle;
  env_h_field=h_field;
  env_k_field=k_field;
  env_k_pos_z=kurama_pos_z;

  G4cout<<"Kurama_gap:"<<env_k_gap<<G4endl;  
  G4cout<<"with KURAMA:"<<env_with_kurama<<G4endl;  
  G4cout<<"Kurama shift to X axis:"<<env_k_move<<G4endl;  

  G4cout<<"B field turn on/off of Helmholtz coil:"<<env_helmholtz_on_off<<G4endl;  
  G4cout<<"B field turn on/off of Kurama:"<<env_kurama_on_off<<G4endl;  
  G4cout<<"B field strength of Helmholtz coil:"<<env_h_field<<G4endl;  
  G4cout<<"B field strength of Kurama magent:"<<env_k_field<<G4endl;  
  G4cout<<"Angle of the K+ spectrometer:"<<env_spec_angle/deg<<G4endl;  
  G4cout<<"Position of the K+ spectrometer:"<<env_k_pos_z<<G4endl;  
  G4cout<<"-------------T____T-------------"<<G4endl;
  G4cout<<"!!!!!!!!!!!!!!!E42!!!!!!!!!!!!!!"<<G4endl;
}

void TPCField::GetFieldValue(const G4double Point[3], G4double Bfield[3]) const
{


  G4double zp, xp, yp;
  //  G4cout<<"test:"<<TPCField::env_k_gap<<G4endl;
  //  G4cout<<"test2:"<<TPCField::env_k_gap<<G4endl;
  //  bool flagInKurama = false;
  //  zp = cos(-TPCAngle*Deg2Rad)*Point[2]-sin(-TPCAngle*Deg2Rad)*Point[0];
  //  xp = sin(-TPCAngle*Deg2Rad)*Point[2]+cos(-TPCAngle*Deg2Rad)*Point[0];
  xp=Point[0];
  yp=Point[1];
  zp=Point[2];
  //  if ( sqrt(xp*xp+zp*zp) < 250.*mm){
  G4int on_off=TPCField::env_helmholtz_on_off;
  G4int on_off2=TPCField::env_kurama_on_off;
  G4double spectrometer_angle=TPCField::env_kurama_on_off;
  G4double k_move_x=TPCField::env_k_move;
  G4double k_gap=TPCField::env_k_gap;
  G4double h_field=TPCField::env_h_field;
  G4double k_field=TPCField::env_k_field;
  //  G4cout<<TPCField::env_Generator<<G4endl;
  //  G4double spectrometer_angle=0;
  //  G4double k_gap=400.;
  //  G4double k_move_x=200.;
  //  G4int on_off=0;
  //  G4cout<<"test:"<<on_off<<G4endl;
  G4double rot_xp=zp*sin(-spectrometer_angle)+xp*cos(-spectrometer_angle);
  G4double rot_zp=zp*cos(-spectrometer_angle)-xp*sin(-spectrometer_angle);

  if(on_off==0){
    //        G4String kurama_gap=getenv("Kurama_gap");
    //        G4double k_gap=atof( kurama_gap.c_str())*mm;
    //    
    //    G4String kurama_move_x=getenv("Kurama_move_x");
    //    G4double k_move_x=atof( kurama_move_x.c_str());

    if( zp > -310.*mm && zp < 310.*mm && xp > -310.*mm && xp < 310.*mm
	&& fabs(yp) < 300.*mm){
      //      G4String helmholtz=getenv("Helmholtz_field");
      //      G4double h_field=atof( helmholtz.c_str()  );
      Bfield[0]=0.*tesla;
      Bfield[1]=h_field*tesla;
      Bfield[2]=0.*tesla;
      
      // G4cout<<"Bfield(y) ="<<Bfield[1]<<", "<<h_field<<G4endl;      
      // G4cout<<"x, y, z ="<<xp*mm<<", "<<yp*mm<<", "<<zp*mm<<std::endl;


    }
    else if( env_with_kurama==1 && (rot_zp >= (env_k_pos_z-400.)*mm &&  rot_zp <= (env_k_pos_z+400.)*mm) &&
	     (rot_xp >= (k_move_x-500.)*mm && rot_xp <= (k_move_x+500.)*mm) && 
	     (yp < k_gap && yp > -k_gap)  ){
      //    G4cout<<k_gap<<G4endl;
      //            G4cout<<"------------------------------------"<<G4endl;
      //    G4cout<<"test0:"<<rot_xp<<":"<<rot_zp<<G4endl;
      //    G4cout<<"test1:"<<xp<<":"<<zp<<G4endl;
      //    G4String kurama=getenv("Kurama_field");
      //      G4double k_field=atof( kurama.c_str()  );
      Bfield[0]=0.*tesla;
      Bfield[1]=k_field*tesla;
      Bfield[2]=0.*tesla;

      //getchar();
      //      G4cout<<"test2:"<<yp<<":"<<k_field<<G4endl;      
    }
    else {
      Bfield[0]=0.*tesla;
      Bfield[1]=0.*tesla;
      Bfield[2]=0.*tesla;
    }

  }else if(on_off==1) { /// by using field map from OPERA-3D
    const G4double xmax = 250.0*mm;
    const G4double ymax = 250.0*mm;
    const G4double zmax = 250.0*mm;
    //    const G4double zmin = -300.0*mm;
    G4int iX,iY,iZ;
    
    //    Bfield[0] = 0.;
    //    Bfield[1] = 0.;
    //    Bfield[2] = 1.0*tesla; 
    //    return;
    
    iX = (int) ((Point[0]-xOPERA3D[0])/BFIELD_GRID_X);
    iY = (int) ((Point[1]-yOPERA3D[0])/BFIELD_GRID_Y);
    iZ = (int) ((Point[2]-zOPERA3D[0])/BFIELD_GRID_Z);
    G4cout<<"iX,iY,iZ"<<iX<<","<<iY<<","<<iZ<<G4endl;
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    if(fabs(Point[2])<zmax && fabs(Point[1])<ymax && fabs(Point[0])<xmax ) {
      G4double s,t,u;
      G4double sp, tp, up;
    /*
     *    bOPERA3D[0] : BX
     *    bOPERA3D[1] : BY
     *    bOPERA3D[2] : BZ
     *    bOut[0]   : Bx (at interior point)
     *    bOut[1]   : By (at interior point)
     *    bOut[2]   : Bz (at interior point)
     */
      s = (Point[0] - xOPERA3D[iX]) / BFIELD_GRID_X;
      t = (Point[1] - yOPERA3D[iY]) / BFIELD_GRID_Y;
      u = (Point[2] - zOPERA3D[iZ]) / BFIELD_GRID_Z;
      sp = 1.0-s; tp = 1.0-t;  up = 1.0-u; 
      
      for (int i=0; i<3; i++){
	G4double b[8]= {
	  bOPERA3D[i][iX][iY][iZ],       bOPERA3D[i][iX+1][iY][iZ],
	  bOPERA3D[i][iX+1][iY][iZ+1],   bOPERA3D[i][iX][iY][iZ+1],
	  bOPERA3D[i][iX][iY+1][iZ],     bOPERA3D[i][iX+1][iY+1][iZ],
	  bOPERA3D[i][iX+1][iY+1][iZ+1], bOPERA3D[i][iX][iY+1][iZ+1] };
	
	G4double bTmp  = up*b[0] + u*b[3]; /* up*(iX,  iY,  iZ) + u*(iX,  iY,  iZ+1) */
	G4double bTmp1 = up*b[1] + u*b[2]; /* up*(iX+1,iY,  iZ) + u*(iX+1,iY,  iZ+1) */
	G4double bTmp2 = up*b[4] + u*b[7]; /* up*(iX,  iY+1,iZ) + u*(iX,  iY+1,iZ+1) */
	G4double bTmp3 = up*b[5] + u*b[6]; /* up*(iX+1,iY+1,iZ) + u*(iX+1,iY+1,iZ+1) */
	
	Bfield[i] =  sp * (tp * bTmp + t * bTmp2) + s * (tp * bTmp1 + t * bTmp3) ;
	G4cout<<Bfield[i]<<G4endl;
      }
    } else {  //othere case
      Bfield[0] = 0.;
      Bfield[1] = 0.;
      Bfield[2] = 0.; 
    }
  }  
  return;
}


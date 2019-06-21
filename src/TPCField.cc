// -*- C++ -*-

#include "TPCField.hh"

#include <globals.hh>
#include <stdlib.h>
#include <fstream>

#include <CLHEP/Units/PhysicalConstants.h>

#include "ConfMan.hh"

namespace
{
  using CLHEP::deg;
  using CLHEP::mm;
  using CLHEP::tesla;
  const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
TPCField::TPCField( const G4String& file1, const G4String& file2 )
{
  if( gConf.Get<G4int>( "ShsFieldMap" ) == 1 ){
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
  if( gConf.Get<G4int>( "ConstructKurama" ) == 1 ){
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

//_____________________________________________________________________________
TPCField::~TPCField( void )
{
}

//_____________________________________________________________________________
void
TPCField::GetFieldValue( const G4double Point[3], G4double Bfield[3] ) const
{
  static const G4int shs_fieldmap = gConf.Get<G4int>( "ShsFieldMap" );
  static const G4int construct_kurama = gConf.Get<G4int>( "ConstructKurama" );
  static const G4double spectrometer_angle = gConf.Get<G4double>( "KuramaAngle" );
  static const G4double k_move_x = gConf.Get<G4double>( "KuramaShift" );
  static const G4double k_gap = 800*mm/2;
  static const G4double h_field = gConf.Get<G4double>( "ShsField" ) * tesla;
  static const G4double k_field = gConf.Get<G4double>( "KuramaField" ) * tesla;

  const G4double xp = Point[0];
  const G4double yp = Point[1];
  const G4double zp = Point[2];
  const G4double rot_xp = zp*sin(-spectrometer_angle)+xp*cos(-spectrometer_angle);
  const G4double rot_zp = zp*cos(-spectrometer_angle)-xp*sin(-spectrometer_angle);

  if( shs_fieldmap == 0 ){
    if( zp > -310.*mm && zp < 310.*mm && xp > -310.*mm && xp < 310.*mm
	&& fabs(yp) < 300.*mm){
      Bfield[0] = 0.*tesla;
      Bfield[1] = h_field;
      Bfield[2] = 0.*tesla;
    }
    else if( construct_kurama == 1 &&
	     (rot_zp >= -400.*mm && rot_zp <= 400.*mm) &&
	     (rot_xp >= (k_move_x-500.)*mm && rot_xp <= (k_move_x+500.)*mm) &&
	     (yp < k_gap && yp > -k_gap)  ){
      Bfield[0]=0.*tesla;
      Bfield[1]=k_field*tesla;
      Bfield[2]=0.*tesla;
    }
    else {
      Bfield[0]=0.*tesla;
      Bfield[1]=0.*tesla;
      Bfield[2]=0.*tesla;
    }

  } else { /// by using field map from OPERA-3D
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

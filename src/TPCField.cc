// -*- C++ -*-

#include "TPCField.hh"

#include <fstream>
#include <iomanip>

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4ThreeVector.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"

namespace
{
  using CLHEP::deg;
  using CLHEP::mm;
  using CLHEP::tesla;
  const auto& gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
}

//_____________________________________________________________________________
TPCField::TPCField( const G4String& file1, const G4String& file2 )
{
  G4cout << FUNC_NAME << G4endl;
  if( gConf.Get<G4int>( "ShsFieldMap" ) == 1 ){
    std::ifstream ifs( file1 );
    G4cout << "   Reading " << file1 << G4endl;
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
    G4cout << "   Finish reading OPERA3D file" << G4endl;
  }
  if( gConf.Get<G4int>( "KuramaFieldMap" ) == 1 ){
    std::ifstream ifs( file2 );
    G4cout << "   Reading " << file2 << G4endl;
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
    G4cout << "   Finish reading OPERA3D file" << G4endl;
  }
  G4cout << "   Initialized" << G4endl;
}

//_____________________________________________________________________________
TPCField::~TPCField( void )
{
}

//_____________________________________________________________________________
void
TPCField::GetFieldValue( const G4double Point[4], G4double* Bfield ) const
{
  static const G4int shs_fieldmap = gConf.Get<G4int>( "ShsFieldMap" );
  static const G4int construct_kurama = gConf.Get<G4int>( "ConstructKurama" );
  static const G4int kurama_fieldmap = gConf.Get<G4int>( "KuramaFieldMap" );
  static const G4double angle = gConf.Get<G4double>( "KuramaAngle" ) * deg;
  static const G4double h_field = gConf.Get<G4double>( "ShsField" ) * tesla;
  static const G4double k_field = gConf.Get<G4double>( "KuramaField" ) * tesla;
  static const auto kurama_pos = gGeom.GetGlobalPosition( "KURAMA" ) * mm;
  static const auto kurama_size = gSize.GetSize( "KuramaField" ) * mm;
  const G4double xp = Point[0];
  const G4double yp = Point[1];
  const G4double zp = Point[2];
  const G4double rot_xp = zp*std::sin(-angle) + xp*std::cos(-angle);
  const G4double rot_zp = zp*std::cos(-angle) - xp*std::sin(-angle);
  const G4ThreeVector pos( rot_xp, yp, rot_zp );

  if( shs_fieldmap == 0 ){
    if( zp > -310.*mm && zp < 310.*mm && xp > -310.*mm && xp < 310.*mm
	&& std::abs(yp) < 300.*mm){
      Bfield[0] = 0.*tesla;
      Bfield[1] = h_field;
      Bfield[2] = 0.*tesla;
    }
    else {
      Bfield[0] = 0.*tesla;
      Bfield[1] = 0.*tesla;
      Bfield[2] = 0.*tesla;
    }
  }
  if( construct_kurama == 1 && kurama_fieldmap == 0 &&
      std::abs( pos.x() - kurama_pos.x() ) <= kurama_size.x()/2 &&
      std::abs( pos.y() - kurama_pos.y() ) <= kurama_size.y()/2 &&
      std::abs( pos.z() - kurama_pos.z() ) <= kurama_size.z()/2 ){
    Bfield[0] = 0.*tesla;
    Bfield[1] = k_field;
    Bfield[2] = 0.*tesla;
  }
  // } else { /// by using field map from OPERA-3D
  //   const G4double xmax = 250.0*mm;
  //   const G4double ymax = 250.0*mm;
  //   const G4double zmax = 250.0*mm;
  //   //    const G4double zmin = -300.0*mm;
  //   G4int iX,iY,iZ;
  //   //    Bfield[0] = 0.;
  //   //    Bfield[1] = 0.;
  //   //    Bfield[2] = 1.0*tesla;
  //   //    return;
  //   iX = (int) ((Point[0]-xOPERA3D[0])/BFIELD_GRID_X);
  //   iY = (int) ((Point[1]-yOPERA3D[0])/BFIELD_GRID_Y);
  //   iZ = (int) ((Point[2]-zOPERA3D[0])/BFIELD_GRID_Z);
  //   G4cout<<"iX,iY,iZ"<<iX<<","<<iY<<","<<iZ<<G4endl;
  //   if( std::abs(Point[2])<zmax && std::abs(Point[1])<ymax &&
  // 	std::abs(Point[0])<xmax ){
  //     G4double s,t,u;
  //     G4double sp, tp, up;
  //     /*
  //      *    bOPERA3D[0] : BX
  //      *    bOPERA3D[1] : BY
  //      *    bOPERA3D[2] : BZ
  //      *    bOut[0]   : Bx (at interior point)
  //      *    bOut[1]   : By (at interior point)
  //      *    bOut[2]   : Bz (at interior point)
  //      */
  //     s = (Point[0] - xOPERA3D[iX]) / BFIELD_GRID_X;
  //     t = (Point[1] - yOPERA3D[iY]) / BFIELD_GRID_Y;
  //     u = (Point[2] - zOPERA3D[iZ]) / BFIELD_GRID_Z;
  //     sp = 1.0-s; tp = 1.0-t;  up = 1.0-u;
  //     for( G4int i=0; i<3; ++i ){
  // 	G4double b[8]= { bOPERA3D[i][iX][iY][iZ],
  // 			 bOPERA3D[i][iX+1][iY][iZ],
  // 			 bOPERA3D[i][iX+1][iY][iZ+1],
  // 			 bOPERA3D[i][iX][iY][iZ+1],
  // 			 bOPERA3D[i][iX][iY+1][iZ],
  // 			 bOPERA3D[i][iX+1][iY+1][iZ],
  // 			 bOPERA3D[i][iX+1][iY+1][iZ+1],
  // 			 bOPERA3D[i][iX][iY+1][iZ+1] };
  // 	G4double bTmp  = up*b[0] + u*b[3]; /* up*(iX,  iY,  iZ) + u*(iX,  iY,  iZ+1) */
  // 	G4double bTmp1 = up*b[1] + u*b[2]; /* up*(iX+1,iY,  iZ) + u*(iX+1,iY,  iZ+1) */
  // 	G4double bTmp2 = up*b[4] + u*b[7]; /* up*(iX,  iY+1,iZ) + u*(iX,  iY+1,iZ+1) */
  // 	G4double bTmp3 = up*b[5] + u*b[6]; /* up*(iX+1,iY+1,iZ) + u*(iX+1,iY+1,iZ+1) */
  // 	Bfield[i] = sp * (tp * bTmp + t * bTmp2) + s * (tp * bTmp1 + t * bTmp3) ;
  // 	G4cout<<Bfield[i]<<G4endl;
  //     }
  //   } else {  //othere case
  //     Bfield[0] = 0.;
  //     Bfield[1] = 0.;
  //     Bfield[2] = 0.;
  //   }
  // }
#if 0
  G4ThreeVector b( Bfield[0], Bfield[1], Bfield[2] );
  if( b.mag() > 0.01*tesla ){
    PrintHelper helper( 4, std::ios::fixed, G4cout );
    G4cout << FUNC_NAME << " X=( "
	   << std::setw(10) << Point[0] << " "
	   << std::setw(10) << Point[1] << " "
	   << std::setw(10) << Point[2] << " "
	   << std::setw(10) << Point[3] << " ), B=( "
	   << std::setw(10) << Bfield[0]/tesla << " "
	   << std::setw(10) << Bfield[1]/tesla << " "
	   << std::setw(10) << Bfield[2]/tesla << " )" << G4endl;
  }
#endif
  return;
}

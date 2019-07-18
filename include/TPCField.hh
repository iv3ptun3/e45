// -*- C++ -*-

#ifndef TPC_FIELD_HH
#define TPC_FIELD_HH

#include <globals.hh>
#include <G4MagneticField.hh>
#include <G4String.hh>

// const double solenoidOffset = 1617.0;

const G4int MAX_DIM_X_OPERA3D = 121;
const G4int MAX_DIM_Y_OPERA3D = 121;
const G4int MAX_DIM_Z_OPERA3D = 121;

const G4int MAX_KURAMA_X_OPERA3D = 121;
const G4int MAX_KURAMA_Y_OPERA3D = 121;
const G4int MAX_KURAMA_Z_OPERA3D = 121;


#define BFIELD_GRID_X 0.5
#define BFIELD_GRID_Y 0.5
#define BFIELD_GRID_Z 0.5

//_____________________________________________________________________________
class TPCField : public G4MagneticField
{
public:
  static G4String ClassName( void );
  TPCField( const G4String& fname, const G4String& fname2 );
  // TPCField( char* fname, char* fname2, G4int on_off, G4int on_off2,
  // 	    G4double spec_angle, G4double k_gap, G4double k_move,
  // 	    G4double h_field, G4double k_field );
  ~TPCField( void );

private:
  TPCField( void );
  TPCField( const TPCField& );
  TPCField& operator =( const TPCField& );

private:
  G4double xOPERA3D[MAX_DIM_X_OPERA3D];
  G4double yOPERA3D[MAX_DIM_Y_OPERA3D];
  G4double zOPERA3D[MAX_DIM_Z_OPERA3D];
  G4double bOPERA3D[3][MAX_DIM_X_OPERA3D][MAX_DIM_Y_OPERA3D][MAX_DIM_Z_OPERA3D];

  G4double xKuOPERA3D[MAX_KURAMA_X_OPERA3D];
  G4double yKuOPERA3D[MAX_KURAMA_Y_OPERA3D];
  G4double zKuOPERA3D[MAX_KURAMA_Z_OPERA3D];
  G4double bKuOPERA3D[3][MAX_KURAMA_X_OPERA3D][MAX_KURAMA_Y_OPERA3D][MAX_KURAMA_Z_OPERA3D];

  // G4double env_k_gap;
  // G4double env_helmholtz_on_off;
  // G4double env_kurama_on_off;
  // G4double env_k_move;
  // G4double env_with_kurama;
  // G4double env_spec_angle;
  // G4double env_h_field;
  // G4double env_k_field;
  // G4double env_k_pos_z;

public:
  virtual void GetFieldValue( const G4double Point[4], G4double* Bfield ) const;
};

//_____________________________________________________________________________
inline G4String
TPCField::ClassName( void )
{
  static G4String s_name("TPCField");
  return s_name;
}

#endif

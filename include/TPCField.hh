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
struct MagnetInfo
{
  enum EMagnetType { kDipole, kQuadrople };
  G4int    type;
  G4String name;
  G4double ra1;
};

//_____________________________________________________________________________
class TPCField : public G4MagneticField
{
public:
  static G4String ClassName( void );
  TPCField( void );
  ~TPCField( void );

private:
  TPCField( const TPCField& );
  TPCField& operator =( const TPCField& );

private:
  G4bool   m_k18_status;
  G4bool   m_kurama_status;
  G4bool   m_shs_status;
  G4String m_kurama_field_map;
  G4String m_shs_field_map;

  G4double xOPERA3D[MAX_DIM_X_OPERA3D];
  G4double yOPERA3D[MAX_DIM_Y_OPERA3D];
  G4double zOPERA3D[MAX_DIM_Z_OPERA3D];
  G4double bOPERA3D[3][MAX_DIM_X_OPERA3D][MAX_DIM_Y_OPERA3D][MAX_DIM_Z_OPERA3D];

  G4double xKuOPERA3D[MAX_KURAMA_X_OPERA3D];
  G4double yKuOPERA3D[MAX_KURAMA_Y_OPERA3D];
  G4double zKuOPERA3D[MAX_KURAMA_Z_OPERA3D];
  G4double bKuOPERA3D[3][MAX_KURAMA_X_OPERA3D][MAX_KURAMA_Y_OPERA3D][MAX_KURAMA_Z_OPERA3D];

public:
  virtual void GetFieldValue( const G4double Point[4], G4double* Bfield ) const;

public:
  G4bool GetStatusK18Field( void ) const { return m_k18_status; }
  G4bool GetStatusKuramaField( void ) const { return m_kurama_status; }
  G4bool GetStatusShsField( void ) const { return m_shs_status; }
  G4bool Initialize( void );
  void   SetKuramaFieldMap( G4String map ){ m_kurama_field_map = map; }
  void   SetShsFieldMap( G4String map ){ m_shs_field_map = map; }
  void   SetStatusK18Field( G4bool flag=true ){ m_k18_status = flag; }
  void   SetStatusKuramaField( G4bool flag=true ){ m_kurama_status = flag; }
  void   SetStatusShsField( G4bool flag=true ){ m_shs_status = flag; }

};

//_____________________________________________________________________________
inline G4String
TPCField::ClassName( void )
{
  static G4String s_name("TPCField");
  return s_name;
}

#endif

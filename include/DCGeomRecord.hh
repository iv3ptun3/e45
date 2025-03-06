// -*- C++ -*-

#ifndef DC_GEOM_RECORD_HH
#define DC_GEOM_RECORD_HH

#include <G4ThreeVector.hh>

#include <TString.h>

typedef G4ThreeVector ThreeVector;

//_____________________________________________________________________________
class DCGeomRecord
{ 
public:
  static TString ClassName( void );
  DCGeomRecord( Int_t id, const TString& name,
                Double_t x, Double_t y, Double_t z, Double_t ta,
                Double_t ra1, Double_t ra2, Double_t length, Double_t resol,
		Double_t w0, Double_t dd, Double_t ofs );
  DCGeomRecord( Int_t id, const TString& name,
		const ThreeVector pos, Double_t ta,
		Double_t ra1, Double_t ra2, Double_t length, Double_t resol,
		Double_t w0, Double_t dd, Double_t ofs );
  ~DCGeomRecord( void );

private:
  DCGeomRecord( const DCGeomRecord& );
  DCGeomRecord& operator =( const DCGeomRecord );

private:
  Int_t       m_id;
  TString     m_name;
  ThreeVector m_pos;
  Double_t    m_tilt_angle;
  Double_t    m_rot_angle1;
  Double_t    m_rot_angle2;
  Double_t    m_length;
  Double_t    m_resolution;
  Double_t    m_w0;
  Double_t    m_dd;
  Double_t    m_offset;

  Double_t m_dxds, m_dxdt, m_dxdu;
  Double_t m_dyds, m_dydt, m_dydu;
  Double_t m_dzds, m_dzdt, m_dzdu;

  Double_t m_dsdx, m_dsdy, m_dsdz;
  Double_t m_dtdx, m_dtdy, m_dtdz;
  Double_t m_dudx, m_dudy, m_dudz;

public:
  const ThreeVector& Position( void )       const { return m_pos; }
  ThreeVector        NormalVector( void )   const;
  ThreeVector        UnitVector( void )     const;
  Int_t              Id( void )             const { return m_id;         }
  TString            Name( void )           const { return m_name;       }
  const ThreeVector& Pos( void )            const { return m_pos;        }
  Double_t           TiltAngle( void )      const { return m_tilt_angle; }
  Double_t           RotationAngle1( void ) const { return m_rot_angle1; }
  Double_t           RotationAngle2( void ) const { return m_rot_angle2; }
  Double_t           Length( void )         const { return m_length;     }
  Double_t           Resolution( void )     const { return m_resolution; }
  void               SetResolution( Double_t res ) { m_resolution = res; }
  Double_t           WirePitch( void )      const { return m_dd;         }

  Double_t dsdx( void ) const { return m_dsdx; }
  Double_t dsdy( void ) const { return m_dsdy; }
  Double_t dsdz( void ) const { return m_dsdz; }
  Double_t dtdx( void ) const { return m_dtdx; }
  Double_t dtdy( void ) const { return m_dtdy; }
  Double_t dtdz( void ) const { return m_dtdz; }
  Double_t dudx( void ) const { return m_dudx; }
  Double_t dudy( void ) const { return m_dudy; }
  Double_t dudz( void ) const { return m_dudz; }

  Double_t dxds( void ) const { return m_dxds; }
  Double_t dxdt( void ) const { return m_dxdt; }
  Double_t dxdu( void ) const { return m_dxdu; }
  Double_t dyds( void ) const { return m_dyds; }
  Double_t dydt( void ) const { return m_dydt; }
  Double_t dydu( void ) const { return m_dydu; }
  Double_t dzds( void ) const { return m_dzds; }
  Double_t dzdt( void ) const { return m_dzdt; }
  Double_t dzdu( void ) const { return m_dzdu; }

  Double_t WirePos( Double_t wire )   const;
  Int_t    WireNumber( Double_t pos ) const;
  void     Print( void ) const;

private:
  void CalcVectors( void );
};

//_____________________________________________________________________________
inline TString
DCGeomRecord::ClassName( void )
{
  static const TString s_name("DCGeomRecord");
  return s_name;
}

//______________________________________________________________________________
struct DCGeomRecordComp
  : public std::__binary_function <DCGeomRecord*, DCGeomRecord*, Bool_t>
{
  Bool_t operator()( const DCGeomRecord* const p1,
		     const DCGeomRecord* const p2 ) const
  { return p1->Id() < p2->Id(); }
};

#endif

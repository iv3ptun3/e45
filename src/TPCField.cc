// -*- C++ -*-

#include "TPCField.hh"

#include <fstream>
#include <iomanip>

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4ThreeVector.hh>
#include <G4TwoVector.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FieldMap.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"

namespace
{
  using CLHEP::deg;
  using CLHEP::cm;
  using CLHEP::mm;
  using CLHEP::tesla;
  const auto& gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
}

//_____________________________________________________________________________
G4bool
MagnetInfo::CalcK18Field( const G4ThreeVector& point, G4double* bfield ) const
{
  PrintHelper helper( 4, std::ios::scientific, G4cout );
  G4ThreeVector dist = point - pos;
  switch( type ){
  case kDipole: {
    G4double phi = G4TwoVector( dist.x(), dist.z() ).phi()*math::Rad2Deg();
    if( std::abs( dist.mag() - rho ) < size.x() &&
	std::abs( dist.y() ) < size.y() &&
	( -180. < phi ) && ( phi < -180.+bend ) ){
      bfield[0] = 0.*tesla;
      bfield[1] = b0;
      bfield[2] = 0.*tesla;
      return true;
    } else {
      return false;
    }
  }
  case kQuadrupole: {
    dist.rotateY( -ra1 );
    if( dist.perp() < 2.*a0 &&
	std::abs( dist.z() ) < size.z() ){
      G4ThreeVector b( b0*dist.y()/a0, b0*dist.x()/a0, 0.*tesla );
#ifdef DEBUG
      G4cout << name << " " << dist << " " << b*(1/tesla) << G4endl;
#endif
      b.rotateY( ra1 );
      bfield[0] = b.x();
      bfield[1] = b.y();
      bfield[2] = b.z();
      return true;
    }
    return false;
  }
  default:
    G4Exception( FUNC_NAME, "", RunMustBeAborted, "" );
    return false;
  }
}

//_____________________________________________________________________________
TPCField::TPCField( void )
  : m_k18_status( false ),
    m_kurama_status( false ),
    m_shs_status( false ),
    m_kurama_field_map(),
    m_shs_field_map(),
    m_magnet_map()
{
}

//_____________________________________________________________________________
TPCField::~TPCField( void )
{
}

//_____________________________________________________________________________
void
TPCField::AddMagnetInfo( const MagnetInfo& mag )
{
  G4cout << "   Add magnet info : " << mag.name << G4endl;
  m_magnet_map[mag.name] = mag;
}

//_____________________________________________________________________________
G4bool
TPCField::Initialize( void )
{
  G4cout << FUNC_NAME << G4endl;

  if( m_shs_status &&
      gConf.Get<G4int>( "ShsFieldMap" ) == 1 ){
		double nmr =  gConf.Get<G4double>( "SHSFLDNMR" ) ;
		double scale_factor = gConf.Get<G4double>("HSFLDCALIB");
		nmr*=scale_factor;
    m_shs_field_map->SetValueCalc( gConf.Get<G4double>( "SHSFLDCALC" ) );
    m_shs_field_map->SetValueNMR( nmr);
    m_shs_field_map->Initialize();
  }
  if( m_kurama_status &&
      gConf.Get<G4int>( "KuramaFieldMap" ) == 1 ){
    m_kurama_field_map->SetValueCalc( gConf.Get<G4double>( "KURAMAFLDCALC" ) );
    m_kurama_field_map->SetValueNMR( gConf.Get<G4double>( "KURAMAFLDNMR" ) );
    m_kurama_field_map->Initialize();
  }

  return true;
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
  static const auto kurama_size = gSize.GetSize( "KuramaField" ) * 0.5 * mm;
  static const auto shs_pos = gGeom.GetGlobalPosition( "HypTPC" ) * mm;
  static const auto shs_size = gSize.GetSize( "ShsField" ) * 0.5 * mm;
  static const auto shs_field_ofset = gSize.GetSize( "ShsFieldOffset" ) * mm;
  static const auto shs_field_ofset_conf = gConf.Get<G4int>( "ShsFieldOffset" );
  const G4double xp = Point[0];
  const G4double yp = Point[1];
  const G4double zp = Point[2];
  const G4double rot_xp = zp*std::sin(-angle) + xp*std::cos(-angle);
  const G4double rot_zp = zp*std::cos(-angle) - xp*std::sin(-angle);
  const G4ThreeVector pos( rot_xp, yp, rot_zp );
  const G4ThreeVector kurama_coord = pos - kurama_pos;

  Bfield[0] = 0.*tesla;
  Bfield[1] = 0.*tesla;
  Bfield[2] = 0.*tesla;

  if( shs_fieldmap == 0 ){
    if( std::abs( pos.x() ) < shs_size.x() &&
	std::abs( pos.y() ) < shs_size.y() &&
	std::abs( pos.z() ) < shs_size.z() ){
      Bfield[0] += 0.*tesla;
      Bfield[1] += h_field;
      Bfield[2] += 0.*tesla;
    }
  } else {
    G4double shs_point[3] = { kurama_coord.x()/cm, kurama_coord.y()/cm, kurama_coord.z()/cm };
    
    if(shs_field_ofset_conf==1){
      for(int shsi=0;shsi<3;shsi++){
	shs_point[shsi]+=shs_field_ofset[shsi]/cm;
      }
    }
    if( m_shs_field_map->IsInsideField( shs_point ) ){
      m_shs_field_map->GetFieldValue( shs_point, Bfield );
    }
  }
  if( construct_kurama == 1 ){
    if( kurama_fieldmap == 0 ){
      if( std::abs( kurama_coord.x() ) < kurama_size.x() &&
	  std::abs( kurama_coord.y() ) < kurama_size.y() &&
	  std::abs( kurama_coord.z() ) < kurama_size.z() ){
	Bfield[0] += 0.*tesla;
	Bfield[1] += k_field;
	Bfield[2] += 0.*tesla;
      }
    } else {
      G4double kurama_point[3] =
	{ kurama_coord.x()/cm, kurama_coord.y()/cm, kurama_coord.z()/cm };
      //      if( m_kurama_field_map->IsInsideField( kurama_coord ) ){
      if( m_kurama_field_map->IsInsideField( kurama_point ) ){
	m_kurama_field_map->GetFieldValue( kurama_point, Bfield );
      }
    }
  }

  // K18Beamline
  if( m_k18_status ){
    for( auto& m : m_magnet_map ){
      if( m.second.CalcK18Field( G4ThreeVector( xp, yp, zp ), Bfield ) )
	return;
    }
  }

#if 0
  G4ThreeVector b( Bfield[0], Bfield[1], Bfield[2] );
  if( b.mag() > 0.001*tesla
      // || true
      ){
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

//_____________________________________________________________________________
void
TPCField::SetKuramaFieldMap( G4String map )
{
  if( m_kurama_field_map )
    delete m_kurama_field_map;
  m_kurama_field_map = new FieldMap( map );
}

//_____________________________________________________________________________
void
TPCField::SetShsFieldMap( G4String map )
{
	G4cout<<"Field "<<map<<G4endl;
  if( m_shs_field_map )
    delete m_shs_field_map;
  m_shs_field_map = new FieldMap( map );
}

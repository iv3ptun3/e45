// -*- C++ -*-

#include "ConfMan.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <libgen.h>
#include <sstream>
#include <vector>

#include <CLHEP/Units/SystemOfUnits.h>

#include "BeamMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "JamMan.hh"

//_____________________________________________________________________________
ConfMan::ConfMan( void )
  : m_conf_key("CONF"),
    m_conf_dir(),
    m_is_ready(false),
    m_file(),
    m_string(),
    m_double(),
    m_int(),
    m_bool()
{
}

//_____________________________________________________________________________
ConfMan::~ConfMan( void )
{
}

//_____________________________________________________________________________
G4bool
ConfMan::Initialize( void )
{
  if( m_is_ready ){
    G4cerr << FUNC_NAME
	   << " already initialied" << G4endl;
    return false;
  }

  std::ifstream ifs( m_file[m_conf_key] );
  if( !ifs.is_open() ){
    G4cerr << FUNC_NAME
	   << " cannot open file : " << m_file[m_conf_key] << G4endl;
    return false;
  }

  G4cout << FUNC_NAME << G4endl
	 << " open file : " << m_file[m_conf_key] << G4endl;

  m_conf_dir = ::dirname( const_cast<char*>( m_file[m_conf_key].data() ) );

  G4String line;
  while( ifs.good() && line.readLine( ifs ) ){
    if( line[0]=='#' ) continue;
    std::istringstream iss( line );
    G4String key, val;
    iss >> key >> val;
    if( key.isNull() || val.isNull() )
      continue;
    G4cout << " key = "   << std::setw(20) << std::left << key
	   << " value = " << std::setw(30) << std::left << val
	   << G4endl;
    m_file[key]   = FilePath(val);
    m_string[key] = val;
    m_double[key] = std::strtod( val, nullptr );
    m_int[key]    = std::strtol( val, nullptr, 10 );
    m_bool[key]   = static_cast<G4bool>( std::strtol( val, nullptr, 10 ) );
  }

  if ( !InitializeParameterFiles() || !InitializeHistograms() )
    return false;

  // if( gUser.IsReady() )
  //   gUser.Print();

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
ConfMan::Initialize( const G4String& file_name )
{
  m_file[m_conf_key] = file_name;
  return Initialize();
}

//_____________________________________________________________________________
G4bool
ConfMan::InitializeHistograms( void )
{
  return true;
}

//_____________________________________________________________________________
G4bool
ConfMan::InitializeParameterFiles( void )
{
  return ( InitializeParameter<DCGeomMan>("DCGEO") &&
	   InitializeParameter<BeamMan>("BEAM") &&
	   InitializeParameter<DetSizeMan>("DSIZE") &&
	   InitializeParameter<JamMan>("JAM") );
}

//_____________________________________________________________________________
// G4bool
// ConfMan::Finalize( void )
// {
//   return FinalizeProcess();
// }

//_____________________________________________________________________________
G4String
ConfMan::FilePath( const G4String& src ) const
{
  std::ifstream tmp( src );
  if ( tmp.good() )
    return src;
  else
    return m_conf_dir + "/" + src;
}

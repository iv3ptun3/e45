// -*- C++ -*-

#include "ConfMan.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

#include <TString.h>
#include <TSystem.h>

#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"

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
Bool_t
ConfMan::Initialize( void )
{
  if( m_is_ready ){
    std::cerr << "#W " << FUNC_NAME
	      << " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs( m_file[m_conf_key] );
  if( !ifs.is_open() ){
    std::cerr << "#E " << FUNC_NAME
	      << " cannot open file : " << m_file[m_conf_key] << std::endl;
    return false;
  }

  std::cout << "#D " << FUNC_NAME << std::endl
	      << " open file : " << m_file[m_conf_key] << std::endl;

  m_conf_dir = gSystem->DirName( m_file[m_conf_key] );

  TString line;
  while( ifs.good() && line.ReadLine( ifs ) ){
    if( line[0]=='#' ) continue;
    line.ReplaceAll(",",  ""); // remove ,
    line.ReplaceAll(":",  ""); // remove :
    line.ReplaceAll("\"", ""); // remove "

    std::istringstream iss( line.Data() );
    TString key, val;
    iss >> key >> val;
    if( key.IsNull() || val.IsNull() )
      continue;

    std::cout << " key = "   << std::setw(20) << std::left << key
	      << " value = " << std::setw(30) << std::left << val
	      << std::endl;

    m_file[key]   = FilePath(val);
    m_string[key] = val;
    m_double[key] = val.Atof();
    m_int[key]    = val.Atoi();
    m_bool[key]   = (Bool_t)val.Atoi();
  }

  if ( !InitializeParameterFiles() || !InitializeHistograms() )
    return false;

  // if( gUser.IsReady() )
  //   gUser.Print();

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::Initialize( const TString& file_name )
{
  m_file[m_conf_key] = file_name;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeHistograms( void )
{
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeParameterFiles( void )
{
  return ( InitializeParameter<DCGeomMan>("DCGEO") &&
	   InitializeParameter<DetSizeMan>("DSIZE") );
}

//_____________________________________________________________________________
// Bool_t
// ConfMan::Finalize( void )
// {
//   return FinalizeProcess();
// }

//_____________________________________________________________________________
TString
ConfMan::FilePath( const TString& src ) const
{
  std::ifstream tmp( src );
  if ( tmp.good() )
    return src;
  else
    return m_conf_dir + "/" + src;
}

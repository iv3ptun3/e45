// -*- C++ -*-

#include "DetSizeMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

#include <G4ThreeVector.hh>

#include "FuncName.hh"

//_____________________________________________________________________________
DetSizeMan::DetSizeMan( void )
  : m_is_ready(false),
    m_file_name(),
    m_param_map()
{
}

//_____________________________________________________________________________
DetSizeMan::~DetSizeMan( void )
{
}

//_____________________________________________________________________________
G4bool
DetSizeMan::Initialize( void )
{
  std::ifstream ifs( m_file_name );
  if( !ifs.is_open() ){
    std::cerr << "#E " << FUNC_NAME << " "
	      << "No such parameter file : " << m_file_name << std::endl;
    return false;
  }

  G4String line;
  while( ifs.good() && line.readLine(ifs) ){
    if( line[0]=='#' ) continue;
    std::istringstream input_line( line );

    G4String first_param;
    input_line >> first_param;

    G4String key = first_param;
    ParamArray param_array;
    G4double   param;
    while( input_line >> param ){
      param_array.push_back( param );
    }
    m_param_map[key] = param_array;
  }

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
DetSizeMan::Initialize( const G4String& filename )
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
G4double
DetSizeMan::Get( const G4String& key, G4int i ) const
{
  std::stringstream param;
  param << key << "(" << i << ")";

  PIterator itr = m_param_map.find(key);

  if( itr==m_param_map.end() ||
      i+1 > (G4int)itr->second.size() ){
    throw std::invalid_argument( std::string(FUNC_NAME+" No such key : "+key) );
  }

  return itr->second.at(i);
}

//_____________________________________________________________________________
G4ThreeVector
DetSizeMan::GetSize( const G4String& key ) const
{
  return G4ThreeVector( Get( key, G4ThreeVector::X ),
			Get( key, G4ThreeVector::Y ),
			Get( key, G4ThreeVector::Z ) );
}

//_____________________________________________________________________________
void
DetSizeMan::Print( void ) const
{
  std::cout << "#D " << FUNC_NAME << std::endl;

  const G4int w = 20;
  PIterator itr, end=m_param_map.end();
  for( itr=m_param_map.begin(); itr!=end; ++itr){
    std::cout << " key = " << std::setw(w) << std::left
		<< itr->first << itr->second.size() << " : ";
    for( G4int i=0, n=itr->second.size(); i<n; ++i ){
      std::cout << std::setw(5) << std::right
		  << itr->second.at(i) << " ";
    }
    std::cout << std::endl;
  }
}

// -*- C++ -*-

#include "DetSizeMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

#include "FuncName.hh"

//_____________________________________________________________________________
DetSizeMan::DetSizeMan( void )
  : m_is_ready(false),
    m_file_name()
{
}

//_____________________________________________________________________________
DetSizeMan::~DetSizeMan( void )
{
}

//_____________________________________________________________________________
Bool_t
DetSizeMan::Initialize( void )
{
  std::ifstream ifs( m_file_name );
  if( !ifs.is_open() ){
    std::cerr << "#E " << FUNC_NAME << " "
	      << "No such parameter file : " << m_file_name << std::endl;
    return false;
  }

  TString line;
  while( ifs.good() && line.ReadLine(ifs) ){
    if( line[0]=='#' ) continue;
    std::istringstream input_line( line.Data() );

    TString first_param;
    input_line >> first_param;

    TString key = first_param;
    ParamArray param_array;
    Double_t   param;
    while( input_line >> param ){
      param_array.push_back( param );
    }
    m_param_map[key] = param_array;
  }

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
DetSizeMan::Initialize( const TString& filename )
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
Int_t
DetSizeMan::GetSize( const TString& key ) const
{
  PIterator itr = m_param_map.find(key);
  if( itr==m_param_map.end() ){
    Print(m_file_name);
    std::cerr << "#E " << FUNC_NAME << " "
		<< "No such key : " << key << std::endl;
    return 0;
  }

  return itr->second.size();
}

//_____________________________________________________________________________
Double_t
DetSizeMan::Get( const TString& key, Int_t i ) const
{
  std::stringstream param;
  param << key << "(" << i << ")";

  PIterator itr = m_param_map.find(key);

  if( itr==m_param_map.end() ||
      i+1 > (Int_t)itr->second.size() ){
    throw std::invalid_argument( std::string(FUNC_NAME+" No such key : "+key) );
  }

  return itr->second.at(i);
}

//_____________________________________________________________________________
void
DetSizeMan::Print( Option_t* ) const
{
  std::cout << "#D " << FUNC_NAME << std::endl;

  const Int_t w = 20;
  PIterator itr, end=m_param_map.end();
  for( itr=m_param_map.begin(); itr!=end; ++itr){
    std::cout << " key = " << std::setw(w) << std::left
		<< itr->first << itr->second.size() << " : ";
    for( Int_t i=0, n=itr->second.size(); i<n; ++i ){
      std::cout << std::setw(5) << std::right
		  << itr->second.at(i) << " ";
    }
    std::cout << std::endl;
  }
}

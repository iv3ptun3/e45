// -*- C++ -*-

#include "JamMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <TFile.h>
#include <TTree.h>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"

#define SKIP_NP0 1

//_____________________________________________________________________________
void
JamInfo::Print( void ) const
{
  PrintHelper helper( 4, std::ios::fixed, G4cout );
  const G4int w = 8;
  G4cout << "   np=" << np << G4endl;
  for( G4int i=0; i<np; ++i ){
    G4cout << "   " << i << " "
	   << "pid=" << std::setw(w) << pid[i] << " "
	   << "px=" << std::setw(w) << px[i] << " "
	   << "py=" << std::setw(w) << py[i] << " "
	   << "pz=" << std::setw(w) << pz[i] << G4endl;
  }
}

//_____________________________________________________________________________
JamMan::JamMan( void )
  : m_is_ready( false ),
    m_file_name(),
    m_file(),
    m_param_array(),
    m_n_param()
{
}

//_____________________________________________________________________________
JamMan::~JamMan( void )
{
}

//_____________________________________________________________________________
G4bool
JamMan::Initialize( void )
{
  if( m_file_name.isNull() )
    return true;

  m_file = new TFile( m_file_name );
  TTree* tree = dynamic_cast<TTree*>( m_file->Get( "tree" ) );

  if( !m_file->IsOpen() || !tree )
    return false;

  m_param_array.clear();
  JamInfo jam;
  tree->SetBranchAddress( "np", &jam.np );
  tree->SetBranchAddress( "pid", jam.pid );
  tree->SetBranchAddress( "px", jam.px );
  tree->SetBranchAddress( "py", jam.py );
  tree->SetBranchAddress( "pz", jam.pz );

  for( Long64_t i=0, n=tree->GetEntries(); i<n; ++i ){
    tree->GetEntry( i );
#if SKIP_NP0
    if( jam.np == 0 )
      continue;
#endif
    m_param_array.push_back( jam );
  }

  m_file->Close();
  m_n_param = m_param_array.size();
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
JamMan::Initialize( const G4String& filename )
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
const JamInfo&
JamMan::Get( void ) const
{
  return m_param_array.at( G4RandFlat::shootInt( m_n_param ) );
}

//_____________________________________________________________________________
void
JamMan::Print( void ) const
{
  PrintHelper helper( 4, std::ios::fixed, G4cout );

  G4cout << FUNC_NAME << G4endl;
  for( const auto& j : m_param_array ){
    j.Print();
  }
  G4cout << "   nparam = " << m_param_array.size() << G4endl;
}

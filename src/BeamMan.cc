// -*- C++ -*-

#include "BeamMan.hh"

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

#include <ConfMan.hh>
#include "FuncName.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
void
BeamInfo::Print( void ) const
{
  PrintHelper helper( 4, std::ios::fixed, G4cout );
  const G4int w = 8;
  G4cout << "   "
	 << "x=" << std::setw(w) << x << " "
	 << "y=" << std::setw(w) << y << " "
	 << "u=" << std::setw(w) << u << " "
	 << "v=" << std::setw(w) << v << " "
	 << "p=(" << std::setw(w) << p.x() << ", "
	 << std::setw(w) << p.y() << ", "
	 << std::setw(w) << p.z() << ")" << G4endl;
}

//_____________________________________________________________________________
BeamMan::BeamMan( void )
  : m_is_ready(false),
    m_file_name(),
    m_file(),
    m_param_array(),
    m_n_param()
{
}

//_____________________________________________________________________________
BeamMan::~BeamMan( void )
{
}

//_____________________________________________________________________________
G4bool
BeamMan::Initialize( void )
{
  const G4double p0 = ConfMan::GetInstance().Get<G4double>( "BeamMom" );

  m_file = new TFile( m_file_name );
  TTree* tree = dynamic_cast<TTree*>( m_file->Get( "tree" ) );

  if( !m_file->IsOpen() || !tree )
    return false;

  BeamInfo beam;
  tree->SetBranchAddress( "x", &beam.x );
  tree->SetBranchAddress( "y", &beam.y );
  tree->SetBranchAddress( "u", &beam.u );
  tree->SetBranchAddress( "v", &beam.v );
  tree->SetBranchAddress( "p", &beam.dp );
  for( Long64_t i=0, n=tree->GetEntries(); i<n; ++i ){
    tree->GetEntry( i );
    G4double dxdz = std::tan( beam.u*CLHEP::mrad );
    G4double dydz = std::tan( beam.v*CLHEP::mrad );
    G4double pp = p0 * ( 1. + 0.01*beam.dp );
    G4double pz = pp / std::sqrt( dxdz*dxdz + dydz*dydz + 1. );
    beam.p.set( pz*dxdz, pz*dydz, pz );
    m_param_array.push_back( beam );
  }

  m_file->Close();
  m_n_param = m_param_array.size();
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
BeamMan::Initialize( const G4String& filename )
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
const BeamInfo&
BeamMan::Get( void ) const
{
  return m_param_array.at( G4RandFlat::shootInt( m_n_param ) );
}

//_____________________________________________________________________________
void
BeamMan::Print( void ) const
{
  PrintHelper helper( 4, std::ios::fixed, G4cout );
  const G4int w = 8;

  G4cout << FUNC_NAME << G4endl;
  for( const auto& b : m_param_array ){
    G4cout << "   "
	   << "x=" << std::setw(w) << b.x << " "
	   << "y=" << std::setw(w) << b.y << " "
	   << "u=" << std::setw(w) << b.u << " "
	   << "v=" << std::setw(w) << b.v << " "
	   << "p=(" << std::setw(w) << b.p.x() << ", "
	   << std::setw(w) << b.p.y() << ", "
	   << std::setw(w) << b.p.z() << ")" << G4endl;
  }
  G4cout << "   nparam = " << m_param_array.size() << G4endl;
}

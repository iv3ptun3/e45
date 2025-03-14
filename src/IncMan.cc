// -*- C++ -*-

#include "IncMan.hh"

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
IncInfo::Print( void ) const
{
  PrintHelper helper( 4, std::ios::fixed, G4cout );
  const G4int w = 8;
  G4cout << " * Reaction * " << G4endl;
  G4cout << "ich=" << std::setw(w) << ich << G4endl;
  G4cout << " * Beam * " << G4endl;
  G4cout << "px="  << std::setw(w) << bpx  << " "
         << "py="  << std::setw(w) << bpy  << " "
         << "pz="  << std::setw(w) << bpz  << G4endl;
  G4cout << " * Fs * " << G4endl;
  G4cout << "   np=" << np << G4endl;
  for( G4int i=0; i<np; ++i ){
    G4cout << "   - " << i << " "
	   << "pid=" << std::setw(w) << pid[i] << " "
	   << "px=" << std::setw(w)  << px[i] << " "
	   << "py=" << std::setw(w)  << py[i] << " "
	   << "pz=" << std::setw(w)  << pz[i] << G4endl;
  }
}

//_____________________________________________________________________________
IncMan::IncMan( void )
  : m_is_ready( false ),
    m_file_name(),
    m_file(),
    m_tree(),
    m_event( new IncInfo ),
    m_n_event()
{
}

//_____________________________________________________________________________
IncMan::~IncMan( void )
{
  if( m_file && m_file->IsOpen() )
    m_file->Close();
}

//_____________________________________________________________________________
G4bool
IncMan::Initialize( void )
{
  if( m_file_name.empty() )
    return true;

  m_file = new TFile( m_file_name );
  m_tree = dynamic_cast<TTree*>( m_file->Get( "tree" ) );

  if( !m_file->IsOpen() || !m_tree )
    return false;

  m_event = new IncInfo;

  m_file = new TFile( m_file_name );
  m_tree = dynamic_cast<TTree*>( m_file->Get( "tree" ) );
  // reaction
  m_tree->SetBranchAddress( "ich", &m_event->ich );
  // beam
  m_tree->SetBranchAddress( "bpx", &m_event->bpx );
  m_tree->SetBranchAddress( "bpy", &m_event->bpy );
  m_tree->SetBranchAddress( "bpz", &m_event->bpz );
  // beam
  m_tree->SetBranchAddress( "kppx", &m_event->kppx );
  m_tree->SetBranchAddress( "kppy", &m_event->kppy );
  m_tree->SetBranchAddress( "kppz", &m_event->kppz );
  // event
  m_tree->SetBranchAddress( "np", &m_event->np );
  m_tree->SetBranchAddress( "pid", m_event->pid );
  m_tree->SetBranchAddress( "px", m_event->px );
  m_tree->SetBranchAddress( "py", m_event->py );
  m_tree->SetBranchAddress( "pz", m_event->pz );

  m_tree->SetBranchStatus( "*", false );
  m_tree->SetBranchStatus( "ich", true );

  m_tree->SetBranchStatus( "bpx", true );
  m_tree->SetBranchStatus( "bpy", true );
  m_tree->SetBranchStatus( "bpz", true );
  m_tree->SetBranchStatus( "kppx", true );
  m_tree->SetBranchStatus( "kppy", true );
  m_tree->SetBranchStatus( "kppz", true );

  m_tree->SetBranchStatus( "np", true );
  m_tree->SetBranchStatus( "pid", true );
  m_tree->SetBranchStatus( "px", true );
  m_tree->SetBranchStatus( "py", true );
  m_tree->SetBranchStatus( "pz", true );

  m_n_event = m_tree->GetEntries();
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
G4bool
IncMan::Initialize( const G4String& filename )
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
IncInfo*
IncMan::Get( void ) const
{
  return Get( G4RandFlat::shootInt( m_n_event ) );
}

//_____________________________________________________________________________
IncInfo*
IncMan::Get( Int_t i ) const
{
  auto curr_file = gFile;
  m_file->cd();
  while(1){
    bool good = 1;
    if(i>=m_n_event)i=0;
    m_tree->GetEntry( i );
    for(int ip=0;ip<m_event->np;++ip){
      if(isnan(m_event->px[ip]))good = 0;
      if(isnan(m_event->py[ip]))good = 0;
      if(isnan(m_event->pz[ip]))good = 0;
    }
    if(good)break;
    G4cout<<"Warning! broken events! "<<i<<" -> "<<i+1<<G4endl;
    ++i;
  }
  if( curr_file )
    curr_file->cd();
  return m_event;
}

//_____________________________________________________________________________
void
IncMan::Print( void ) const
{
  PrintHelper helper( 4, std::ios::fixed, G4cout );

  G4cout << FUNC_NAME << G4endl;
  if( m_event )
    m_event->Print();
  G4cout << "   n_event = " << m_n_event << G4endl;
}

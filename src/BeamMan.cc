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

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
G4double
BeamInfo::GetX( G4double offset ) const
{
  return x + std::tan( -1.*u*CLHEP::mrad )*offset;
}

//_____________________________________________________________________________
G4double
BeamInfo::GetY( G4double offset ) const
{
  return y + std::tan( -1.*v*CLHEP::mrad )*offset;
}

//_____________________________________________________________________________
G4int 
BeamInfo::GetTrigPat(G4int flag) const
{
	return trigpat[flag];
}

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
  : m_is_ready( false ),
    m_file_name(),
    m_file(),
    m_param_array(),
    m_n_param(),
    m_is_vi( false ),
    m_primary_z( 0. )
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
  const auto& gConf = ConfMan::GetInstance();
  const auto& gGeom = DCGeomMan::GetInstance();
  const G4double p0 = gConf.Get<G4double>( "BeamMom" );

  if( m_file_name.isNull() )
    return true;
  
	m_param_array.clear();
  G4int generator = gConf.Get<G4int>( "Generator" );
	m_is_vi = ( gConf.Get<G4int>( "Generator" ) == 10 );
	if( abs(generator) == 135 or abs(generator) == 493 or abs(generator) ==938 ){
		m_is_k18 = 1;
	}
	if( abs(generator) == 100 ){
		m_is_kurama = 1;
	}

  m_primary_z = gGeom.GetLocalZ( "Vertex" );
  m_target_z = gGeom.GetLocalZ( "SHSTarget" );
//  double bh2_z = gGeom.GetLocalZ( "BH2" );
  if( !m_is_vi )
    m_primary_z -= 1318.9*CLHEP::mm; // from VO
	TTree* tree = nullptr;
	m_file = new TFile( m_file_name );
	if(m_is_k18){
		tree = dynamic_cast<TTree*>( m_file->Get( "k18track" ) );
	}
	else if (m_is_kurama){
		tree = dynamic_cast<TTree*>( m_file->Get( "kurama" ) );
	}
	else{
		tree = dynamic_cast<TTree*>( m_file->Get( "tree" ) );
	}


  if( !m_file->IsOpen() || !tree )
    return false;
	int ntBeam,evnum,runnum;
	double xout[5];
	double yout[5];
	double uout[5];
	double vout[5];
	double pBeam[5];
	double qBeam[5];
	double m2Beam[5];
	int trigpat[32];

  BeamInfo beam;
	if(m_is_k18){
		tree->SetBranchAddress( "ntK18",&ntBeam);
		tree->SetBranchAddress( "evnum",&evnum);
		tree->SetBranchAddress( "runnum",&runnum);
		tree->SetBranchAddress( "xout",xout);
		tree->SetBranchAddress( "yout",yout);
		tree->SetBranchAddress( "uout",uout);
		tree->SetBranchAddress( "vout",vout);
		tree->SetBranchAddress( "pHS",pBeam);
		tree->SetBranchAddress( "trigpat",trigpat);
	}
	else if (m_is_kurama){
		tree->SetBranchAddress( "ntKurama",&ntBeam);
		tree->SetBranchAddress( "evnum",&evnum);
		tree->SetBranchAddress( "runnum",&runnum);
		tree->SetBranchAddress( "xtgtKurama",xout);
		tree->SetBranchAddress( "ytgtKurama",yout);
		tree->SetBranchAddress( "utgtKurama",uout);
		tree->SetBranchAddress( "vtgtKurama",vout);
		tree->SetBranchAddress( "pKurama",pBeam);
		tree->SetBranchAddress( "qKurama",qBeam);
		tree->SetBranchAddress( "m2",m2Beam);
		tree->SetBranchAddress( "trigpat",trigpat);
	}
	else{
		tree->SetBranchAddress( "x", &beam.x );
		tree->SetBranchAddress( "y", &beam.y );
		tree->SetBranchAddress( "u", &beam.u );
		tree->SetBranchAddress( "v", &beam.v );
		tree->SetBranchAddress( "p", &beam.dp );
	}
	std::cout<<"BeamEvents = "<<tree->GetEntries()<<std::endl;
  for( Long64_t i=0, n=tree->GetEntries(); i<n; ++i ){
    tree->GetEntry( i );
		if(m_is_k18 or m_is_kurama){
			beam.x=0;
			beam.y=0;
			beam.z=0;
			beam.u=0;
			beam.v=0;
			beam.p.set(0,0,0);
			beam.evnum = evnum;
			beam.runnum = runnum;
			beam.ntBeam = ntBeam;
			for(int it=0;it<ntBeam;++it){
				beam.x = xout[0];
				beam.y = yout[0];
				beam.u = uout[0];// u,v definition = dxdz,dydz, not mrad. 
				beam.v = yout[0];
				double pz = pBeam[0] / sqrt(uout[0]*uout[0]+vout[0]*vout[0]+1);
				beam.p.set(pz * uout[0],pz* vout[0], pz);
				if(m_is_kurama){
					beam.z = m_target_z;
					beam.m2 = m2Beam[0];
					beam.q = qBeam[0];
				}
				else{
					beam.z = m_primary_z;//VO
				}
			}
			for(int itrg=0;itrg<32;++itrg){
				beam.trigpat[itrg] = trigpat[itrg];
			}
		}
		else{
			beam.x *= -1.*CLHEP::cm; // -cm -> mm
			beam.y *= -1.*CLHEP::cm; // -cm -> mm
			G4double dxdz = std::tan( -1.*beam.u*CLHEP::mrad ); // -mrad -> tan
			G4double dydz = std::tan( -1.*beam.v*CLHEP::mrad ); // -mrad -> tan
			G4double pp = p0 * ( 1. + beam.dp*CLHEP::perCent ); // dp/p[%] -> GeV/c
			G4double pz = pp / std::sqrt( dxdz*dxdz + dydz*dydz + 1. );
			beam.x += dxdz * m_primary_z;
			beam.y += dydz * m_primary_z;
			beam.z = m_primary_z;
			beam.p.set( pz*dxdz, pz*dydz, pz );
		}
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

const BeamInfo&
BeamMan::Get( G4int iev ) const
{
  return m_param_array.at( iev%m_n_param );
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

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
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
namespace{
	int ntBeam,evnum,runnum;
	double xout[5];
	double yout[5];
	double uout[5];
	double vout[5];
	double pBeam[5];
	double qBeam[5];
	double m2Beam[5];
	int trigpat[32];
	int trigflag[32];
	
	int ntK18,ntKurama;
  TTreeReaderValue<int>* ntTPCK18;
  TTreeReaderValue<int>* ntTPCKurama;

  TTreeReaderValue<vector<int>> *isgoodTPCK18 = nullptr;
  TTreeReaderValue<vector<int>> *isgoodTPCKurama = nullptr;

  TTreeReaderValue<vector<int>> *inside = nullptr;
  TTreeReaderValue<vector<double>> *pHS = nullptr;
  TTreeReaderValue<vector<double>>* xbTPC=nullptr;
  TTreeReaderValue<vector<double>>* ybTPC=nullptr;
  TTreeReaderValue<vector<double>>* ubTPC=nullptr;
  TTreeReaderValue<vector<double>>* vbTPC=nullptr;
  
  TTreeReaderValue<vector<double>>* xsTPC=nullptr;
  TTreeReaderValue<vector<double>>* ysTPC=nullptr;
  TTreeReaderValue<vector<double>>* usTPC=nullptr;
  TTreeReaderValue<vector<double>>* vsTPC=nullptr;

  TTreeReaderValue<vector<double>>* vtxTPC=nullptr;
  TTreeReaderValue<vector<double>>* vtyTPC=nullptr;
  TTreeReaderValue<vector<double>>* vtzTPC=nullptr;

  TTreeReaderValue<vector<double>>* xtgtHS=nullptr;
  TTreeReaderValue<vector<double>>* ytgtHS=nullptr;
  TTreeReaderValue<vector<double>>* ztgtHS=nullptr;
  TTreeReaderValue<vector<double>>* utgtHS=nullptr;
  TTreeReaderValue<vector<double>>* vtgtHS=nullptr;
  TTreeReaderValue<vector<double>>* pTPCKurama=nullptr;
  TTreeReaderValue<vector<double>>* qTPCKurama=nullptr;
  TTreeReaderValue<vector<double>>* m2TPCKurama=nullptr;
  TTreeReaderValue<vector<double>>* MissMass=nullptr;
  TTreeReaderValue<vector<double>>* utgtTPCKurama=nullptr;
  TTreeReaderValue<vector<double>>* vtgtTPCKurama=nullptr;
  TTreeReaderValue<vector<double>>* xtgtKurama=nullptr;
  TTreeReaderValue<vector<double>>* ytgtKurama=nullptr;
  
// For kurama data  
  TTreeReaderValue<vector<vector<double>>> *xvpHS = nullptr;
  TTreeReaderValue<vector<vector<double>>> *yvpHS = nullptr;
  TTreeReaderValue<vector<vector<double>>> *zvpHS = nullptr;
  TTreeReaderValue<vector<vector<double>>> *xvpKurama = nullptr;
  TTreeReaderValue<vector<vector<double>>> *yvpKurama = nullptr;
  TTreeReaderValue<vector<vector<double>>> *zvpKurama = nullptr;

  TTreeReaderValue<vector<double>> *xoutK18 = nullptr;
  TTreeReaderValue<vector<double>> *youtK18 = nullptr;
  TTreeReaderValue<vector<double>> *uoutK18 = nullptr;
  TTreeReaderValue<vector<double>> *voutK18 = nullptr;
  TTreeReaderValue<vector<double>> *p_3rd = nullptr;
  TTreeReaderValue<vector<vector<double>>> *layerK18 = nullptr;
  TTreeReaderValue<vector<vector<double>>> *wireK18 = nullptr;
  TTreeReaderValue<vector<vector<double>>> *localhitposK18 = nullptr;

  TTreeReaderValue<vector<double>> *xout_ = nullptr;
  TTreeReaderValue<vector<double>> *yout_ = nullptr;
  TTreeReaderValue<vector<double>> *zout_ = nullptr;
  TTreeReaderValue<vector<double>> *pxout = nullptr;
  TTreeReaderValue<vector<double>> *pyout = nullptr;
  TTreeReaderValue<vector<double>> *pzout = nullptr;
  TTreeReaderValue<vector<vector<double>>> *layer = nullptr;
  TTreeReaderValue<vector<vector<double>>> *wire = nullptr;
  TTreeReaderValue<vector<vector<double>>> *localhitpos = nullptr;
  TTreeReaderValue<vector<double>> *m2Kurama = nullptr;

  TTreeReaderValue<double> *KFpvalXi = nullptr;
  TTreeReaderValue<double> * KFXiPx = nullptr;
  TTreeReaderValue<double> * KFXiPy = nullptr;
  TTreeReaderValue<double> * KFXiPz = nullptr;

  TTreeReaderValue<bool> *XiFlight = nullptr;
  TTreeReaderValue<bool> *XiProd = nullptr;

  TTreeReaderValue<double>* xiprodvtx_x = nullptr;
  TTreeReaderValue<double>* xiprodvtx_y = nullptr;
  TTreeReaderValue<double>* xiprodvtx_z = nullptr;
  TTreeReaderValue<double>* xiprodmom_x = nullptr;
  TTreeReaderValue<double>* xiprodmom_y = nullptr;
  TTreeReaderValue<double>* xiprodmom_z = nullptr;


}
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
	m_mm_array.clear();
  G4int generator = gConf.Get<G4int>( "Generator" );
	m_is_vi = ( gConf.Get<G4int>( "Generator" ) == 10 );
	if( abs(generator) == 135 or abs(generator) == 493 or abs(generator) ==938 ){
		m_is_k18 = 1;
		G4cout<<"Generating K18 Beam"<<G4endl;
	}
	if( abs(generator) == 100 ){
		m_is_kurama = 1;
	}
	if(	generator == 181321 ){
		m_is_missmassXi = 1;
	}
	if( generator == 1001321){
		m_is_TPCXi = 1;
	}
  m_primary_z = gGeom.GetLocalZ( "Vertex" );
  m_target_z = gGeom.GetLocalZ( "SHSTarget" );
  if( !m_is_vi )
    m_primary_z -= 1318.9*CLHEP::mm; // from VO
	TTree* tree = nullptr;
	TTreeReader* reader = nullptr;
	m_file = new TFile( m_file_name );
	if(m_is_k18){
		tree = (TTree*) m_file->Get( "k18track" );
	}
	else if (m_is_kurama){
		tree = (TTree*) m_file->Get( "kurama" );
	}
	else if(m_is_missmassXi or m_is_TPCXi){
		tree = (TTree*)  m_file->Get( "tpc"  );
	}
	else{
		tree = (TTree*) m_file->Get("tree");
	}
	G4cout<<tree->GetEntries()<<G4endl;

  if( !m_file->IsOpen() || !tree )
    return false;
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
		tree->SetBranchAddress( "trigflag",trigflag);
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
	else if (m_is_missmassXi){
		reader = new TTreeReader("tpc",m_file);
//		tree->SetBranchAddress( "evnum",&evnum);
//		tree->SetBranchAddress( "runnum",&runnum);
		ntTPCK18 = new TTreeReaderValue<int>(*reader,"ntK18");		
		isgoodTPCK18 = new TTreeReaderValue<vector<int>>(*reader,"isgoodTPCK18");		
		inside = new TTreeReaderValue<vector<int>>(*reader,"inside"); 
		MissMass = new TTreeReaderValue<vector<double>>(*reader,"MissMassCorrDETPC"); 
		pHS = new TTreeReaderValue<vector<double>>(*reader,"pK18"); 
		xbTPC = new TTreeReaderValue<vector<double>>(*reader,"xbTPC"); 
		ybTPC = new TTreeReaderValue<vector<double>>(*reader,"ybTPC"); 
		ubTPC = new TTreeReaderValue<vector<double>>(*reader,"ubTPC"); 
		vbTPC = new TTreeReaderValue<vector<double>>(*reader,"vbTPC"); 
		vtxTPC = new TTreeReaderValue<vector<double>>(*reader,"vtxTPC");
		vtyTPC = new TTreeReaderValue<vector<double>>(*reader,"vtyTPC");
		vtzTPC = new TTreeReaderValue<vector<double>>(*reader,"vtzTPC");


		
//		tree->SetBranchAddress("ntKurama",&ntKurama);
		ntTPCKurama = new TTreeReaderValue<int>(*reader,"ntKurama");		
		isgoodTPCKurama = new TTreeReaderValue<vector<int>>(*reader,"isgoodTPCKurama");		
		pTPCKurama = new TTreeReaderValue<vector<double>>(*reader,"pCorrDETPC"); 
		qTPCKurama = new TTreeReaderValue<vector<double>>(*reader,"qTPCKurama"); 
		xsTPC = new TTreeReaderValue<vector<double>>(*reader,"xsTPC"); 
		ysTPC = new TTreeReaderValue<vector<double>>(*reader,"ysTPC"); 
		usTPC = new TTreeReaderValue<vector<double>>(*reader,"usTPC"); 
		vsTPC = new TTreeReaderValue<vector<double>>(*reader,"vsTPC"); 
		m2TPCKurama = new TTreeReaderValue<vector<double>>(*reader,"m2TPCKurama"); 
		
// For kurama data  
		p_3rd = new TTreeReaderValue<vector<double>>(*reader,"p_3rd");
		xoutK18 = new TTreeReaderValue<vector<double>>(*reader,"xoutK18");
		youtK18 = new TTreeReaderValue<vector<double>>(*reader,"youtK18");
		uoutK18 = new TTreeReaderValue<vector<double>>(*reader,"uoutK18");
		voutK18 = new TTreeReaderValue<vector<double>>(*reader,"voutK18");
		xvpHS = new TTreeReaderValue<vector<vector<double>>>(*reader,"xvpHS"); 
		yvpHS = new TTreeReaderValue<vector<vector<double>>>(*reader,"yvpHS"); 
		zvpHS = new TTreeReaderValue<vector<vector<double>>>(*reader,"zvpHS"); 
		xtgtHS = new TTreeReaderValue<vector<double>>(*reader,"xtgtHS");
		ytgtHS = new TTreeReaderValue<vector<double>>(*reader,"ytgtHS");
		ztgtHS = new TTreeReaderValue<vector<double>>(*reader,"ztgtHS");
		layerK18 = new TTreeReaderValue<vector<vector<double>>>(*reader,"layerK18");
		wireK18 = new TTreeReaderValue<vector<vector<double>>>(*reader,"wireK18");
		localhitposK18 = new TTreeReaderValue<vector<vector<double>>>(*reader,"localhitposK18");

		xvpKurama = new TTreeReaderValue<vector<vector<double>>>(*reader,"xvpKurama"); 
		yvpKurama = new TTreeReaderValue<vector<vector<double>>>(*reader,"yvpKurama"); 
		zvpKurama = new TTreeReaderValue<vector<vector<double>>>(*reader,"zvpKurama"); 
		xtgtKurama = new TTreeReaderValue<vector<double>>(*reader,"xtgtKurama");
		ytgtKurama = new TTreeReaderValue<vector<double>>(*reader,"ytgtKurama");
		xout_ = new TTreeReaderValue<vector<double>>(*reader,"xout"); 
		yout_ = new TTreeReaderValue<vector<double>>(*reader,"yout"); 
		zout_ = new TTreeReaderValue<vector<double>>(*reader,"zout"); 
		pxout = new TTreeReaderValue<vector<double>>(*reader,"pxout"); 
		pyout = new TTreeReaderValue<vector<double>>(*reader,"pyout"); 
		pzout = new TTreeReaderValue<vector<double>>(*reader,"pzout"); 
		layer = new TTreeReaderValue<vector<vector<double>>>(*reader,"layer"); 
		wire = new TTreeReaderValue<vector<vector<double>>>(*reader,"wire"); 
		localhitpos = new TTreeReaderValue<vector<vector<double>>>(*reader,"localhitpos"); 
		
	}
	else if(m_is_TPCXi){
		reader = new TTreeReader("tpc",m_file);

		ntTPCK18 = new TTreeReaderValue<int>(*reader,"ntK18");		
		pHS = new TTreeReaderValue<vector<double>>(*reader,"pK18");
		xbTPC = new TTreeReaderValue<vector<double>>(*reader,"xtgtK18");	
		ybTPC = new TTreeReaderValue<vector<double>>(*reader,"ytgtK18");
		ubTPC = new TTreeReaderValue<vector<double>>(*reader,"utgtK18");
		vbTPC = new TTreeReaderValue<vector<double>>(*reader,"vtgtK18");

		ntTPCKurama = new TTreeReaderValue<int>(*reader,"ntKurama");		
		isgoodTPCKurama = new TTreeReaderValue<vector<int>>(*reader,"isgoodTPCKurama");		
		pTPCKurama = new TTreeReaderValue<vector<double>>(*reader,"pCorrDETPC"); 
		qTPCKurama = new TTreeReaderValue<vector<double>>(*reader,"qTPCKurama"); 
		xsTPC = new TTreeReaderValue<vector<double>>(*reader,"xtgtTPCKurama"); 
		ysTPC = new TTreeReaderValue<vector<double>>(*reader,"ytgtTPCKurama");	
		usTPC = new TTreeReaderValue<vector<double>>(*reader,"utgtTPCKurama");	
		vsTPC = new TTreeReaderValue<vector<double>>(*reader,"vtgtTPCKurama");	
		m2TPCKurama = new TTreeReaderValue<vector<double>>(*reader,"m2TPCKurama"); 
		
		
		inside = new TTreeReaderValue<vector<int>>(*reader,"insideTPC"); 
		MissMass = new TTreeReaderValue<vector<double>>(*reader,"MissMassCorrDETPC"); 
		vtxTPC = new TTreeReaderValue<vector<double>>(*reader,"vtxTPC");
		vtyTPC = new TTreeReaderValue<vector<double>>(*reader,"vtyTPC");
		vtzTPC = new TTreeReaderValue<vector<double>>(*reader,"vtzTPC");

		KFpvalXi = new TTreeReaderValue<double>(*reader,"KFpvalXi");
		KFXiPx = new TTreeReaderValue<double>(*reader,"KFXimom_x");
		KFXiPy = new TTreeReaderValue<double>(*reader,"KFXimom_y");
		KFXiPz = new TTreeReaderValue<double>(*reader,"KFXimom_z");

		XiFlight = new TTreeReaderValue<bool>(*reader,"XiFlight");
		XiProd = new TTreeReaderValue<bool>(*reader,"XiProd");

		xiprodvtx_x = new TTreeReaderValue<double>(*reader,"xiprodvtx_x");
		xiprodvtx_y = new TTreeReaderValue<double>(*reader,"xiprodvtx_y");
		xiprodvtx_z = new TTreeReaderValue<double>(*reader,"xiprodvtx_z");
		xiprodmom_x = new TTreeReaderValue<double>(*reader,"xiprodmom_x");
		xiprodmom_y = new TTreeReaderValue<double>(*reader,"xiprodmom_y");
		xiprodmom_z = new TTreeReaderValue<double>(*reader,"xiprodmom_z");
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
		if(i%100000==0)G4cout<<Form("Event %lld/%lld",i,tree->GetEntries())<<G4endl;
    if(reader)reader->Next();
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
			if(m_is_k18 and trigflag[14]<0 and trigflag[23]<0) continue;
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
			m_nBeam++;
			m_param_array.push_back( beam );
		}
		else if(m_is_missmassXi){
			int in = (*inside)->at(0);
			double mm = (*MissMass)->at(0);
			int ntK18 = **ntTPCK18;
			int ntKurama = **ntTPCKurama;
			int goodTPCK18 = (*isgoodTPCK18)->at(0);
			int goodTPCKurama = (*isgoodTPCKurama)->at(0);
			for(int itk18=0;itk18<ntK18;++itk18){
				if(ntK18!= 1 or ntKurama !=1) continue;
				if(!goodTPCK18 or !goodTPCKurama) continue;
				double ub = (*ubTPC)->at(itk18);
				double vb = (*vbTPC)->at(itk18);
				double nb = hypot(hypot(1,ub),vb);
				double pb = (*pHS)->at(itk18);
				double pzb = pb/nb;
				G4ThreeVector TVKm(pzb*ub,pzb*vb,pzb);
			for(int itkurama=0;itkurama<ntKurama;++itkurama){
				double qKp = (*qTPCKurama)->at(itkurama);
				double m2Kp = (*m2TPCKurama)->at(itkurama);
				double pKp = (*pTPCKurama)->at(itkurama);
				
				double us = (*usTPC)->at(itkurama);
				double vs = (*vsTPC)->at(itkurama);
				double ns = hypot(hypot(1,us),vs);
				double pzs = pKp/ns;
				G4ThreeVector TVKp(pzs*us,pzs*vs,pzs);
				if((*vtzTPC)->at(0)==0 )continue;

				if(pKp < 1.4 and qKp > 0 and m2Kp > 0.12 and m2Kp < 0.4 and in and abs(mm-1.321)<0.1){
					MMVertex MMVert;
					MMVert.x = ((*vtxTPC)->at(0));
					MMVert.y = ((*vtyTPC)->at(0));
					MMVert.ntK18 = ntK18;
					MMVert.ntKurama = ntKurama;
					MMVert.Moms.push_back(TVKm);
					MMVert.Moms.push_back(TVKp);
					MMVert.Moms.push_back(TVKm-TVKp);

					//					vector<double>& xoutk18 = *(xoutK18->Get()); 
					MMVert.p_3rd = *(p_3rd->Get());
					MMVert.xoutK18 = *(xoutK18->Get());
					MMVert.youtK18 = *(youtK18->Get());
					MMVert.uoutK18 = *(uoutK18->Get());
					MMVert.voutK18 = *(voutK18->Get());	
					MMVert.xvpHS = *(xvpHS->Get());	
					MMVert.yvpHS = *(yvpHS->Get());	
					MMVert.zvpHS = *(zvpHS->Get());	
					MMVert.xtgtHS = *(xtgtHS->Get());	
					MMVert.ytgtHS = *(ytgtHS->Get());	
					MMVert.ztgtHS = *(ztgtHS->Get());	
					MMVert.layerK18 = *(layerK18->Get());	
					MMVert.wireK18 = *(wireK18->Get());	
					MMVert.localhitposK18 = *(localhitposK18->Get());	
					
					
					MMVert.xvpKurama = *(xvpKurama->Get());	
					MMVert.yvpKurama = *(yvpKurama->Get());	
					MMVert.zvpKurama = *(zvpKurama->Get());	
					MMVert.xtgtKurama = *(xtgtKurama->Get());	
					MMVert.ytgtKurama = *(ytgtKurama->Get());	
					
					MMVert.xoutK18 = *(xoutK18->Get());	
					MMVert.youtK18 = *(youtK18->Get());	
					MMVert.uoutK18 = *(uoutK18->Get());	
					MMVert.voutK18 = *(voutK18->Get());	
					MMVert.xout = *(xout_->Get());	
					MMVert.yout = *(yout_->Get());	
					MMVert.zout = *(zout_->Get());	
					MMVert.pxout = *(pxout->Get());	
					MMVert.pyout = *(pyout->Get());	
					MMVert.pzout = *(pzout->Get());	
					MMVert.layer = *(layer->Get());	
					MMVert.wire = *(wire->Get());	
					MMVert.localhitpos = *(localhitpos->Get());	
					
					m_mm_array.push_back(MMVert);		
				}
			}
			}
		}
		else if(m_is_TPCXi){
			int ntK18 = **ntTPCK18;
			int ntKurama = **ntTPCKurama;
			for(int itk18=0;itk18<ntK18;++itk18){
				if(ntK18!= 1 or ntKurama !=1) continue;
				int in = (*inside)->at(0);
				double mm = (*MissMass)->at(0);
				double ub = (*ubTPC)->at(itk18);
				double vb = (*vbTPC)->at(itk18);
				double nb = hypot(hypot(1,ub),vb);
				double pb = (*pHS)->at(itk18);
				double pzb = pb/nb;
				G4ThreeVector TVKm(pzb*ub,pzb*vb,pzb);
			for(int itkurama=0;itkurama<ntKurama;++itkurama){


				double qKp = (*qTPCKurama)->at(itkurama);
				double m2Kp = (*m2TPCKurama)->at(itkurama);
				double pKp = (*pTPCKurama)->at(itkurama);
				
				double us = (*usTPC)->at(itkurama);
				double vs = (*vsTPC)->at(itkurama);
				double ns = hypot(hypot(1,us),vs);
				double pzs = pKp/ns;
				G4ThreeVector TVKp(pzs*us,pzs*vs,pzs);
				if((*vtzTPC)->at(0)==0 )continue;
				double KFpval = **KFpvalXi;
				if(!**XiProd) continue;
				if(pKp < 1.4 and qKp > 0 and m2Kp > 0.12
				 and m2Kp < 0.4 and in
				 and abs(mm-1.321)<0.1 and KFpval>0.005){
					MMVertex MMVert;
					double PxXi = **xiprodmom_x;
					double PyXi = **xiprodmom_y;
					double PzXi = **xiprodmom_z;
					MMVert.x = **xiprodvtx_x;
					MMVert.y = **xiprodvtx_y;
					MMVert.z = **xiprodvtx_z+6;

					if(abs(MMVert.x)> 15 or abs(MMVert.y)>10 or abs(MMVert.z+143)>10) continue;	
					MMVert.xtgtHS = (*xbTPC->Get());
					MMVert.ytgtHS = (*ybTPC->Get());
					MMVert.xtgtKurama = (*xsTPC->Get());
					MMVert.ytgtKurama = (*ysTPC->Get());
					MMVert.ntK18 = ntK18;
					MMVert.ntKurama = ntKurama;
					MMVert.Moms.push_back(TVKm);
					MMVert.Moms.push_back(TVKp);
					G4ThreeVector TVXi(PxXi,PyXi,PzXi);	
					MMVert.Moms.push_back(TVXi);
					m_mm_array.push_back(MMVert);		
				}
			}
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
			m_param_array.push_back( beam );
		}
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
	auto b =	m_param_array.at(iev);
	return m_param_array.at( iev%m_param_array.size() );
}

const MMVertex&
BeamMan::GetVertex( void ) const
{
 	int nev = m_mm_array.size(); 
	return m_mm_array.at( G4RandFlat::shootInt( nev ) );
}

const MMVertex&
BeamMan::GetVertex( G4int iev ) const
{
 	int nev = m_mm_array.size(); 
	return m_mm_array.at( iev%nev );
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

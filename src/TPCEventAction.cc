// -*- C++ -*-

#include "TPCEventAction.hh"

#include <fstream>

#include <G4RunManager.hh>
#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4VVisManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4VTrajectory.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>

#include "ConfMan.hh"
#include "FuncName.hh"
#include "TPCAnaManager.hh"
#include "TPCScintSD.hh"
#include "TPCACSD.hh"
#include "TPCDCSD.hh"
#include "TPCFTOFSD.hh"
#include "TPCHTOFSD.hh"
#include "TPCNBARSD.hh"
#include "TPCPadSD.hh"
#include "TPCSCHSD.hh"
#include "TPCTargetSD.hh"

namespace
{
  using CLHEP::MeV;
  auto& gAnaMan = TPCAnaManager::GetInstance();
  const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
TPCEventAction::TPCEventAction( void )
  : G4UserEventAction()
{
}

//_____________________________________________________________________________
TPCEventAction::~TPCEventAction( void )
{
}

//_____________________________________________________________________________
void
TPCEventAction::BeginOfEventAction( const G4Event* )
{
  gAnaMan.BeginOfEventAction();
}

//_____________________________________________________________________________
void
TPCEventAction::EndOfEventAction( const G4Event* anEvent )
{
  G4int eventID = anEvent-> GetEventID();
  if( eventID % 100 == 0 ){
    G4cout << FUNC_NAME << G4endl
	   << "   Event number = " << eventID << G4endl;
  }
  G4SDManager* SDManager= G4SDManager::GetSDMpointer();

  // get "Hit Collection of This Event"
  G4HCofThisEvent* HCTE= anEvent-> GetHCofThisEvent();
  if(! HCTE) return;  // no hits in this events. nothing to do!

  static const G4int idcounter = SDManager->GetCollectionID("TPC/hit");
  if( idcounter > 0 ){
    auto padHC = (G4THitsCollection<TPCPadHit>*)HCTE->GetHC( idcounter );
    G4int nhits= padHC -> entries();
    G4int pidtr[MaxHitsTPC]={0};
    G4int ptidtpc[MaxHitsTPC]={0};
    G4int ptidtpc_pid[MaxHitsTPC]={0};
    G4double pmtpc[MaxHitsTPC]={0};
    G4int qqtpc[MaxHitsTPC]={0};
    G4double pxtpc[MaxHitsTPC]={0};
    G4double pytpc[MaxHitsTPC]={0};
    G4double pztpc[MaxHitsTPC]={0};
    G4double pptpc[MaxHitsTPC]={0};
    G4double vtxxtpc[MaxHitsTPC]={0};
    G4double vtxytpc[MaxHitsTPC]={0};
    G4double vtxztpc[MaxHitsTPC]={0};
    G4double vtxpxtpc[MaxHitsTPC]={0};
    G4double vtxpytpc[MaxHitsTPC]={0};
    G4double vtxpztpc[MaxHitsTPC]={0};
    // G4double vtxpptpc[MaxHitsTPC]={0};
    G4double vtxenetpc[MaxHitsTPC]={0};
    G4double detpc[MaxHitsTPC]={0};
    G4int laytpc[MaxHitsTPC]={0};
    G4double lentpc[MaxHitsTPC]={0};
    G4int nparticle=0;
    // G4cout << "TPC  " << nhits << G4endl;
    for( G4int i=0; i<nhits; ++i ){
      G4ThreeVector vtxpos = (*padHC)[i]-> GetVtxPosition();
      G4ThreeVector vtxmom = (*padHC)[i]-> GetVtxMomentum();
      G4double vtxene =(*padHC)[i]-> GetVtxEnergy();
      G4ThreeVector xyz = (*padHC)[i]-> GetPosition();
      G4ThreeVector mom = (*padHC)[i]-> GetMomentum();
      G4double tof= (*padHC)[i]-> GetTOF();
      G4int tid = (*padHC)[i]-> GetTrackID();
      G4int ptid = (*padHC)[i]-> GetParentID();
      G4int ptid_pid = (*padHC)[i]-> GetParentID_pid();
      G4int pid = (*padHC)[i]-> GetParticleID();
      G4double mass = (*padHC)[i]-> GetMass();
      G4int charge = (*padHC)[i]-> GetCharge();
      // std::cout<<"pid="<<pid<<", mass="<<mass<<", charge="<<charge<<std::endl;
      // getchar();
      G4int ilay = (*padHC)[i]-> GetPadLay();
      //      G4double mass = (*padHC)[i]-> GetPDGMass(); //mass(GeV)
      G4int parentid = (*padHC)[i]-> GetParentID();
      G4double tlength = (*padHC)[i]-> GettLength();
      G4int irow=0.;
      G4double beta = (*padHC)[i]-> GetBeta();
      G4double edep = (*padHC)[i]-> GetEdep();
      G4double slength = (*padHC)[i]-> GetsLength();
      //      G4VTrajectoryPoint *tp -> padHC->GetPoint(i);
      //    G4int nhits= padHC -> entries();
      if( nparticle==0 ){
	qqtpc[nparticle]=charge;
	pmtpc[nparticle]=mass;
	detpc[nparticle]=detpc[nparticle]+edep;
	// vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));
	if(ilay>-1){
	  laytpc[nparticle]=laytpc[nparticle]+1;
	}
	pidtr[nparticle]=pid;
	pxtpc[nparticle]=mom[0];
	pytpc[nparticle]=mom[1];
	pztpc[nparticle]=mom[2];

	//////////////////////vertex information /////////////////////////
	vtxpxtpc[nparticle]=vtxmom[0];
	vtxpytpc[nparticle]=vtxmom[1];
	vtxpztpc[nparticle]=vtxmom[2];
	// vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));

	vtxxtpc[nparticle]=vtxpos[0];
	vtxytpc[nparticle]=vtxpos[1];
	vtxztpc[nparticle]=vtxpos[2];

	vtxenetpc[nparticle]=vtxene;

	pptpc[nparticle]=sqrt(pow(mom[0],2.)+pow(mom[1],2.)+pow(mom[2],2.));
	lentpc[nparticle]=tlength;
	ptidtpc[nparticle]=ptid;
	ptidtpc_pid[nparticle]=ptid_pid;
	nparticle=nparticle+1;

      }else if( nparticle>0 ){
	//	G4cout<<nparticle<<G4endl;
      }
      if( (pidtr[nparticle-1] != pid) || (pidtr[nparticle-1] == pid && vtxpxtpc[nparticle-1] != vtxmom[0] && vtxpytpc[nparticle-1] != vtxmom[1] && vtxpztpc[nparticle-1] != vtxmom[2])){
	  qqtpc[nparticle]=charge;
	  pmtpc[nparticle]=mass;
	  detpc[nparticle]=detpc[nparticle]+edep;
	  if(ilay>-1){
	    laytpc[nparticle]=laytpc[nparticle]+1;
	  }
	  pidtr[nparticle]=pid;
	  pxtpc[nparticle]=mom[0];
	  pytpc[nparticle]=mom[1];
	  pztpc[nparticle]=mom[2];

	vtxpxtpc[nparticle]=vtxmom[0];
	vtxpytpc[nparticle]=vtxmom[1];
	vtxpztpc[nparticle]=vtxmom[2];
	// vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));
	vtxenetpc[nparticle]=vtxene;

	vtxxtpc[nparticle]=vtxpos[0];
	vtxytpc[nparticle]=vtxpos[1];
	vtxztpc[nparticle]=vtxpos[2];

	  pptpc[nparticle]=sqrt(pow(mom[0],2.)+pow(mom[1],2.)+pow(mom[2],2.));
	  lentpc[nparticle]=tlength;
	  ptidtpc[nparticle]=ptid;
	  ptidtpc_pid[nparticle]=ptid_pid;
	  nparticle=nparticle+1;
	}else if (pidtr[nparticle-1] == pid && vtxpxtpc[nparticle-1] == vtxmom[0] && vtxpytpc[nparticle-1] == vtxmom[1] && vtxpztpc[nparticle-1] == vtxmom[2]){
	  if( ptidtpc[nparticle-1] != ptid){
	    qqtpc[nparticle]=charge;
	    pmtpc[nparticle]=mass;
	    detpc[nparticle]=detpc[nparticle]+edep;
	    if(ilay>-1){
	      laytpc[nparticle]=laytpc[nparticle]+1;
	    }
	    pidtr[nparticle]=pid;
	    pxtpc[nparticle]=mom[0];
	    pytpc[nparticle]=mom[1];
	    pztpc[nparticle]=mom[2];

	vtxpxtpc[nparticle]=vtxmom[0];
	vtxpytpc[nparticle]=vtxmom[1];
	vtxpztpc[nparticle]=vtxmom[2];
	// vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));

	vtxenetpc[nparticle]=vtxene;

	vtxxtpc[nparticle]=vtxpos[0];
	vtxytpc[nparticle]=vtxpos[1];
	vtxztpc[nparticle]=vtxpos[2];

	    pptpc[nparticle]=sqrt(pow(mom[0],2.)+pow(mom[1],2.)+pow(mom[2],2.));
	    lentpc[nparticle]=tlength;
	    ptidtpc[nparticle]=ptid;
	    ptidtpc_pid[nparticle]=ptid_pid;
	    nparticle=nparticle+1;
	  }
	  else{
	    detpc[nparticle-1]=detpc[nparticle-1]+edep;
	    if(ilay>-1){
	      laytpc[nparticle-1]=laytpc[nparticle-1]+1;
	    }
	    lentpc[nparticle-1]=tlength;
	  }
	}
      //      if(ilay>-1){ //--> ilay 1 : target ilay 0 : TPC, layer is from 2 to 38.
      if(ilay>-1){ //-->  -1 : TPC, layer is from 0 to 38. 2012.10.30
	gAnaMan.SetCounterData(nparticle-1,tof, xyz, mom, tid, pid,ilay,irow,beta,edep/MeV,parentid,tlength,slength);
      }
    }
    for(G4int i=0;i<nparticle;i++){
      gAnaMan.SetTPCData(i, pidtr[i], ptidtpc[i], ptidtpc_pid[i], pxtpc[i],pytpc[i],pztpc[i],pptpc[i], qqtpc[i], pmtpc[i], detpc[i], lentpc[i],laytpc[i],
			     vtxpxtpc[i],vtxpytpc[i],vtxpztpc[i],
			     vtxxtpc[i],vtxytpc[i],vtxztpc[i], vtxenetpc[i]);
    }
  }

  static const G4int id_scint = SDManager-> GetCollectionID("HTOF/hit");
  if( id_scint > 0 ){
    const auto HC = (G4THitsCollection<TPCHTOFHit>*)HCTE->GetHC( id_scint );
    G4int nhits = HC->entries();
    for( G4int i=0; i<nhits; ++i ){
      auto hit = (*HC)[i];
      gAnaMan.SetHTOFData( hit );
    }
  }

  if( gConf.Get<G4int>( "UseAC" ) == 1. ){
    static const G4int id_ac = SDManager->GetCollectionID("AC/hit");
    if( id_ac > 0 ){
      auto acHC = (G4THitsCollection<TPCACHit>*)HCTE->GetHC( id_ac );
      G4int nhits= acHC -> entries();
      for(G4int i=0; i< nhits; i++) {
	G4ThreeVector vtxpos = (*acHC)[i]-> GetVtxPosition();
	G4ThreeVector vtxmom = (*acHC)[i]-> GetVtxMomentum();
	G4double vtxene =(*acHC)[i]-> GetVtxEnergy();
	G4ThreeVector xyz = (*acHC)[i]-> GetPosition();
	G4ThreeVector mom = (*acHC)[i]-> GetMomentum();
	G4int ptid = (*acHC)[i]-> GetParentID();
	G4double tof= (*acHC)[i]-> GetTOF();
	G4int tid = (*acHC)[i]-> GetTrackID();
	G4int pid = (*acHC)[i]-> GetParticleID();
	G4int did = (*acHC)[i]-> GetDetectorID();
	G4double mass = (*acHC)[i]-> GetParticleMassID();
	G4int qq = (*acHC)[i]-> GetParticleQqID();
	G4double tlength = (*acHC)[i]-> GetLength();
	gAnaMan.SetACData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
      }
    }
  }

  static const G4int id_nbar = SDManager-> GetCollectionID("NBAR/hit");
  if( id_nbar > 0 && gConf.Get<G4int>( "UseNBar" ) == 1. ){
    auto nbarHC = (G4THitsCollection<TPCNBARHit>*)HCTE->GetHC( id_nbar );
    G4int nhits = nbarHC->entries();
    for(G4int i=0; i< nhits; i++) {
      G4ThreeVector vtxpos = (*nbarHC)[i]-> GetVtxPosition();
      G4ThreeVector vtxmom = (*nbarHC)[i]-> GetVtxMomentum();
      G4double vtxene =(*nbarHC)[i]-> GetVtxEnergy();
      G4ThreeVector xyz = (*nbarHC)[i]-> GetPosition();
      G4ThreeVector mom = (*nbarHC)[i]-> GetMomentum();
      G4int ptid = (*nbarHC)[i]-> GetParentID();
      G4double tof= (*nbarHC)[i]-> GetTOF();
      G4int tid = (*nbarHC)[i]-> GetTrackID();
      G4int pid = (*nbarHC)[i]-> GetParticleID();
      G4int did = (*nbarHC)[i]-> GetDetectorID();
      G4double mass = (*nbarHC)[i]-> GetParticleMassID();
      G4int qq = (*nbarHC)[i]-> GetParticleQqID();
      G4double tlength = (*nbarHC)[i]-> GetLength();
      gAnaMan.SetNBARData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
    }
  }

  static const G4int id_target = SDManager->GetCollectionID("TAR/hit");
  if( id_target > 0 ){
    auto targetHC = (G4THitsCollection<TPCTargetHit>*)HCTE->GetHC( id_target );
    G4int nhits = targetHC->entries();
    for(G4int i=0; i< nhits; i++) {
      G4double kinene =(*targetHC)[i]-> GetKinEnergy();;
      G4ThreeVector mom = (*targetHC)[i]-> GetMomentum();
      G4double test= sqrt(pow(mom.getX(),2)+pow(mom.getY(),2)+pow(mom.getZ(),2));
      if(kinene==0.00000000000000000000 && test != 0.  ){
	G4ThreeVector vtxpos = (*targetHC)[i]-> GetVtxPosition();
	G4ThreeVector vtxmom = (*targetHC)[i]-> GetVtxMomentum();
	G4double vtxene =(*targetHC)[i]-> GetVtxEnergy();
	G4ThreeVector xyz = (*targetHC)[i]-> GetPosition();
	// G4double tof= (*targetHC)[i]-> GetTOF();
	G4int tid = (*targetHC)[i]-> GetTrackID();
	G4int ptid = (*targetHC)[i]-> GetParentID();
	G4int pid = (*targetHC)[i]-> GetParticleID();
	// G4double mass = (*targetHC)[i]-> GetMass();
	// G4int charge = (*targetHC)[i]-> GetCharge();
	// // G4double mass = (*targetHC)[i]-> GetPDGMass(); //mass(GeV)
	// G4int parentid = (*targetHC)[i]-> GetParentID();
	// G4double tlength = (*targetHC)[i]-> GettLength();
	// G4int irow=0.;
	// G4double beta = (*targetHC)[i]-> GetBeta();
	// G4double edep = (*targetHC)[i]-> GetEdep();
	gAnaMan.SetTargetData(i,xyz, mom, tid, pid,ptid,vtxpos,vtxmom,vtxene);
      }
    }
  }

  static const G4int id_dc = SDManager->GetCollectionID("DC/hit");
  if( id_dc > 0 ){
    auto dcHC = (G4THitsCollection<TPCDCHit>*)HCTE->GetHC( id_dc );
    G4int nhits= dcHC -> entries();
    for(G4int i=0; i< nhits; i++) {
      G4ThreeVector vtxpos = (*dcHC)[i]-> GetVtxPosition();
      G4ThreeVector vtxmom = (*dcHC)[i]-> GetVtxMomentum();
      G4double vtxene =(*dcHC)[i]-> GetVtxEnergy();
      G4ThreeVector xyz = (*dcHC)[i]-> GetPosition();
      G4ThreeVector mom = (*dcHC)[i]-> GetMomentum();
      G4int ptid = (*dcHC)[i]-> GetParentID();
      G4double tof= (*dcHC)[i]-> GetTOF();
      G4int tid = (*dcHC)[i]-> GetTrackID();
      G4int pid = (*dcHC)[i]-> GetParticleID();
      G4int did = (*dcHC)[i]-> GetDetectorID();
      G4double mass = (*dcHC)[i]-> GetParticleMassID();
      G4int qq = (*dcHC)[i]-> GetParticleQqID();
      G4double tlength = (*dcHC)[i]-> GetLength();
      gAnaMan.SetDCData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
    }
  }

  static const G4int id_ch = SDManager->GetCollectionID("SCH/hit");
  if( id_ch > 0 ){
    auto chHC = (G4THitsCollection<TPCSCHHit>*)HCTE->GetHC( id_ch );
    G4int nhits= chHC -> entries();
    for(G4int i=0; i< nhits; i++) {
      G4ThreeVector vtxpos = (*chHC)[i]-> GetVtxPosition();
      G4ThreeVector vtxmom = (*chHC)[i]-> GetVtxMomentum();
      G4double vtxene =(*chHC)[i]-> GetVtxEnergy();
      G4ThreeVector xyz = (*chHC)[i]-> GetPosition();
      G4ThreeVector mom = (*chHC)[i]-> GetMomentum();
      G4int ptid = (*chHC)[i]-> GetParentID();
      G4double tof= (*chHC)[i]-> GetTOF();
      G4int tid = (*chHC)[i]-> GetTrackID();
      G4int pid = (*chHC)[i]-> GetParticleID();
      G4int did = (*chHC)[i]-> GetDetectorID();
      G4double mass = (*chHC)[i]-> GetParticleMassID();
      G4int qq = (*chHC)[i]-> GetParticleQqID();
      G4double tlength = (*chHC)[i]-> GetLength();
      gAnaMan.SetSCHData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
    }
  }

  static const G4int id_ftof = SDManager->GetCollectionID("FTOF/hit");
  if( id_ftof > 0 ){
    auto ftofHC = (G4THitsCollection<TPCFTOFHit>*)HCTE->GetHC( id_ftof );
    G4int nhits = ftofHC->entries();
    for(G4int i=0; i< nhits; i++) {
      G4ThreeVector vtxpos = (*ftofHC)[i]-> GetVtxPosition();
      G4ThreeVector vtxmom = (*ftofHC)[i]-> GetVtxMomentum();
      G4double vtxene =(*ftofHC)[i]-> GetVtxEnergy();
      G4ThreeVector xyz = (*ftofHC)[i]-> GetPosition();
      G4ThreeVector mom = (*ftofHC)[i]-> GetMomentum();
      G4int ptid = (*ftofHC)[i]-> GetParentID();
      G4double tof= (*ftofHC)[i]-> GetTOF();
      G4int tid = (*ftofHC)[i]-> GetTrackID();
      G4int pid = (*ftofHC)[i]-> GetParticleID();
      G4int did = (*ftofHC)[i]-> GetDetectorID();
      G4double mass = (*ftofHC)[i]-> GetParticleMassID();
      G4int qq = (*ftofHC)[i]-> GetParticleQqID();
      G4double tlength = (*ftofHC)[i]-> GetLength();
      gAnaMan.SetFTOFData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
    }
  }

  auto trajectoryContainer = anEvent->GetTrajectoryContainer();
  if( trajectoryContainer && G4VVisManager::GetConcreteInstance() ){
    G4int n_trajectories = trajectoryContainer->entries();
    for( G4int i=0; i<n_trajectories; ++i ){
      auto trj = (G4Trajectory*)( (*(anEvent->GetTrajectoryContainer()))[i] );
      trj->DrawTrajectory();
    }
  }
  gAnaMan.EndOfEventAction();
}

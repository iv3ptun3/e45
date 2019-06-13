// ====================================================================
//   TPCEventAction.cc
//
// ====================================================================
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"

#include "TPCEventAction.hh"
#include "TPCTargetSD.hh"
#include "TPCPadSD.hh"
#include "TPCScintSD.hh"
#include "TPCACSD.hh"
#include "TPCDCSD.hh"
#include "TPCCHSD.hh"
#include "TPCNBARSD.hh"
#include "TPCFTOFSD.hh"

#include "TPCAnaManager.hh"

#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VTrajectory.hh"
#include <fstream>
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"


extern std::ofstream ofs;

////////////////////////////////
TPCEventAction::TPCEventAction(TPCAnaManager* ana)
  :AnaManager(ana)
////////////////////////////////
{}

/////////////////////////////////
TPCEventAction::~TPCEventAction()
/////////////////////////////////
{}

///////////////////////////////////////////////////////////////
void TPCEventAction::BeginOfEventAction(const G4Event* anEvent)
///////////////////////////////////////////////////////////////
{
  //  G4int nVtx= anEvent-> GetNumberOfPrimaryVertex();
  /*  for( G4int i=0; i< nVtx; i++) {
    const G4PrimaryVertex* primaryVertex= anEvent-> GetPrimaryVertex(i);
    }*/

  G4String env_ac_use = getenv("AC_Use");
  ac_use = atoi( env_ac_use.c_str() );

  G4String env_n_bar_use = getenv("N_Bar_Use");
  n_bar_use = atoi( env_n_bar_use.c_str() );

  G4String env_experiment_num = getenv("Experiment_NUM");
  experiment_num = atoi( env_experiment_num.c_str() );

  G4String env_with_kurama = getenv("With_KURAMA");
  with_kurama = atoi( env_with_kurama.c_str() );


  AnaManager->BeginOfEventAction();
}

/////////////////////////////////////////////////////////////
void TPCEventAction::EndOfEventAction(const G4Event* anEvent)
/////////////////////////////////////////////////////////////
{
  G4int eventID = anEvent-> GetEventID();
  //  G4cout<<"event #:"<<eventID<<G4endl;
  if(eventID % 100 == 0){
    G4cout << ">>> Event ID: " << eventID << G4endl;
  }
  G4SDManager* SDManager= G4SDManager::GetSDMpointer();

  // get "Hit Collection of This Event"
  G4HCofThisEvent* HCTE= anEvent-> GetHCofThisEvent();
  if(! HCTE) return;  // no hits in this events. nothing to do!


  static G4int idcounter= -1;
  if(idcounter<0) idcounter= SDManager-> GetCollectionID("TPC/hit");
  G4THitsCollection<TPCPadHit>* padHC = 0;
  padHC = (G4THitsCollection<TPCPadHit>*)HCTE-> GetHC(idcounter);

  if(padHC){
    G4int nhits= padHC -> entries();
    G4int pidtr[MaxTrackTPC]={0};
    G4int ptidtpc[MaxTrackTPC]={0};
    G4int ptidtpc_pid[MaxTrackTPC]={0};
    G4double pmtpc[MaxTrackTPC]={0};
    G4int qqtpc[MaxTrackTPC]={0};
    G4double pxtpc[MaxTrackTPC]={0};
    G4double pytpc[MaxTrackTPC]={0};
    G4double pztpc[MaxTrackTPC]={0};
    G4double pptpc[MaxTrackTPC]={0};

    G4double vtxxtpc[MaxTrackTPC]={0};
    G4double vtxytpc[MaxTrackTPC]={0};
    G4double vtxztpc[MaxTrackTPC]={0};
    G4double vtxpxtpc[MaxTrackTPC]={0};
    G4double vtxpytpc[MaxTrackTPC]={0};
    G4double vtxpztpc[MaxTrackTPC]={0};
    G4double vtxpptpc[MaxTrackTPC]={0};
    G4double vtxenetpc[MaxTrackTPC]={0};

    G4double detpc[MaxTrackTPC]={0};
    G4int laytpc[MaxTrackTPC]={0};
    G4double lentpc[MaxTrackTPC]={0};


    G4int nparticle=0;
    //            G4cout<<"pad nhits:"<<nhits<<G4endl;
    for(G4int i=0; i< nhits; i++) {
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
	vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));
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
	vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));

	vtxxtpc[nparticle]=vtxpos[0];
	vtxytpc[nparticle]=vtxpos[1];
	vtxztpc[nparticle]=vtxpos[2];

	vtxenetpc[nparticle]=vtxene;

	pptpc[nparticle]=sqrt(pow(mom[0],2.)+pow(mom[1],2.)+pow(mom[2],2.));
	lentpc[nparticle]=tlength;
	ptidtpc[nparticle]=ptid;
	ptidtpc_pid[nparticle]=ptid_pid;
	nparticle=nparticle+1;
	
      }else if( nparticle>0 )
	//	G4cout<<nparticle<<G4endl;
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
	vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));
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
	vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));

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
	AnaManager->SetCounterData(nparticle-1,tof, xyz, mom, tid, pid,ilay,irow,beta,edep/MeV,parentid,tlength,slength);
      }
    }
    for(G4int i=0;i<nparticle;i++){
      AnaManager->SetTPCData(i, pidtr[i], ptidtpc[i], ptidtpc_pid[i], pxtpc[i],pytpc[i],pztpc[i],pptpc[i], qqtpc[i], pmtpc[i], detpc[i], lentpc[i],laytpc[i], 
			     vtxpxtpc[i],vtxpytpc[i],vtxpztpc[i],
			     vtxxtpc[i],vtxytpc[i],vtxztpc[i], vtxenetpc[i]);
    }
  }
  //  G4cout<<"33333333333333333333333333"<<G4endl;
  static G4int idscint = -1;
  if(idscint < 0) idscint = SDManager-> GetCollectionID("SCINT/hit");
  G4THitsCollection<TPCScintHit>* scintHC = 0;
  scintHC = (G4THitsCollection<TPCScintHit>*)HCTE-> GetHC(idscint);

  if(scintHC){
    G4int nhits= scintHC -> entries();
    //    G4cout<<nhits<<G4endl;
    for(G4int i=0; i< nhits; i++) {
      G4ThreeVector vtxpos = (*scintHC)[i]-> GetVtxPosition();
      G4ThreeVector vtxmom = (*scintHC)[i]-> GetVtxMomentum();
      G4double vtxene =(*scintHC)[i]-> GetVtxEnergy();

      G4ThreeVector xyz = (*scintHC)[i]-> GetPosition();
      G4ThreeVector mom = (*scintHC)[i]-> GetMomentum();
      G4int ptid = (*scintHC)[i]-> GetParentID();
      G4double tof= (*scintHC)[i]-> GetTOF();
      G4int tid = (*scintHC)[i]-> GetTrackID();
      G4int pid = (*scintHC)[i]-> GetParticleID();
      G4int did = (*scintHC)[i]-> GetDetectorID();
      G4double mass = (*scintHC)[i]-> GetParticleMassID();
      G4int qq = (*scintHC)[i]-> GetParticleQqID();
      G4double tlength = (*scintHC)[i]-> GetLength();
      //      printf("didsc : %d \n", did);
      AnaManager->SetScintData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
      //      printf("didsc : %d \n", did);
    }
  }

  if(ac_use==1.){
    static G4int idac = -1;
    if(idac < 0) idac = SDManager-> GetCollectionID("AC/hit");
    G4THitsCollection<TPCACHit>* acHC = 0;
    acHC = (G4THitsCollection<TPCACHit>*)HCTE-> GetHC(idac);

    if(acHC){
      G4int nhits= acHC -> entries();
      //    G4cout<<nhits<<G4endl;
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
	//      printf("didsc : %d \n", did);
	AnaManager->SetACData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
	//      printf("didsc : %d \n", did);
      }
    }
  }



  if(n_bar_use==1.){
    static G4int idnbar = -1;
    if(idnbar < 0) idnbar = SDManager-> GetCollectionID("NBAR/hit");
    G4THitsCollection<TPCNBARHit>* nbarHC = 0;
    nbarHC = (G4THitsCollection<TPCNBARHit>*)HCTE-> GetHC(idnbar);

    if(nbarHC){
      G4int nhits= nbarHC -> entries();
      //    G4cout<<nhits<<G4endl;
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
	//      printf("didsc : %d \n", did);
	AnaManager->SetNBARData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
	//      printf("didsc : %d \n", did);
      }
    }
  }

  static G4int idtarget = -1;
  if(idtarget < 0) idtarget = SDManager-> GetCollectionID("TAR/hit");
  //  G4cout<<"idtarget:"<<idtarget<<G4endl;
  G4THitsCollection<TPCTargetHit>* targetHC = 0;
  targetHC = (G4THitsCollection<TPCTargetHit>*)HCTE-> GetHC(idtarget);
  //  G4cout<<targetHC<<G4endl;

  if(targetHC){
    G4int nhits= targetHC -> entries();
    //    G4cout<<"target nhits:"<<nhits<<G4endl;
    for(G4int i=0; i< nhits; i++) {
      G4double kinene =(*targetHC)[i]-> GetKinEnergy();;
      G4ThreeVector mom = (*targetHC)[i]-> GetMomentum();
      G4double test= sqrt(pow(mom.getX(),2)+pow(mom.getY(),2)+pow(mom.getZ(),2));
      //      G4cout<<i<<"th hits,"<<mom<<": kine energy:"<<kinene<<G4endl;
      if(kinene==0.00000000000000000000 && test != 0.  ){

	G4ThreeVector vtxpos = (*targetHC)[i]-> GetVtxPosition();
	G4ThreeVector vtxmom = (*targetHC)[i]-> GetVtxMomentum();
	G4double vtxene =(*targetHC)[i]-> GetVtxEnergy();
	
	G4ThreeVector xyz = (*targetHC)[i]-> GetPosition();

	
	G4double tof= (*targetHC)[i]-> GetTOF();
	G4int tid = (*targetHC)[i]-> GetTrackID();
	G4int ptid = (*targetHC)[i]-> GetParentID();
	G4int pid = (*targetHC)[i]-> GetParticleID();
	G4double mass = (*targetHC)[i]-> GetMass();
	G4int charge = (*targetHC)[i]-> GetCharge();
	//      G4double mass = (*targetHC)[i]-> GetPDGMass(); //mass(GeV)
	G4int parentid = (*targetHC)[i]-> GetParentID();
	G4double tlength = (*targetHC)[i]-> GettLength();
	
	G4int irow=0.;
	G4double beta = (*targetHC)[i]-> GetBeta();
	G4double edep = (*targetHC)[i]-> GetEdep();
	AnaManager->SetTargetData(i,xyz, mom, tid, pid,ptid,vtxpos,vtxmom,vtxene);
	//  void SetTargetData(G4ThreeVector pos, G4ThreeVector mom, 
	//		     G4int track, G4int particle,
	//		     G4int parentid, G4ThreeVector vtxpos, 
	//		     G4ThreeVector vtxmom, G4double vtxene);
      }
    //    G4cout<<target
    }
  }


  if(experiment_num==42||experiment_num==27){
  ////////////////DCs
  static G4int iddc = -1;
  if(iddc < 0) iddc = SDManager-> GetCollectionID("DC/hit");
  G4THitsCollection<TPCDCHit>* dcHC = 0;
  dcHC = (G4THitsCollection<TPCDCHit>*)HCTE-> GetHC(iddc);

  if(dcHC){
    G4int nhits= dcHC -> entries();
    //    G4cout<<nhits<<G4endl;
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
      //      printf("didsc : %d \n", did);
      AnaManager->SetDCData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
      //      printf("didsc : %d \n", did);
    }
  }


  ////////////////CH
  static G4int idch = -1;
  if(idch < 0) idch = SDManager-> GetCollectionID("CH/hit");
  G4THitsCollection<TPCCHHit>* chHC = 0;
  chHC = (G4THitsCollection<TPCCHHit>*)HCTE-> GetHC(idch);

  if(chHC){
    G4int nhits= chHC -> entries();
    //    G4cout<<nhits<<G4endl;
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
      //      printf("didsc : %d \n", did);
      AnaManager->SetCHData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
      //      printf("didsc : %d \n", did);
    }
  }



  ////////////////FTOF
  static G4int idftof = -1;
  if(idftof < 0) idftof = SDManager-> GetCollectionID("FTOF/hit");
  G4THitsCollection<TPCFTOFHit>* ftofHC = 0;
  ftofHC = (G4THitsCollection<TPCFTOFHit>*)HCTE-> GetHC(idftof);
  if(ftofHC){
    G4int nhits= ftofHC -> entries();
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
      //      printf("ftof : %d \n", did);
      //      printf("ftof id: %d \n", tid);
      AnaManager->SetFTOFData(tof, xyz, mom, tid, pid,did,mass,qq,ptid,vtxpos,vtxmom,vtxene,tlength);
      //      printf("didsc : %d \n", did);
    }
  }

  }



  G4TrajectoryContainer* trajectoryContainer = anEvent->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  if (1) {
    if (G4VVisManager::GetConcreteInstance())
      {
        for (G4int i=0; i<n_trajectories; i++)
          {
            G4Trajectory* trj = (G4Trajectory*)
              ((*(anEvent->GetTrajectoryContainer()))[i]);
            G4String particleName = trj->GetParticleDefinition()->GetParticleName();
	  trj->DrawTrajectory(50);      }
      }
  // Hits                                                                       
    //      if( HCTE ){                                                              
    //        if( padHC )  padHC->DrawAllHits();                                    
    //      }                                                                       
  }


  AnaManager->EndOfEventAction();
}


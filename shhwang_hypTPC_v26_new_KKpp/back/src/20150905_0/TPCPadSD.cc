// ====================================================================
//   TPCPadSD.cc
//
// ====================================================================
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "TPCPadSD.hh"
#include "TPCPadHit.hh"
#include "Randomize.hh"

////////////////////////////////////////////////
TPCPadSD::TPCPadSD(const G4String& name)
  : G4VSensitiveDetector(name)
////////////////////////////////////////////////
{
  collectionName.insert("hit");
}

/////////////////////////////////
TPCPadSD::~TPCPadSD()
/////////////////////////////////
{
}

////////////////////////////////////////////////
void TPCPadSD::Initialize(G4HCofThisEvent* HCTE)
////////////////////////////////////////////////
{
  ntrk=0;
  // create hit collection(s)
  hitsCollection = new G4THitsCollection<TPCPadHit>( SensitiveDetectorName,
						       collectionName[0]);

  // push H.C. to "Hit Collection of This Event"
  G4int hcid = GetCollectionID(0);
  HCTE-> AddHitsCollection(hcid, hitsCollection);


  if(env_gemdischarge>0){  
    deadarea = env_deadarea;
    if(env_gemfixdead==0){  
      select_plane=CLHEP::RandFlat::shoot(0.,1.);
      num_plane=0;
      if(select_plane<0.25){
	num_plane=1;
      }else if(select_plane>=0.25 && select_plane<0.5){
	num_plane=2;
      }else if(select_plane>=0.5 && select_plane<0.75){
	num_plane=3;
      }else if(select_plane>=0.75){
	num_plane=4;
      }
      num_deadarea = 250/deadarea;
      if(deadarea*num_deadarea < 250){
	num_deadarea=num_deadarea+1.;
      }
      select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);	
    }else if(env_gemfixdead==1){  
      num_plane=env_gemdeadplane;
      select_dead=env_gemdeadplanedivision;
    }
  } 
}

///////////////////////////////////////////////////////////
G4bool TPCPadSD::ProcessHits(G4Step* aStep, 
				G4TouchableHistory* ROhist)
///////////////////////////////////////////////////////////
{
  // get step information from "PreStepPoint"
  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();
  const G4Track* aTrack = aStep->GetTrack();
  G4int copyNo = preStepPoint -> GetPhysicalVolume()->GetCopyNo();

  if(preStepPoint-> GetStepStatus() != fGeomBoundary) return false;
  //  if(preStepPoint-> GetStepStatus() == fGeomBoundary){
  G4String particleName;
  if(aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge() == 0.)
    return false;
  particleName = aStep-> GetTrack()-> GetDefinition()-> GetParticleName();

  G4String particleType;
  particleType = aTrack->GetDefinition()->GetParticleType();

  ///lepton rejection
  if(particleType == "lepton")
    return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4ThreeVector pos= preStepPoint-> GetPosition();
  G4double hitx=pos.getX();
  G4double hitz=pos.getZ();

  //dead layer cross x and z
  /*
  if(abs(pos.getX())<5.)
    return false;
  if(abs(pos.getZ())<5.)
    return false;
  */
  //dead layer 45 deg 
  G4double deadlayer=12.;
  if(hitx < (hitz+deadlayer/2.*sqrt(2)) && hitx > (hitz-deadlayer/2.*sqrt(2))  )    return false;
  if(hitx < (-hitz+deadlayer/2.*sqrt(2)) && hitx > (-hitz-deadlayer/2.*sqrt(2))  )    return false;
  //////////end dead layer

  if(env_gemdischarge>0){
    if(env_gemdischarge==1){
      if(num_plane==1){
	if(hitz>=0 && fabs(hitx)<=hitz){	  ////select discharge area
	  if(hitz>deadarea*(select_dead-1) && hitz<deadarea*select_dead){
	    return false;
	  }
	}
      }else if(num_plane==2){
	if(hitx>=0 && hitx>fabs(hitz)){
	  if(hitx>deadarea*(select_dead-1) && hitx<deadarea*select_dead){
	    //	    G4cout<<"num_plane:select_dead--->"<<num_plane<<":"<<select_dead<<G4endl;
	    return false;
	  }
	}
      }else if(num_plane==3){
	if(hitz<0 && fabs(hitx)<fabs(hitz)){
	  //	  G4int select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);	
	  //	  G4cout<<"num_plane:select_dead--->"<<num_plane<<":"<<select_dead<<G4endl;
	  if(fabs(hitz)>deadarea*(select_dead-1) && fabs(hitz)<deadarea*select_dead){
	    return false;
	  }
	}
      }else if(num_plane==4){
	if(hitx<0 && fabs(hitx)>fabs(hitz)){
	  //	  G4int select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);	
	  //	  G4cout<<"num_plane:select_dead--->"<<num_plane<<":"<<select_dead<<G4endl;
	  if(fabs(hitx)>deadarea*(select_dead-1) && fabs(hitx)<deadarea*select_dead){
	    return false;
	  }
	}
      }
    }  ////end design #1

    ///////////////design 2
    if(env_gemdischarge==2){
      if(num_plane==1){
	if(hitz>=0 && fabs(hitx)<=hitz){	  ////select discharge area
	  if( hitx < -6*sqrt(2) + hitz - deadarea*sqrt(2)*(select_dead-1) 
		     && hitx > -6*sqrt(2) + hitz-deadarea*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==2){
	if(hitx>=0 && hitx>fabs(hitz)){
	  if( hitx > 6*sqrt(2)+ -hitz + deadarea*sqrt(2)*(select_dead-1) 
		     && hitx < 6*sqrt(2) + -hitz + deadarea*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==3){
	if(hitz<0 && fabs(hitx)<fabs(hitz)){
	  //	  G4int select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);	
	  //	  select_dead=3;
	  if( hitx > 6*sqrt(2)+ hitz + deadarea*sqrt(2)*(select_dead-1) 
		     && hitx < 6*sqrt(2)+ hitz + deadarea*sqrt(2)*(select_dead) ){
	    //	    G4cout<<"ok"<<G4endl;
	    return false;
	  }
	}
      }else if(num_plane==4){
	if(hitx<0 && fabs(hitx)>fabs(hitz)){
	  //	  G4int select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);	
	  //	  select_dead=3;
	  if( hitx < -6*sqrt(2) -hitz - deadarea*sqrt(2)*(select_dead-1) 
		     && hitx > -6*sqrt(2) -hitz - deadarea*sqrt(2)*(select_dead) ){
	    //	    G4cout<<"ok"<<G4endl;
	    return false;
	  }
	}
      }
    }  //end design #2

    ///////////////design 5, designed shape
    if(env_gemdischarge==5){
      G4double divi_width[6]={44.,43.,43.,43.,43.,43.0};
      //      G4double divi_width3[15]={25.,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,25.,25.,25.,25.};
      G4double divi_width3[20]={13.};
      //      for(G4int ii=0;ii<20;ii++){
      //	divi_width3[ii]=13.;
      //	//	G4cout<<divi_width3[ii]<<G4endl;
      //      }
      if(num_plane==1){
	if(select_dead>6){
	  G4cout<<"exceed maximum division"<<G4endl;
	  exit(-1);
	}
	if(hitz>=0 && fabs(hitx)<=hitz){	  ////select discharge area
	  if( hitx < -6.5*sqrt(2) + hitz - divi_width[select_dead-1]*sqrt(2)*(select_dead-1) 
		     && hitx > -6.5*sqrt(2) + hitz-divi_width[select_dead-1]*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==2){
	if(select_dead>6.){
	  G4cout<<"exceed maximum division"<<G4endl;
	  exit(-1);
	}
	if(hitx>=0 && hitx>fabs(hitz)){
	  if( hitx > 6.5*sqrt(2)+ -hitz + divi_width[select_dead-1]*sqrt(2)*(select_dead-1) 
		     && hitx < 6.5*sqrt(2) + -hitz + divi_width[select_dead-1]*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==3){
	//	if(select_dead>15.){
	if(select_dead>20.){
	  G4cout<<"exceed maximum division"<<G4endl;
	  exit(-1);
	}
	if(hitz<0 && fabs(hitx)<fabs(hitz)){
	  if( hitx > 6.5*sqrt(2)+ hitz + divi_width3[select_dead-1]*sqrt(2)*(select_dead-1) 
		     && hitx < 6.5*sqrt(2)+ hitz + divi_width3[select_dead-1]*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==4){
	if(select_dead>6.){
	  G4cout<<"exceed maximum division"<<G4endl;
	  exit(-1);
	}
	if(hitx<0 && fabs(hitx)>fabs(hitz)){
	  if( hitx < -6.5*sqrt(2) -hitz - divi_width[select_dead-1]*sqrt(2)*(select_dead-1) 
		     && hitx > -6.5*sqrt(2) -hitz - divi_width[select_dead-1]*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }
    }  //end design #5


    ///////////////design 4 anticlockwise
    //    G4cout<<env_gemdischarge<<G4endl;
    if(env_gemdischarge==4){
      if(num_plane==1){
	if(hitz>=0 && fabs(hitx)<=hitz){	  ////select discharge area
	  if( hitx > +6*sqrt(2) - hitz + deadarea*sqrt(2)*(select_dead-1) 
		     && hitx < +6*sqrt(2) - hitz + deadarea*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==2){
	if(hitx>=0 && hitx>fabs(hitz)){
	  if( hitx > 6*sqrt(2)+ +hitz + deadarea*sqrt(2)*(select_dead-1) 
	      && hitx < 6*sqrt(2) + +hitz + deadarea*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==3){
	if(hitz<0 && fabs(hitx)<fabs(hitz)){
	  //	    G4cout<<"n"<<G4endl;
	  if( hitx < -6*sqrt(2)- hitz - deadarea*sqrt(2)*(select_dead-1) 
		     && hitx > -6*sqrt(2) - hitz - deadarea*sqrt(2)*(select_dead) ){
	    //	    G4cout<<"ok"<<G4endl;
	    return false;
	  }
	}
      }else if(num_plane==4){
	if(hitx<0 && fabs(hitx)>fabs(hitz)){
	  if( hitx < -6*sqrt(2) +hitz - deadarea*sqrt(2)*(select_dead-1) 
		     && hitx > -6*sqrt(2) +hitz - deadarea*sqrt(2)*(select_dead) ){
	    //	    G4cout<<"ok"<<G4endl;
	    return false;
	  }
	}
      }
    }  //end design #4 anticlockwise


///////////////design 3, test on plane #3 
    if(env_gemdischarge==3){
      if(num_plane==1){
	if(hitz>=0 && fabs(hitx)<=hitz){	  ////select discharge area
	  if( hitx < hitz - deadarea*sqrt(2)*(select_dead-1) 
		     && hitx > hitz-deadarea*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==2){
	if(hitx>=0 && hitx>fabs(hitz)){
	  if( hitx > -hitz + deadarea*sqrt(2)*(select_dead-1) 
		     && hitx < -hitz + deadarea*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==3){
	if(hitz<0 && fabs(hitx)<fabs(hitz)){
	  
	  if( hitx > hitz + deadarea*sqrt(2)*(select_dead-1) 
		     && hitx < hitz + deadarea*sqrt(2)*(select_dead) && hitx<0.){
		//	    G4cout<<num_plane<<":"<<deadarea<<":"<<select_dead<<G4endl;
	    return false;
	  }
	}
      }else if(num_plane==4){
	if(hitx<0 && fabs(hitx)>fabs(hitz)){
	  //	  G4int select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);	
	  //	  select_dead=3;
	  if( hitx < -hitz - deadarea*sqrt(2)*(select_dead-1) 
		     && hitx > -hitz - deadarea*sqrt(2)*(select_dead) ){
	    //	    G4cout<<"ok"<<G4endl;
	    return false;
	  }
	}
      }
    }  //end design #3 test



      //  if(hitx < (hitz+5*sqrt(2)) && hitx > (hitz-5*sqrt(2))  )
      //    return false;
      //  if(hitx < (-hitz+5*sqrt(2)) && hitx > (-hitz-5*sqrt(2))  )
      //    return false;
    }//end discharge
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep<0.0005){
    return false;
  }
  G4ThreeVector VertexPosition = aTrack->GetVertexPosition();
  G4ThreeVector VertexMomentum = aTrack->GetVertexMomentumDirection();
  G4double VertexEnergy = aTrack -> GetVertexKineticEnergy(); // Ek = sqrt(p^2+m^2)-m
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 

  G4ThreeVector mom= preStepPoint-> GetMomentum();
  G4double tof= preStepPoint-> GetGlobalTime();
  G4double beta= preStepPoint-> GetBeta();
  G4int tid =  aStep-> GetTrack()-> GetTrackID();
  G4int pid =  aStep-> GetTrack()-> GetDefinition() -> GetPDGEncoding();
  G4double mass = aStep -> GetTrack()->GetDynamicParticle()->GetMass();
  //  G4cout<<mass<<G4endl;
  G4int charge = aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge();
  //  G4double edep = aStep->GetTotalEnergyDeposit();
  //  G4double edep = aStep->GetNonIonizingEnergyDeposit();
  //  G4double deltaep = aStep->GetDeltaEnergy();
  G4double tlength = aStep->GetTrack()-> GetTrackLength();
  G4double slength = aStep->GetTrack()-> GetStepLength();
  G4double length = aStep-> GetStepLength();
  G4int parentID =  aStep-> GetTrack()-> GetParentID();

  //  G4int test =  aStep-> GetTrack()-> GetParentID() -> GetPDGCharge();
  //  G4cout<<"-------------------"<<G4endl;
  //  G4cout<<"pid :"<<pid<<G4endl;
  //  G4cout<<"test :"<<tlength<<G4endl;
  //  G4cout<<"steplength :"<<slength<<G4endl;
  //  G4cout<<"length :"<<length<<G4endl;
  //  G4cout<<"dedx :"<<edep<<std::setw(2)<<":"<<deltaep<<G4endl;
  //  G4cout<<"copyNo :"<<copyNo<<G4endl;
  //  G4cout<<"tlength"<<tlength<<G4endl; 

  G4int iLay=copyNo;
  G4int iRow=0;
  G4String name = physVol->GetName();

  if(name=="TPC_PV"){
    name="PadPV-1";
  }

  //  if(name=="Target_PV"){
  //    name="PadPV-1";
  //  }

  sscanf(name,"PadPV%d",&iLay);
  TPCPadHit* ahit= new TPCPadHit(pos, mom, tof, tid, pid, iLay, iRow, beta, edep,parentID,tlength, mass, charge, VertexPosition, VertexMomentum, VertexEnergy,slength );

  hitsCollection-> insert(ahit);
  return true;
  //}End of fboundary

}

////////////////////////////////////////////////////
void TPCPadSD::EndOfEvent(G4HCofThisEvent* HCTE)
////////////////////////////////////////////////////
{
}

////////////////////////////
void TPCPadSD::DrawAll()
////////////////////////////
{
}

/////////////////////////////
void TPCPadSD::PrintAll()
/////////////////////////////
{
  hitsCollection-> PrintAllHits();
}


////////////////////////////
void TPCPadSD::TPCPadSD_Set()
////////////////////////////
{

  G4String GEMDischarge = getenv("GEMDischarge");
  env_gemdischarge=atoi(GEMDischarge.c_str());
  G4String DeadArea = getenv("DeadArea");
  env_deadarea = atof(DeadArea.c_str());
  G4String GEMFixDead = getenv("GEMFixDead");
  env_gemfixdead = atoi(GEMFixDead.c_str());
  G4String GEMDeadPlane = getenv("GEMDeadPlane");
  env_gemdeadplane=atoi(GEMDeadPlane.c_str());
  G4String GEMDeadPlaneDivision = getenv("GEMDeadPlaneDivision");
  env_gemdeadplanedivision=atoi(GEMDeadPlaneDivision.c_str());
  G4cout<<"Study on GEM discharge:"<<env_gemdischarge<<G4endl;
  if(env_gemdischarge==5){
    G4cout<<"Designed a configuration of GEM discharge area"<<G4endl;
    G4cout<<"Width of GEM electrodes are fixed!!!!"<<G4endl;
  }    
}

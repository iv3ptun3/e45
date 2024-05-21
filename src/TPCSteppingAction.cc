// -*- C++ -*-

#include "TPCSteppingAction.hh"

#include <G4Material.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4SteppingManager.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>
#include <G4TrackStatus.hh>
#include <G4VPhysicalVolume.hh>
#include <TString.h>
#include "ConfMan.hh"
namespace
{
  const auto& gConf = ConfMan::GetInstance();
	auto& gTrackBuffer = TPCTrackBuffer::GetInstance();
}

//_____________________________________________________________________________
TPCSteppingAction::TPCSteppingAction( void )
{
}

//_____________________________________________________________________________
TPCSteppingAction::~TPCSteppingAction( void )
{
}

//_____________________________________________________________________________
void
TPCSteppingAction::UserSteppingAction( const G4Step* theStep )
{
  static const G4bool KillStepInIron = gConf.Get<G4bool>( "KillStepInIron" );

  auto theTrack = theStep->GetTrack();
	int TrackID = theTrack->GetTrackID();
  auto theParticle = theTrack->GetParticleDefinition();
  int ParticleID = theParticle->GetPDGEncoding();
	auto particleName = theParticle->GetParticleName();
	auto prePoint = theStep->GetPreStepPoint();
	auto prePosition = prePoint->GetPosition();
	auto preMomentum = prePoint->GetMomentum();
  auto prePV = prePoint->GetPhysicalVolume();
  auto prePVName = prePV->GetName();
  auto postPoint = theStep->GetPostStepPoint();
	auto postMomentum = postPoint->GetMomentum();
  auto theProcess = postPoint->GetProcessDefinedStep()->GetProcessName();
  auto NStep = theTrack->GetCurrentStepNumber();
	int parentID = theTrack->GetParentID();
  // check if it is alive
  //  if( theTrack->GetTrackStatus() != fAlive ) { return; }
	//std::cout<<Form(" Processing track %d, Step %d" ,TrackID,NStep)<<std::endl;
  // check if it is primary
	if(particleName == "xi-"){
		if(NStep==1){
			gTrackBuffer.SetVertexMomentum(preMomentum/CLHEP::GeV,0);		
			gTrackBuffer.SetVertexLV(G4LorentzVector(preMomentum/CLHEP::GeV,hypot(theParticle->GetPDGMass()/CLHEP::GeV,preMomentum.mag()/CLHEP::GeV)),0);	
		}
		gTrackBuffer.SetMomentum(postMomentum/CLHEP::GeV,0);		
		gTrackBuffer.SetLV(G4LorentzVector(postMomentum/CLHEP::GeV,hypot(theParticle->GetPDGMass()/CLHEP::GeV,postMomentum.mag()/CLHEP::GeV)),0);	
//		G4cout<<"pre :"<<Form("(%g,%g,%g)",(preMomentum).x(),preMomentum.y(),preMomentum.z())<<G4endl;
//		G4cout<<"post :"<<Form("(%g,%g,%g)",(postMomentum).x(),postMomentum.y(),postMomentum.z())<<G4endl;
//		G4cout<<"dif :"<<(preMomentum-postMomentum).mag()<<G4endl;
	}
	if(particleName == "lambda" ){
		if(NStep==1){
			gTrackBuffer.SetVertexMomentum(preMomentum/CLHEP::GeV,1);		
			gTrackBuffer.SetVertexLV(G4LorentzVector(preMomentum/CLHEP::GeV,hypot(theParticle->GetPDGMass()/CLHEP::GeV,preMomentum.mag()/CLHEP::GeV)),1);	
		}
		gTrackBuffer.SetMomentum(postMomentum/CLHEP::GeV,1);		
		gTrackBuffer.SetLV(G4LorentzVector(postMomentum/CLHEP::GeV,hypot(theParticle->GetPDGMass()/CLHEP::GeV,postMomentum.mag()/CLHEP::GeV)),1);	
		gTrackBuffer.SetLambdaID(TrackID);
	}
	if(particleName == "proton" and  parentID == gTrackBuffer.GetLambdaID()){
		if(NStep==1){
			gTrackBuffer.SetVertexMomentum(preMomentum/CLHEP::GeV,2);		
			gTrackBuffer.SetVertexLV(G4LorentzVector(preMomentum/CLHEP::GeV,hypot(theParticle->GetPDGMass()/CLHEP::GeV,preMomentum.mag()/CLHEP::GeV)),2);	
		}
		gTrackBuffer.SetMomentum(postMomentum/CLHEP::GeV,2);
		gTrackBuffer.SetLV(G4LorentzVector(postMomentum/CLHEP::GeV,hypot(theParticle->GetPDGMass()/CLHEP::GeV,postMomentum.mag()/CLHEP::GeV)),2);
	}
	if(TrackID<1000){	
		if(NStep==1){
			gTrackBuffer.SetTrack(TrackID,parentID,ParticleID,prePosition,preMomentum);
		}
	}
	if( theTrack->GetParentID() != 0 ){
		auto CreationProcess = theTrack->GetCreatorProcess()->GetProcessName();
//		G4cout<<CreationProcess<<G4endl;
		if(CreationProcess != "Decay"){
			theTrack->SetTrackStatus( fStopAndKill );
			return; 
		}
		if(NStep > 2000){
			theTrack->SetTrackStatus( fStopAndKill );
			return; 
		}
	}
	if( particleName == "gamma" ) { 
		theTrack->SetTrackStatus( fStopAndKill );
		return; 
	}
	/*
	if(hypot(prePosition.x(),prePosition.z()) > 400 and theTrack->GetParentID() != 0){
		theTrack->SetTrackStatus( fStopAndKill );
		return; 
	}
	*/
  // check if it is NOT muon
  //  auto definition = theTrack->GetDefinition();
  //  if( ( definition == G4MuonPlus::MuonPlusDefinition() ) ||
  //      ( definition == G4MuonMinus::MuonMinusDefinition() ) )
  //  { return; }

  //  G4cout<<"start stepping action:"<<prePVName<<G4endl;

#ifdef DEBUG
  if( theProcess != "eIoni" &&
      theProcess != "hIoni" &&
      theProcess != "msc" &&
      theProcess != "eBeam" &&
      theProcess != "Transportation" ){
    G4cout << particleName << " " << theProcess << G4endl;
  }
#endif

  if( KillStepInIron ){
    auto preMaterial = prePoint->GetMaterial();
    if( preMaterial->GetName() == "Iron" ){
      theTrack->SetTrackStatus( fStopAndKill );
      return;
    }
  }

  // if( prePVName.contains( "Coil" ) || prePVName.contains( "Guard" ) ){
  //   theTrack->SetTrackStatus( fStopAndKill );
  //   return;
  // }
}

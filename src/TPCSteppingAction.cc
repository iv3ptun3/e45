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
  auto theTrack = theStep->GetTrack();
  auto theParticle = theTrack->GetParticleDefinition();
  auto particleName = theParticle->GetParticleName();
  auto prePoint = theStep->GetPreStepPoint();
  auto prePV = prePoint->GetPhysicalVolume();
  auto prePVName = prePV->GetName();
  auto preMaterial = prePoint->GetMaterial();
  auto postPoint = theStep->GetPostStepPoint();
  auto theProcess = postPoint->GetProcessDefinedStep()->GetProcessName();

  // check if it is alive
  //  if( theTrack->GetTrackStatus() != fAlive ) { return; }

  // check if it is primary
  //  if( theTrack->GetParentID() != 0 ) { return; }

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

  if( preMaterial->GetName() == "Iron" ){
    theTrack->SetTrackStatus( fStopAndKill );
    return;
  }

  if( prePVName.contains( "Coil" ) ){
    theTrack->SetTrackStatus( fStopAndKill );
    return;
  }

  // if( prePVName.contains( "Target" ) )
  //   G4cout<<"end stepping action = "<<prePVName<<G4endl;
  //  G4StepPoint * thePostPoint = theStep->GetPostStepPoint();
  //  G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume();
  //  G4String thePostPVname = thePostPV->GetName();
  //  if(thePostPVname(0,4)!="calo") { return; }
  // then suspend the track
  //  theTrack->SetTrackStatus(fSuspend);
}

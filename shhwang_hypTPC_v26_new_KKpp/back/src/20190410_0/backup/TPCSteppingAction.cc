/*
  SteppingAction.cc
  2007/4  K.Shirotori

  2010/4/24 T.Takahashi
*/

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "ConfMan.hh"

//const double TargetCenter = 40.0;

SteppingAction::SteppingAction( DetectorConstruction* det,
				EventAction* evt )
  : detector(det), eventaction(evt)
{ }

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* theStep )
{ 
  ConfMan *confMan = ConfMan::GetConfManager();
  int StepFlag = confMan->StepFlag();

  if(!StepFlag) return;

  G4Track * theTrack = theStep->GetTrack();

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(theStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 

#if 0
  G4ThreeVector pos = theTrack->GetPosition();

  G4cerr << physVol->GetName() << " pos=" << pos << G4endl;
 

  //  G4StepPoint *preStepPoint = theStep->GetPreStepPoint();
  //  G4TouchableHandle Touchable = preStepPoint->GetTouchableHandle();
  //  G4ThreeVector posl = Touchable->GetHistory()->

  //    GetTopTransform().TransformPoint( pos );

#endif

  // All particle which hits SKS is killed.
  if( GetFStopSKS(StepFlag) && 
      detector->IsVolumeStopper(physVol) ){
    theTrack->SetTrackStatus( fStopAndKill );
    return;
  }

  G4String particle = theTrack->GetDefinition()->GetParticleName();

  // Gamma ray which hits SKS is killed.
  if( GetFStopSKSGam(StepFlag) && particle=="gamma" &&
      detector->IsVolumeStopper(physVol) ){
    theTrack->SetTrackStatus( fStopAndKill );
    return;
  }

  // All neutrino are killed.
  if( GetFStopNu(StepFlag) &&
      ( particle =="anti_nu_e" || particle =="anti_nu_mu" 
	|| particle =="nu_e" || particle =="nu_mu" ) ){
    theTrack->SetTrackStatus(fStopAndKill);
    return;
  }
  // Gamma is killed.
  if( GetFStopGam(StepFlag) && particle =="gamma" ){
    theTrack->SetTrackStatus(fStopAndKill);
    return;
  }
  // e+/e- are killed.
  if( GetFStopE(StepFlag) && ( particle =="e-" || particle =="e+" ) ){
    theTrack->SetTrackStatus(fStopAndKill);
    return;
  }
  
  return;
}


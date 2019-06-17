//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "TPCSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Material.hh"

TPCSteppingAction::TPCSteppingAction()
{;}

TPCSteppingAction::~TPCSteppingAction()
{;}

void TPCSteppingAction::UserSteppingAction(const G4Step * theStep)
{
  auto theTrack = theStep->GetTrack();
  auto particle = theTrack->GetParticleDefinition();
  auto particleName = particle->GetParticleName();


  // check if it is alive
  //  if(theTrack->GetTrackStatus()!=fAlive) { return; }

  // check if it is primary
  //  if(theTrack->GetParentID()!=0) { return; }

  // check if it is NOT muon
  //  G4ParticleDefinition * particleType = theTrack->GetDefinition();
  //  if((particleType==G4MuonPlus::MuonPlusDefinition())
  //   ||(particleType==G4MuonMinus::MuonMinusDefinition()))
  //  { return; }

  // check if particles enters to the calorimeter volume, the process will stop.
  G4StepPoint * thePrePoint = theStep->GetPreStepPoint();
  G4VPhysicalVolume * thePrePV = thePrePoint->GetPhysicalVolume();
  G4String thePrePVname = thePrePV->GetName();
  //  G4cout<<"start stepping action:"<<thePrePVname<<G4endl;

  G4Material * material =  thePrePoint -> GetMaterial();
  auto postStepPoint = theStep->GetPostStepPoint();
  auto process = postStepPoint->GetProcessDefinedStep()->GetProcessName();

#if 0
  if( process != "eIoni" &&
      process != "hIoni" &&
      process != "msc" &&
      process != "eBeam" &&
      process != "Transportation" ){
    G4cout << particleName << " " << process << G4endl;
  }
#endif

  G4String m_name=material->GetName();
  if(m_name=="Iron" ){
    theTrack->SetTrackStatus( fStopAndKill );
    return;
  }

  if(thePrePVname=="HelmPV"){
    theTrack->SetTrackStatus( fStopAndKill );
    return;
  }
  //  G4cout<<"end stepping action:"<<thePrePVname<<G4endl;
  //  G4StepPoint * thePostPoint = theStep->GetPostStepPoint();
  //  G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume();
  //  G4String thePostPVname = thePostPV->GetName();
  //  if(thePostPVname(0,4)!="calo") { return; }
  // then suspend the track
  //  theTrack->SetTrackStatus(fSuspend);
}

// -*- C++ -*-

#include "TPCTargetVPSD.hh"

#include <G4Step.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VTouchable.hh>
#include "TString.h"
#include "FuncName.hh"
#include "TPCTargetVPHit.hh"

//_____________________________________________________________________________
TPCTargetVPSD::TPCTargetVPSD( const G4String& name )
  : G4VSensitiveDetector( name ),
    m_hits_collection()
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
TPCTargetVPSD::~TPCTargetVPSD( void )
{
}

//_____________________________________________________________________________
void
TPCTargetVPSD::Initialize( G4HCofThisEvent* HCTE )
{
  m_hits_collection = new G4THitsCollection<TPCTargetVPHit>( SensitiveDetectorName,
							 collectionName[0] );
  HCTE->AddHitsCollection( GetCollectionID(0), m_hits_collection );
}

//_____________________________________________________________________________
G4bool
TPCTargetVPSD::ProcessHits( G4Step* aStep, G4TouchableHistory* /* ROhist */ )
{
  const auto preStepPoint = aStep->GetPreStepPoint();
  const auto aTrack = aStep->GetTrack();
  const auto Definition = aTrack->GetDefinition();
  const G4String particleName = Definition->GetParticleName();
  const G4String particleType = Definition->GetParticleType();
	auto pos = preStepPoint->GetPosition(); 	
  const auto postStepPoint = aStep->GetPostStepPoint();
	auto postpos = postStepPoint->GetPosition(); 	
//	G4cout<<Form("VPHit! (%g,%g,%g) ",pos.x(),pos.y(),pos.z())<<G4endl;
  if( preStepPoint->GetStepStatus() != fGeomBoundary )
    return false;
  if( Definition->GetPDGCharge() == 0. )
    return false;

  // if( particleName == "e-" )
  //   return false;
  // if( particleName == "e+" )
  //   return false;
  // if( particleName != "kaon+" )
  //   return false;
  // if( particleName != "pi-" && particleName != "pi+" )
  //   return false;
  // if( particleName != "pi+" && particleName != "pi-" &&
  //     particleName != "proton" )
  //   return false;
  // if( particleType == "lepton" )
  //   return false;

  m_hits_collection->insert( new TPCTargetVPHit( SensitiveDetectorName, aStep ) );

  return true;
}

//_____________________________________________________________________________
void
TPCTargetVPSD::EndOfEvent( G4HCofThisEvent* /* HCTE */ )
{
}

//_____________________________________________________________________________
void
TPCTargetVPSD::DrawAll( void )
{
}

//_____________________________________________________________________________
void
TPCTargetVPSD::PrintAll( void )
{
  m_hits_collection->PrintAllHits();
}

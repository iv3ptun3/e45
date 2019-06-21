// -*- C++ -*-

#include "TPCCHSD.hh"

#include <G4VPhysicalVolume.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VTouchable.hh>
#include <G4TouchableHistory.hh>
#include "TPCCHHit.hh"

//_____________________________________________________________________________
TPCCHSD::TPCCHSD( const G4String& name )
  : G4VSensitiveDetector( name )
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
TPCCHSD::~TPCCHSD( void )
{
}

//_____________________________________________________________________________
void
TPCCHSD::Initialize( G4HCofThisEvent* HCTE )
{
  hitsCollection = new G4THitsCollection<TPCCHHit>( SensitiveDetectorName,
						    collectionName[0] );
  G4int hcid = GetCollectionID(0);
  HCTE->AddHitsCollection( hcid, hitsCollection );
}

//_____________________________________________________________________________
G4bool
TPCCHSD::ProcessHits( G4Step* aStep, G4TouchableHistory* /* ROhist */ )
{
  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();
  G4String particleName;

  if(aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge() == 0.)
    return false;

  particleName = aStep-> GetTrack()-> GetDefinition()-> GetParticleName();
  /*
  // e+/e- rejection
  if( particleName == "e-")
    return false;
  if( particleName == "e+")
    return false;
  */
  const G4Track* aTrack = aStep->GetTrack();
  G4String particleType;
  particleType = aTrack->GetDefinition()->GetParticleType();

    if(particleType == "lepton")
      return false;
  //
  //  if( (particleName != "kaon+")){
  //    return false;
  //  }

  //  if( (particleName != "pi-") && (particleName != "pi+")){
  //    return false;
  //  }
  //  if((particleName != "pi+")&& (particleName != "pi-")
  //     && (particleName != "proton")){
  //    return false;
  //  }
  if(preStepPoint-> GetStepStatus() != fGeomBoundary) return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();


  G4ThreeVector VertexPosition = aTrack->GetVertexPosition();
  G4ThreeVector VertexMomentum = aTrack->GetVertexMomentumDirection();
  G4double VertexEnergy = aTrack -> GetVertexKineticEnergy(); // Ek = sqrt(p^2+m^2)-m

  G4ThreeVector pos= preStepPoint-> GetPosition();
  G4ThreeVector mom= preStepPoint-> GetMomentum();
  G4double tof= preStepPoint-> GetGlobalTime();
  G4int tid =  aStep-> GetTrack()-> GetTrackID();
  G4int pid =  aStep-> GetTrack()-> GetDefinition() -> GetPDGEncoding();
  G4double mass =  aStep-> GetTrack()-> GetDynamicParticle() -> GetMass();
  G4int qq =  aStep-> GetTrack()-> GetDynamicParticle() -> GetCharge();
  G4double tlength = aStep->GetTrack()-> GetTrackLength();
  //  G4double slength = aStep->GetTrack()-> GetStepLength();
  //  G4cout<<mass<<":"<<tlength<<":"<<slength<<G4endl;
  G4int parentID =  aStep-> GetTrack()-> GetParentID();
  //    G4cout <<"test: "<<parentID << G4endl;
  // Get Pad number
  G4int iDet;
  G4String name = physVol->GetName();
  //  G4cout << name << G4endl;

  //  sscanf(name,"CHPV%d",&iDet);
  G4int copyNo = preStepPoint -> GetPhysicalVolume()->GetCopyNo();
  iDet=copyNo;

  // create a new hit and push them to "Hit Coleltion"
  TPCCHHit* ahit= new TPCCHHit(pos, mom, tof, tid, pid, iDet, mass, qq,
			       parentID, VertexPosition, VertexMomentum,
			       VertexEnergy, tlength );
  hitsCollection-> insert(ahit);

  return true;
}

//_____________________________________________________________________________
void
TPCCHSD::EndOfEvent( G4HCofThisEvent* /* HCTE */ )
{
}

//_____________________________________________________________________________
void
TPCCHSD::DrawAll( void )
{
}

//_____________________________________________________________________________
void
TPCCHSD::PrintAll( void )
{
  hitsCollection-> PrintAllHits();
}

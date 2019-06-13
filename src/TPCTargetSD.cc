// ====================================================================
//   TPCTargetSD.cc
//
// ====================================================================
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "TPCTargetSD.hh"
#include "TPCTargetHit.hh"

////////////////////////////////////////////////
TPCTargetSD::TPCTargetSD(const G4String& name)
  : G4VSensitiveDetector(name)
////////////////////////////////////////////////
{
  collectionName.insert("hit");
}

/////////////////////////////////
TPCTargetSD::~TPCTargetSD()
/////////////////////////////////
{
}


////////////////////////////////////////////////
void TPCTargetSD::Initialize(G4HCofThisEvent* HCTE)
////////////////////////////////////////////////
{
  // create hit collection(s)
  hitsCollection = new G4THitsCollection<TPCTargetHit>( SensitiveDetectorName,
						       collectionName[0]);

  // push H.C. to "Hit Collection of This Event"
  G4int hcid = GetCollectionID(0);
  HCTE-> AddHitsCollection(hcid, hitsCollection);
}

///////////////////////////////////////////////////////////
G4bool TPCTargetSD::ProcessHits(G4Step* aStep, 
				G4TouchableHistory* ROhist)
///////////////////////////////////////////////////////////
{

  // get step information from "PreStepPoint"
  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();
  const G4Track* aTrack = aStep->GetTrack();

  //  if(preStepPoint-> GetStepStatus() != fGeomBoundary) return false;

  G4String particleName;
  if(aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge() == 0.)
    return false;
  particleName = aStep-> GetTrack()-> GetDefinition()-> GetParticleName();

  G4String particleType;
  particleType = aTrack->GetDefinition()->GetParticleType();

  ///lepton rejection
  if(particleType == "lepton")
    return false;


  //  G4cout<<"test2"<<G4endl;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 


  //  G4cout<<"test3"<<G4endl;
  G4ThreeVector pos= preStepPoint-> GetPosition();
  G4ThreeVector VertexPosition = aTrack->GetVertexPosition();

  G4ThreeVector VertexMomentum = aTrack->GetVertexMomentumDirection();
  G4double VertexEnergy = aTrack -> GetVertexKineticEnergy(); // Ek = sqrt(p^2+m^2)-m

  //  G4cout<<"test4"<<G4endl;

  G4ThreeVector mom= preStepPoint-> GetMomentum();
  G4double tof= preStepPoint-> GetGlobalTime();
  G4double beta= preStepPoint-> GetBeta();
  G4int tid =  aStep-> GetTrack()-> GetTrackID();
  G4int pid =  aStep-> GetTrack()-> GetDefinition() -> GetPDGEncoding();
  G4double mass = aStep -> GetTrack()->GetDynamicParticle()->GetMass();
  G4int charge = aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge();
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double tlength = aStep->GetTrack()-> GetTrackLength();
  G4int parentID =  aStep-> GetTrack()-> GetParentID();

  //  G4cout<<"test1"<<G4endl;

  ////kinetic energy (2012. 6. 30.)
  //  G4ThreeVector momentum = aTrack->GetMomentum();
  //  G4double kinEnergy     = aTrack->GetKineticEnergy();
  //  G4double kinEnergy     = aTrack->GetKinEnergy();
  G4double kinEnergy     = aTrack->GetKineticEnergy();
  //  G4double globalTime    = track->GetGlobalTime();
  
  //  G4cout<<pid<<":"<<kinEnergy<<G4endl;
  G4int iLay=0;
  G4int iRow=0;
  G4String name = physVol->GetName();

  TPCTargetHit* ahit= new TPCTargetHit(pos, mom, tid, pid, parentID, mass, charge, VertexPosition, VertexMomentum, VertexEnergy, kinEnergy);

  hitsCollection-> insert(ahit);
  return true;

}

////////////////////////////////////////////////////
void TPCTargetSD::EndOfEvent(G4HCofThisEvent* HCTE)
////////////////////////////////////////////////////
{
}

////////////////////////////
void TPCTargetSD::DrawAll()
////////////////////////////
{
}

/////////////////////////////
void TPCTargetSD::PrintAll()
/////////////////////////////
{
  hitsCollection-> PrintAllHits();
}

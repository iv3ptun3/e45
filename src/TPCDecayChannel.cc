#ifndef TPCDecayChannel_h
#define TPCDecayChannel
#include "TPCDecayChannel.hh"
#include <G4ParticleTable.hh>
#include "ConfMan.hh"
#include "TString.h"
namespace
{
  const auto& gConf = ConfMan::GetInstance();
	auto& gTrackBuffer = TPCTrackBuffer::GetInstance();
}
TPCPolarizedDecayChannel
::TPCPolarizedDecayChannel(
				G4String ParentName,
				G4double BR,
				G4double Alpha,
				G4int Order,
				G4String Daughter1Name,
				G4String Daughter2Name){
	order = Order;
	Initialize(ParentName,BR,Alpha,
			Daughter1Name,Daughter2Name);
}

void
TPCPolarizedDecayChannel
::Initialize(
				G4String ParentName,
				G4double BR,
				G4double Alpha,
				G4String Daughter1Name,
				G4String Daughter2Name){
	
	rbranch = BR;
	parent_name = new G4String(ParentName);
	daughters_name = new G4String*[2];
	daughters_name[0] = new G4String(Daughter1Name);
	daughters_name[1] = new G4String(Daughter2Name);
	numberOfDaughters = 2;
	PDFCTheta = TF1("PDF","1 + [0]*[1]*x",-1,1);	
	PDFCTheta.SetParameter(0,Alpha);
	PDFCTheta.SetParameter(1,1);

	Daughters[0] = G4ParticleTable::GetParticleTable()->FindParticle(Daughter1Name);
	Daughters[1] = G4ParticleTable::GetParticleTable()->FindParticle(Daughter2Name);
	
	DaughterMass[0] = Daughters[0]->GetPDGMass();
	DaughterMass[1] = Daughters[1]->GetPDGMass();
	parent = G4ParticleTable::GetParticleTable()->FindParticle(ParentName);

	ParentMass = parent->GetPDGMass();
	

}
G4RotationMatrix
TPCPolarizedDecayChannel::SpinToMomentum(){
	double th_m = Polarity.theta();
	double ph_m = Polarity.phi();
	G4RotationMatrix Rot;
	Rot.rotateZ(-ph_m);
	Rot.rotateY(th_m);
	auto mom_rot = Rot*MomVector;
	auto ph_v = mom_rot.phi();
	Rot.rotateZ(-ph_v);//p,e,s to x,y,z
	auto RollBack = Rot;
	G4RotationMatrix SpinToMom;
	SpinToMom.rotateY(-acos(-1)/2);
	return Rot*SpinToMom;
}

G4DecayProducts*
TPCPolarizedDecayChannel::DecayIt(G4double ParentMas){
	
	LoadPolarityMomentum();	
	G4double P_d = Pmx(ParentMass,DaughterMass[0],DaughterMass[1]);
	double PI = acos(-1);
	double Phi = 2*PI*G4UniformRand();
	double Theta = acos(PDFCTheta.GetRandom());
//	Theta = MomVector.theta();
//	Phi = MomVector.phi();
	double px = P_d*sin(Theta)*cos(Phi);
	double py = P_d*sin(Theta)*sin(Phi);
	double pz = P_d*cos(Theta);
	G4ThreeVector DaughterMom1(px,py,pz);
	G4ThreeVector DaughterMom2= -DaughterMom1;

	auto STM = SpinToMomentum();
	DaughterMom1 = STM*DaughterMom1;
	DaughterMom2 = STM*DaughterMom2;
	SavePolarityMomentum(DaughterMom1);

	DaughterMom[0] = DaughterMom1;
	DaughterMom[1] = DaughterMom2;
	G4ThreeVector dummy(0,0,0);
	G4DynamicParticle* parent_dy = new G4DynamicParticle(parent,dummy,0);
	G4DecayProducts* products = new G4DecayProducts(*parent_dy); 
	delete parent_dy;
	G4DynamicParticle* Daughter1 = new G4DynamicParticle(Daughters[0],DaughterMom[0]); 
	G4DynamicParticle* Daughter2 = new G4DynamicParticle(Daughters[1],DaughterMom[1]); 
	products->PushProducts(Daughter1);
	products->PushProducts(Daughter2);
	return products;
}
void TPCPolarizedDecayChannel::LoadPolarityMomentum(){
	Polarity = gTrackBuffer.GetPolarity(order);
	MomVector = gTrackBuffer.GetMomentum(order);
	if(Polarity.mag()==0)G4cout<<"Warning! "<<*parent_name<<" Polarity not set!"<<G4endl;
	if(MomVector.mag()==0)G4cout<<"Warning! "<<*parent_name<<" Momentum cannot be tracked!"<<G4endl;
}
void TPCPolarizedDecayChannel::SavePolarityMomentum(G4ThreeVector MomD){
	auto LVParent = G4LorentzVector(MomVector,hypot(ParentMass/CLHEP::GeV,MomVector.mag()));
	auto ParentFrame = LVParent.boostVector();
	auto LVDaughter = G4LorentzVector(MomD,hypot(DaughterMass[0],MomD.mag()))/CLHEP::GeV;
	LVDaughter.boost(ParentFrame);
	auto TVDaughter = LVDaughter.vect();
	gTrackBuffer.SetMomentum(TVDaughter,order+1);		
	gTrackBuffer.SetLV(LVDaughter,order+1);		
	auto SpinDaughter = TVDaughter.cross(MomVector);
	SpinDaughter = SpinDaughter*(1./SpinDaughter.mag());
	gTrackBuffer.SetPolarity(SpinDaughter,order+1);
}





G4double TPCPolarizedDecayChannel::Pmx(G4double e,G4double p1, G4double p2){
   G4double ppp = (e+p1+p2)*(e+p1-p2)*(e-p1+p2)*(e-p1-p2)/(4.0*e*e);
	 if (ppp>0) return std::sqrt(ppp);
	 else       return -1.;
}
/*
XiToLdPiChannel
::XiToLdPiChannel(
				G4String ParentName,
				G4double BR,
				G4double Alpha,
				G4String Daughter1Name,
				G4String Daughter2Name){
	Initialize(ParentName,BR,Alpha,
			Daughter1Name,Daughter2Name);
}
void XiToLdPiChannel::LoadPolarityMomentum(){
	Polarity = gTrackBuffer.GetPolarity(order);
	MomVector = gTrackBuffer.GetMomentum(order);
	if(Polarity.mag()==0)G4cout<<"Warning! "<<*parent_name<<" Polarity not set!"<<G4endl;
	if(MomVector.mag()==0)G4cout<<"Warning! "<<*parent_name<<" Momentum cannot be tracked!"<<G4endl;
}
void XiToLdPiChannel::SavePolarityMomentum(G4ThreeVector MomM,G4ThreeVector MomD){
	auto LVParent = G4LorentzVector(MomM,hypot(ParentMass,MomM.mag()));
	auto ParentFrame = LVParent.boostVector();
	auto LVDaughter = G4LorentzVector(MomD,hypot(DaughterMass[0],MomD.mag()));
	LVDaughter.boost(ParentFrame);
	auto TVDaughter = LVDaughter.vect();
	gTrackBuffer.SetMomentum(TVDaughter,order+1);		
	gTrackBuffer.SetLV(LVDaughter,order+1);		
	auto SpinDaughter = TVDaughter.cross(MomM);
	SpinDaughter = SpinDaughter*(1./SpinDaughter.mag());
	gTrackBuffer.SetPolarity(SpinDaughter,order+1);
}

LdToPPiChannel
::LdToPPiChannel(
				G4String ParentName,
				G4double BR,
				G4double Alpha,
				G4String Daughter1Name,
				G4String Daughter2Name){
	Initialize(ParentName,BR,Alpha,
			Daughter1Name,Daughter2Name);
}

void LdToPPiChannel::LoadPolarityMomentum(){
	Polarity = gTrackBuffer.GetPolarity(1);
	MomVector = gTrackBuffer.GetMomentum(1);
	if(Polarity.mag()==0)G4cout<<"Warning! "<<*parent_name<<" Polarity not set!"<<G4endl;
	if(MomVector.mag()==0)G4cout<<"Warning! "<<*parent_name<<" Momentum cannot be tracked!"<<G4endl;
}

void LdToPPiChannel::SavePolarityMomentum(G4ThreeVector MomM,G4ThreeVector MomD){
	auto LVParent = G4LorentzVector(MomM,hypot(ParentMass,MomM.mag()));
	auto ParentFrame = LVParent.boostVector();
	auto LVDaughter = G4LorentzVector(MomD,hypot(DaughterMass[0],MomD.mag()));
	LVDaughter.boost(ParentFrame);
	auto TVDaughter = LVDaughter.vect();
	gTrackBuffer.SetMomentum(TVDaughter,2);		
	gTrackBuffer.SetLV(LVDaughter,2);		
	auto SpinDaughter = TVDaughter.cross(MomM);
	SpinDaughter = SpinDaughter*(1./SpinDaughter.mag());
	gTrackBuffer.SetPolarity(SpinDaughter,2);
}

*/
#endif

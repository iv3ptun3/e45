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
				G4double* DecayParameter,
				G4int Order,
				G4String Daughter1Name,
				G4String Daughter2Name
				){
	order = Order;
	Alpha = DecayParameter[0];
	double Angle = DecayParameter[1]*TMath::DegToRad();
	Beta = sqrt(1-Alpha*Alpha)*sin(Angle);
	Gamma = sqrt(1-Alpha*Alpha)*cos(Angle);
	Initialize(ParentName,BR,Alpha,
			Daughter1Name,Daughter2Name);
}

void
TPCPolarizedDecayChannel
::Initialize(
				G4String ParentName,
				G4double BR,
				G4double alpha,
				G4String Daughter1Name,
				G4String Daughter2Name
				){
	rbranch = BR;
	parent_name = new G4String(ParentName);
	daughters_name = new G4String*[2];
	daughters_name[0] = new G4String(Daughter1Name);
	daughters_name[1] = new G4String(Daughter2Name);
	numberOfDaughters = 2;
	PDFCTheta = TF1("PDF","1 + [0]*[1]*x",-1,1);	
	PDFCTheta.SetParameter(0,alpha);

	Daughters[0] = G4ParticleTable::GetParticleTable()->FindParticle(Daughter1Name);
	Daughters[1] = G4ParticleTable::GetParticleTable()->FindParticle(Daughter2Name);
	
	DaughterMass[0] = Daughters[0]->GetPDGMass();
	DaughterMass[1] = Daughters[1]->GetPDGMass();
	parent = G4ParticleTable::GetParticleTable()->FindParticle(ParentName);

	ParentMass = parent->GetPDGMass();
	

}

G4DecayProducts*
TPCPolarizedDecayChannel::DecayIt(G4double ParentMas){
	
	LoadPolarityMomentum();	
	PDFCTheta.SetParameter(1,Polarization);
	G4double P_d = Pmx(ParentMass,DaughterMass[0],DaughterMass[1]);
	double PI = acos(-1);
	double Phi = 2*PI*G4UniformRand();
	Theta = acos(PDFCTheta.GetRandom());
	double px = P_d*sin(Theta)*cos(Phi);
	double py = P_d*sin(Theta)*sin(Phi);
	double pz = P_d*cos(Theta);
	auto Zaxis = Polarity*(1./Polarity.mag()); 
	auto Xaxis = MomVector.cross(Polarity);
	Xaxis = Xaxis*(1./Xaxis.mag());
	auto Yaxis = Zaxis.cross(Xaxis);
	auto DaughterMom1 = px*Xaxis+py*Yaxis+pz*Zaxis;
	auto DaughterMom2 =-DaughterMom1;
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
	Polarization = gTrackBuffer.GetPolarization(order); 
	MomVector = gTrackBuffer.GetMomentum(order);
//	G4cout<<*parent_name<<" Polarization P = "<<Polarization<<G4endl;
	if(abs(Polarization) >1.01 or isnan(Polarization)){
		G4cout<<"Warning! "<<*parent_name<<" Polarization unphysical! P = "<<Polarization<<G4endl;
		Polarization = 0;
	}
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
	gTrackBuffer.SetVertexMomentum(TVDaughter,order+1);		
	gTrackBuffer.SetLV(LVDaughter,order+1);		
	gTrackBuffer.SetVertexLV(LVDaughter,order+1);		
	MomD = MomD*(1./MomD.mag());
	auto PolM = Polarity*(1./Polarity.mag());
	auto PxM = PolM.cross(MomD); 
	auto MxPxM =MomD.cross(PxM);
	auto SpinDaughter =( (Alpha+Polarization*cos(Theta))*MomD + Polarization*Beta*PxM + Polarization*Gamma * MxPxM) * 1./(1 + Alpha * Polarization * cos(Theta));
	auto PolDaughter = SpinDaughter.mag();
	SpinDaughter = SpinDaughter*(1./SpinDaughter.mag());
	gTrackBuffer.SetPolarity(SpinDaughter,order+1);
	gTrackBuffer.SetPolarization(PolDaughter,order+1);
}
void TPCPolarizedDecayChannel::SetPolarization(G4double pol){
	gTrackBuffer.SetPolarization(pol,order);
}




G4double TPCPolarizedDecayChannel::Pmx(G4double e,G4double p1, G4double p2){
   G4double ppp = (e+p1+p2)*(e+p1-p2)*(e-p1+p2)*(e-p1-p2)/(4.0*e*e);
	 if (ppp>0) return std::sqrt(ppp);
	 else       return -1.;
}
#endif

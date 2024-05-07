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
	auto LVDBK = LVDaughter;
	LVDaughter.boost(ParentFrame);
	auto TVDaughter = LVDaughter.vect();
	auto TVVertParent = gTrackBuffer.GetVertexMomentum(order);		
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
	if(abs(PolDaughter) >1.01 or isnan(PolDaughter)){
		auto Px=gTrackBuffer.GetMomentumOfTrack_x();
		auto Py=gTrackBuffer.GetMomentumOfTrack_y();
		auto Pz=gTrackBuffer.GetMomentumOfTrack_z();
		auto Vx=gTrackBuffer.GetVertexOfTrack_x();
		auto Vy=gTrackBuffer.GetVertexOfTrack_y();
		auto Vz=gTrackBuffer.GetVertexOfTrack_z();
		auto Kmx =	Px[1];
		auto Kmy =	Py[1];
		auto Kmz =	Pz[1];
		auto Kpx =	Px[2];
		auto Kpy =	Py[2];
		auto Kpz =	Pz[2];
		auto Xpx =	Px[3];
		auto Xpy =	Py[3];
		auto Xpz =	Pz[3];
		auto Vpx =	Vx[1];
		auto Vpy =	Vy[1];
		auto Vpz =	Vz[1];
		
		double VpXx,VpXy,VpXz;
		double VpLx,VpLy,VpLz;
		auto pid = gTrackBuffer.GetPIDOfTrack();
		auto parid = gTrackBuffer.GetParentIDOfTrack();
		G4ThreeVector VtxXiDecay;	
		G4ThreeVector VtxLdDecay;	
		for(int in=0;in<1000;++in){
			if(parid[in] == 3312){
				VpXx =	Vx[in];
				VpXy =	Vy[in];
				VpXz =	Vz[in];
				VtxXiDecay=G4ThreeVector(VpXx,VpXy,VpXz);
			}
			if(parid[in] == 3112){
				VpXx =	Vx[in];
				VpXy =	Vy[in];
				VpXz =	Vz[in];
				VtxLdDecay=G4ThreeVector(VpXx,VpXy,VpXz);
			}
		}
		G4ThreeVector TVKm(Kmx,Kmy,Kmz);
		G4ThreeVector TVKp(Kpx,Kpy,Kpz);
		G4ThreeVector Vtx(Vpx,Vpy,Vpz);

		G4cout<<"Warning! Unphysical Polarization "<<*parent_name<<G4endl;
		G4cout<<Form("PKm=(%g,%g,%g)",TVKm.x(),TVKm.y(),TVKm.z())<<G4endl;
		G4cout<<Form("PKp=(%g,%g,%g)",TVKp.x(),TVKp.y(),TVKp.z())<<G4endl;
		G4cout<<Form("ProdVert=(%g,%g,%g)",Vtx.x(),Vtx.y(),Vtx.z())<<G4endl;
		G4cout<<Form("XiDVert=(%g,%g,%g)",VtxXiDecay.x(),VtxXiDecay.y(),VtxXiDecay.z())<<G4endl;
		G4cout<<Form("LdDVert=(%g,%g,%g)",VtxLdDecay.x(),VtxLdDecay.y(),VtxLdDecay.z())<<G4endl;
		
		
		
		G4cout<<Form("Parent Mom = (%g,%g,%g) ",LVParent.x(),LVParent.y(),LVParent.z())<<G4endl;
		G4cout<<Form("Parent VertMom = (%g,%g,%g) ",TVVertParent.x(),TVVertParent.y(),TVVertParent.z())<<G4endl;
		G4cout<<Form("Daughter Mom = (%g,%g,%g) ",LVDaughter.x(),LVDaughter.y(),LVDaughter.z())<<G4endl;
		G4cout<<Form("Daughter MomCM = (%g,%g,%g) ",LVDBK.x(),LVDBK.y(),LVDBK.z())<<G4endl;
	}
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

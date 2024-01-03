#ifndef TPCDecayChannel_h
#define TPCDecayChannel
#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"
#include "G4Decay.hh"
#include "TPCTrackBuffer.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4RotationMatrix.hh"
#include "Randomize.hh"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVector3.h"
using namespace std;
class TPCPolarizedDecayChannel : public G4VDecayChannel
{
	public:
		TPCPolarizedDecayChannel(){};
		TPCPolarizedDecayChannel(
				G4String Parent,
				G4double BR,
				G4double Alpha,
				G4int Order,
				G4String Daughter1,
				G4String Daughter2);

	public:
		 ~TPCPolarizedDecayChannel(){};
	public:
		virtual	G4DecayProducts* DecayIt(G4double);
	public:
		G4double Pmx(G4double e,G4double p1, G4double p2);
		G4ThreeVector SetPolarity();
		G4ThreeVector GetPolarity();
		void SetAlpha(double Alpha){alpha = Alpha;}
	protected:
		G4double ParentMass;
		G4ThreeVector Polarity;//Z axis on mother frame
		G4ThreeVector MomVector;//X axis on mother frame
		TF1 PDFCTheta;
		int order;
		double alpha;
		double GetTheta();
		G4ParticleDefinition* Daughters[2];
		G4ParticleDefinition* parent;
		G4ThreeVector DaughterMom[2];
		G4RotationMatrix SpinToMomentum();
		void Initialize(
				G4String Parent,
				G4double BR,
				G4double Alpha,
				G4String Daughter1,
				G4String Daughter2);
		void LoadPolarityMomentum();
		void SavePolarityMomentum(G4ThreeVector MomDauthger);
		G4double DaughterMass[2];
		//		TLorentzVector ToCM(TLorentzVector LV);
//		TLorentzVector ToLab(TLorentzVector LV);
};
/*
class XiToLdPiChannel : public TPCPolarizedDecayChannel{
	public:
		XiToLdPiChannel(){};
	
		XiToLdPiChannel(
				G4String Parent,
				G4double BR,
				G4double Alpha,
				G4String Daughter1,
				G4String Daughter2);

		~XiToLdPiChannel(){};
	protected:	
		virtual void LoadPolarityMomentum();
		virtual void SavePolarityMomentum(G4ThreeVector MomM, G4ThreeVector MomD);
};
class LdToPPiChannel : public TPCPolarizedDecayChannel{
	public:
		LdToPPiChannel(){};
	
		LdToPPiChannel(
				G4String Parent,
				G4double BR,
				G4double Alpha,
				G4String Daughter1,
				G4String Daughter2);

		~LdToPPiChannel(){};
	protected:	
		virtual void LoadPolarityMomentum();
		virtual void SavePolarityMomentum(G4ThreeVector MomM, G4ThreeVector MomD);
};

*/
#endif

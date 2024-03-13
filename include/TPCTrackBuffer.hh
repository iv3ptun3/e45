#ifndef TPC_TRACK_BUFFER_HH
#define TPC_TRACK_BUFFER_HH 
#include <TMath.h>
#include <stdio.h>
#include <iostream>
#include "globals.hh"
#include <G4ThreeVector.hh>
#include "G4LorentzVector.hh"
namespace{
	const int MaxBufferSize = 1000;
	const auto qnan = TMath::QuietNaN();
}
class TPCTrackBuffer{
	private:
		int NumberOfTracks;
		int PIDOfTrack[MaxBufferSize];
		int ParentIDOfTrack[MaxBufferSize];
		double VertexOfTrack_x[MaxBufferSize];
		double VertexOfTrack_y[MaxBufferSize];
		double VertexOfTrack_z[MaxBufferSize];
		double MomentumOfTrack[MaxBufferSize];
		double MomentumOfTrack_x[MaxBufferSize];
		double MomentumOfTrack_y[MaxBufferSize];
		double MomentumOfTrack_z[MaxBufferSize];
		G4ThreeVector Polarity[10];
		G4double Polarization[10];
		G4ThreeVector Momentum[10];
		G4ThreeVector VertexMomentum[10];
		G4LorentzVector LVCM;	
		G4LorentzVector LV[10];
		G4LorentzVector VertexLV[10];
		int LambdaID;
	public:
		TPCTrackBuffer(){
			ClearTrackBuffer();
		}
		~TPCTrackBuffer(){
		}
  	static TPCTrackBuffer& GetInstance( void );
		G4String ClassName();
		void ClearTrackBuffer(){
			//		std::cout<<"ClearBuffer"<<std::endl;
			NumberOfTracks = -1;
			LambdaID = -1;
			LVCM = G4LorentzVector(0,0,0,0);
			for(int i=0;i<10;++i){
//				Polarization[i] = 0;
				Polarity[i] = G4ThreeVector(0,0,0); 
				Momentum[i] = G4ThreeVector(0,0,0);
				VertexMomentum[i] = G4ThreeVector(0,0,0);
				LV[i] = G4LorentzVector(0,0,0,0);
				VertexLV[i] = G4LorentzVector(0,0,0,0);
			}
			for(int i=0;i<MaxBufferSize;++i){
				PIDOfTrack[i]=qnan;
				ParentIDOfTrack[i]=qnan;
				VertexOfTrack_x[i]=qnan;
				VertexOfTrack_y[i]=qnan;
				VertexOfTrack_z[i]=qnan;
				MomentumOfTrack[i]=qnan;
				MomentumOfTrack_x[i]=qnan;
				MomentumOfTrack_y[i]=qnan;
				MomentumOfTrack_z[i]=qnan;
			}	
		}
		void SetTrack(int ID,int ParentID,int PID,G4ThreeVector vert, G4ThreeVector mom){
			if(ID > NumberOfTracks)NumberOfTracks = ID; 
			PIDOfTrack[ID] = PID;
			ParentIDOfTrack[ID] = ParentID;
			VertexOfTrack_x[ID]=vert.x();
			VertexOfTrack_y[ID]=vert.y();
			VertexOfTrack_z[ID]=vert.z();
			MomentumOfTrack[ID]=mom.mag();
			MomentumOfTrack_x[ID]=mom.x();
			MomentumOfTrack_y[ID]=mom.y();
			MomentumOfTrack_z[ID]=mom.z();
		}
		int GetNumberOfTracks(){
			return NumberOfTracks;
		}
		int* GetPIDOfTrack(){
			return PIDOfTrack;
		}
		int* GetParentIDOfTrack(){
			return ParentIDOfTrack;
		}
		double* GetVertexOfTrack_x(){
			return VertexOfTrack_x;
		}
		double* GetVertexOfTrack_y(){
			return VertexOfTrack_y;
		}
		double* GetVertexOfTrack_z(){
			return VertexOfTrack_z;
		}
		double* GetMomentumOfTrack(){
			return MomentumOfTrack;
		}
		double* GetMomentumOfTrack_x(){
			return MomentumOfTrack_x;
		}
		double* GetMomentumOfTrack_y(){
			return MomentumOfTrack_y;
		}
		double* GetMomentumOfTrack_z(){
			return MomentumOfTrack_z;
		}
		void SetPolarity(G4ThreeVector pol,int i){
			Polarity[i] = 	pol;
		}
		void SetPolarization(G4double pol,int i){
			Polarization[i] = pol;
		}
		void SetMomentum(G4ThreeVector mom,int i){
			Momentum[i] = 	mom;
		}
		void SetVertexMomentum(G4ThreeVector mom,int i){
			VertexMomentum[i] = 	mom;
		}
		void SetLV(G4LorentzVector lv,int i){
			LV[i] = 	lv;
		}
		void SetVertexLV(G4LorentzVector lv,int i){
			VertexLV[i] = 	lv;
		}
		void SetCMLV(G4LorentzVector lv){
			LVCM = 	lv;
		}
		G4ThreeVector GetPolarity(int i){
			return Polarity[i];
		}
		G4double GetPolarization(int i){
			return Polarization[i];
		}
		G4ThreeVector GetMomentum(int i){
			return Momentum[i];
		}
		G4ThreeVector GetVertexMomentum(int i){
			return VertexMomentum[i];
		}
		G4LorentzVector GetLV(int i){
			return LV[i];
		}
		G4LorentzVector GetVertexLV(int i){
			return VertexLV[i];
		}
		G4LorentzVector GetCMLV(){
			return LVCM;
		}
		void SetLambdaID(int id ){
			LambdaID = id;
		}
		int GetLambdaID(){
			return LambdaID;
		}
};
inline G4String
TPCTrackBuffer::ClassName( void )
{
  static G4String s_name("TPCTrackBuffer");
  return s_name;
}
 inline TPCTrackBuffer& 
 TPCTrackBuffer::GetInstance( void ){
  static TPCTrackBuffer s_instance;
  return s_instance;
 }
#endif

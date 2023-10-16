#ifndef TPC_TRACK_BUFFER_HH
#define TPC_TRACK_BUFFER_HH 
#include <stdio.h>
#include <iostream>
#include "globals.hh"
#include <G4ThreeVector.hh>
namespace{
	const int MaxBufferSize = 1000;
}
class TPCTrackBuffer{
	private:
		int NumberOfTracks;
		int PIDOfTrack[MaxBufferSize];
		double VertexOfTrack_x[MaxBufferSize];
		double VertexOfTrack_y[MaxBufferSize];
		double VertexOfTrack_z[MaxBufferSize];
		double MomentumOfTrack_x[MaxBufferSize];
		double MomentumOfTrack_y[MaxBufferSize];
		double MomentumOfTrack_z[MaxBufferSize];
	public:
		TPCTrackBuffer(){
		}
		~TPCTrackBuffer(){
		}
  	static TPCTrackBuffer& GetInstance( void );
		G4String ClassName();
		void ClearTrackBuffer(){
			//		std::cout<<"ClearBuffer"<<std::endl;
			NumberOfTracks = -1;
			for(int i=0;i<MaxBufferSize;++i){
				PIDOfTrack[i]=0;
				VertexOfTrack_x[i]=0;
				VertexOfTrack_y[i]=0;
				VertexOfTrack_z[i]=0;
				MomentumOfTrack_x[i]=0;
				MomentumOfTrack_y[i]=0;
				MomentumOfTrack_z[i]=0;
			}	
		}
		void SetTrack(int ID,int PID,G4ThreeVector vert, G4ThreeVector mom){
			if(ID > NumberOfTracks)NumberOfTracks = ID; 
			PIDOfTrack[ID] = PID;
			VertexOfTrack_x[ID]=vert.x();
			VertexOfTrack_y[ID]=vert.y();
			VertexOfTrack_z[ID]=vert.z();
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
		double* GetVertexOfTrack_x(){
			return VertexOfTrack_x;
		}
		double* GetVertexOfTrack_y(){
			return VertexOfTrack_y;
		}
		double* GetVertexOfTrack_z(){
			return VertexOfTrack_z;
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

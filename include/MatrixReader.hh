#ifndef MatrixReader
#define MatrixReader
#include <fstream>
namespace{
	using namespace std;
	int ReadParam(ifstream& file, int* buf, TString& buf_l){
		if(!file.is_open()){
			cout<<"file not open"<<endl;
			return false;
		}
		TString line;
		if(file.good()&&line.ReadLine(file)){
			buf_l = line;
			if(line.IsNull()){
				return 0;
			}   
			line.ReplaceAll(" ","");
			line.ReplaceAll(",","");
			line.ReplaceAll("\"","");
			if(line.First("TOFsegment")+1){
				line.ReplaceAll("#TOFsegment","");
				buf[0] = line.Atoi(); 
				return 1;
			}else if(line.First("SCH")+1){
				line.ReplaceAll("SCH","");
				auto line_ch =(TString)line(0,2);
				int sch_ch = line_ch.Atoi();
				buf[0] = sch_ch;
				int nl = line.Length();
				buf[1] = nl;
				
				for(int il = 2;il < nl; ++il){
					buf[il  ] = line[il]-48;
				}
				return 2;
			}
			return -1;
		}
		else{
			return 0;
		}
	}
	TString mat2d;
	TString mat3d;
	int Mat2D[24][64];//ToF,SCH
	static int Mat3D[28][64][8];//ToF+BVH,SCH,BH2
	void ReadMatrix(){
		ifstream matf2,matf3;
		int buf[200];
		TString hh;
		matf2.open(mat2d);
//		matf3.open(mat3d);
		int flag = ReadParam(matf2,buf,hh);	
		int ToF_seg = -1,SCH_seg = -1,SCH_nl = -1;
		while(flag){
			if(flag == 1){
				ToF_seg = buf[0];
			}
			if(flag == 2){
				SCH_seg = buf[0];
				SCH_nl = buf[1];
				Mat2D[ToF_seg][SCH_seg] = buf[2];
			}
			flag = ReadParam(matf2,buf,hh);	
		}
/*
		flag = ReadParam(matf3,buf,hh);	
		while(flag){
			if(flag == 1){
				ToF_seg = buf[0];
			}
			if(flag == 2){
				SCH_seg = buf[0];
				SCH_nl = buf[1];
				for(int il = 2; il < SCH_nl;++il){
					Mat3D[ToF_seg][SCH_seg][il-2] = buf[il];
				}
			}
			flag = ReadParam(matf3,buf,hh);	
		}
		*/
	};

}

#endif

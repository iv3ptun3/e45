#include "TPCParamMan.hh"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <CLHEP/Units/SystemOfUnits.h>

#include "DeleteUtility.hh"

#include "FuncName.hh"
namespace
{
const Int_t RowMask = 0xff;
const Int_t LayerMask = 0xff;
const Int_t RowShift = 0;
const Int_t LayerShift = 8;

inline Int_t
MakeKey(Int_t layer, Int_t row)
{
  return (((layer & LayerMask) << LayerShift) |
          ((row & RowMask) << RowShift));
}
}
TPCParamMan::TPCParamMan()
  : m_is_ready(false),
    m_file_name()
{
}
TPCParamMan::~TPCParamMan()
{
  ClearResCont();
}
void
TPCParamMan::ClearResCont()
{
  del::ClearMap(m_ResContainer);
}

Bool_t
TPCParamMan::Initialize()
{
  if(m_is_ready){
    G4cerr << FUNC_NAME << " already initialied" << G4endl;
    return false;
  }

  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    G4cerr << FUNC_NAME << " file open fail : "
		<< m_file_name <<" /" <<G4endl;
    return false;
  }

  ClearResCont();

  Int_t line_number = 0;
  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    ++line_number;
    if(line.IsNull() || line[0]=='#') continue;
    std::istringstream input_line(line.Data());
    Int_t layer, row, aty;
    Double_t p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11;
    if(input_line >> layer >> row >> aty >> p0 >> p1){
      Int_t key = MakeKey(layer, row);
      switch(aty){
      case kAdc: {
									 break;
								 }
      case kTdc: {
									 break;
								 }
      case kY: {
									 break;
								 }
      case kCobo: {
									 break;
								 }
      case kRes: {
	if(input_line >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8 >> p9 >> p10 >> p11){
	  TPCResParam *pre_param = m_ResContainer[key];
	  std::vector<Double_t> params{ p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11 }; 
	  TPCResParam *param = new TPCResParam(params);
	  m_ResContainer[key] = param;

	  Int_t HS = layer; Int_t InOut = row;
    if(HS==0&&InOut==0) m_Res_HSOFF_Inner = params;
    else if(HS==0&&InOut==1) m_Res_HSOFF_Outer = params;
    else if(HS==0&&InOut==2) m_ResCl_HSOFF_Inner = params;
    else if(HS==0&&InOut==3) m_ResCl_HSOFF_Outer = params;
    else if(HS==1&&InOut==0) m_Res_HSON_Inner = params;
    else if(HS==1&&InOut==1) m_Res_HSON_Outer = params;
    else if(HS==1&&InOut==2) m_ResCl_HSON_Inner = params;
    else if(HS==1&&InOut==3){
      m_ResCl_HSON_Outer = params;
    } 
    else G4cerr << FUNC_NAME << ": Invalid TPC resolution parameter"
			   << G4endl << " p0, 1: HS On 0 : HS Off"
			   << G4endl << " p1, 1: Outer layers 0 : inner layers "
			   << std::endl;

	  if(pre_param){
	    G4cerr << FUNC_NAME << ": duplicated key "
			<< " following record is deleted." << std::endl
			<< " layer = " << layer << ","
			<< " row = " << row
			<< " aty = " << aty
			<< G4endl;
	    delete pre_param;
	  }
	  break;
	}
      }
      default:
        G4cerr << FUNC_NAME << ": Invalid Input" << std::endl
                    << " ===> L" << line_number << " " << line << std::endl;
        break;
      } // switch
    } else {
      G4cerr << FUNC_NAME << ": Invalid Input" << std::endl
		  << " ===> L" << line_number << " " << line << std::endl;
    }
  } // while

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________


Bool_t
TPCParamMan::Initialize(const G4String& file_name)
{
  m_file_name = file_name;
  return Initialize();
}


TPCResParam*
TPCParamMan::GetResmap(Int_t B, Int_t InnerOrOuter) const
{
  Int_t key = MakeKey(B, InnerOrOuter);
  ResIterator itr = m_ResContainer.find(key);
  if(itr != m_ResContainer.end())
    return itr->second;
  else
    return nullptr;
}

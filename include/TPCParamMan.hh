#ifndef TPCParamMan_hh
#define TPCParamMan_hh
#include <map>
#include <TMath.h>
#include <globals.hh>
#include <TString.h>
class TPCResParam
{
public:
  TPCResParam(const std::vector<Double_t> params)
    : m_params(params)
    {}
  ~TPCResParam()
    {}

private:
  TPCResParam();
  TPCResParam(const TPCResParam&);
  TPCResParam& operator =(const TPCResParam&);

private:
  std::vector<Double_t> m_params;

public:
  const std::vector<Double_t>& Params() const { return m_params; }
};

//_____________________________________________________________________________
class TPCParamMan
{
public:
  static const G4String& ClassName();
  static TPCParamMan&   GetInstance();
  ~TPCParamMan();

private:
  TPCParamMan();
  TPCParamMan(const TPCParamMan&);
  TPCParamMan& operator =(const TPCParamMan&);

private:
  //  enum eAorT { kAdc, kTdc, kY, kATY };
  enum eAorT { kAdc, kTdc, kY, kCobo, kRes };
  typedef std::map<Int_t, TPCResParam*> ResContainer;
  typedef ResContainer::const_iterator ResIterator;
  Bool_t        m_is_ready;
  G4String       m_file_name;
  ResContainer  m_ResContainer;
  std::vector<Double_t> m_Res_HSON_Inner;
  std::vector<Double_t> m_Res_HSON_Outer;
  std::vector<Double_t> m_Res_HSOFF_Inner;
  std::vector<Double_t> m_Res_HSOFF_Outer;

public:
  Bool_t Initialize();
  Bool_t Initialize(const G4String& file_name);
  Bool_t IsReady() const { return m_is_ready; }
  void   SetFileName(const G4String& file_name) { m_file_name = file_name; }

private:
  void          ClearResCont();
  TPCResParam*  GetResmap(Int_t B, Int_t InnerOrOuter) const;

public:
  static const std::vector<Double_t>& TPCResolutionParams(Bool_t HSOn, Bool_t Inner);

};

//_____________________________________________________________________________
inline const G4String&
TPCParamMan::ClassName()
{
  static G4String g_name("TPCParamMan");
  return g_name;
}

//_____________________________________________________________________________
inline TPCParamMan&
TPCParamMan::GetInstance()
{
  static TPCParamMan s_instance;
  return s_instance;
}

//______________________________________________________________________________
inline const std::vector<Double_t>&
TPCParamMan::TPCResolutionParams(Bool_t HSOn, Bool_t Inner)
{
  if(!HSOn&&Inner) return GetInstance().m_Res_HSOFF_Inner;
  if(!HSOn&&!Inner) return GetInstance().m_Res_HSOFF_Outer;
  if(HSOn&&Inner) return GetInstance().m_Res_HSON_Inner;
  //if(HSOn&&!Inner) return GetInstance().m_Res_HSON_Outer;
  return GetInstance().m_Res_HSON_Outer;
}

#endif

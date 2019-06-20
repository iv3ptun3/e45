// -*- C++ -*-

#ifndef DET_SIZE_MAN_HH
#define DET_SIZE_MAN_HH

#include <string>
#include <map>
#include <vector>

#include <Rtypes.h>
#include <TString.h>

//_____________________________________________________________________________
class DetSizeMan
{
public:
  static TString     ClassName( void );
  static DetSizeMan& GetInstance( void );
  ~DetSizeMan( void );

private:
  DetSizeMan( void );
  DetSizeMan( const DetSizeMan&  );
  DetSizeMan& operator =( const DetSizeMan& );

private:
  typedef std::vector<Double_t>         ParamArray;
  typedef std::map<TString, ParamArray> ParamMap;
  typedef ParamMap::const_iterator      PIterator;
  Bool_t   m_is_ready;
  TString  m_file_name;
  ParamMap m_param_map;

public:
  Bool_t   Initialize( void );
  Bool_t   Initialize( const TString& filename );
  Bool_t   IsReady( void ) const { return m_is_ready; }
  Int_t    GetSize( const TString& key ) const;
  Double_t Get( const TString& key, Int_t i=0 ) const;
  void     Print( Option_t* option=nullptr ) const;
};

//_____________________________________________________________________________
inline TString
DetSizeMan::ClassName( void )
{
  static TString s_name("DetSizeMan");
  return s_name;
}

//_____________________________________________________________________________
inline DetSizeMan&
DetSizeMan::GetInstance( void )
{
  static DetSizeMan s_instance;
  return s_instance;
}

#endif

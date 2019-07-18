// -*- C++ -*-

#ifndef BEAM_MAN_HH
#define BEAM_MAN_HH

#include <string>
#include <map>
#include <vector>

#include <globals.hh>
#include <G4ThreeVector.hh>

#include <Rtypes.h>

class TFile;

//_____________________________________________________________________________
struct BeamInfo
{
  // beam information at FF position (VO+1200).
  Double_t x; // [mm]
  Double_t y; // [mm]
  Double_t u; // [mrad]
  Double_t v; // [mrad]
  Double_t dp; // [%]
  G4ThreeVector p;
  void Print( void ) const;
};

//_____________________________________________________________________________
class BeamMan
{
public:
  static G4String ClassName( void );
  static BeamMan& GetInstance( void );
  ~BeamMan( void );

private:
  BeamMan( void );
  BeamMan( const BeamMan&  );
  BeamMan& operator =( const BeamMan& );

private:
  typedef std::vector<BeamInfo> ParamArray;
  G4bool     m_is_ready;
  G4String   m_file_name;
  TFile*     m_file;
  ParamArray m_param_array;
  G4int      m_n_param;

public:
  G4bool          Initialize( void );
  G4bool          Initialize( const G4String& filename );
  G4bool          IsReady( void ) const { return m_is_ready; }
  const BeamInfo& Get( void ) const;
  void            Print( void ) const;
};

//_____________________________________________________________________________
inline G4String
BeamMan::ClassName( void )
{
  static G4String s_name("BeamMan");
  return s_name;
}

//_____________________________________________________________________________
inline BeamMan&
BeamMan::GetInstance( void )
{
  static BeamMan s_instance;
  return s_instance;
}

#endif

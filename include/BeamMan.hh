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
  Double_t x; // [mm]
  Double_t y; // [mm]
  Double_t u; // [mrad]
  Double_t v; // [mrad]
  Double_t dp; // [%]
  G4ThreeVector p; // [GeV/c]
  Double_t z; // [mm]
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
  G4bool        m_is_ready;
  G4String      m_file_name;
  TFile*        m_file;
  ParamArray    m_param_array;
  G4int         m_n_param;
  G4bool        m_is_vi; // true:VI or false:VO
  G4double      m_primary_z; // from VI or VO
  G4ThreeVector m_vi_pos;

public:
  const BeamInfo&      Get( void ) const;
  G4double             GetPrimaryZ( void ) const { return m_primary_z; }
  const G4ThreeVector& GetVIPosition( void ) const { return m_vi_pos; }
  G4bool               Initialize( void );
  G4bool               Initialize( const G4String& filename );
  G4bool               IsReady( void ) const { return m_is_ready; }
  void                 Print( void ) const;
  void                 SetPrimaryZ( G4double z ){ m_primary_z = z; }
  void                 SetVIPosition( G4ThreeVector pos ){ m_vi_pos = pos; }
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

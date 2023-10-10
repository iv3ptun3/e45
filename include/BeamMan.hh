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
  G4double      x; // [mm]
  G4double      y; // [mm]
  G4double      u; // [mrad] -> dxdz in k18-analyzer
  G4double      v; // [mrad]
  G4double      dp; // [%]
  G4ThreeVector p; // [GeV/c]
  G4double      z; // [mm]
  G4double      m2; // [GeV / c2]
  G4double      q; 
  G4int					trigpat[32];
	G4double GetX( G4double offset=0. ) const;
  G4double GetY( G4double offset=0. ) const;
	G4int GetTrigPat(G4int flag) const;
  void     Print( void ) const;
	G4int evnum;
	G4int runnum;
	G4int ntBeam;
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
  G4bool        m_is_k18;
  G4bool        m_is_kurama;
  G4double      m_primary_z; // from VI or VO
  G4double      m_target_z;
  G4ThreeVector m_vi_pos;

public:
  const BeamInfo&      Get( void ) const;
  const BeamInfo&      Get( G4int iev ) const;
  G4double             GetPrimaryZ( void ) const { return m_primary_z; }
  const G4ThreeVector& GetVIPosition( void ) const { return m_vi_pos; }
  G4bool               Initialize( void );
  G4bool               Initialize( const G4String& filename );
  G4bool               IsReady( void ) const { return m_is_ready; }
  G4bool               IsK18( void ) const { return m_is_k18; }
  G4bool               IsKurama( void ) const { return m_is_kurama; }
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

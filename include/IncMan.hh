// -*- C++ -*-

#ifndef INC_MAN_HH
#define INC_MAN_HH

#include <string>
#include <map>
#include <vector>

#include <globals.hh>
#include <G4ThreeVector.hh>

#include <Rtypes.h>

#include "TPCAnaManager.hh"

class TFile;
class TTree;

//_____________________________________________________________________________
struct IncInfo
{
  // reaction
  Int_t ich;
  // beam
  Double_t bpx; // [GeV/c]
  Double_t bpy; // [GeV/c]
  Double_t bpz; // [GeV/c]
  // Kp
  Double_t kppx; // [GeV/c]
  Double_t kppy; // [GeV/c]
  Double_t kppz; // [GeV/c]
  // fs
  Int_t np;
  Int_t pid[MaxHits];
  Double_t px[MaxHits]; // [GeV/c]
  Double_t py[MaxHits]; // [GeV/c]
  Double_t pz[MaxHits]; // [GeV/c]
  void Print( void ) const;
};

//_____________________________________________________________________________
class IncMan
{
public:
  static G4String ClassName( void );
  static IncMan& GetInstance( void );
  ~IncMan( void );

private:
  IncMan( void );
  IncMan( const IncMan&  );
  IncMan& operator =( const IncMan& );

private:
  G4bool     m_is_ready;
  G4String   m_file_name;
  TFile*     m_file;
  TTree*     m_tree;
  IncInfo*   m_event;
  G4int      m_n_event;

public:
  IncInfo*   Get( void ) const;
  IncInfo*   Get( Int_t i ) const;
  G4bool     Initialize( void );
  G4bool     Initialize( const G4String& filename );
  G4bool     IsReady( void ) const { return m_is_ready; }
  void       Print( void ) const;
};

//_____________________________________________________________________________
inline G4String
IncMan::ClassName( void )
{
  static G4String s_name("IncMan");
  return s_name;
}

//_____________________________________________________________________________
inline IncMan&
IncMan::GetInstance( void )
{
  static IncMan s_instance;
  return s_instance;
}

#endif

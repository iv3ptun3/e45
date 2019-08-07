// -*- C++ -*-

#ifndef TPC_WC_SD_HH
#define TPC_WC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCWCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCWCSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName( void );
  TPCWCSD( const G4String& name );
  virtual ~TPCWCSD( void );

private:
  G4THitsCollection<TPCWCHit>* m_hits_collection;
  G4double                     m_refractive_index;

public:
  G4double GetRefractiveIndex( void ) const { return m_refractive_index; }
  void     SetRefractiveIndex( G4double index ){ m_refractive_index = index; }

public:
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

//_____________________________________________________________________________
inline G4String
TPCWCSD::ClassName( void )
{
  static G4String s_name("TPCWCSD");
  return s_name;
}

#endif

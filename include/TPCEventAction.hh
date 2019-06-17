// -*- C++ -*-

#ifndef TPC_EVENT_ACTION_HH
#define TPC_EVENT_ACTION_HH

#include <G4UserEventAction.hh>
#include <G4Types.hh>

class G4Event;

//_____________________________________________________________________________
class TPCEventAction : public G4UserEventAction
{
public:
  TPCEventAction( void );
  virtual ~TPCEventAction( void );
  virtual void BeginOfEventAction( const G4Event* anEvent );
  virtual void EndOfEventAction( const G4Event* anEvent );

private:
  G4int ac_use;
  G4int n_bar_use;
  G4int experiment_num;
  G4int with_kurama;
  //  G4int pidtr[50];
  //  G4int nparticle;
};

#endif

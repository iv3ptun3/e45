// -*- C++ -*-

#ifndef TPC_EVENT_ACTION_HH
#define TPC_EVENT_ACTION_HH

#include <G4String.hh>
#include <G4UserEventAction.hh>
#include <G4Types.hh>

class G4Event;

//_____________________________________________________________________________
class TPCEventAction : public G4UserEventAction
{
public:
  static G4String ClassName( void );
  TPCEventAction( void );
  virtual ~TPCEventAction( void );
  virtual void BeginOfEventAction( const G4Event* anEvent );
  virtual void EndOfEventAction( const G4Event* anEvent );
};

//_____________________________________________________________________________
inline G4String
TPCEventAction::ClassName( void )
{
  static G4String s_name("TPCEventAction");
  return s_name;
}

#endif

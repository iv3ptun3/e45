// -*- C++ -*-

#ifndef TPC_RUN_ACTION_HH
#define TPC_RUN_ACTION_HH

#include <G4UserRunAction.hh>
#include <G4String.hh>

class G4Run;
class TPCAnaManager;

//_____________________________________________________________________________
class TPCRunAction : public G4UserRunAction
{
public:
  static G4String ClassName( void );
  TPCRunAction( void );
  virtual ~TPCRunAction( void );
  virtual void BeginOfRunAction( const G4Run* aRun );
  virtual void EndOfRunAction( const G4Run* aRun );
};

//_____________________________________________________________________________
inline G4String
TPCRunAction::ClassName( void )
{
  static G4String s_name("TPCRunAction");
  return s_name;
}

#endif

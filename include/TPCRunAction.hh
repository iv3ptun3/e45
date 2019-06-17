// -*- C++ -*-

#ifndef TPC_RUN_ACTION_HH
#define TPC_RUN_ACTION_HH

#include <G4UserRunAction.hh>

class G4Run;
class TPCAnaManager;

class TPCRunAction : public G4UserRunAction
{
public:
  TPCRunAction( void );
  virtual ~TPCRunAction( void );
  virtual void BeginOfRunAction( const G4Run* aRun );
  virtual void EndOfRunAction( const G4Run* aRun );
};

#endif

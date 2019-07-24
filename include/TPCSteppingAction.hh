// -*- C++ -*-

#ifndef TPC_STEPPING_ACTION_HH
#define TPC_STEPPING_ACTION_HH

#include <G4UserSteppingAction.hh>

//_____________________________________________________________________________
class TPCSteppingAction : public G4UserSteppingAction
{
public:
  TPCSteppingAction( void );
  virtual ~TPCSteppingAction( void );
  virtual void UserSteppingAction( const G4Step* theStep );
};

#endif

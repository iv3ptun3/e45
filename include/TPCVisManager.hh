// -*- C++ -*-

#ifndef TPC_VIS_MANAGER_HH
#define TPC_VIS_MANAGER_HH

#include <G4VisManager.hh>

//_____________________________________________________________________________
class TPCVisManager : public G4VisManager
{
public:
  TPCVisManager( void );
  virtual ~TPCVisManager( void );

private:
  virtual void RegisterGraphicsSystems( void );
};

#endif

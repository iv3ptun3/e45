// -*- C++ -*-

#ifndef TPC_PHYSICS_LIST_H
#define TPC_PHYSICS_LIST_H

#include <G4VModularPhysicsList.hh>

//_____________________________________________________________________________
class TPCPhysicsList : public G4VModularPhysicsList
{
public:
  TPCPhysicsList( void );
  virtual ~TPCPhysicsList( void );

private:
  G4VPhysicsConstructor*              m_em_physics_list;
  std::vector<G4VPhysicsConstructor*> m_hadron_physics_list;

protected:
  void ConstructBaryons( void );
  void ConstructBosons( void );
  void ConstructEM( void );
  void ConstructGeneral( void );
  void ConstructHadron( void );
  void ConstructIons( void );
  void ConstructLeptons( void );
  void ConstructMesons( void );
  void ConstructParticle( void );
  void ConstructProcess( void );
  void ConstructShortLived( void );
  void ConstructStableHyperons( void );
  void SetCuts( void );
};

#endif

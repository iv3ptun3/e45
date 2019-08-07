// -*- C++ -*-

#ifndef VHIT_INFO_HH
#define VHIT_INFO_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>

class G4Step;

//_____________________________________________________________________________
class VHitInfo
{
public:
  static G4String ClassName( void );
  VHitInfo( const G4String& name, G4Step* step=nullptr );
  ~VHitInfo( void );

private:
  G4String      m_detector_name;
  G4String      m_particle_name;
  G4ThreeVector m_position;
  G4ThreeVector m_momentum;
  G4double      m_time;
  G4double      m_energy_deposit;
  G4int         m_track_id;
  G4int         m_particle_id;
  G4int         m_detector_id;
  G4double      m_step_length;
  G4double      m_mass;
  G4double      m_charge;
  G4int         m_parent_id;
  G4double      m_track_length;
  G4ThreeVector m_vertex_position;
  G4ThreeVector m_vertex_momentum;
  G4double      m_vertex_kinetic_energy; // sqrt(p^2+m^2)-m

public:
  void     AddEnergyDeposit( G4double de ){ m_energy_deposit += de; }
  G4String GetDetectorName( void ) const { return m_detector_name; }
  G4String GetParticleName( void ) const { return m_particle_name; }
  const G4ThreeVector& GetPosition( void ) const { return m_position; }
  const G4ThreeVector& GetMomentum( void ) const { return m_momentum; }
  G4double GetTime( void ) const { return m_time; }
  G4double GetEnergyDeposit( void ) const { return m_energy_deposit; }
  G4int    GetTrackID( void ) const { return m_track_id; }
  G4int    GetParticleID( void ) const { return m_particle_id; }
  G4int    GetDetectorID( void ) const { return m_detector_id; }
  G4double GetMass( void ) const { return m_mass; }
  G4double GetCharge( void ) const { return m_charge; }
  G4int    GetParentID( void ) const { return m_parent_id; }
  G4double GetTrackLength( void ) const { return m_track_length; }
  const G4ThreeVector& GetVertexPosition( void ) const
  { return m_vertex_position; }
  const G4ThreeVector& GetVertexMomentum( void ) const
  { return m_vertex_momentum; }
  G4double      GetVertexKineticEnergy( void ) const
  { return m_vertex_kinetic_energy; }
  void          Print( void ) const;
};

//_____________________________________________________________________________
inline G4String
VHitInfo::ClassName( void )
{
  static G4String s_name("VHitInfo");
  return s_name;
}

#endif

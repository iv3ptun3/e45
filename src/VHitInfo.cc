// -*- C++ -*-

#include "VHitInfo.hh"

#include <G4Step.hh>
#include <G4Track.hh>

#include "FuncName.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
VHitInfo::VHitInfo( const G4String& name, G4Step* step )
  : m_detector_name( name ),
    m_position(),
    m_momentum(),
    m_time(),
    m_energy_deposit(),
    m_track_id(),
    m_particle_id(),
    m_detector_id(),
    m_mass(),
    m_charge(),
    m_parent_id(),
    m_track_length(),
    m_vertex_position(),
    m_vertex_momentum(),
    m_vertex_kinetic_energy()
{
  if( step ){
    const auto track = step->GetTrack();
    const auto point = step->GetPreStepPoint();
    m_position = point->GetPosition();
    m_momentum = point->GetMomentum();
    m_time = point->GetGlobalTime();
    m_energy_deposit = step->GetTotalEnergyDeposit();
    m_track_id = track->GetTrackID();
    m_particle_id = track->GetDefinition()->GetPDGEncoding();
    m_detector_id = point->GetPhysicalVolume()->GetCopyNo();
    m_mass = track->GetDynamicParticle()->GetMass();
    m_charge = track->GetDynamicParticle()->GetCharge();
    m_parent_id = track->GetParentID();
    m_track_length = track->GetTrackLength();
    m_vertex_position = track->GetVertexPosition();
    m_vertex_momentum = track->GetVertexMomentumDirection();
    m_vertex_kinetic_energy = track->GetVertexKineticEnergy();
  }
}

//_____________________________________________________________________________
VHitInfo::~VHitInfo( void )
{
}

//_____________________________________________________________________________
void
VHitInfo::Print( void ) const
{
  PrintHelper helper( 3, std::ios::fixed, G4cout );
  G4cout << FUNC_NAME << " " << m_detector_name << G4endl
	 << "   Position            = " << m_position*(1/CLHEP::mm)
	 << " mm" << G4endl
	 << "   Momentum            = " << m_momentum*(1/CLHEP::GeV)
	 << " GeV/c" << G4endl
	 << "   Time                = " << m_time/CLHEP::ns << " ns" << G4endl
	 << "   EnergyDeposit       = " << m_energy_deposit/CLHEP::MeV
	 << " MeV" << G4endl
	 << "   TrackID             = " << m_track_id << G4endl
	 << "   ParticleID          = " << m_particle_id << G4endl
	 << "   DetectorID          = " << m_detector_id << G4endl
	 << "   Mass                = " << m_mass/CLHEP::MeV << " MeV" << G4endl
	 << "   Charge              = " << m_charge << G4endl
	 << "   ParentID            = " << m_parent_id << G4endl
	 << "   TrackLength         = " << m_track_length/CLHEP::mm
	 << " mm" << G4endl
	 << "   VertexPosition      = " << m_vertex_position*(1/CLHEP::mm)
	 << " mm" << G4endl
	 << "   VertexMomentum      = " << m_vertex_momentum
	 << " GeV/c" << G4endl
	 << "   VertexKineticEnergy = " << m_vertex_kinetic_energy/CLHEP::MeV
	 << " MeV" << G4endl;
}

// -*- C++ -*-

#include "TPCPadSD.hh"

#include <G4VPhysicalVolume.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VTouchable.hh>
#include <G4TouchableHistory.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4DynamicParticle.hh>
#include <G4DecayProducts.hh>
#include <G4PhysicsLogVector.hh>
#include <G4ParticleChangeForDecay.hh>
#include <G4DecayProcessType.hh>
#include <Randomize.hh>

#include "ConfMan.hh"
#include "FuncName.hh"
#include "TPCPadHit.hh"
#include "padHelper.hh"

namespace
{
  const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
TPCPadSD::TPCPadSD( const G4String& name )
  : G4VSensitiveDetector( name ),
    m_gem_discharge( gConf.Get<G4int>( "GemDischarge" ) ),
    m_gem_fix_dead( gConf.Get<G4int>( "GemFixDead" ) ),
    m_gem_dead_plane( gConf.Get<G4int>( "GemDeadPlane" ) ),
    m_gem_dead_plane_division( gConf.Get<G4int>( "GemDeadPlaneDivision" ) ),
    m_dead_area( gConf.Get<G4double>( "DeadArea" ) )
{
  collectionName.insert("hit");
  G4cout << FUNC_NAME << G4endl
	 << "   Study on GEM discharge = " << m_gem_discharge << G4endl;
  if( m_gem_discharge == 5 ){
    G4cout << "   Designed a configuration of GEM discharge area" << G4endl;
    G4cout << "   Width of GEM electrodes are fixed!!!!" << G4endl;
  }
}

//_____________________________________________________________________________
TPCPadSD::~TPCPadSD( void )
{
}

//_____________________________________________________________________________
void
TPCPadSD::Initialize( G4HCofThisEvent* HCTE )
{
  ntrk=0;
  hitsCollection = new G4THitsCollection<TPCPadHit>( SensitiveDetectorName,
						     collectionName[0] );
  // push H.C. to "Hit Collection of This Event"
  G4int hcid = GetCollectionID(0);
  HCTE->AddHitsCollection( hcid, hitsCollection );

  if( m_gem_discharge > 0 ){
    if( m_gem_fix_dead == 0 ){
      select_plane=CLHEP::RandFlat::shoot(0.,1.);
      num_plane=0;
      if(select_plane<0.25){
	num_plane=1;
      }else if(select_plane>=0.25 && select_plane<0.5){
	num_plane=2;
      }else if(select_plane>=0.5 && select_plane<0.75){
	num_plane=3;
      }else if(select_plane>=0.75){
	num_plane=4;
      }
      num_deadarea = 250/m_dead_area;
      if(m_dead_area*num_deadarea < 250){
	num_deadarea=num_deadarea+1.;
      }
      select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);
    }else if( m_gem_fix_dead == 1 ){
      num_plane = m_gem_dead_plane;
      select_dead = m_gem_dead_plane_division;
    }
  }
}

//_____________________________________________________________________________
G4bool
TPCPadSD::ProcessHits( G4Step* aStep, G4TouchableHistory* /* ROhist */ )
{
  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();
  const G4Track* aTrack = aStep->GetTrack();
  G4int copyNo = preStepPoint -> GetPhysicalVolume()->GetCopyNo();

  if(preStepPoint-> GetStepStatus() != fGeomBoundary) return false;
  //  if(preStepPoint-> GetStepStatus() == fGeomBoundary){
  G4String particleName;
  
	if(aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge() == 0.)
    return false;
	
  particleName = aStep-> GetTrack()-> GetDefinition()-> GetParticleName();

  G4String particleType;
  particleType = aTrack->GetDefinition()->GetParticleType();

  ///lepton rejection
  if(particleType == "lepton")
    return false;
  if(particleType == "gamma")
    return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4ThreeVector pos= preStepPoint-> GetPosition();
  G4double hitx=pos.getX();
  G4double hitz=pos.getZ();

  //dead layer cross x and z
  /*
  if(abs(pos.getX())<5.)
    return false;
  if(abs(pos.getZ())<5.)
    return false;
  */
  //dead layer 45 deg
  G4double deadlayer=12.;
  if(hitx < (hitz+deadlayer/2.*sqrt(2)) && hitx > (hitz-deadlayer/2.*sqrt(2))  )    return false;
  if(hitx < (-hitz+deadlayer/2.*sqrt(2)) && hitx > (-hitz-deadlayer/2.*sqrt(2))  )    return false;
  //////////end dead layer

  if( m_gem_discharge > 0 ){
    if( m_gem_discharge == 1 ){
      if(num_plane==1){
	if(hitz>=0 && fabs(hitx)<=hitz){	  ////select discharge area
	  if(hitz>m_dead_area*(select_dead-1) && hitz<m_dead_area*select_dead){
	    return false;
	  }
	}
      }else if(num_plane==2){
	if(hitx>=0 && hitx>fabs(hitz)){
	  if(hitx>m_dead_area*(select_dead-1) && hitx<m_dead_area*select_dead){
	    //	    G4cout<<"num_plane:select_dead--->"<<num_plane<<":"<<select_dead<<G4endl;
	    return false;
	  }
	}
      }else if(num_plane==3){
	if(hitz<0 && fabs(hitx)<fabs(hitz)){
	  //	  G4int select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);
	  //	  G4cout<<"num_plane:select_dead--->"<<num_plane<<":"<<select_dead<<G4endl;
	  if(fabs(hitz)>m_dead_area*(select_dead-1) && fabs(hitz)<m_dead_area*select_dead){
	    return false;
	  }
	}
      }else if(num_plane==4){
	if(hitx<0 && fabs(hitx)>fabs(hitz)){
	  //	  G4int select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);
	  //	  G4cout<<"num_plane:select_dead--->"<<num_plane<<":"<<select_dead<<G4endl;
	  if(fabs(hitx)>m_dead_area*(select_dead-1) && fabs(hitx)<m_dead_area*select_dead){
	    return false;
	  }
	}
      }
    }  ////end design #1

    ///////////////design 2
    if( m_gem_discharge == 2 ){
      if(num_plane==1){
	if(hitz>=0 && fabs(hitx)<=hitz){	  ////select discharge area
	  if( hitx < -6*sqrt(2) + hitz - m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx > -6*sqrt(2) + hitz-m_dead_area*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==2){
	if(hitx>=0 && hitx>fabs(hitz)){
	  if( hitx > 6*sqrt(2)+ -hitz + m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx < 6*sqrt(2) + -hitz + m_dead_area*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==3){
	if(hitz<0 && fabs(hitx)<fabs(hitz)){
	  //	  G4int select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);
	  //	  select_dead=3;
	  if( hitx > 6*sqrt(2)+ hitz + m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx < 6*sqrt(2)+ hitz + m_dead_area*sqrt(2)*(select_dead) ){
	    //	    G4cout<<"ok"<<G4endl;
	    return false;
	  }
	}
      }else if(num_plane==4){
	if(hitx<0 && fabs(hitx)>fabs(hitz)){
	  //	  G4int select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);
	  //	  select_dead=3;
	  if( hitx < -6*sqrt(2) -hitz - m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx > -6*sqrt(2) -hitz - m_dead_area*sqrt(2)*(select_dead) ){
	    //	    G4cout<<"ok"<<G4endl;
	    return false;
	  }
	}
      }
    }  //end design #2

    ///////////////design 5, designed shape
    if( m_gem_discharge == 5 ){
      G4double divi_width[6]={44.,43.,43.,43.,43.,43.0};
      //      G4double divi_width3[15]={25.,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,25.,25.,25.,25.};
      G4double divi_width3[20]={13.};
      //      for(G4int ii=0;ii<20;ii++){
      //	divi_width3[ii]=13.;
      //	//	G4cout<<divi_width3[ii]<<G4endl;
      //      }
      if(num_plane==1){
	if(select_dead>6){
	  G4cout<<"exceed maximum division"<<G4endl;
	  exit(-1);
	}
	if(hitz>=0 && fabs(hitx)<=hitz){	  ////select discharge area
	  if( hitx < -6.5*sqrt(2) + hitz - divi_width[select_dead-1]*sqrt(2)*(select_dead-1)
		     && hitx > -6.5*sqrt(2) + hitz-divi_width[select_dead-1]*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==2){
	if(select_dead>6.){
	  G4cout<<"exceed maximum division"<<G4endl;
	  exit(-1);
	}
	if(hitx>=0 && hitx>fabs(hitz)){
	  if( hitx > 6.5*sqrt(2)+ -hitz + divi_width[select_dead-1]*sqrt(2)*(select_dead-1)
		     && hitx < 6.5*sqrt(2) + -hitz + divi_width[select_dead-1]*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==3){
	//	if(select_dead>15.){
	if(select_dead>20.){
	  G4cout<<"exceed maximum division"<<G4endl;
	  exit(-1);
	}
	if(hitz<0 && fabs(hitx)<fabs(hitz)){
	  if( hitx > 6.5*sqrt(2)+ hitz + divi_width3[select_dead-1]*sqrt(2)*(select_dead-1)
		     && hitx < 6.5*sqrt(2)+ hitz + divi_width3[select_dead-1]*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==4){
	if(select_dead>6.){
	  G4cout<<"exceed maximum division"<<G4endl;
	  exit(-1);
	}
	if(hitx<0 && fabs(hitx)>fabs(hitz)){
	  if( hitx < -6.5*sqrt(2) -hitz - divi_width[select_dead-1]*sqrt(2)*(select_dead-1)
		     && hitx > -6.5*sqrt(2) -hitz - divi_width[select_dead-1]*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }
    }  //end design #5


    ///////////////design 4 anticlockwise
    if( m_gem_discharge == 4 ){
      if(num_plane==1){
	if(hitz>=0 && fabs(hitx)<=hitz){	  ////select discharge area
	  if( hitx > +6*sqrt(2) - hitz + m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx < +6*sqrt(2) - hitz + m_dead_area*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==2){
	if(hitx>=0 && hitx>fabs(hitz)){
	  if( hitx > 6*sqrt(2)+ +hitz + m_dead_area*sqrt(2)*(select_dead-1)
	      && hitx < 6*sqrt(2) + +hitz + m_dead_area*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==3){
	if(hitz<0 && fabs(hitx)<fabs(hitz)){
	  //	    G4cout<<"n"<<G4endl;
	  if( hitx < -6*sqrt(2)- hitz - m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx > -6*sqrt(2) - hitz - m_dead_area*sqrt(2)*(select_dead) ){
	    //	    G4cout<<"ok"<<G4endl;
	    return false;
	  }
	}
      }else if(num_plane==4){
	if(hitx<0 && fabs(hitx)>fabs(hitz)){
	  if( hitx < -6*sqrt(2) +hitz - m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx > -6*sqrt(2) +hitz - m_dead_area*sqrt(2)*(select_dead) ){
	    //	    G4cout<<"ok"<<G4endl;
	    return false;
	  }
	}
      }
    }  //end design #4 anticlockwise


///////////////design 3, test on plane #3
    if( m_gem_discharge == 3 ){
      if(num_plane==1){
	if(hitz>=0 && fabs(hitx)<=hitz){	  ////select discharge area
	  if( hitx < hitz - m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx > hitz-m_dead_area*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==2){
	if(hitx>=0 && hitx>fabs(hitz)){
	  if( hitx > -hitz + m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx < -hitz + m_dead_area*sqrt(2)*(select_dead) ){
	    return false;
	  }
	}
      }else if(num_plane==3){
	if(hitz<0 && fabs(hitx)<fabs(hitz)){

	  if( hitx > hitz + m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx < hitz + m_dead_area*sqrt(2)*(select_dead) && hitx<0.){
		//	    G4cout<<num_plane<<":"<<m_dead_area<<":"<<select_dead<<G4endl;
	    return false;
	  }
	}
      }else if(num_plane==4){
	if(hitx<0 && fabs(hitx)>fabs(hitz)){
	  //	  G4int select_dead=CLHEP::RandFlat::shoot(1.,num_deadarea+1.);
	  //	  select_dead=3;
	  if( hitx < -hitz - m_dead_area*sqrt(2)*(select_dead-1)
		     && hitx > -hitz - m_dead_area*sqrt(2)*(select_dead) ){
	    //	    G4cout<<"ok"<<G4endl;
	    return false;
	  }
	}
      }
    }  //end design #3 test



      //  if(hitx < (hitz+5*sqrt(2)) && hitx > (hitz-5*sqrt(2))  )
      //    return false;
      //  if(hitx < (-hitz+5*sqrt(2)) && hitx > (-hitz-5*sqrt(2))  )
      //    return false;
    }//end discharge
//  G4double edep = aStep->GetTotalEnergyDeposit();
 

	// if (edep<0.0005){
  //   return false;
  // }
  G4ThreeVector VertexPosition = aTrack->GetVertexPosition();
  G4ThreeVector VertexMomentum = aTrack->GetVertexMomentumDirection();
  G4double VertexEnergy = aTrack -> GetVertexKineticEnergy(); // Ek = sqrt(p^2+m^2)-m
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();

  G4ThreeVector mom= preStepPoint-> GetMomentum();
  G4double tof= preStepPoint-> GetGlobalTime();
  G4double beta= preStepPoint-> GetBeta();
  G4int tid =  aStep-> GetTrack()-> GetTrackID();
  G4int pid =  aStep-> GetTrack()-> GetDefinition() -> GetPDGEncoding();
  G4double mass = aStep -> GetTrack()->GetDynamicParticle()->GetMass();
  //  G4cout<<mass<<G4endl;
  G4int charge = aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge();
  //  G4double edep = aStep->GetTotalEnergyDeposit();
  //  G4double edep = aStep->GetNonIonizingEnergyDeposit();
  //  G4double deltaep = aStep->GetDeltaEnergy();
  G4double tlength = aStep->GetTrack()-> GetTrackLength();
  G4double slength = aStep->GetTrack()-> GetStepLength();
  // G4double length = aStep-> GetStepLength();
  G4int parentID =  aStep-> GetTrack()-> GetParentID();

  G4int parentID_pid = aStep-> GetTrack()->GetDynamicParticle()->GetDefinition()->GetPDGEncoding();

  //  G4int test =  aStep-> GetTrack()-> GetParentID() -> GetPDGCharge();
  //  G4cout<<"-------------------"<<G4endl;
  //  G4cout<<"pid :"<<pid<<G4endl;
  //  G4cout<<"test :"<<tlength<<G4endl;
  //  G4cout<<"steplength :"<<slength<<G4endl;
  //  G4cout<<"length :"<<length<<G4endl;
  //  G4cout<<"dedx :"<<edep<<std::setw(2)<<":"<<deltaep<<G4endl;
  //  G4cout<<"copyNo :"<<copyNo<<G4endl;
  //  G4cout<<"tlength"<<tlength<<G4endl;

  G4int iLay_copyNo=copyNo;
  G4int iPad = padHelper::findPadID(hitz, hitx);
  G4int iLay= padHelper::getLayerID(iPad);
  G4int iRow= padHelper::getRowID(iPad);
  
	G4double PadLen = padHelper::getLength(iLay);
	G4ThreeVector PadPos(hitx,0,hitz + 143); 
	G4ThreeVector MomT(mom.x(),0,mom.z());	
	G4double alpha = PadPos.theta()-MomT.theta();
	G4double PathT = PadLen * 1./cos(alpha);
	G4double Pitch = mom.y()/MomT.mag();
	G4double Path = PathT * sqrt(1+Pitch*Pitch);
	slength = Path;

	
	G4double edepMean =TPCdEdx(mass,beta)*Path; 
	//	G4double IonEn = (0.9 * 188 + 0.1 * 41.7)*eV;
//	G4double edepSig = sqrt(edepMean / IonEn ) * IonEn;
	G4double edepSig =TPCdEdxSig(mass,mom.mag())*Path;
	if(edepSig / edepMean < 0.01 or edepSig/edepMean > 0.5)edepSig = 0.2*edepMean;
	G4double edep = G4RandGauss::shoot(edepMean,edepSig);
	if(edep <0.1* edepMean)edep = 0.1*edepMean;
#ifdef DEBUG
  
  //for test 
  G4double radius = sqrt( hitx*hitx + (hitz+143.)*(hitz+143.));
  TVector3 Point = padHelper::getPoint(iPad);
  G4int iPad_re = padHelper::findPadID(Point.z(), Point.x());
  G4cout<<"hitx = "<< hitx
   	<<", pointx "<< Point.x()
   	<<", hitz =" << hitz
  	<<", pointz "<< Point.z()
   	<<", radius = "<< radius
   	<<", iPad ="<<iPad
	<<", iPad_re ="<<iPad_re
   	<<", iLay_copyNo = " <<iLay_copyNo
   	<< ", iLay = "<<iLay<<G4endl;

  G4cout<<"dx ="<<hitx-Point.x()
	<<", dz ="<<hitz-Point.z()<<std::endl;
#endif
  
  



  // if(iLay_copyNo!=iLay)
  //   G4cerr << FUNC_NAME 
  // 	   << " iLay_copyNo = " <<iLay_copyNo
  // 	   << ", iLay = "<<iLay<<G4endl;
  

  G4String name = physVol->GetName();

  



  if(name=="TPC_PV"){
    name="PadPV-1";
  }

  //  if(name=="Target_PV"){
  //    name="PadPV-1";
  //  }

  //  sscanf(name,"PadPV%d",&iLay);
  sscanf(name,"PadPV%d",&iLay_copyNo);
  TPCPadHit* ahit= new TPCPadHit(pos, mom, tof, tid, pid, iLay, iRow, beta, edep,parentID,tlength, mass, charge, VertexPosition, VertexMomentum, VertexEnergy,slength, parentID_pid );

  hitsCollection-> insert(ahit);
  return true;
  //}End of fboundary
}

//_____________________________________________________________________________
void
TPCPadSD::EndOfEvent(G4HCofThisEvent* /* HCTE */ )
{
}

//_____________________________________________________________________________
void
TPCPadSD::DrawAll( void )
{
}
G4double
TPCPadSD::TPCdEdx(Double_t mass/*MeV/c2*/, Double_t beta){

  Double_t rho=0.; //[g cm-3]
  Double_t ZoverA=0.; //[mol g-1]
  Double_t I=0.; //[eV]
  Double_t density_effect_par[6]={0.}; //Sternheimer’s parameterization
  //P10  
	rho = TMath::Power(10.,-3)*(0.9*1.662 + 0.1*0.6672);
	ZoverA = 17.2/37.6;
	I = 0.9*188.0 + 0.1*41.7;
	density_effect_par[0] = 0.9*0.19714 + 0.1*0.09253;
	density_effect_par[1] = 0.9*2.9618 + 0.1*3.6257;
	density_effect_par[2] = 0.9*1.7635 + 0.1*1.6263;
	density_effect_par[3] = 0.9*4.4855 + 0.1*3.9716;
	density_effect_par[4] = 0.9*11.9480 + 0.1*9.5243;
	density_effect_par[5] = 0.;

  Double_t Z = 1.;
  Double_t me = 0.5109989461; //[MeV]
  Double_t K = 0.307075; //[MeV cm2 mol-1]
  Double_t constant = rho*K*ZoverA; //[MeV cm-1]
	constant = constant;
  Double_t I2 = I*I; //Mean excitaion energy [eV]
  Double_t beta2 = beta*beta;
  Double_t gamma2 = 1./(1.-beta2);
  Double_t MeVToeV = TMath::Power(10.,6);
  Double_t Wmax = 2*me*beta2*gamma2/((me/mass+1.)*(me/mass+1.)+2*(me/mass)*(TMath::Sqrt(gamma2)-1));
  Double_t delta = DensityEffectCorrection(TMath::Sqrt(beta2*gamma2), density_effect_par);
  Double_t dedx = constant*Z*Z/beta2*(0.5*TMath::Log(2*me*beta2*gamma2*Wmax*MeVToeV*MeVToeV/I2) - beta2 - 0.5*delta);
	
	G4double conversion_factor = 11073.3;
  return conversion_factor*dedx;

}
double
TPCPadSD::TPCdEdxSig(Double_t mass/*MeV/c2*/, Double_t mom){
	double par[3]={0,0,0} ;
	if(mass < 0.2){//pion;
		par[0] = 7.792;
		par[1] =	-8.704;
		par[2] = 4.477;
	}
	else if(mass > 0.7){//proton or heavier;
		par[0] = 33.92;
		par[1] = -26.24;
		par[2] = 6.259; 
	}
	else {//Kaon -> temp value
		par[0] = 7.792;
		par[1] =	-8.704;
		par[2] = 4.477;
	}
	double value = par[0]+par[1]*mom+par[2]*mom*mom;
	return value;
}
G4double
TPCPadSD::DensityEffectCorrection(Double_t betagamma, Double_t *par){

    //reference : Sternheimer’s parameterizatio(PDG)
    //notation : par[0] : a, par[1] : k, par[2] : x0, par[3] : x1, par[4] : _C, par[5] : delta0
    Double_t constant = 2*TMath::Log(10);
    Double_t delta = 0.;
    Double_t X = log10(betagamma);
    if(X<=par[2]) delta = par[5]*TMath::Power(10., 2*(X - par[2]));
    else if(par[2]<X && X<par[3]) delta = constant*X - par[4] + par[0]*pow((par[3] - X), par[1]);
    else if(X>=par[3]) delta = constant*X - par[4];

  return delta;

}
//_____________________________________________________________________________
void
TPCPadSD::PrintAll( void )
{
  hitsCollection-> PrintAllHits();
}

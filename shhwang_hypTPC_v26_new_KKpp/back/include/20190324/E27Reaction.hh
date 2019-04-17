// ====================================================================
//   E27Reaction.hh
//
// ====================================================================

#ifndef E27REACTION_H 
#define E27REACTION_H 
 
#include "TPCAnaManager.hh"
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>

class TPCPrimaryGeneratorAction;
class G4Event;

class E27Reaction
{
public:
  E27Reaction(TPCPrimaryGeneratorAction * PGAction)
    :pGen(PGAction)
  {}
  ~E27Reaction(){}

  void E27_beamthrough(G4Event* anEvent);//#2701
  void E27_Kptest(G4Event* anEvent);//#2702
  void E27_Kpp_F_LambdaP(G4Event* anEvent);//#2703
  void E27_Kpp_F_SigmaZP(G4Event* anEvent);//#2704
  void E27_Kpp_F_LambdaPizP(G4Event* anEvent);//#2705
  void E27_Kpp_F_SigmaZPizP(G4Event* anEvent);//#2706
  void E27_Kpp_F_SigmaPPimP(G4Event* anEvent);//#2707
  void E27_K11B_Lambda10Be(G4Event* anEvent);//#2708
  void E27_Kptest2(G4Event* anEvent);//#2709  

private:
  TPCPrimaryGeneratorAction *pGen;
  
};


#endif

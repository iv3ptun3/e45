// ====================================================================
//   KKPpReaction.hh
//
// ====================================================================

#ifndef KKppREACTION_H 
#define KKppREACTION_H 
 
#include "TPCAnaManager.hh"
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>

class TPCPrimaryGeneratorAction;
class G4Event;

class KKppReaction
{
public:
  KKppReaction(TPCPrimaryGeneratorAction * PGAction)
    :pGen(PGAction)
  {}
  ~KKppReaction(){}

  void KKpp_LL1(G4Event* anEvent);//#3001
  void KKpp_LL2(G4Event* anEvent);//#3002
  void KKpp_LSmPip(G4Event* anEvent);//#3003
  void KKpp_LSpPim(G4Event* anEvent);//#3004

  void JAMInput(G4Event* anEvent, TTree *t1);//#3101
  void KKpp_BeamThrough1(G4Event* anEvent);//#3102
  void JAMInput_K0(G4Event* anEvent, TTree* t1);//#3103
  void JAMInput_K0bar(G4Event* anEvent, TTree* t1);//#3104
private:
  TPCPrimaryGeneratorAction *pGen;
  
};


#endif

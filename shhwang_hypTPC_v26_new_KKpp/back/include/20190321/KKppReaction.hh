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

private:
  TPCPrimaryGeneratorAction *pGen;
  
};


#endif

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

  void E27_beamthrough(G4Event* anEvent);
  void E27_Kptest(G4Event* anEvent);

private:
  TPCPrimaryGeneratorAction *pGen;
  
};


#endif

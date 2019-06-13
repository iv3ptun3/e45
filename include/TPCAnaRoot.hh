#ifndef TPCANAROOT_H
#define TPCANAROOT_H 1

#include "G4ThreeVector.hh"
#include "globals.hh"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TBranch.h"
#include "common.hh"

class TH1F;
class TH2F;
class TFile;
class TTree;

//const int MaxTrack = 1560*4;
//const int MaxTrack = 78*4;
const G4int MaxTrack = 54*20;

const G4int MaxTrackTPC = 50;

const G4int MaxNthLay = 40;
const G4int MaxNthPad = 250;

const G4int MaxTrig = 500;
const G4int MaxDCHit = 500;
const G4int MaxPrimaryParticle = 10;

struct Tree1Ev {

  G4int ev;        // Event number  
  G4double pg[4];  // 4-momentum for inncident beam
  G4int gen;        // generator number
  G4int mode;        // mode number

  G4double mm_d; //missing-mass of d(pi,K) reaction 
  G4double mm_p;  //missing-mass for proton target kinematic
  G4double mm; //missing-mass (should be same as mm_d)
  G4double theta; //theta of scat K
  G4double theta_scat; //theta of (pi, K) reaction
  G4double theta_CM; //theta of (pi, K) reaction (CM flame)

  G4int HitNum_K;
  G4int HitNumAC_K;
  G4int HitNumNBAR_K;
  G4int HitNumDC_K;
  G4int HitNumFTOF_K;
  G4int HitNumCH_K;
  G4int HitNumScint_K;
  //  int tpctrNum_K;
  G4int HitNumTarget_K;


  G4int HitNum_p;
  G4int HitNumAC_p;
  G4int HitNumNBAR_p;
  G4int HitNumDC_p;
  G4int HitNumFTOF_p;
  G4int HitNumCH_p;
  G4int HitNumScint_p;
  //  int tpctrNum_p;
  G4int HitNumTarget_p;


  G4int npid;
  G4double x0[MaxPrimaryParticle][3];   
  G4double p0[MaxPrimaryParticle][5];  
  G4double pt0[MaxPrimaryParticle];  
  G4double mass0[MaxPrimaryParticle];  
  G4int pid0[MaxPrimaryParticle];  
  G4double theta0[MaxPrimaryParticle];  

  /* number of ntrks in TPC by shhwang*/
  G4int ntrtpc;
  G4double trpptpc[MaxTrackTPC];
  G4double trpxtpc[MaxTrackTPC];
  G4double trpytpc[MaxTrackTPC];
  G4double trpztpc[MaxTrackTPC];
  G4double trpttpc[MaxTrackTPC];

  G4double trpptpcfit[MaxTrackTPC];
  G4double trpxtpcfit[MaxTrackTPC];
  G4double trpytpcfit[MaxTrackTPC];
  G4double trpztpcfit[MaxTrackTPC];
  G4double trpttpcfit[MaxTrackTPC];

  G4int trqqtpc[MaxTrackTPC];
  G4int trpidtpc[MaxTrackTPC];
  G4int trparentidtpc[MaxTrackTPC];
  G4int trparentid_pid_tpc[MaxTrackTPC];
  G4double trpmtpc[MaxTrackTPC];
  G4double trdetpc[MaxTrackTPC];
  G4double trlentpc[MaxTrackTPC];
  G4double trdedxtpc[MaxTrackTPC];
  G4double trdedxtrtpc[MaxTrackTPC]; //trancated mean, but now just mean

  G4int trlaytpc[MaxTrackTPC];


  G4double trvtxpxtpc[MaxTrackTPC];
  G4double trvtxpytpc[MaxTrackTPC];
  G4double trvtxpztpc[MaxTrackTPC];
  G4double trvtxpptpc[MaxTrackTPC];

  G4double trvtxxtpc[MaxTrackTPC];
  G4double trvtxytpc[MaxTrackTPC];
  G4double trvtxztpc[MaxTrackTPC];

  G4double trvtxxtpcfit[MaxTrackTPC];
  G4double trvtxytpcfit[MaxTrackTPC];
  G4double trvtxztpcfit[MaxTrackTPC];

  /////PAD multiplicity & ASAD multiplicy
  G4int nthlay[MaxTrack];
  G4int nthpad[MaxTrack];
  G4int laypad[MaxTrack][MaxNthLay][MaxNthPad]; //[layer][pad number]


  ///////////////
  G4int nttpc;                 // Number of Hit in Pads
  G4int ntrk[MaxTrack];        // Number of Track

  G4int ititpc[MaxTrack];      // Track ID
  G4int idtpc[MaxTrack];       // Particle ID
  G4double xtpc[MaxTrack];     // coordinates
  G4double ytpc[MaxTrack];     // coordinates
  G4double ztpc[MaxTrack];     // coordinates

  G4double x0tpc[MaxTrack];    // coordinates
  G4double y0tpc[MaxTrack];    // coordinates
  G4double z0tpc[MaxTrack];    // coordinates
  G4double resoX[MaxTrack];    // coordinates

  G4double pxtpc[MaxTrack];    // momentum
  G4double pytpc[MaxTrack];    // momentum
  G4double pztpc[MaxTrack];    // momentum
  G4double pptpc[MaxTrack];    // momentum
  G4double masstpc[MaxTrack];    // mass

  G4double betatpc[MaxTrack];    // beta

  G4double edeptpc[MaxTrack];    // Energy deposit
  G4double dedxtpc[MaxTrack];    // Energy deposit/dx
  G4double slengthtpc[MaxTrack];    // Energy deposit/dx

  G4int laytpc[MaxTrack];      // number of pad layer
  G4int rowtpc[MaxTrack];      // number of pad raw  
  G4double toftpc[MaxTrack];   // tof
  G4int parentID[MaxTrack];      // parent id
  G4double cir_r[MaxTrack];   // fit radius
  G4double cir_x[MaxTrack];   // fit center x
  G4double cir_z[MaxTrack];   // fit center z
  G4double cir_fit[MaxTrack];   // fit center fit
  G4int vtx_flag[MaxTrack]; // flag, how to estimate vtx 
  G4double a_fory[MaxTrack]; // co-efficient a for linear track (y, theta) 
  G4double b_fory[MaxTrack]; // co-efficient b for linear track (y, theta) 


  G4int ntsc;                 // Number of Hit in Scint.
  G4int tidsc[MaxTrig];       // Track ID   
  G4int pidsc[MaxTrig];	    // Particle ID
  G4int didsc[MaxTrig];	    // detector ID
  G4double masssc[MaxTrig];	    // particle mass ID
  G4int qqsc[MaxTrig];	    // particle mass ID
  G4double xsc[MaxTrig];      // coordinates
  G4double ysc[MaxTrig];      // coordinates
  G4double zsc[MaxTrig];      // coordinates
  G4double pxsc[MaxTrig];     // momentum
  G4double pysc[MaxTrig];     // momentum
  G4double pzsc[MaxTrig];     // momentum
  G4double ppsc[MaxTrig];     // momentum
  G4double tofsc[MaxTrig];    // tof
  G4int scpID[MaxTrig];    //parent id

  G4double trvtxxscint[MaxTrig];
  G4double trvtxyscint[MaxTrig];
  G4double trvtxzscint[MaxTrig];

  G4double trvtxpxscint[MaxTrig];
  G4double trvtxpyscint[MaxTrig];
  G4double trvtxpzscint[MaxTrig];
  G4double trvtxppscint[MaxTrig];
  G4double lengthsc[MaxTrig];


  ///ac
  G4int ntac;                 // Number of Hit in Acint.
  G4int tidac[MaxTrig];       // Track ID   
  G4int pidac[MaxTrig];	    // Particle ID
  G4int didac[MaxTrig];	    // detector ID
  G4double massac[MaxTrig];	    // particle mass ID
  G4int qqac[MaxTrig];	    // particle mass ID
  G4double xac[MaxTrig];      // coordinates
  G4double yac[MaxTrig];      // coordinates
  G4double zac[MaxTrig];      // coordinates
  G4double pxac[MaxTrig];     // momentum
  G4double pyac[MaxTrig];     // momentum
  G4double pzac[MaxTrig];     // momentum
  G4double ppac[MaxTrig];     // momentum
  G4double tofac[MaxTrig];    // tof
  G4int acpID[MaxTrig];    //parent id

  G4double trvtxxac[MaxTrig];
  G4double trvtxyac[MaxTrig];
  G4double trvtxzac[MaxTrig];

  G4double trvtxpxac[MaxTrig];
  G4double trvtxpyac[MaxTrig];
  G4double trvtxpzac[MaxTrig];
  G4double trvtxppac[MaxTrig];
  G4double lengthac[MaxTrig];


  ////////nbar
  G4int ntnbar;                 // Number of Hit in Nbarint.
  G4int tidnbar[MaxTrig];       // Trnbark ID   
  G4int pidnbar[MaxTrig];	    // Particle ID
  G4int didnbar[MaxTrig];	    // detector ID
  G4double massnbar[MaxTrig];	    // particle mass ID
  G4int qqnbar[MaxTrig];	    // particle mass ID
  G4double xnbar[MaxTrig];      // coordinates
  G4double ynbar[MaxTrig];      // coordinates
  G4double znbar[MaxTrig];      // coordinates
  G4double pxnbar[MaxTrig];     // momentum
  G4double pynbar[MaxTrig];     // momentum
  G4double pznbar[MaxTrig];     // momentum
  G4double ppnbar[MaxTrig];     // momentum
  G4double tofnbar[MaxTrig];    // tof
  G4int nbarpID[MaxTrig];    //parent id

  G4double trvtxxnbar[MaxTrig];
  G4double trvtxynbar[MaxTrig];
  G4double trvtxznbar[MaxTrig];

  G4double trvtxpxnbar[MaxTrig];
  G4double trvtxpynbar[MaxTrig];
  G4double trvtxpznbar[MaxTrig];
  G4double trvtxppnbar[MaxTrig];
  G4double lengthnbar[MaxTrig];




  ///dc
  G4int ntdc;                 // Number of Hit in Dcint.
  G4int tiddc[MaxTrig];       // Trdck ID   
  G4int piddc[MaxTrig];	    // Particle ID
  G4int diddc[MaxTrig];	    // detector ID
  G4double massdc[MaxTrig];	    // particle mass ID
  G4int qqdc[MaxTrig];	    // particle mass ID
  G4double xdc[MaxTrig];      // coordinates
  G4double ydc[MaxTrig];      // coordinates
  G4double zdc[MaxTrig];      // coordinates
  G4double pxdc[MaxTrig];     // momentum
  G4double pydc[MaxTrig];     // momentum
  G4double pzdc[MaxTrig];     // momentum
  G4double ppdc[MaxTrig];     // momentum
  G4double tofdc[MaxTrig];    // tof
  G4int dcpID[MaxTrig];    //parent id

  G4double trvtxxdc[MaxTrig];
  G4double trvtxydc[MaxTrig];
  G4double trvtxzdc[MaxTrig];

  G4double trvtxpxdc[MaxTrig];
  G4double trvtxpydc[MaxTrig];
  G4double trvtxpzdc[MaxTrig];
  G4double trvtxppdc[MaxTrig];
  G4double lengthdc[MaxTrig];


  ///ch
  G4int ntch;                 // Number of Hit in Chint.
  G4int tidch[MaxTrig];       // Trchk ID   
  G4int pidch[MaxTrig];	    // Particle ID
  G4int didch[MaxTrig];	    // detector ID
  G4double massch[MaxTrig];	    // particle mass ID
  G4int qqch[MaxTrig];	    // particle mass ID
  G4double xch[MaxTrig];      // coordinates
  G4double ych[MaxTrig];      // coordinates
  G4double zch[MaxTrig];      // coordinates
  G4double pxch[MaxTrig];     // momentum
  G4double pych[MaxTrig];     // momentum
  G4double pzch[MaxTrig];     // momentum
  G4double ppch[MaxTrig];     // momentum
  G4double tofch[MaxTrig];    // tof
  G4int chpID[MaxTrig];    //parent id

  G4double trvtxxch[MaxTrig];
  G4double trvtxych[MaxTrig];
  G4double trvtxzch[MaxTrig];

  G4double trvtxpxch[MaxTrig];
  G4double trvtxpych[MaxTrig];
  G4double trvtxpzch[MaxTrig];
  G4double trvtxppch[MaxTrig];
  G4double lengthch[MaxTrig];



  ///ftof
  G4int ntftof;                 // Number of Hit in Ftofint.
  G4int tidftof[MaxTrig];       // Trftofk ID   
  G4int pidftof[MaxTrig];	    // Particle ID
  G4int didftof[MaxTrig];	    // detector ID
  G4double massftof[MaxTrig];	    // particle mass ID
  G4int qqftof[MaxTrig];	    // particle mass ID
  G4double xftof[MaxTrig];      // coordinates
  G4double yftof[MaxTrig];      // coordinates
  G4double zftof[MaxTrig];      // coordinates
  G4double pxftof[MaxTrig];     // momentum
  G4double pyftof[MaxTrig];     // momentum
  G4double pzftof[MaxTrig];     // momentum
  G4double ppftof[MaxTrig];     // momentum
  G4double tofftof[MaxTrig];    // tof
  G4int ftofpID[MaxTrig];    //parent id

  G4double trvtxxftof[MaxTrig];
  G4double trvtxyftof[MaxTrig];
  G4double trvtxzftof[MaxTrig];

  G4double trvtxpxftof[MaxTrig];
  G4double trvtxpyftof[MaxTrig];
  G4double trvtxpzftof[MaxTrig];
  G4double trvtxppftof[MaxTrig];
  G4double lengthftof[MaxTrig];

  /////////////It does not use in current G4
  /*  G4int ntfdc;                     // Number of Hit in FDC
  G4int tidfdc[MaxTrackFDC];       // Track ID   
  G4int pidfdc[MaxTrackFDC];	 // Particle ID
  G4int didfdc[MaxTrackFDC];	 // detector ID
  G4double xfdc[MaxTrackFDC];      // coordinates
  G4double yfdc[MaxTrackFDC];      // coordinates
  G4double zfdc[MaxTrackFDC];      // coordinates
  G4double pxfdc[MaxTrackFDC];     // momentum
  G4double pyfdc[MaxTrackFDC];     // momentum
  G4double pzfdc[MaxTrackFDC];     // momentum
  G4double toffdc[MaxTrackFDC];    // tof
  */

  G4int targethits;                     // Number of Hit in TARGET
  G4int targetpid[MaxTrack];       // Track ID   
  G4int targetparentid[MaxTrack];	 // Particle ID
  G4int targettid[MaxTrack];	 // detector ID
  G4double targetpos[MaxTrack][3];      // coordinates
  G4double targetvtx[MaxTrack][3];      // coordinates

};


class TPCAnaRoot
{
private:
  TFile *rootfile;
  TH1F *hpos[3];
  TH1F *hmom[3];
  TH1F *htime;

  TTree *tree;
  Tree1Ev tree1ev;

public:
  TPCAnaRoot();
  ~TPCAnaRoot();
  void BeginOfRunAction(int runnum);
  void EndOfRunAction();
  void BeginOfEventAction();
  void BeginOfTrackingAction();
  G4int EndOfTrackingAction();
  void FillMom(G4double *mom);
  void FillPos(G4double *Pos);
  void FillNtrk(G4int ntrk);
  void FillPos0(G4double *Pos,G4double Reso);
  void FillTime(G4double time);
  void FillBeta(G4double beta);
  void FillEdep(G4double edep);
  void FilldEdx(G4double dedx);
  void FillsLength(G4double slength);

  void FillTrackID(G4int tid){
    tree1ev.ititpc[tree1ev.nttpc] = tid;
  }
  void FillParticleID(G4int pid){
    tree1ev.idtpc[tree1ev.nttpc] = pid;
  }

  void FillPadLay(G4int ilay){
    tree1ev.laytpc[tree1ev.nttpc] = ilay;
  }

  void FillPadRow(G4int irow){
    tree1ev.rowtpc[tree1ev.nttpc] = irow;
  }

  void incHit(){
    tree1ev.nttpc += 1;
  }

  //  void tpcHit(){
  //    tree1ev.ntrtpc += 1;
  //  }
  
  void FillPrimaryParticle(int id, double* x0, double* p0, int pid);
  void FillNumOfK(int HitNum_K, int HitNumAC_K, int HitNumNBAR_K, 
		  int HitNumDC_K, int HitNumCH_K, int HitNumFTOF_K, 
		  int HitNumScint_K, int HitNumTarget_K);

  void FillNumOfp(int HitNum_p, int HitNumAC_p, int HitNumNBAR_p, 
		  int HitNumDC_p, int HitNumCH_p, int HitNumFTOF_p, 
		  int HitNumScint_p, int HitNumTarget_p);

  void FillBeam(G4double px, G4double py, G4double pz);
  void FillPrimaryInfo(G4double mm_d, G4double mm_p, G4double theta, 
		       G4double theta_scat, G4double theta_CM);

  void FillLayerPad(G4int nlay,G4int npad);

  void FillGenMode(G4int gen,G4int mode);

  void FillTPCData(G4double tpcpx1,G4double tpcpy1,
		   G4double tpcpz1,G4double tpcpp1, 
		   G4int tpcpid1, G4int tpcparentid1, G4int tpcparentid_pid1, 
		   G4int tpcqq1, 
		   G4double tpcpm1, G4double tpcde1, G4double tpclen1, 
		   G4double tpcdedx1,G4double tpcdedxtr1, G4int tpclay1,
		   G4double tpcvtxpx1,G4double tpcvtxpy1, G4double tpcvtxqz1, 
		   G4double tpcvtxx1,G4double tpcvtxy1, G4double tpcvtxz1,
		   G4double tpcvtxx1fit,G4double tpcvtxy1fit, G4double tpcvtxz1fit,
		   G4double tpcpxfit1,G4double tpcpyfit1,G4double tpcpzfit1,G4double tpcptfit1,
		   G4double cir_r, G4double cir_x, G4double cir_z, G4double cir_fit,
		   G4int vtx_flag, G4double a_fory, G4double b_fory
		   );

  void FillScintData(G4double time, G4double* pos, G4double* mom,
		     G4int tid, G4int pid, G4int did, G4double mass,G4int qq,G4int parentid,
		     G4double scintvtxpx1,G4double scintvtxpy1, G4double scintvtxqz1, 
		     G4double scintvtxx1,G4double scintvtxy1, G4double scintvtxz1, G4double tlength1);  

  void FillACData(G4double time, G4double* pos, G4double* mom,
		     G4int tid, G4int pid, G4int did, G4double mass,G4int qq,G4int parentid,
		     G4double acvtxpx1,G4double acvtxpy1, G4double acvtxqz1, 
		     G4double acvtxx1,G4double acvtxy1, G4double acvtxz1, G4double tlength1);  


  void FillNBARData(G4double time, G4double* pos, G4double* mom,
		     G4int tid, G4int pid, G4int did, G4double mass,G4int qq,G4int parentid,
		     G4double nbarvtxpx1,G4double nbarvtxpy1, G4double nbarvtxqz1, 
		     G4double nbarvtxx1,G4double nbarvtxy1, G4double nbarvtxz1, G4double tlength1);  


  void FillDCData(G4double time, G4double* pos, G4double* mom,
		     G4int tid, G4int pid, G4int did, G4double mass,G4int qq,G4int parentid,
		     G4double vtxpx1,G4double vtxpy1, G4double vtxqz1, 
		     G4double vtxx1,G4double vtxy1, G4double vtxz1, G4double tlength1);  

  void FillCHData(G4double time, G4double* pos, G4double* mom,
		     G4int tid, G4int pid, G4int did, G4double mass,G4int qq,G4int parentid,
		     G4double vtxpx1,G4double vtxpy1, G4double vtxqz1, 
		     G4double vtxx1,G4double vtxy1, G4double vtxz1, G4double tlength1);  

  void FillFTOFData(G4double time, G4double* pos, G4double* mom,
		     G4int tid, G4int pid, G4int did, G4double mass,G4int qq,G4int parentid,
		     G4double vtxpx1,G4double vtxpy1, G4double vtxqz1, 
		     G4double vtxx1,G4double vtxy1, G4double vtxz1, G4double tlength1);  



  /*  void FillFDCData(G4double time, G4double* pos, G4double* mom,
			  G4int tid, G4int pid, G4int did);  
  */

  void FillTargetData(G4int nhit, G4int particleid, G4int parentid, G4int trackid,
		      G4ThreeVector pos, G4ThreeVector vtx);  


  void FillTree();
};

#endif

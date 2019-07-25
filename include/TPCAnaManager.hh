// -*- C++ -*-

#ifndef TPC_ANA_MANAGER_HH
#define TPC_ANA_MANAGER_HH

#include <G4ThreeVector.hh>

#include <TVector3.h>

struct Track;

class VHitInfo;

//_____________________________________________________________________________
static const G4int MaxHits    = 500;
static const G4int MaxHitsTPC = 500;

//const int MaxTrack = 1560*4;
//const int MaxTrack = 78*4;
const G4int MaxTrack = 54*20;

const G4int MaxNthLay = 40;
const G4int MaxNthPad = 250;

const G4int MaxPrimaryParticle = 10;

void initTrack(Track* aTrack);
void initTrack_ku(Track* aTrack);
int setInitialPara(Track* aTrack, double* initPara);
int setVirtualPlane(Track* aTrack);
void minuitInit(double printLevel);

static const int MAXtpctrNum=30;
static const int MAXtpctrhitNum=500;

//_____________________________________________________________________________
struct CounterData
{
  G4int ntrk;
  G4double resoX;
  G4int trackID;
  G4int particleID;
  G4double time;
  G4double beta;
  G4double edep;
  G4double dedx;
  G4double slength;
  G4double mass;
  G4double pos0[3];
  G4double pos[3];
  G4double mom[4];
  G4int iLay;
  G4int iPad;
  G4int iRow;
  G4int parentID;
};

//_____________________________________________________________________________
struct TPCData
{
  G4int tpctr;
  G4int tpcpid;
  G4int tpcparentid;
  G4int tpcparentid_pid;
  G4double tpcpx;
  G4double tpcpy;
  G4double tpcpz;
  G4double tpcpp;

  G4double tpcpxfit;
  G4double tpcpyfit;
  G4double tpcpzfit;
  G4double tpcppfit;
  G4double tpcptfit;

  G4int tpcqq;
  G4double tpcpm;
  G4double tpcde;
  G4double tpclen;
  G4double tpcdedx;
  G4int tpclay;
  G4double tpcvtxpx;
  G4double tpcvtxpy;
  G4double tpcvtxpz;
  G4double tpcvtxx;
  G4double tpcvtxy;
  G4double tpcvtxz;
  //  G4double tpcene2
};

struct TargetData
{
  G4int targettr;
  G4double targetpx;
  G4double targetpy;
  G4double targetpz;
  G4double targetpp;
  G4int targetqq;

  G4double targetpm;
  G4double targetparticleid;
  G4double targetparentid;
  G4double targettrackid;

  G4double targetde;
  G4double targetlen;
  G4double targetdedx;
  G4ThreeVector targetvtxmom;
  //  G4double targetvtxpy;
  //  G4double targetvtxpz;
  //  G4double targetvtxx;
  //  G4double targetvtxy;
  G4ThreeVector targetvtx;

  //  G4double targetposx;
  //  G4double targetposy;
  //  G4double targetposz;
  G4ThreeVector targetpos;

};

//_____________________________________________________________________________
struct ScintData
{
  G4int trackID;
  G4double massSH;
  G4int qqSH;
  G4int particleID;
  G4int detectorID;
  G4int parentID;
  G4double time;
  G4double length;
  G4double pos[3];
  G4double mom[4];
  G4double scintvtxpx;
  G4double scintvtxpy;
  G4double scintvtxpz;
  G4double scintvtxx;
  G4double scintvtxy;
  G4double scintvtxz;
};

//_____________________________________________________________________________
struct ACData
{
  G4int trackID;
  G4double massSH;
  G4int qqSH;
  G4int particleID;
  G4int detectorID;
  G4int parentID;
  G4double time;
  G4double length;
  G4double pos[3];
  G4double mom[4];
  G4double acvtxpx;
  G4double acvtxpy;
  G4double acvtxpz;
  G4double acvtxx;
  G4double acvtxy;
  G4double acvtxz;
};

//_____________________________________________________________________________
struct NBARData
{
  G4int trackID;
  G4double massSH;
  G4int qqSH;
  G4int particleID;
  G4int detectorID;
  G4int parentID;
  G4double time;
  G4double length;
  G4double pos[3];
  G4double mom[4];
  G4double nbarvtxpx;
  G4double nbarvtxpy;
  G4double nbarvtxpz;
  G4double nbarvtxx;
  G4double nbarvtxy;
  G4double nbarvtxz;
};

//_____________________________________________________________________________
struct DCData
{
  G4int trackID;
  G4double massSH;
  G4int qqSH;
  G4int particleID;
  G4int detectorID;
  G4int parentID;
  G4double time;
  G4double length;
  G4double pos[3];
  G4double mom[4];
  G4double vtxpx;
  G4double vtxpy;
  G4double vtxpz;
  G4double vtxx;
  G4double vtxy;
  G4double vtxz;
};

//_____________________________________________________________________________
struct SCHData
{
  G4int trackID;
  G4double massSH;
  G4int qqSH;
  G4int particleID;
  G4int detectorID;
  G4int parentID;
  G4double time;
  G4double length;
  G4double pos[3];
  G4double mom[4];
  G4double vtxpx;
  G4double vtxpy;
  G4double vtxpz;
  G4double vtxx;
  G4double vtxy;
  G4double vtxz;
};

//_____________________________________________________________________________
struct FTOFData
{
  G4int trackID;
  G4double massSH;
  G4int qqSH;
  G4int particleID;
  G4int detectorID;
  G4int parentID;
  G4double time;
  G4double length;
  G4double pos[3];
  G4double mom[4];
  G4double vtxpx;
  G4double vtxpy;
  G4double vtxpz;
  G4double vtxx;
  G4double vtxy;
  G4double vtxz;
};

//_____________________________________________________________________________
struct FDCData
{
  G4int trackID;
  G4int particleID;
  G4int detectorID;
  G4double time;
  G4double pos[3];
  G4double mom[4];
};

//_____________________________________________________________________________
struct PrimaryBeam
{
  G4double pg[4];                        // 4-mom 0: px, 1: py, 2: pz, 3: ene
  G4int gen;
  G4int mode;
};

//_____________________________________________________________________________
struct PrimaryParticle
{
  G4int NumOfParticle;                   // Number of Primary particle
  G4double x0[MaxPrimaryParticle][3];    // Vertex position  0: x, 1: y, 2: z
  G4double p0[MaxPrimaryParticle][5];    // 4-mom, mass 0: px, 1: py, 2: pz, 3: ene, 4: mass
  G4int pid0[MaxPrimaryParticle];
};

//_____________________________________________________________________________
struct PrimaryInfo
{
  G4double mm_d;
  G4double mm_p;
  G4double theta;
  G4double theta_scat;
  G4double theta_CM;
};

//_____________________________________________________________________________
struct Event
{
  Int_t evnum; // Event number
  TVector3* pb; // momentum of inncident beam
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
  G4int HitNumSCH_K;
  G4int HitNumScint_K;
  //  int tpctrNum_K;
  G4int HitNumTarget_K;


  G4int HitNum_p;
  G4int HitNumAC_p;
  G4int HitNumNBAR_p;
  G4int HitNumDC_p;
  G4int HitNumFTOF_p;
  G4int HitNumSCH_p;
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
  G4double trpptpc[MaxHitsTPC];
  G4double trpxtpc[MaxHitsTPC];
  G4double trpytpc[MaxHitsTPC];
  G4double trpztpc[MaxHitsTPC];
  G4double trpttpc[MaxHitsTPC];

  G4double trpptpcfit[MaxHitsTPC];
  G4double trpxtpcfit[MaxHitsTPC];
  G4double trpytpcfit[MaxHitsTPC];
  G4double trpztpcfit[MaxHitsTPC];
  G4double trpttpcfit[MaxHitsTPC];

  G4int trqqtpc[MaxHitsTPC];
  G4int trpidtpc[MaxHitsTPC];
  G4int trparentidtpc[MaxHitsTPC];
  G4int trparentid_pid_tpc[MaxHitsTPC];
  G4double trpmtpc[MaxHitsTPC];
  G4double trdetpc[MaxHitsTPC];
  G4double trlentpc[MaxHitsTPC];
  G4double trdedxtpc[MaxHitsTPC];
  G4double trdedxtrtpc[MaxHitsTPC]; //trancated mean, but now just mean

  G4int trlaytpc[MaxHitsTPC];


  G4double trvtxpxtpc[MaxHitsTPC];
  G4double trvtxpytpc[MaxHitsTPC];
  G4double trvtxpztpc[MaxHitsTPC];
  G4double trvtxpptpc[MaxHitsTPC];

  G4double trvtxxtpc[MaxHitsTPC];
  G4double trvtxytpc[MaxHitsTPC];
  G4double trvtxztpc[MaxHitsTPC];

  G4double trvtxxtpcfit[MaxHitsTPC];
  G4double trvtxytpcfit[MaxHitsTPC];
  G4double trvtxztpcfit[MaxHitsTPC];

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
  G4int tidsc[MaxHits];       // Track ID
  G4int pidsc[MaxHits];	    // Particle ID
  G4int didsc[MaxHits];	    // detector ID
  G4double masssc[MaxHits];	    // particle mass ID
  G4int qqsc[MaxHits];	    // particle mass ID
  G4double xsc[MaxHits];      // coordinates
  G4double ysc[MaxHits];      // coordinates
  G4double zsc[MaxHits];      // coordinates
  G4double pxsc[MaxHits];     // momentum
  G4double pysc[MaxHits];     // momentum
  G4double pzsc[MaxHits];     // momentum
  G4double ppsc[MaxHits];     // momentum
  G4double tofsc[MaxHits];    // tof
  G4int scpID[MaxHits];    //parent id

  G4double trvtxxscint[MaxHits];
  G4double trvtxyscint[MaxHits];
  G4double trvtxzscint[MaxHits];

  G4double trvtxpxscint[MaxHits];
  G4double trvtxpyscint[MaxHits];
  G4double trvtxpzscint[MaxHits];
  G4double trvtxppscint[MaxHits];
  G4double lengthsc[MaxHits];

  // HTOF
  G4int nhHtof;
  G4int tidHtof[MaxHits];       // Track ID
  G4int pidHtof[MaxHits];	    // Particle ID
  G4int didHtof[MaxHits];	    // detector ID
  G4double massHtof[MaxHits];	    // particle mass ID
  G4int qqHtof[MaxHits];	    // particle mass ID
  G4double xHtof[MaxHits];      // coordinates
  G4double yHtof[MaxHits];      // coordinates
  G4double zHtof[MaxHits];      // coordinates
  G4double pxHtof[MaxHits];     // momentum
  G4double pyHtof[MaxHits];     // momentum
  G4double pzHtof[MaxHits];     // momentum
  G4double ppHtof[MaxHits];     // momentum
  G4double tofHtof[MaxHits];    // tof
  G4int HtofpID[MaxHits];    //parent id
  G4double trvtxxHtof[MaxHits];
  G4double trvtxyHtof[MaxHits];
  G4double trvtxzHtof[MaxHits];
  G4double trvtxpxHtof[MaxHits];
  G4double trvtxpyHtof[MaxHits];
  G4double trvtxpzHtof[MaxHits];
  G4double trvtxppHtof[MaxHits];
  G4double lengthHtof[MaxHits];

  ///ac
  G4int ntac;                 // Number of Hit in Acint.
  G4int tidac[MaxHits];       // Track ID
  G4int pidac[MaxHits];	    // Particle ID
  G4int didac[MaxHits];	    // detector ID
  G4double massac[MaxHits];	    // particle mass ID
  G4int qqac[MaxHits];	    // particle mass ID
  G4double xac[MaxHits];      // coordinates
  G4double yac[MaxHits];      // coordinates
  G4double zac[MaxHits];      // coordinates
  G4double pxac[MaxHits];     // momentum
  G4double pyac[MaxHits];     // momentum
  G4double pzac[MaxHits];     // momentum
  G4double ppac[MaxHits];     // momentum
  G4double tofac[MaxHits];    // tof
  G4int acpID[MaxHits];    //parent id

  G4double trvtxxac[MaxHits];
  G4double trvtxyac[MaxHits];
  G4double trvtxzac[MaxHits];

  G4double trvtxpxac[MaxHits];
  G4double trvtxpyac[MaxHits];
  G4double trvtxpzac[MaxHits];
  G4double trvtxppac[MaxHits];
  G4double lengthac[MaxHits];


  ////////nbar
  G4int ntnbar;                 // Number of Hit in Nbarint.
  G4int tidnbar[MaxHits];       // Trnbark ID
  G4int pidnbar[MaxHits];	    // Particle ID
  G4int didnbar[MaxHits];	    // detector ID
  G4double massnbar[MaxHits];	    // particle mass ID
  G4int qqnbar[MaxHits];	    // particle mass ID
  G4double xnbar[MaxHits];      // coordinates
  G4double ynbar[MaxHits];      // coordinates
  G4double znbar[MaxHits];      // coordinates
  G4double pxnbar[MaxHits];     // momentum
  G4double pynbar[MaxHits];     // momentum
  G4double pznbar[MaxHits];     // momentum
  G4double ppnbar[MaxHits];     // momentum
  G4double tofnbar[MaxHits];    // tof
  G4int nbarpID[MaxHits];    //parent id

  G4double trvtxxnbar[MaxHits];
  G4double trvtxynbar[MaxHits];
  G4double trvtxznbar[MaxHits];

  G4double trvtxpxnbar[MaxHits];
  G4double trvtxpynbar[MaxHits];
  G4double trvtxpznbar[MaxHits];
  G4double trvtxppnbar[MaxHits];
  G4double lengthnbar[MaxHits];




  ///dc
  G4int ntdc;                 // Number of Hit in Dcint.
  G4int tiddc[MaxHits];       // Trdck ID
  G4int piddc[MaxHits];	    // Particle ID
  G4int diddc[MaxHits];	    // detector ID
  G4double massdc[MaxHits];	    // particle mass ID
  G4int qqdc[MaxHits];	    // particle mass ID
  G4double xdc[MaxHits];      // coordinates
  G4double ydc[MaxHits];      // coordinates
  G4double zdc[MaxHits];      // coordinates
  G4double pxdc[MaxHits];     // momentum
  G4double pydc[MaxHits];     // momentum
  G4double pzdc[MaxHits];     // momentum
  G4double ppdc[MaxHits];     // momentum
  G4double tofdc[MaxHits];    // tof
  G4int dcpID[MaxHits];    //parent id

  G4double trvtxxdc[MaxHits];
  G4double trvtxydc[MaxHits];
  G4double trvtxzdc[MaxHits];

  G4double trvtxpxdc[MaxHits];
  G4double trvtxpydc[MaxHits];
  G4double trvtxpzdc[MaxHits];
  G4double trvtxppdc[MaxHits];
  G4double lengthdc[MaxHits];

  // SCH
  G4int ntsch;                 // Number of Hit in Chint.
  G4int tidsch[MaxHits];       // Trchk ID
  G4int pidsch[MaxHits];	    // Particle ID
  G4int didsch[MaxHits];	    // detector ID
  G4double masssch[MaxHits];	    // particle mass ID
  G4int qqsch[MaxHits];	    // particle mass ID
  G4double xsch[MaxHits];      // coordinates
  G4double ysch[MaxHits];      // coordinates
  G4double zsch[MaxHits];      // coordinates
  G4double pxsch[MaxHits];     // momentum
  G4double pysch[MaxHits];     // momentum
  G4double pzsch[MaxHits];     // momentum
  G4double ppsch[MaxHits];     // momentum
  G4double tofsch[MaxHits];    // tof
  G4int schpID[MaxHits];    //parent id

  G4double trvtxxsch[MaxHits];
  G4double trvtxysch[MaxHits];
  G4double trvtxzsch[MaxHits];

  G4double trvtxpxsch[MaxHits];
  G4double trvtxpysch[MaxHits];
  G4double trvtxpzsch[MaxHits];
  G4double trvtxppsch[MaxHits];
  G4double lengthsch[MaxHits];



  ///ftof
  G4int ntftof;                 // Number of Hit in Ftofint.
  G4int tidftof[MaxHits];       // Trftofk ID
  G4int pidftof[MaxHits];	    // Particle ID
  G4int didftof[MaxHits];	    // detector ID
  G4double massftof[MaxHits];	    // particle mass ID
  G4int qqftof[MaxHits];	    // particle mass ID
  G4double xftof[MaxHits];      // coordinates
  G4double yftof[MaxHits];      // coordinates
  G4double zftof[MaxHits];      // coordinates
  G4double pxftof[MaxHits];     // momentum
  G4double pyftof[MaxHits];     // momentum
  G4double pzftof[MaxHits];     // momentum
  G4double ppftof[MaxHits];     // momentum
  G4double tofftof[MaxHits];    // tof
  G4int ftofpID[MaxHits];    //parent id

  G4double trvtxxftof[MaxHits];
  G4double trvtxyftof[MaxHits];
  G4double trvtxzftof[MaxHits];

  G4double trvtxpxftof[MaxHits];
  G4double trvtxpyftof[MaxHits];
  G4double trvtxpzftof[MaxHits];
  G4double trvtxppftof[MaxHits];
  G4double lengthftof[MaxHits];

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

//_____________________________________________________________________________
class TPCAnaManager
{
public:
  static G4String ClassName( void );
  static TPCAnaManager& GetInstance( void );
  ~TPCAnaManager( void );

private:
  TPCAnaManager( void );
  TPCAnaManager( const TPCAnaManager& );
  TPCAnaManager& operator=( const TPCAnaManager& );

private:
  std::vector<VHitInfo*> m_htof_hc;

  TargetData targetData[MaxTrack];
  CounterData counterData[MaxTrack];
  TPCData tpcData[MAXtpctrNum];
  ScintData scintData[MaxHits];
  ACData acData[MaxHits];
  NBARData nbarData[MaxHits];

  DCData dcData[MaxHits];
  SCHData schData[MaxHits];
  FTOFData ftofData[MaxHits];

  PrimaryBeam primaryBeam;
  PrimaryParticle primaryParticle;
  PrimaryInfo primaryInfo;
  int HitNum;
  int HitNumAC;
  int HitNumNBAR;
  int HitNumDC;
  int HitNumFTOF;
  int HitNumSCH;
  int HitNumScint;
  int tpctrNum;
  int HitNumTarget;

  int HitNum_K;
  int HitNumAC_K;
  int HitNumNBAR_K;
  int HitNumDC_K;
  int HitNumFTOF_K;
  int HitNumSCH_K;
  int HitNumScint_K;
  //  int tpctrNum_K;
  int HitNumTarget_K;

  int HitNum_p;
  int HitNumAC_p;
  int HitNumNBAR_p;
  int HitNumDC_p;
  int HitNumFTOF_p;
  int HitNumSCH_p;
  int HitNumScint_p;  //  int tpctrNum_K;
  int HitNumTarget_p;

  // G4double vtxxfit[MAXtpctrNum];//read fit parameters
  // G4double vtxyfit[MAXtpctrNum];//read fit parameters
  // G4double vtxzfit[MAXtpctrNum];//read fit parameters

  // G4double vtxpxfit[MAXtpctrNum];//read fit parameters
  // G4double vtxpyfit[MAXtpctrNum];//read fit parameters
  // G4double vtxpzfit[MAXtpctrNum];//read fit parameters


  G4double mean[MAXtpctrNum];//read fit parameters
  G4double trmean[MAXtpctrNum];//read fit parameters
  G4double cir_r[MAXtpctrNum];//read fit parameters
  G4double error[MAXtpctrNum];//read fit parameters
  G4double chi2[MAXtpctrNum];//read fit parameters
  G4double ndf[MAXtpctrNum];//read fit parameters
  G4double Pz[MAXtpctrNum];//read fit parameters

  ////////////////////getenv parameters
  G4double pad_length_in;
  G4double pad_length_out;
  G4double pad_gap;
  G4double pad_in_width;
  G4double pad_out_width;
  G4double pad_in_num;
  G4double pad_out_num;
  G4double truncated_mean_cut;
  G4int trigger_env;
  G4double target_pos_z;
  G4int env_on_off_helm;
  G4double env_helm_field;
  G4int env_pad_config;
  G4int env_Experiment_num;

  /*
    G4double angle[40]={0};
    G4double seg_angle[40]={0};
    G4int numpads[40]={0};

    G4double pad_in[40]={0};
    G4double pad_out[40]={0};
    G4double tpc_rad=250;
  */
  G4double angle[40];
  G4double seg_angle[40];
  G4double seg_width[40];
  G4int numpads[40];

  G4double pad_in[40];
  G4double pad_out[40];
  G4double tpc_rad;

public:
  void BeginOfRunAction(int runnum);
  void EndOfRunAction();
  void BeginOfEventAction();
  int  EndOfEventAction();

  void SetTPCData(G4int tpctr, G4int tpcpid, G4int tpcparentid, G4int tpcparentid_pid, G4double tpcpx, G4double tpcpy,G4double tpcpz,G4double tpcpp,  G4int tpcqq, G4double tpcpm, G4double tpcde, G4double tpclen, G4int tpclay,
		  G4double vtxpxtpc2,G4double vtxpytpc2,G4double vtxpztpc2,
		  G4double vtxxtpc2,G4double vtxytpc2,G4double vtxztpc2, G4double vtxenetpc2);

  void SetCounterData(G4int ntrk, G4double time, G4ThreeVector pos, G4ThreeVector mom,
		      G4int track, G4int particle,
		      G4int iLay, G4int iRow, G4double beta, G4double edep,
		      G4int parentid, G4double tlength, G4double slength);
  void SetHTOFData( VHitInfo* hit );
  void SetTargetData( G4int nhits, G4ThreeVector pos, G4ThreeVector mom,
		      G4int track, G4int particle,
		      G4int parentid, G4ThreeVector vtxpos,
		      G4ThreeVector vtxmom, G4double vtxene );

  //  void SetCounterData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
  //		      G4int track, G4int particle, G4int ilay, G4int iRaw);
  void SetScintData( G4double time, G4ThreeVector pos, G4ThreeVector mom,
		     G4int track, G4int particle, G4int detector,
		     G4double mass, G4int qq,G4int parentid,
		     G4ThreeVector vtxpos, G4ThreeVector vtxmom,
		     G4double vtxene, G4double tlength );
  void SetACData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
		 G4int track, G4int particle, G4int detector,G4double mass, G4int qq,G4int parentid, G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength);


  void SetNBARData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
		   G4int track, G4int particle, G4int detector,G4double mass, G4int qq,G4int parentid, G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength);


  void SetDCData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
		 G4int track, G4int particle, G4int detector,G4double mass, G4int qq,
		 G4int parentid, G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength);

  void SetSCHData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
		  G4int track, G4int particle, G4int detector,G4double mass, G4int qq,
		  G4int parentid, G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength);

  void SetFTOFData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
		   G4int track, G4int particle, G4int detector,G4double mass, G4int qq,
		   G4int parentid, G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength);

  //  void SetFDCData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
  //				 G4int track, G4int particle, G4int detector);
  void SetPrimaryBeam( const G4ThreeVector& p );
  void SetPrimaryBeam( G4double px, G4double py, G4double pz );

  void SetGeneratorID(G4int gen);
  void SetModeID(G4int mode);

  void SetNumberOfPrimaryParticle(G4int num);
  void SetPrimaryParticle(G4double px, G4double py, G4double pz);
  void SetPrimaryParticle(G4int id, G4double px, G4double py, G4double pz, G4double mass);
  void SetPrimaryParticle(G4int id, G4double px, G4double py, G4double pz, G4double mass, G4int pid0);
  void SetPrimaryVertex(G4int id, G4double x, G4double y, G4double z);
  void SetPrimaryInfo(G4double mm_d, G4double mm_p, G4double theta, G4double theta_scat, G4double theta_cm);

  //int CircleIntersect(double x1, double y1, double r1, double x2, double y2, double r2, double inter1[2], double inter2[2])
  int CircleIntersect(double x1, double y1, double r1, double x2, double y2, double r2,
		      double ca1, double cb1, double ct01, int qq1,
		      double ca2, double cb2, double ct02, int qq2,
		      double inter1[3], double inter2[3])
  {
    // function inputs: x1, y1, r1, x2, y2, r2
    // function output: inter1, inter2 = coordinates of intersections

    double d,e,f,g,a,b,c;
    double x,y,discrim;

    if( x1 == x2 && y1 == y2 ){
      G4cout << x1 << " " << y1 << " " << r1 << G4endl;
      return 0;
    }
    //    G4cout << x1 << " " << y1 << " " << r1 << G4endl;
    //    G4cout << x2 << " " << y2 << " " << r2 << G4endl;

    d = -0.5*( r1*r1 - r2*r2 - x1*x1 - y1*y1 + x2*x2 + y2*y2 );
    e =  0.5*( r1*r1 + r2*r2 - x1*x1 - y1*y1 - x2*x2 - y2*y2 );
    if( fabs(y1-y2) < 1.0e-20 ) {
      x = d/(x1 - x2);
      a = 1.0;
      b = 0.0;
      c = x*x - x*(x1 + x2) - e;
      discrim = -4*a*c;
      if( discrim < 0 ) return 0;
      y = sqrt(discrim) / (2*a);
      inter1[0] = x;
      inter1[1] = y;
      inter2[0] = x;
      inter2[1] = -y;
      return 1;
    }
    f = (x1 - x2) / (y1 - y2);
    g = d / (y1 - y2);
    // cout << "d=" << d << " e=" << e << " f=" << f << " g=" << g << endl;

    a = 1. + f*f;
    b = f*(y1 + y2) - 2*f*g - (x1 + x2);
    c = g*g - g*(y1 + y2) - e;
    //    G4cout << "a=" << a << " b=" << b << " c=" << c << G4endl;

    discrim = b*b - 4*a*c;
    //    G4cout << "discrim = " << discrim << G4endl;
    if( discrim < 0 ) return 0;
    inter1[0] = (-b + sqrt(discrim) ) / (2*a);
    inter1[1] = g - f*inter1[0];
    inter2[0] = (-b - sqrt(discrim) ) / (2*a);
    inter2[1] = g - f*inter2[0];


    double theta11 = atan2(inter1[1]-y1, inter1[0]-x1);
    double theta21 = atan2(inter1[1]-y2, inter1[0]-x2);
    double theta12 = atan2(inter2[1]-y1, inter2[0]-x1);
    double theta22 = atan2(inter2[1]-y2, inter2[0]-x2);

    double tmp_y11 = -1.*(double)qq1*ca1*r1*(theta11-ct01)+cb1;
    double tmp_y21 = -1.*(double)qq2*ca2*r2*(theta21-ct02)+cb2;
    double tmp_y12 = -1.*(double)qq1*ca1*r1*(theta12-ct01)+cb1;
    double tmp_y22 = -1.*(double)qq2*ca2*r2*(theta22-ct02)+cb2;

    // std::cout<<"theta11="<<theta11<<", theta21="<<theta21
    // 	     <<", theta12="<<theta12<<", theta22="<<theta22<<std::endl;

    // std::cout<<"tmp_y11="<<tmp_y11<<", tmp_y21="<<tmp_y21<<std::endl;
    // std::cout<<"tmp_y12="<<tmp_y12<<", tmp_y22="<<tmp_y22<<std::endl;
    //getchar();

    inter1[2] = (tmp_y11+tmp_y21)/2.;
    inter2[2] = (tmp_y12+tmp_y22)/2.;


    return 1;
  }


  double linearFitter(const int np,
		      const double *x,
		      const double *y, double *er,
		      double *a, double *b){

    int i;

    double alpha=0.;
    double beta=0.;
    double gamma=0.;
    double AA=0;
    double BB=0;

    for(i=0;i<np;i++){
      alpha+=x[i]*x[i]/er[i]/er[i];
      beta+=x[i]/er[i]/er[i];
      gamma+=1./er[i]/er[i];
      AA+=y[i]*x[i]/er[i]/er[i];
      BB+=y[i]/er[i]/er[i];

      //  G4cout<<"x test: "<<x[i]<<G4endl;
      //  G4cout<<"y test: "<<y[i]<<G4endl;
    }

    //  G4cout<<"beta test: "<<beta<<G4endl;
    //  G4cout<<"alpha test: "<<alpha<<G4endl;
    //  G4cout<<"gamma test: "<<gamma<<G4endl;

    *a=( gamma*AA -  beta*BB )/(alpha*gamma-beta*beta);
    *b=(-beta *AA + alpha*BB )/(alpha*gamma-beta*beta);

    //  G4cout<<"a test: "<<(*a)<<G4endl;
    //  G4cout<<"b test: "<<(*b)<<G4endl;

    return 1.;
  }


  double circleFit(const double *mX,const double *mY,const double *mZ,
		   const int npoints, double* mXCenter, double* mYCenter,
		   double* mRadius, double* Pz_, double* a_forz,
		   double* b_forz, double* theta0_fory)
  {
    double xx, yy, xx2, yy2;
    double f, g, h, p, q, t, g0, g02, a, b, c, d;
    double xroot, ff, fp, xd, yd, g1;
    double dx, dy, dradius2, xnom;

    double xgravity = 0.0;
    double ygravity = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    double xy = 0.0;
    double xx2y2 = 0.0;
    double yx2y2 = 0.0;
    double x2y22 = 0.0;
    double radius2 = 0.0;

    double mVariance = 0.0;

    if (npoints <= 3){
      fprintf(stderr,"circleFit: npoints %d <= 3\n",npoints);
      return -1;
    }else  if (npoints > 499){
      fprintf(stderr,"circleFit: npoints %d > 499\n",npoints);
      return -1;
    }

    for (int i=0; i<npoints; i++) {
      xgravity += mX[i];
      ygravity += mY[i];
    }
    xgravity /= npoints;
    ygravity /= npoints;

    for (int i=0; i<npoints; i++) {
      xx  = mX[i]-xgravity;
      yy  = mY[i]-ygravity;
      xx2 = xx*xx;
      yy2 = yy*yy;
      x2  += xx2;
      y2  += yy2;
      xy  += xx*yy;
      xx2y2 += xx*(xx2+yy2);
      yx2y2 += yy*(xx2+yy2);
      x2y22 += (xx2+yy2)*(xx2+yy2);
    }
    if (xy == 0.){
      fprintf(stderr,"circleFit: xy = %f,    grav=%f, %f\n",xy,xgravity,ygravity);
      return -1;
    }

    f = (3.*x2+y2)/npoints;
    g = (x2+3.*y2)/npoints;
    h = 2*xy/npoints;
    p = xx2y2/npoints;
    q = yx2y2/npoints;
    t = x2y22/npoints;
    g0 = (x2+y2)/npoints;
    g02 = g0*g0;
    a = -4.0;
    b = (f*g-t-h*h)/g02;
    c = (t*(f+g)-2.*(p*p+q*q))/(g02*g0);
    d = (t*(h*h-f*g)+2.*(p*p*g+q*q*f)-4.*p*q*h)/(g02*g02);
    xroot = 1.0;
    for (int i=0; i<5; i++) {
      ff = (((xroot+a)*xroot+b)*xroot+c)*xroot+d;
      fp = ((4.*xroot+3.*a)*xroot+2.*b)*xroot+c;
      xroot -= ff/fp;
    }
    g1 = xroot*g0;
    xnom = (g-g1)*(f-g1)-h*h;
    if (xnom == 0.){
      fprintf(stderr,"circleFit: xnom1 = %f\n",xnom);
      return -1;
    }


    yd = (q*(f-g1)-h*p)/xnom;
    xnom = f-g1;
    if (xnom == 0.){
      fprintf(stderr,"circleFit: xnom2 = %f\n",xnom);
      return -1;
    }

    xd = (p-h*yd )/xnom;

    radius2 = xd*xd+yd*yd+g1;
    *mXCenter = xd+xgravity;
    *mYCenter = yd+ygravity;
    for (int i=0; i<npoints; i++) {
      dx = mX[i]-(*mXCenter);
      dy = mY[i]-(*mYCenter);
      dradius2 = dx*dx+dy*dy;
      mVariance += dradius2+radius2-2.*sqrt(dradius2*radius2);
    }

    *mRadius  = (double) sqrt(radius2);
    double RadiusMes=(double) sqrt(radius2);

    /////linear fit for Pz
    double rr[500],zer[500];
    //Linear fit but exact
    //  rr[0]=0;

    rr[0]=0;
    zer[0]=1.;

    for(int i=1; i<npoints; i++){
      G4double aa = (sqrt((pow(mX[i-1]-(*mXCenter),2)+pow(mY[i-1]-(*mYCenter),2))*(pow(mX[i]-(*mXCenter),2)+pow(mY[i]-(*mYCenter),2))));
      G4double diff=acos(((mX[i-1]-(*mXCenter))*(mX[i]-(*mXCenter))+(mY[i-1]-(*mYCenter))*(mY[i]-(*mYCenter)))/aa);
      rr[i]=rr[i-1]+RadiusMes*diff;

      zer[i]=1.0;//must be corrected
    }

    double aa,bb;
    linearFitter(npoints,rr,mZ,zer,&aa,&bb);

    //  *Pz=aa*Pt
    *Pz_=aa*(RadiusMes*(0.299792458)*fabs(env_helm_field));
    *a_forz=aa;
    *b_forz=bb;
    *theta0_fory = atan2(mY[0]-(*mYCenter),mX[0]-(*mXCenter));


    //  *Pz=aa*(RadiusMes*0.299792458);
    return  mVariance;
  }

};

//_____________________________________________________________________________
inline G4String
TPCAnaManager::ClassName( void )
{
  static G4String s_name("TPCAnaManager");
  return s_name;
}

//_____________________________________________________________________________
inline TPCAnaManager&
TPCAnaManager::GetInstance( void )
{
  static TPCAnaManager s_instance;
  return s_instance;
}

#endif

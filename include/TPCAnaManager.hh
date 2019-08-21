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
static const G4int MaxPrimaryParticle = 10;

//const int MaxTrack = 1560*4;
//const int MaxTrack = 78*4;
const G4int MaxTrack = 54*20;

const G4int MaxNthLay = 40;
const G4int MaxNthPad = 250;


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
  Int_t gen;        // generator number
  Int_t mode;        // mode number

  Double_t mm_d; //missing-mass of d(pi,K) reaction
  Double_t mm_p;  //missing-mass for proton target kinematic
  Double_t mm; //missing-mass (should be same as mm_d)
  Double_t theta; //theta of scat K
  Double_t theta_scat; //theta of (pi, K) reaction
  Double_t theta_CM; //theta of (pi, K) reaction (CM flame)

  Int_t HitNum_K;
  //  int tpctrNum_K;

  Int_t HitNum_p;
  //  int tpctrNum_p;

  Int_t nhPrm;
  Int_t pidPrm[MaxPrimaryParticle];
  Double_t xPrm[MaxPrimaryParticle];
  Double_t yPrm[MaxPrimaryParticle];
  Double_t zPrm[MaxPrimaryParticle];
  Double_t pxPrm[MaxPrimaryParticle];
  Double_t pyPrm[MaxPrimaryParticle];
  Double_t pzPrm[MaxPrimaryParticle];
  Double_t ptPrm[MaxPrimaryParticle];
  Double_t mPrm[MaxPrimaryParticle];
  Double_t thetaPrm[MaxPrimaryParticle];

  // BH2
  Int_t nhBh2;
  Int_t tidBh2[MaxHits];
  Int_t pidBh2[MaxHits];
  Int_t didBh2[MaxHits];
  Int_t prtBh2[MaxHits];
  Int_t qBh2[MaxHits];
  Double_t massBh2[MaxHits];
  Double_t xBh2[MaxHits];
  Double_t yBh2[MaxHits];
  Double_t zBh2[MaxHits];
  Double_t pxBh2[MaxHits];
  Double_t pyBh2[MaxHits];
  Double_t pzBh2[MaxHits];
  Double_t ppBh2[MaxHits];
  Double_t deBh2[MaxHits];
  Double_t tBh2[MaxHits];
  Double_t vtxBh2[MaxHits];
  Double_t vtyBh2[MaxHits];
  Double_t vtzBh2[MaxHits];
  Double_t vtpxBh2[MaxHits];
  Double_t vtpyBh2[MaxHits];
  Double_t vtpzBh2[MaxHits];
  Double_t vtppBh2[MaxHits];
  Double_t lengthBh2[MaxHits];

  /* number of ntrks in TPC by shhwang*/
  Int_t ntrtpc;
  Double_t trpptpc[MaxHitsTPC];
  Double_t trpxtpc[MaxHitsTPC];
  Double_t trpytpc[MaxHitsTPC];
  Double_t trpztpc[MaxHitsTPC];
  Double_t trpttpc[MaxHitsTPC];

  Double_t trpptpcfit[MaxHitsTPC];
  Double_t trpxtpcfit[MaxHitsTPC];
  Double_t trpytpcfit[MaxHitsTPC];
  Double_t trpztpcfit[MaxHitsTPC];
  Double_t trpttpcfit[MaxHitsTPC];

  Int_t trqqtpc[MaxHitsTPC];
  Int_t trpidtpc[MaxHitsTPC];
  Int_t trparentidtpc[MaxHitsTPC];
  Int_t trparentid_pid_tpc[MaxHitsTPC];
  Double_t trpmtpc[MaxHitsTPC];
  Double_t trdetpc[MaxHitsTPC];
  Double_t trlentpc[MaxHitsTPC];
  Double_t trdedxtpc[MaxHitsTPC];
  Double_t trdedxtrtpc[MaxHitsTPC]; //trancated mean, but now just mean

  Int_t trlaytpc[MaxHitsTPC];

  Double_t vtpxtpc[MaxHitsTPC];
  Double_t vtpytpc[MaxHitsTPC];
  Double_t vtpztpc[MaxHitsTPC];
  Double_t vtpptpc[MaxHitsTPC];

  Double_t vtxtpc[MaxHitsTPC];
  Double_t vtytpc[MaxHitsTPC];
  Double_t vtztpc[MaxHitsTPC];

  Double_t vtxtpcfit[MaxHitsTPC];
  Double_t vtytpcfit[MaxHitsTPC];
  Double_t vtztpcfit[MaxHitsTPC];

  /////PAD multiplicity & ASAD multiplicy
  Int_t nthlay[MaxTrack];
  Int_t nthpad[MaxTrack];
  Int_t laypad[MaxTrack][MaxNthLay][MaxNthPad]; //[layer][pad number]


  ///////////////
  Int_t nttpc;                 // Number of Hit in Pads
  Int_t ntrk[MaxTrack];        // Number of Track

  Int_t ititpc[MaxTrack];      // Track ID
  Int_t idtpc[MaxTrack];       // Particle ID
  Double_t xtpc[MaxTrack];     // coordinates
  Double_t ytpc[MaxTrack];     // coordinates
  Double_t ztpc[MaxTrack];     // coordinates

  Double_t x0tpc[MaxTrack];    // coordinates
  Double_t y0tpc[MaxTrack];    // coordinates
  Double_t z0tpc[MaxTrack];    // coordinates
  Double_t resoX[MaxTrack];    // coordinates

  Double_t pxtpc[MaxTrack];    // momentum
  Double_t pytpc[MaxTrack];    // momentum
  Double_t pztpc[MaxTrack];    // momentum
  Double_t pptpc[MaxTrack];    // momentum
  Double_t masstpc[MaxTrack];    // mass

  Double_t betatpc[MaxTrack];    // beta

  Double_t edeptpc[MaxTrack];    // Energy deposit
  Double_t dedxtpc[MaxTrack];    // Energy deposit/dx
  Double_t slengthtpc[MaxTrack];    // Energy deposit/dx

  Int_t laytpc[MaxTrack];      // number of pad layer
  Int_t rowtpc[MaxTrack];      // number of pad raw
  Double_t toftpc[MaxTrack];   // tof
  Int_t parentID[MaxTrack];      // parent id
  Double_t cir_r[MaxTrack];   // fit radius
  Double_t cir_x[MaxTrack];   // fit center x
  Double_t cir_z[MaxTrack];   // fit center z
  Double_t cir_fit[MaxTrack];   // fit center fit
  Int_t vtx_flag[MaxTrack]; // flag, how to estimate vtx
  Double_t a_fory[MaxTrack]; // co-efficient a for linear track (y, theta)
  Double_t b_fory[MaxTrack]; // co-efficient b for linear track (y, theta)

  // TARGET
  Int_t nhTgt;
  Int_t tidTgt[MaxHits];
  Int_t pidTgt[MaxHits];
  Int_t prtTgt[MaxHits];
  Double_t xTgt[MaxHits];
  Double_t yTgt[MaxHits];
  Double_t zTgt[MaxHits];
  Double_t vtxTgt[MaxHits];
  Double_t vtyTgt[MaxHits];
  Double_t vtzTgt[MaxHits];
  // HTOF
  Int_t nhHtof;
  Int_t tidHtof[MaxHits];
  Int_t pidHtof[MaxHits];
  Int_t didHtof[MaxHits];
  Int_t prtHtof[MaxHits];
  Int_t qHtof[MaxHits];
  Double_t massHtof[MaxHits];
  Double_t xHtof[MaxHits];
  Double_t yHtof[MaxHits];
  Double_t zHtof[MaxHits];
  Double_t pxHtof[MaxHits];
  Double_t pyHtof[MaxHits];
  Double_t pzHtof[MaxHits];
  Double_t ppHtof[MaxHits];
  Double_t deHtof[MaxHits];
  Double_t tHtof[MaxHits];
  Double_t vtxHtof[MaxHits];
  Double_t vtyHtof[MaxHits];
  Double_t vtzHtof[MaxHits];
  Double_t vtpxHtof[MaxHits];
  Double_t vtpyHtof[MaxHits];
  Double_t vtpzHtof[MaxHits];
  Double_t vtppHtof[MaxHits];
  Double_t lengthHtof[MaxHits];
  // SDC
  Int_t nhSdc;
  Int_t tidSdc[MaxHits];
  Int_t pidSdc[MaxHits];
  Int_t didSdc[MaxHits];
  Int_t prtSdc[MaxHits];
  Int_t qSdc[MaxHits];
  Double_t massSdc[MaxHits];
  Double_t xSdc[MaxHits];
  Double_t ySdc[MaxHits];
  Double_t zSdc[MaxHits];
  Double_t pxSdc[MaxHits];
  Double_t pySdc[MaxHits];
  Double_t pzSdc[MaxHits];
  Double_t ppSdc[MaxHits];
  Double_t deSdc[MaxHits];
  Double_t tSdc[MaxHits];
  Double_t vtxSdc[MaxHits];
  Double_t vtySdc[MaxHits];
  Double_t vtzSdc[MaxHits];
  Double_t vtpxSdc[MaxHits];
  Double_t vtpySdc[MaxHits];
  Double_t vtpzSdc[MaxHits];
  Double_t vtppSdc[MaxHits];
  Double_t lengthSdc[MaxHits];
  // SCH
  Int_t nhSch;
  Int_t tidSch[MaxHits];
  Int_t pidSch[MaxHits];
  Int_t didSch[MaxHits];
  Int_t prtSch[MaxHits];
  Int_t qSch[MaxHits];
  Double_t massSch[MaxHits];
  Double_t xSch[MaxHits];
  Double_t ySch[MaxHits];
  Double_t zSch[MaxHits];
  Double_t pxSch[MaxHits];
  Double_t pySch[MaxHits];
  Double_t pzSch[MaxHits];
  Double_t ppSch[MaxHits];
  Double_t deSch[MaxHits];
  Double_t tSch[MaxHits];
  Double_t vtxSch[MaxHits];
  Double_t vtySch[MaxHits];
  Double_t vtzSch[MaxHits];
  Double_t vtpxSch[MaxHits];
  Double_t vtpySch[MaxHits];
  Double_t vtpzSch[MaxHits];
  Double_t vtppSch[MaxHits];
  Double_t lengthSch[MaxHits];
  // FTOF
  Int_t nhFtof;
  Int_t tidFtof[MaxHits];
  Int_t pidFtof[MaxHits];
  Int_t didFtof[MaxHits];
  Int_t prtFtof[MaxHits];
  Int_t qFtof[MaxHits];
  Double_t massFtof[MaxHits];
  Double_t xFtof[MaxHits];
  Double_t yFtof[MaxHits];
  Double_t zFtof[MaxHits];
  Double_t pxFtof[MaxHits];
  Double_t pyFtof[MaxHits];
  Double_t pzFtof[MaxHits];
  Double_t ppFtof[MaxHits];
  Double_t deFtof[MaxHits];
  Double_t tFtof[MaxHits];
  Double_t vtxFtof[MaxHits];
  Double_t vtyFtof[MaxHits];
  Double_t vtzFtof[MaxHits];
  Double_t vtpxFtof[MaxHits];
  Double_t vtpyFtof[MaxHits];
  Double_t vtpzFtof[MaxHits];
  Double_t vtppFtof[MaxHits];
  Double_t lengthFtof[MaxHits];
  // LAC
  Int_t nhLac;
  Int_t tidLac[MaxHits];
  Int_t pidLac[MaxHits];
  Int_t didLac[MaxHits];
  Int_t prtLac[MaxHits];
  Int_t qLac[MaxHits];
  Double_t massLac[MaxHits];
  Double_t xLac[MaxHits];
  Double_t yLac[MaxHits];
  Double_t zLac[MaxHits];
  Double_t pxLac[MaxHits];
  Double_t pyLac[MaxHits];
  Double_t pzLac[MaxHits];
  Double_t ppLac[MaxHits];
  Double_t deLac[MaxHits];
  Double_t tLac[MaxHits];
  Double_t vtxLac[MaxHits];
  Double_t vtyLac[MaxHits];
  Double_t vtzLac[MaxHits];
  Double_t vtpxLac[MaxHits];
  Double_t vtpyLac[MaxHits];
  Double_t vtpzLac[MaxHits];
  Double_t vtppLac[MaxHits];
  Double_t lengthLac[MaxHits];
  // WC
  Int_t nhWc;
  Int_t tidWc[MaxHits];
  Int_t pidWc[MaxHits];
  Int_t didWc[MaxHits];
  Int_t prtWc[MaxHits];
  Int_t qWc[MaxHits];
  Double_t massWc[MaxHits];
  Double_t xWc[MaxHits];
  Double_t yWc[MaxHits];
  Double_t zWc[MaxHits];
  Double_t pxWc[MaxHits];
  Double_t pyWc[MaxHits];
  Double_t pzWc[MaxHits];
  Double_t ppWc[MaxHits];
  Double_t deWc[MaxHits];
  Double_t tWc[MaxHits];
  Double_t vtxWc[MaxHits];
  Double_t vtyWc[MaxHits];
  Double_t vtzWc[MaxHits];
  Double_t vtpxWc[MaxHits];
  Double_t vtpyWc[MaxHits];
  Double_t vtpzWc[MaxHits];
  Double_t vtppWc[MaxHits];
  Double_t lengthWc[MaxHits];
  // VP
  Int_t nhVp;
  Int_t tidVp[MaxHits];
  Int_t pidVp[MaxHits];
  Int_t didVp[MaxHits];
  Int_t prtVp[MaxHits];
  Int_t qVp[MaxHits];
  Double_t massVp[MaxHits];
  Double_t xVp[MaxHits];
  Double_t yVp[MaxHits];
  Double_t zVp[MaxHits];
  Double_t pxVp[MaxHits];
  Double_t pyVp[MaxHits];
  Double_t pzVp[MaxHits];
  Double_t ppVp[MaxHits];
  Double_t deVp[MaxHits];
  Double_t tVp[MaxHits];
  Double_t vtxVp[MaxHits];
  Double_t vtyVp[MaxHits];
  Double_t vtzVp[MaxHits];
  Double_t vtpxVp[MaxHits];
  Double_t vtpyVp[MaxHits];
  Double_t vtpzVp[MaxHits];
  Double_t vtppVp[MaxHits];
  Double_t lengthVp[MaxHits];
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
  CounterData counterData[MaxTrack];
  TPCData tpcData[MAXtpctrNum];

  PrimaryInfo primaryInfo;
  int HitNum;
  int tpctrNum;
  int HitNum_K;
  int HitNum_p;

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
  G4int env_on_off_helm;
  G4double env_helm_field;
  G4int env_pad_config;
  G4int env_Experiment_num;

  G4double angle[40];
  G4double seg_angle[40];
  G4double seg_width[40];
  G4int numpads[40];
  G4double pad_in[40];
  G4double pad_out[40];
  G4double tpc_rad;

public:
  void BeginOfRunAction( G4int runnum );
  void EndOfRunAction( void );
  void BeginOfEventAction( void );
  int  EndOfEventAction( void );
  void SetTPCData( G4int tpctr, G4int tpcpid, G4int tpcparentid,
		   G4int tpcparentid_pid, G4double tpcpx, G4double tpcpy,
		   G4double tpcpz,G4double tpcpp,  G4int tpcqq, G4double tpcpm,
		   G4double tpcde, G4double tpclen, G4int tpclay,
		   G4double vtxpxtpc2,G4double vtxpytpc2,G4double vtxpztpc2,
		   G4double vtxxtpc2,G4double vtxytpc2,G4double vtxztpc2,
		   G4double vtxenetpc2 );
  void SetBH2Data( const VHitInfo* hit );
  void SetCounterData( G4int ntrk, G4double time, G4ThreeVector pos,
		       G4ThreeVector mom, G4int track, G4int particle,
		       G4int iLay, G4int iRow, G4double beta, G4double edep,
		       G4int parentid, G4double tlength, G4double slength );
  void SetHTOFData( const VHitInfo* hit );
  void SetFTOFData( const VHitInfo* hit );
  void SetLACData( const VHitInfo* hit );
  void SetSCHData( const VHitInfo* hit );
  void SetSDCData( const VHitInfo* hit );
  void SetVPData( const VHitInfo* hit );
  void SetWCData( const VHitInfo* hit );
  void SetGeneratorID(G4int gen);
  void SetModeID(G4int mode);
  void SetNumberOfPrimaryParticle( G4int n );
  void SetTargetData( const VHitInfo* hit );
  void SetPrimaryBeam( const G4ThreeVector& p );
  void SetPrimaryBeam( G4double px, G4double py, G4double pz );
  void SetPrimaryParticle( G4double px, G4double py, G4double pz );
  void SetPrimaryParticle( G4int id, const G4ThreeVector& p, G4double mass,
			   G4int pid=-9999 );
  void SetPrimaryParticle( G4int id, G4double px, G4double py, G4double pz,
			   G4double mass, G4int pid=-9999 );
  void SetPrimaryVertex( G4int id, const G4ThreeVector& x );
  void SetPrimaryVertex( G4int id, G4double x, G4double y, G4double z );
  void SetPrimaryInfo( G4double mm_d, G4double mm_p, G4double theta,
		       G4double theta_scat, G4double theta_cm );
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

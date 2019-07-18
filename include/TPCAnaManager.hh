// -*- C++ -*-

#ifndef TPC_ANA_MANAGER_HH
#define TPC_ANA_MANAGER_HH

#include <G4ThreeVector.hh>

#include "TPCAnaRoot.hh"

struct Track;

void initTrack(Track* aTrack);
void initTrack_ku(Track* aTrack);
int setInitialPara(Track* aTrack, double* initPara);
int setVirtualPlane(Track* aTrack);
void minuitInit(double printLevel);

static const int MAXtpctrNum=30;
static const int MAXtpctrhitNum=500;

class TH1F;
class TH2F;
class TFile;
class TTree;

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
class TPCAnaManager
{
public:
  static TPCAnaManager& GetInstance( void );
  ~TPCAnaManager( void );

private:
  TPCAnaManager( void );
  TPCAnaManager( const TPCAnaManager& );
  TPCAnaManager& operator=( const TPCAnaManager& );

private:
  TargetData targetData[MaxTrack];
  CounterData counterData[MaxTrack];
  TPCData tpcData[MAXtpctrNum];
  ScintData scintData[MaxTrig];
  ACData acData[MaxTrig];
  NBARData nbarData[MaxTrig];

  DCData dcData[MaxTrig];
  SCHData schData[MaxTrig];
  FTOFData ftofData[MaxTrig];

  TPCAnaRoot    anaRoot;
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

  void SetTargetData( G4int nhits, G4ThreeVector pos, G4ThreeVector mom,
		      G4int track, G4int particle,
		      G4int parentid, G4ThreeVector vtxpos,
		      G4ThreeVector vtxmom, G4double vtxene );

  //  void SetCounterData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
  //		      G4int track, G4int particle, G4int ilay, G4int iRaw);
  void SetScintData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
		    G4int track, G4int particle, G4int detector,G4double mass, G4int qq,G4int parentid, G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength);
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
inline TPCAnaManager&
TPCAnaManager::GetInstance( void )
{
  static TPCAnaManager s_instance;
  return s_instance;
}

#endif

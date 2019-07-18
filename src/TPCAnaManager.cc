// -*- C++ -*-

#include "TPCAnaManager.hh"

#include <algorithm>

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include "ConfMan.hh"
#include "ResHypTPC.hh"
#include "RungeKuttaTracker.hh"
#include "minuit2.hh"
#include "switch.h"
#include "track.hh"

#define	MAX_DIM_FOR_RKF 3

namespace
{
  using CLHEP::mm;
  const ConfMan& gConf = ConfMan::GetInstance();
}

//void chi2(int *npar, double *grad, double *fval,
//	  double *xval, int *iflag, void (*futil)() );

//MinuitFCN fcn =  chi2;
//MinuitFCN *fcn
//extern int chi2prb_(float* chi2,int* ndf,float* prb);
// const int MaxFCNCall = 300;
// const double EPS = 1.;
extern double RKChi2[MAX_ITERATION];
extern double RKPara[MAX_ITERATION][NUM_PARA_RK];
//void rungeKuttaFehlberg(int n, double x[],
//			  void diff(int, double *, double, double *),
//			  double *pt, double *ph, double timeEnd,
//			  int kcnt, double hmin, double hmax, double TOL);

//#include <string>

/*ResHypTPC::init()
{
    y_resolution = 0.5;
  sigma_amp = 13.3/53.9;
  neff = 26;//MIP
  neff_nmpv_correction;// = 0.7*1.2;
  diffuse_GEM = 0.1;
  }*/



///////////RK variable
#define SQ(a) ((a)*(a))
enum ROTATION_MODE {G2L,L2G};

/////////////////////////end RK par

//_____________________________________________________________________________
TPCAnaManager::TPCAnaManager( void )
{
}

//_____________________________________________________________________________
TPCAnaManager::~TPCAnaManager( void )
{
}

//_____________________________________________________________________________
void
TPCAnaManager::BeginOfRunAction( int runnum )
{
  target_pos_z=-143.;
  truncated_mean_cut = gConf.Get<G4double>("TruncatedMeanCut");
  env_Experiment_num = gConf.Get<G4int>("Experiment");
  //out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  pad_length_in = gConf.Get<G4double>("PadLengthIn");
  pad_length_out = gConf.Get<G4double>("PadLengthOut");
  pad_gap = gConf.Get<G4double>("PadGap");

  ////pad configure
  env_pad_config = gConf.Get<G4int>("PadConfigure");
  pad_in_num = gConf.Get<G4int>("PadNumIn");
  pad_out_num = gConf.Get<G4int>("PadNumOut");
  pad_in_width = gConf.Get<G4double>("PadWidthOut");
  pad_out_width = gConf.Get<G4double>("PadWidthOut");

  env_on_off_helm = gConf.Get<G4int>("ShsFieldMap");

  if( env_on_off_helm == 0 ){
    env_helm_field = gConf.Get<G4int>("ShsField");
  }else{
    G4cout << "Env of the Helmholt_fieldmap is wrong" << G4endl;
    exit(-1);
  }

  for(G4int i=0.;i<40;i++){
    angle[i]=0;
    seg_angle[i]=0;
    seg_width[i]=0;
    numpads[i]=0;

    pad_in[i]=0;
    pad_out[i]=0;
  }
    tpc_rad=250;
    G4double cen_diff=fabs(target_pos_z);

    if( env_pad_config ==1 ){
      for(G4int i=0;i<pad_in_num+pad_out_num;i++){
	if(i<pad_in_num){
	  pad_in[i]=10.+(pad_length_in+pad_gap)*i;
	  pad_out[i]=10.+(pad_length_in+pad_gap)*i+pad_length_in;
	  angle[i]=360.;
	}else {
	  pad_in[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num);
	  pad_out[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num) + pad_length_out;
	  angle[i]=180.-acos((pow(pad_out[i],2)+pow(cen_diff,2)-pow(tpc_rad,2))/(2*pad_out[i]*cen_diff))*180./acos(-1.);
	}
      //      G4cout<<angle[i]<<G4endl;
      //      G4cout<<pad_in[i]<<G4endl;
      }


    }else if( env_pad_config ==2 ){
      for(G4int i=0;i<pad_in_num+pad_out_num;i++){
	if(i<pad_in_num){
	  pad_in[i]=10.+(pad_length_in+pad_gap)*i;
	  pad_out[i]=10.+(pad_length_in+pad_gap)*i+pad_length_in;
	  angle[i]=360.;
	  if(i==0){
	    numpads[i]=48.;
	  }else if(i<pad_in_num){
	    numpads[i]=24.*2.*(i+1.)/2.;
	  }
	}else {
	  pad_in[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num);
	  pad_out[i]=10.+(pad_length_in+pad_gap)*pad_in_num+(pad_length_out+pad_gap)*(i-pad_in_num) + pad_length_out;
	}
      }
      angle[10]=180.-155.35;
      angle[11]=180.-144.8;
      angle[12]=180.-138.;
      angle[13]=180.-116.73;
      angle[14]=180.-106.;
      angle[15]=180.-98.77;
      angle[16]=180.-94.29;
      angle[17]=180.-89.8;
      angle[18]=180.-87.18;
      angle[19]=180.-84.16;
      angle[20]=180.-81.48;
      angle[21]=180.-73.39;
      angle[22]=180.-65.51011;
      angle[23]=180.-60.19;
      angle[24]=180.-56.35239;
      angle[25]=180.-52.85;
      angle[26]=180.-50.14;
      angle[27]=180.-47.17;
      angle[28]=180.-41.24;
      angle[29]=180.-29.;
      angle[30]=180.-23.23;
      angle[31]=180.-18.69;

      numpads[10]=208.;
      numpads[11]=218.;
      numpads[12]=230.;
      numpads[13]=214.;
      numpads[14]=212.;
      numpads[15]=214.;
      numpads[16]=220.;
      numpads[17]=224.;
      numpads[18]=232.;
      numpads[19]=238.;
      numpads[20]=244.;
      numpads[21]=232.;
      numpads[22]=218.;
      numpads[23]=210.;
      numpads[24]=206.;
      numpads[25]=202.;
      numpads[26]=200.;
      numpads[27]=196.;
      numpads[28]=178.;
      numpads[29]=130.;
      numpads[30]=108.;
      numpads[31]=90.;
      G4int all_channels=0;
      G4int all_channels2=0;
      G4int num_pad_check=0;


      for(G4int i=0;i<pad_in_num+pad_out_num;i++){
	if(i<pad_in_num){
	  seg_angle[i]=360./double(numpads[i]);
	  seg_width[i]=pad_in[i]*(angle[i])*CLHEP::pi/180./numpads[i];

	  num_pad_check=angle[i]/seg_angle[i];
	}else if(i>=pad_in_num){
	  seg_angle[i]=(180.-angle[i])*2/double(numpads[i]);
	  seg_width[i]=pad_in[i]*(180-angle[i])*2.*acos(-1.)/180./numpads[i];
	  num_pad_check=(180.-angle[i])*2/seg_angle[i];
	}

	G4cout<<i<<" degree :"<<seg_angle[i]<<G4endl;
	G4cout<<i<<" width :"<<seg_angle[i]*acos(-1.)/180.*pad_in[i]<<G4endl;

	all_channels=all_channels+numpads[i];
	all_channels2=all_channels2+num_pad_check;
      }
      G4cout<<"------------------------"<<G4endl;
      G4cout<<"Total pads:"<<all_channels<<G4endl;
      G4cout<<"Total pads(check):"<<all_channels<<G4endl;
      G4cout<<"------------------------"<<G4endl;
    }


  anaRoot.BeginOfRunAction(runnum);

}

void TPCAnaManager::EndOfRunAction()
{
  anaRoot.EndOfRunAction();
}


void TPCAnaManager::BeginOfEventAction()
{
  HitNum=0;
  tpctrNum=0;
  HitNumAC=0;
  HitNumNBAR=0;
  HitNumDC=0;
  HitNumSCH=0;
  HitNumFTOF=0;
  HitNumScint=0;
  HitNumTarget=0;

  //for K+
  HitNum_K=0;
  //  tpctrNum_K=0;
  HitNumAC_K=0;
  HitNumNBAR_K=0;
  HitNumDC_K=0;
  HitNumSCH_K=0;
  HitNumFTOF_K=0;
  HitNumScint_K=0;
  HitNumTarget_K=0;

  //for proton
  HitNum_p=0;
  //  tpctrNum_K=0;
  HitNumAC_p=0;
  HitNumNBAR_p=0;
  HitNumDC_p=0;
  HitNumSCH_p=0;
  HitNumFTOF_p=0;
  HitNumScint_p=0;
  HitNumTarget_p=0;

  anaRoot.BeginOfEventAction();
}


int TPCAnaManager::EndOfEventAction()
{
  // Fill kaon minus energy distribution

  anaRoot.FillBeam(primaryBeam.pg[0], primaryBeam.pg[1], primaryBeam.pg[2]);

  anaRoot.FillGenMode(primaryBeam.gen, primaryBeam.mode);

  // anaRoot.FillNumOfK(HitNum_K, HitNumAC_K, HitNumNBAR_K,
  // 		     HitNumDC_K, HitNumSCH_K, HitNumFTOF_K,
  // 		     HitNumScint_K, HitNumTarget_K);
  //  std::cout<<"HitNumDC_K="<<HitNumDC_K<<", HitNumFTOF_K="<<HitNumFTOF_K<<std::endl;
  // anaRoot.FillNumOfp(HitNum_p, HitNumAC_p, HitNumNBAR_p,
  // 		     HitNumDC_p, HitNumSCH_p, HitNumFTOF_p,
  // 		     HitNumScint_p, HitNumTarget_p);
  // std::cout<<"HitNumDC_p="<<HitNumDC_p<<", HitNumFTOF_p="<<HitNumFTOF_p<<std::endl;

  // Fill Primary particle distribution
  for(G4int id=0; id<primaryParticle.NumOfParticle;id++){
    anaRoot.FillPrimaryParticle(id, &primaryParticle.x0[id][0],
    				&primaryParticle.p0[id][0], primaryParticle.pid0[id]);
  }

  //Fill Primary Infomation for E27
  if(env_Experiment_num ==27 ||env_Experiment_num ==45)
    anaRoot.FillPrimaryInfo(primaryInfo.mm_d, primaryInfo.mm_p, primaryInfo.theta, primaryInfo.theta_scat, primaryInfo.theta_CM);
  else
    anaRoot.FillPrimaryInfo(0., 0., 0., 0., 0.);





  //  if( HitNumScint > 0 ){
  ///this "if" is trigger
  //  if( HitNum > 0 && HitNumFTOF>0){
  if( HitNum > 0){
    //  std::cout<<"here"<<std::endl;
    //  if( HitNum > 0 && HitNumFTOF>0 && HitNumScint>0){

    // Fill TPC hits condition
    // G4int detID = 0;

    G4int c[MAX_TRACK] = {};

    for(G4int i=0;i<MAX_TRACK;i++){
      mean[i]=0.;
      trmean[i]=0.;
    }

    G4double vtxxfit[MAX_TRACK];//read fit parameters
    G4double vtxyfit[MAX_TRACK];//read fit parameters
    G4double vtxzfit[MAX_TRACK];//read fit parameters

    G4double vtxpxfit[MAX_TRACK];//read fit parameters
    // G4double vtxpyfit[MAX_TRACK];//read fit parameters
    G4double vtxpzfit[MAX_TRACK];//read fit parameters

    for(G4int i=0;i<MAX_TRACK;i++){
      vtxxfit[i]=-9999.9999;
      vtxyfit[i]=-9999.9999;
      vtxzfit[i]=-9999.9999;
      vtxpxfit[i]=-9999.9999;
      // vtxpyfit[i]=-9999.9999;
      vtxpzfit[i]=-9999.9999;
      Pz[i]=-9999.9999;
    }

    G4double x[MAX_TRACK][MAXtpctrhitNum]={{-9999.9999},{-9999.9999}};;
    G4double z[MAX_TRACK][MAXtpctrhitNum]={{-9999.9999},{-9999.9999}};
    G4double y[MAX_TRACK][MAXtpctrhitNum]={{-9999.9999},{-9999.9999}};
    G4double ede[MAX_TRACK][MAXtpctrhitNum]={{0.},{0.}};

    ////// shhwang position read
    ///shhwang code


    if(tpctrNum>9){
      G4cout<<"Error--> over the number of tracks in the TPC:"<<tpctrNum<<G4endl;
    }

    G4int sh_paID[MAX_TRACK] = {};
    for( G4int i=0; i<HitNum; i++){
      G4int ii=counterData[i].ntrk;
      x[ii][c[ii]]=counterData[i].pos[0];
      z[ii][c[ii]]=counterData[i].pos[2];
      y[ii][c[ii]]=counterData[i].pos[1];
      sh_paID[ii]=counterData[i].parentID;
      ede[ii][c[ii]]=counterData[i].dedx;
      c[ii]=c[ii]+1;
    }

    G4double test[MAX_TRACK]={-1};
    G4double cx[MAX_TRACK]={-9999.9999};
    G4double cz[MAX_TRACK]={-9999.9999};
    G4double cir_x[MAX_TRACK]={-9999.9999};
    G4double cir_z[MAX_TRACK]={-9999.9999};
    G4double rad[MAX_TRACK]={-9999.9999};
    G4double a_fory[MAX_TRACK]={-9999.9999};
    G4double b_fory[MAX_TRACK]={-9999.9999};
    G4double theta0_fory[MAX_TRACK]={-9999.9999};
    G4int vtx_flag[MAX_TRACK]={-1};

    for(G4int i=0;i<MAX_TRACK;i++){
      cir_r[i]=-9999.9999;
      cir_x[i]=-9999.9999;
      cir_z[i]=-9999.9999;
      mean[i]=-9999.9999;
    }



    for( G4int kk=0; kk<tpctrNum; kk++){
      if(c[kk]>3.){
	//	G4cout<<"start circle fit"<<G4endl;
	test[kk]=circleFit(x[kk],z[kk],y[kk],c[kk],&cx[kk],&cz[kk],&rad[kk],&Pz[kk],
			   &a_fory[kk], &b_fory[kk], &theta0_fory[kk]);
	if(test[kk]!=-1.){
	  cir_r[kk]=rad[kk];
	  cir_x[kk]=cx[kk];
	  cir_z[kk]=cz[kk];
	}

      }
    }
    G4double mom_theta[MAX_TRACK]={0.};

    // calcute vtx with production points
    for(G4int i=0;i<MAX_TRACK;i++){
      G4double rho1 = rad[i];
      G4double cx1 = cx[i];
      G4double cz1 = cz[i];
      G4double cx2 = primaryParticle.x0[0][0];
      G4double cz2 = primaryParticle.x0[0][2];
      G4double theta12=atan2(cz2-cz1, cx2-cx1);
      G4double ca1=a_fory[i];
      G4double cb1=b_fory[i];
      G4double ct01=theta0_fory[i];


      G4double cent_dist=sqrt(pow(cx1-cx2,2)+pow(cz1-cz2,2));

      vtxxfit[i]=cos(theta12)*cent_dist+cx1;
      vtxzfit[i]=sin(theta12)*cent_dist+cz1;
      vtxyfit[i]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

      mom_theta[i]=atan2(vtxzfit[i]-cz[i],vtxxfit[i]-cx[i])-acos(-1.)/2;

      vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
      vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);

      vtx_flag[i]=1;
    }

    ////think about parent ID
    ////--> find the track with same parent ID
    //// sh_

    for(G4int i=0;i<MAX_TRACK;i++){
      for(G4int j=i;j<MAX_TRACK;j++){
	if(i!=j && (test[i]>0 && test[j]>0) ){
	  if(sh_paID[i]==sh_paID[j] && sh_paID[i]>0. && sh_paID[j]>0.){
	    //	    G4cout<<"vtx1"<<env_helm_field<<G4endl;
	    G4double rho1=rad[i];
	    G4double rho2=rad[j];

	    G4double cx1=cx[i];
	    G4double cz1=cz[i];
	    G4double ca1=a_fory[i];
	    G4double cb1=b_fory[i];
	    G4double ct01=theta0_fory[i];

	    G4double cx2=cx[j];
	    G4double cz2=cz[j];
	    G4double ca2=a_fory[j];
	    G4double cb2=b_fory[j];
	    G4double ct02=theta0_fory[j];


	    G4double cent_dist=sqrt(pow(cx1-cx2,2)+pow(cz1-cz2,2));

	    double point1[3]={0};
	    double point2[3]={0};
	    G4int k;

	    if((cent_dist-(rho1+rho2))>0.){


	      G4double theta12=atan2(cz2-cz1,cx2-cx1);
	      G4double centr=rho1+(cent_dist-(rho1+rho2))/2;

	      point1[0]=cos(theta12)*centr+cx1;
	      point1[1]=sin(theta12)*centr+cz1;
	      point1[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

	      G4double theta21=atan2(cz1-cz2,cx1-cx2);
	      G4double centr1=rho2+(cent_dist-(rho1+rho2))/2;
	      point2[0]=cos(theta21)*centr1+cx2;
	      point2[1]=sin(theta21)*centr1+cz2;
	      point2[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta21-ct02)+cb2;

	      vtxxfit[i]=point1[0];
	      vtxzfit[i]=point1[1];
	      vtxyfit[i]=(point1[2]+point2[2])/2.;
	      vtxxfit[j]=point2[0];
	      vtxzfit[j]=point2[1];
	      vtxyfit[j]=(point1[2]+point2[2])/2.;
	      vtx_flag[i]=2;
	      vtx_flag[j]=2;
	      //
	    }else  if((cent_dist+fmin(rho1,rho2))<fmax(rho1,rho2)){
	      if(rho1>=rho2){ //rho1>rho2
		G4double theta12=atan2(cz2-cz1,cx2-cx1);
		G4double centr=rho1-(rho1-cent_dist-rho2)/2; //rho1>rho2
		point1[0]=cos(theta12)*centr+cx1;
		point1[1]=sin(theta12)*centr+cz1;
		point1[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

		G4double theta21=atan2(cz2-cz1,cx2-cx1);
		G4double centr1=rho2+(rho1-cent_dist-rho2)/2.; //rho1>rho2
		point2[0]=cos(theta21)*centr1+cx2;
		point2[1]=sin(theta21)*centr1+cz2;
		point2[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta21-ct02)+cb2;
		//		G4cout<<"test1"<<G4endl;

	      }else if(rho2>rho1){ //rho1<rho2
		G4double theta12=atan2(cz1-cz2,cx1-cx2);
		G4double centr=rho2-(rho2-cent_dist-rho1)/2; //rho1<rho2
		point1[0]=cos(theta12)*centr+cx2;
		point1[1]=sin(theta12)*centr+cz2;
		point1[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta12-ct02)+cb2;

		G4double theta21=atan2(cz1-cz2,cx1-cx2);
		G4double centr1=rho1+(rho2-cent_dist-rho1)/2; //rho1<rho2
		point2[0]=cos(theta21)*centr1+cx1;
		point2[1]=sin(theta21)*centr1+cz1;
		point2[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta21-ct01)+cb1;
	      }

	      vtxxfit[i]=point1[0];
	      vtxzfit[i]=point1[1];
	      vtxyfit[i]=(point1[2]+point2[2])/2.;
	      // vtxxfit[j]=point1[0];
	      // vtxzfit[j]=point1[1];
	      // vtxyfit[j]=point1[2];
	      vtxxfit[j]=point2[0];
	      vtxzfit[j]=point2[1];
	      vtxyfit[j]=(point1[2]+point2[2])/2.;

	      vtx_flag[i]=3;
	      vtx_flag[j]=3;
	    } else {

	      //k = CircleIntersect(cx1,cz1,rho1,cx2,cz2,rho2,point1,point2);
	      k = CircleIntersect(cx1,cz1,rho1,cx2,cz2,rho2,ca1,cb1,ct01,tpcData[i].tpcqq,ca2,cb2,ct02,tpcData[j].tpcqq,point1,point2);
	      if(k == 0) {
		G4cout << "no solution" << G4endl;
	      }else if(k>0){


		G4double dist1=sqrt(pow(point1[0]-tpcData[i].tpcvtxx,2)+pow(point1[1]-tpcData[i].tpcvtxz,2));
		G4double dist2=sqrt(pow(point2[0]-tpcData[i].tpcvtxx,2)+pow(point2[1]-tpcData[i].tpcvtxz,2));

		if(dist1<=dist2){//point1 is correct
		  vtxxfit[i]=point1[0];
		  vtxzfit[i]=point1[1];
		  vtxyfit[i]=point1[2];

		  vtxxfit[j]=point1[0];
		  vtxzfit[j]=point1[1];
		  vtxyfit[j]=point1[2];
		}else if(dist1>dist2){//point1 is correct
		  vtxxfit[i]=point2[0];
		  vtxzfit[i]=point2[1];
		  vtxyfit[i]=point2[2];

		  vtxxfit[j]=point2[0];
		  vtxzfit[j]=point2[1];
		  vtxyfit[j]=point2[2];
		}
		vtx_flag[i]=4;
		vtx_flag[j]=4;
	      }

	      mom_theta[i]=atan2(vtxzfit[i]-cz[i],vtxxfit[i]-cx[i])-acos(-1.)/2;
	      mom_theta[j]=atan2(vtxzfit[j]-cz[j],vtxxfit[j]-cx[j])-acos(-1.)/2;


	      // std::cout<<"x01="<<x[i][0]<<", x2="<<x[j][0]<<std::endl;
	      // std::cout<<"y01="<<y[i][0]<<", y2="<<y[j][0]<<std::endl;
	      // std::cout<<"z01="<<z[i][0]<<", z2="<<z[j][0]<<std::endl;

	      // std::cout<<"vtx fit="<<vtxxfit[i]<<", true vtx="<<tpcData[i].tpcvtxx<<std::endl;
	      // std::cout<<"vty fit="<<vtxyfit[i]<<", true vty="<<tpcData[i].tpcvtxy<<std::endl;
	      // std::cout<<"vtz fit="<<vtxzfit[i]<<", true vtz="<<tpcData[i].tpcvtxz<<std::endl;
	      //getchar();


	      vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      vtxpxfit[j]=cos(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
	      vtxpzfit[j]=sin(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
	      //	    G4cout<<"bfield:"<<env_helm_field<<G4endl;


	    }

	    ///from vertex particle, but it need more than 2
	  }else if(sh_paID[i]==sh_paID[j] && sh_paID[i]==0. && sh_paID[j]==0.){
	    //	    G4cout<<"vtx2"<<env_helm_field<<G4endl;
	    G4double rho1=rad[i];
	    G4double rho2=rad[j];

	    G4double cx1=cx[i];
	    G4double cz1=cz[i];
	    G4double ca1=a_fory[i];
	    G4double cb1=b_fory[i];
	    G4double ct01=theta0_fory[i];

	    G4double cx2=cx[j];
	    G4double cz2=cz[j];
	    G4double ca2=a_fory[j];
	    G4double cb2=b_fory[j];
	    G4double ct02=theta0_fory[j];

	    G4double cent_dist=sqrt(pow(cx1-cx2,2)+pow(cz1-cz2,2));

	    double point1[3]={0};
	    double point2[3]={0};
	    G4int k;

	    if((cent_dist-(rho1+rho2))>0.){


	      G4double theta12=atan2(cz2-cz1,cx2-cx1);
	      G4double centr=rho1+(cent_dist-(rho1+rho2))/2;

	      point1[0]=cos(theta12)*centr+cx1;
	      point1[1]=sin(theta12)*centr+cz1;
	      point1[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

	      G4double theta21=atan2(cz1-cz2,cx1-cx2);
	      G4double centr1=rho2+(cent_dist-(rho1+rho2))/2;
	      point2[0]=cos(theta21)*centr1+cx2;
	      point2[1]=sin(theta21)*centr1+cz2;
	      point2[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta21-ct02)+cb2;

	      vtxxfit[i]=point1[0];
	      vtxzfit[i]=point1[1];
	      vtxyfit[i]=(point1[2]+point2[2])/2.;
	      vtxxfit[j]=point1[0];
	      vtxzfit[j]=point1[1];
	      vtxyfit[j]=(point1[2]+point2[2])/2.;

	      vtx_flag[i]=5;
	      vtx_flag[j]=5;
	    }else  if((cent_dist+fmin(rho1,rho2))<fmax(rho1,rho2)){

	      if(rho1>=rho2){ //rho1>rho2
		G4double theta12=atan2(cz2-cz1,cx2-cx1);
		G4double centr=rho1-(rho1-cent_dist-rho2)/2; //rho1>rho2
		point1[0]=cos(theta12)*centr+cx1;
		point1[1]=sin(theta12)*centr+cz1;
		point1[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

		G4double theta21=atan2(cz2-cz1,cx2-cx1);
		G4double centr1=rho2+(rho1-cent_dist-rho2)/2.; //rho1>rho2
		point2[0]=cos(theta21)*centr1+cx2;
		point2[1]=sin(theta21)*centr1+cz2;
		point2[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta21-ct02)+cb2;
		//		G4cout<<"test1"<<G4endl;

	      }else if(rho2>rho1){ //rho1<rho2
		G4double theta12=atan2(cz1-cz2,cx1-cx2);
		G4double centr=rho2-(rho2-cent_dist-rho1)/2; //rho1<rho2
		point1[0]=cos(theta12)*centr+cx2;
		point1[1]=sin(theta12)*centr+cz2;
		point1[2]=-1.*tpcData[j].tpcqq*ca2*rho2*(theta12-ct02)+cb2;

		G4double theta21=atan2(cz1-cz2,cx1-cx2);
		G4double centr1=rho1+(rho2-cent_dist-rho1)/2; //rho1<rho2
		point2[0]=cos(theta21)*centr1+cx1;
		point2[1]=sin(theta21)*centr1+cz1;
		point2[2]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta21-ct01)+cb1;
	      }

	      vtxxfit[i]=point1[0];
	      vtxzfit[i]=point1[1];
	      vtxyfit[i]=(point1[2]+point2[2])/2.;
	      vtxxfit[j]=point2[0];
	      vtxzfit[j]=point2[1];
	      vtxyfit[j]=(point1[2]+point2[2])/2.;
	      vtx_flag[i]=6;
	      vtx_flag[j]=6;
	    } else {

	      //k = CircleIntersect(cx1,cz1,rho1,cx2,cz2,rho2,point1,point2);
	      k = CircleIntersect(cx1,cz1,rho1,cx2,cz2,rho2,ca1,cb1,ct01,tpcData[i].tpcqq,ca2,cb2,ct02,tpcData[j].tpcqq,point1,point2);
	      if(k == 0) {
		G4cout << "no solution" << G4endl;
	      }else if(k>0){

		G4double dist1=sqrt(pow(point1[0]-tpcData[i].tpcvtxx,2)+pow(point1[1]-tpcData[i].tpcvtxz,2));
		G4double dist2=sqrt(pow(point2[0]-tpcData[i].tpcvtxx,2)+pow(point2[1]-tpcData[i].tpcvtxz,2));

		if(dist1<=dist2){//point1 is correct
		  vtxxfit[i]=point1[0];
		  vtxzfit[i]=point1[1];
		  vtxyfit[i]=point1[2];

		  vtxxfit[j]=point1[0];
		  vtxzfit[j]=point1[1];
		  vtxyfit[j]=point1[2];
		}else if(dist1>dist2){//point1 is correct
		  vtxxfit[i]=point2[0];
		  vtxzfit[i]=point2[1];
		  vtxyfit[i]=point2[2];

		  vtxxfit[j]=point2[0];
		  vtxzfit[j]=point2[1];
		  vtxyfit[j]=point2[2];
		}
		vtx_flag[i]=7;
		vtx_flag[j]=7;
	      }

	      mom_theta[i]=atan2(vtxzfit[i]-cz[i],vtxxfit[i]-cx[i])-acos(-1.)/2;
	      mom_theta[j]=atan2(vtxzfit[j]-cz[j],vtxxfit[j]-cx[j])-acos(-1.)/2;

	      vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      vtxpxfit[j]=cos(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
	      vtxpzfit[j]=sin(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
	      //	    G4cout<<"env_helm_field:"<<env_helm_field<<G4endl;
	    }
	  }
	  /////vertex reconstruction with beam
	  /*	  else if(sh_paID[i]==0. && sh_paID[i] != sh_paID[j] ){
		  G4double rho1=rad[i];
		  G4double cx1=cx[i];

		  G4double cent_dist=sqrt(pow(cx1-cx2,2)+pow(cz1-cz2,2));

		  double point1[2]={0};
		  G4int k;

		  if((cent_dist-(rho1+rho2))>0.){


		  G4double theta12=atan2(cz2-cz1,cx2-cx1);
		  G4double centr=rho1+(cent_dist-(rho1+rho2))/2;

		  point1[0]=cos(theta12)*centr+cx1;
		  point1[1]=sin(theta12)*centr+cz1;

		  G4double theta21=atan2(cz1-cz2,cx1-cx2);
		  G4double centr1=rho2+(cent_dist-(rho1+rho2))/2;
		  point2[0]=cos(theta21)*centr1+cx2;
		  point2[1]=sin(theta21)*centr1+cz2;

		  vtxxfit[i]=point1[0];
		  vtxzfit[i]=point1[1];
		  vtxxfit[j]=point1[0];
		  vtxzfit[j]=point1[1];
		  }else  if((cent_dist+fmin(rho1,rho2))<fmax(rho1,rho2)){

		  if(rho1>=rho2){ //rho1>rho2
		  G4double theta12=atan2(cz2-cz1,cx2-cx1);
		  G4double centr=rho1-(rho1-cent_dist-rho2)/2; //rho1>rho2
		  point1[0]=cos(theta12)*centr+cx1;
		  point1[1]=sin(theta12)*centr+cz1;

		  G4double theta21=atan2(cz2-cz1,cx2-cx1);
		  G4double centr1=rho2+(rho1-cent_dist-rho2)/2.; //rho1>rho2
		  point2[0]=cos(theta21)*centr1+cx2;
		  point2[1]=sin(theta21)*centr1+cz2;
		  //		G4cout<<"test1"<<G4endl;

		  }else if(rho2>rho1){ //rho1<rho2
		  G4double theta12=atan2(cz1-cz2,cx1-cx2);
		  G4double centr=rho2-(rho2-cent_dist-rho1)/2; //rho1<rho2
		  point1[0]=cos(theta12)*centr+cx2;
		  point1[1]=sin(theta12)*centr+cz2;

		  G4double theta21=atan2(cz1-cz2,cx1-cx2);
		  G4double centr1=rho1+(rho2-cent_dist-rho1)/2; //rho1<rho2
		  point2[0]=cos(theta21)*centr1+cx1;
		  point2[1]=sin(theta21)*centr1+cz1;
		  }

		  vtxxfit[i]=point1[0];
		  vtxzfit[i]=point1[1];
		  vtxxfit[j]=point1[0];
		  vtxzfit[j]=point1[1];
		  } else {

		  k = CircleIntersect(cx1,cz1,rho1,cx2,cz2,rho2,point1,point2);
		  if(k == 0) {
		  G4cout << "no solution" << G4endl;
		  }else if(k>0){

		  G4double dist1=sqrt(pow(point1[0]-tpcData[i].tpcvtxx,2)+pow(point1[1]-tpcData[i].tpcvtxz,2));
		  G4double dist2=sqrt(pow(point2[0]-tpcData[i].tpcvtxx,2)+pow(point2[1]-tpcData[i].tpcvtxz,2));

		  if(dist1<=dist2){//point1 is correct
		  vtxxfit[i]=point1[0];
		  vtxzfit[i]=point1[1];

		  vtxxfit[j]=point1[0];
		  vtxzfit[j]=point1[1];
		  }else if(dist1>dist2){//point1 is correct
		  vtxxfit[i]=point2[0];
		  vtxzfit[i]=point2[1];

		  vtxxfit[j]=point2[0];
		  vtxzfit[j]=point2[1];
		  }
		  }

		  mom_theta[i]=atan2(vtxzfit[i]-cz[i],vtxxfit[i]-cx[i])-acos(-1.)/2;
		  mom_theta[j]=atan2(vtxzfit[j]-cz[j],vtxxfit[j]-cx[j])-acos(-1.)/2;

		  vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
		  vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
		  vtxpxfit[j]=cos(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
		  vtxpzfit[j]=sin(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);

		  }
		  }
	  */

	}
      }
    }


    /*
    //------------------&^-^-------------------//
    ////rungekutta study
    //------------------&^-^-------------------//
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    //////////// rungekutta tracking based on LEPS TPC analyzer
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    int sector, lay;
    //    cir_r[kk]=rad[kk];
    //    cir_x[kk]=cx[kk];
    //    cir_z[kk]=cz[kk];
    double rkpar[5]={0};
    int iflag=0.;

    Switch sw;
    Track tracks[MAX_TRACK];
    initTrack(tracks);
    for( int kk=0; kk<tpctrNum; kk++){
    ///rkpara: x, y, u(px),v(py), q/p
    tracks[kk].rKInitPara[0]=x[kk][0];
    tracks[kk].rKInitPara[1]=y[kk][0];
    tracks[kk].rKInitPara[2]=vtxpxfit[kk]/sqrt(vtxpxfit[kk]*vtxpxfit[kk]+vtxpyfit[kk]*vtxpyfit[kk]+vtxpzfit[kk]*vtxpzfit[kk]);
    tracks[kk].rKInitPara[3]=vtxpyfit[kk]/sqrt(vtxpxfit[kk]*vtxpxfit[kk]+vtxpyfit[kk]*vtxpyfit[kk]+vtxpzfit[kk]*vtxpzfit[kk]);
    tracks[kk].rKInitPara[4]=tpcData[kk].tpcqq/sqrt(vtxpxfit[kk]*vtxpxfit[kk]+vtxpyfit[kk]*vtxpyfit[kk]+vtxpzfit[kk]*vtxpzfit[kk]);

    for( int j = 0; j < c[kk]; j++){
    tracks[kk].x[j][0] = x[kk][j];
    //	std::cout<<"j:"<<j<<"::"<<tracks[kk].x[j][0]<<std::endl;
    tracks[kk].x[j][1] = y[kk][j];
    tracks[kk].x[j][2] = z[kk][j];
    tracks[kk].numHits ++;
    //	std::cout<<"num_hits:"<<c[kk]<<std::endl;
    }
    }
    for( G4int kk=0; kk<tpctrNum; kk++){
    if(c[kk]>5.){
    RungeKuttaTracker rungekuttatrack(tracks+kk);
    }
    }

    */

    ///////////////////////vertex momentum for P_t
    G4int trn[MAX_TRACK];
    //// trancated mean --> now mean
    for(G4int i=0;i<tpctrNum;i++){
      trn[i]=c[i]*(truncated_mean_cut);
      G4double trtmp[MAX_TRACK]={0.000000000};
      for(G4int iii=0;iii<MAX_TRACK;iii++){
	trtmp[iii]=0.000000000;
      }

      for(G4int l=0;l<trn[i];l++){ //--> loop truncated number
	for(G4int k=0;k<c[i];k++){
	  if(l==0){
	    if(trtmp[l]<ede[i][k]){
	      trtmp[l]=ede[i][k];
	    }
	  }else if(l>0){
	    if(trtmp[l-1]>ede[i][k]){
	      if(trtmp[l]<ede[i][k]){
		trtmp[l]=ede[i][k];
	      }
	    }
	  }

	}//--loop end
      }
      for(G4int j=0;j<c[i];j++){
	if(trn[i]!=0.){
	  G4int sh_ch=1;
	  for(G4int jj=0;jj<trn[i];jj++){
	    if(ede[i][j]==trtmp[jj]){
	      sh_ch=-1;
	    }
	  }
	  if(sh_ch>0){
	    trmean[i]=trmean[i]+ede[i][j]/(c[i]-trn[i]);
	  }
	}else if(trn[i]==0.){
	  trmean[i]=trmean[i]+ede[i][j]/(c[i]-trn[i]);
	}

      }
      //      }
    }
    for( G4int i=0; i<HitNum; i++){

      // std::cout<<"hoge!!!"<<std::endl;
      // getchar();
      anaRoot.FillNtrk(counterData[i].ntrk);
      anaRoot.FillTime(counterData[i].time);
      anaRoot.FillPos(counterData[i].pos);
      anaRoot.FillPos0(counterData[i].pos0,counterData[i].resoX);
      anaRoot.FillMom(counterData[i].mom);
      anaRoot.FillTrackID(counterData[i].trackID);
      anaRoot.FillParticleID(counterData[i].particleID);
      anaRoot.FillPadLay(counterData[i].iLay);
      anaRoot.FillPadRow(counterData[i].iRow);
      anaRoot.FillBeta(counterData[i].beta);
      anaRoot.FillEdep(counterData[i].edep);
      anaRoot.FilldEdx(counterData[i].dedx);
      anaRoot.FillsLength(counterData[i].slength);

      anaRoot.FillLayerPad(counterData[i].iLay,counterData[i].iPad);///multiplicity

      anaRoot.incHit();

    }

    //
    // Scinti.
    //

    for( G4int i=0; i<HitNumScint; i++){
      anaRoot.FillScintData(scintData[i].time, scintData[i].pos,
			    scintData[i].mom,
			    scintData[i].trackID, scintData[i].particleID,
			    scintData[i].detectorID,scintData[i].massSH,scintData[i].qqSH,scintData[i].parentID,
			    scintData[i].scintvtxpx,scintData[i].scintvtxpy,scintData[i].scintvtxpz,
			    scintData[i].scintvtxx,scintData[i].scintvtxy,scintData[i].scintvtxz,
			    scintData[i].length
			    );

    }


    //
    // AC
    //

    for( G4int i=0; i<HitNumAC; i++){
      anaRoot.FillACData(acData[i].time, acData[i].pos,
			 acData[i].mom,
			 acData[i].trackID, acData[i].particleID,
			 acData[i].detectorID,acData[i].massSH,acData[i].qqSH,acData[i].parentID,
			 acData[i].acvtxpx,acData[i].acvtxpy,acData[i].acvtxpz,
			 acData[i].acvtxx,acData[i].acvtxy,acData[i].acvtxz,
			 acData[i].length
			 );

    }



    //
    // nbar
    //

    for( G4int i=0; i<HitNumNBAR; i++){
      anaRoot.FillNBARData(nbarData[i].time, nbarData[i].pos,
			   nbarData[i].mom,
			   nbarData[i].trackID, nbarData[i].particleID,
			   nbarData[i].detectorID,nbarData[i].massSH,nbarData[i].qqSH,nbarData[i].parentID,
			   nbarData[i].nbarvtxpx,nbarData[i].nbarvtxpy,nbarData[i].nbarvtxpz,
			   nbarData[i].nbarvtxx,nbarData[i].nbarvtxy,nbarData[i].nbarvtxz,
			   nbarData[i].length
			   );

    }


    //
    // DC
    //

    for( G4int i=0; i<HitNumDC; i++){
      anaRoot.FillDCData(dcData[i].time, dcData[i].pos,
			 dcData[i].mom,
			 dcData[i].trackID, dcData[i].particleID,
			 dcData[i].detectorID,dcData[i].massSH,dcData[i].qqSH,dcData[i].parentID,
			 dcData[i].vtxpx,dcData[i].vtxpy,dcData[i].vtxpz,
			 dcData[i].vtxx,dcData[i].vtxy,dcData[i].vtxz,
			 dcData[i].length
			 );

    }


    //
    // SCH
    //

    for( G4int i=0; i<HitNumSCH; i++){
      anaRoot.FillSCHData(schData[i].time, schData[i].pos,
			 schData[i].mom,
			 schData[i].trackID, schData[i].particleID,
			 schData[i].detectorID,schData[i].massSH,schData[i].qqSH,schData[i].parentID,
			 schData[i].vtxpx,schData[i].vtxpy,schData[i].vtxpz,
			 schData[i].vtxx,schData[i].vtxy,schData[i].vtxz,
			 schData[i].length
			 );

    }


    //
    // ftof
    //

    //    std::cout<<"tof:"<<HitNumFTOF<<std::endl;
    //    std::cout<<"dc:"<<HitNumDC<<std::endl;
    //    std::cout<<"SCH:"<<HitNumSCH<<std::endl;

    for( G4int i=0; i<HitNumFTOF; i++){
      anaRoot.FillFTOFData(ftofData[i].time, ftofData[i].pos,
			   ftofData[i].mom,
			   ftofData[i].trackID, ftofData[i].particleID,
			   ftofData[i].detectorID,ftofData[i].massSH,ftofData[i].qqSH,ftofData[i].parentID,
			   ftofData[i].vtxpx,ftofData[i].vtxpy,ftofData[i].vtxpz,
			   ftofData[i].vtxx,ftofData[i].vtxy,ftofData[i].vtxz,
			   ftofData[i].length
			   );

    }


    //
    // Target.
    //

    for( G4int i=0; i<HitNumTarget; i++){
      anaRoot.FillTargetData(i,
			     targetData[i].targetparticleid,
			     targetData[i].targetparentid,
			     targetData[i].targettrackid,
			     targetData[i].targetpos,
			     targetData[i].targetvtx
			     );
    }


    //
    // TPC
    //
    for( G4int i=0; i<tpctrNum; i++){
      //      G4cout<<"abs"<<abs(env_helm_field)<<G4endl;
      //      G4cout<<"fabs"<<fabs(env_helm_field)<<G4endl;
      anaRoot.FillTPCData(tpcData[i].tpcpx,
			  tpcData[i].tpcpy,tpcData[i].tpcpz,
			  tpcData[i].tpcpp,
			  tpcData[i].tpcpid, tpcData[i].tpcparentid, tpcData[i].tpcparentid_pid,
			  tpcData[i].tpcqq,
			  tpcData[i].tpcpm,tpcData[i].tpcde,
			  tpcData[i].tpclen,mean[i],trmean[i],
			  tpcData[i].tpclay,
			  tpcData[i].tpcvtxpx,tpcData[i].tpcvtxpy,tpcData[i].tpcvtxpz,
			  tpcData[i].tpcvtxx,tpcData[i].tpcvtxy,tpcData[i].tpcvtxz,
			  vtxxfit[i],vtxyfit[i],vtxzfit[i],
			  vtxpxfit[i],Pz[i],vtxpzfit[i],cir_r[i]*(0.299792458)*fabs(env_helm_field),
			  cir_r[i],cir_x[i],cir_z[i],test[i],
			  vtx_flag[i], a_fory[i], b_fory[i]
			  );
    }

    // detID = 0;

    //------------------&^-^-------------------//
    ////rungekutta study for kurama
    //------------------&^-^-------------------//
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    //////////// rungekutta tracking based on LEPS TPC analyzer
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    if(HitNumDC>0){
      Track kurama_tr[MAX_TRACK];
      //    std::cout<<"start kurama tracking"<<std::endl;
      kurama_tr[0].numHits=0.;
      //    std::cout<<kurama_tr[0].numHits<<std::endl;
      initTrack_ku(kurama_tr);
      int num_cand_track_id[30];
      int num_cand_track=0;
      double dc_res[22]=
	{0.050,0.05,0.05,0.05,//ssd1
	 0.05,0.05,0.05,0.05,//ssd2
	 0.300,0.300,0.300,0.300,0.300,0.300,//dc1
	 0.300,0.300,0.300,0.300,//dc2
	 0.300,0.300,0.300,0.300//dc3
	};
      int dc_idwi[22]=
	{1,2,1,2,//ssd1
	 1,2,1,2,
	 1,1,2,2,1,1,//dc1
	 1,2,1,2,//dc2
	 1,2,1,2//dc3
	};

      int dc_id[22]=
	{0,0,0,0,//ssd1
	 0,0,0,0,//ssd2
	 101,102,103,104,105,106,//dc1
	 121,122,123,124,//dc2
	 131,132,133,134//dc3
	};
      int num_hit_track[30] = {};

      for(int i=0 ; i<HitNumDC ; i++){
	//      std::cout << dcData[i].trackID << ":" <<dcData[i].detectorID << std::endl;
	//      std::cout<<"position:"<<dcData[i].pos[2]<<std::endl;
	if(i==0){
	  num_cand_track_id[0]=dcData[i].trackID;
	  num_cand_track++;
	  num_hit_track[num_cand_track-1]=1;
	  //initial para for rungekutta
	  kurama_tr[num_cand_track-1].rKInitPara[0]=dcData[i].pos[0];
	  kurama_tr[num_cand_track-1].rKInitPara[1]=dcData[i].pos[1];
	  kurama_tr[num_cand_track-1].rKInitPara[2]=dcData[i].mom[0]/sqrt(dcData[i].mom[0]*dcData[i].mom[0]+dcData[i].mom[1]*dcData[i].mom[1]+dcData[i].mom[2]*dcData[i].mom[2]);
	  kurama_tr[num_cand_track-1].rKInitPara[3]=dcData[i].mom[1]/sqrt(dcData[i].mom[0]*dcData[i].mom[0]+dcData[i].mom[1]*dcData[i].mom[1]+dcData[i].mom[2]*dcData[i].mom[2]);
	  kurama_tr[num_cand_track-1].rKInitPara[4]=dcData[i].qqSH/sqrt(dcData[i].mom[0]*dcData[i].mom[0]+dcData[i].mom[1]*dcData[i].mom[1]+dcData[i].mom[2]*dcData[i].mom[2]);


	  //      for( int j = 0; j < c[kk]; j++){
	  //	kurama_tr[0].x[j][0] = x[kk][j];
	  //	kurama_tr[0].x[j][1] = y[kk][j];
	  //	kurama_tr[0].x[j][2] = z[kk][j];
	  //	kurama_tr[0].numHits ++;
	  //	kurama_tr[0].
	}else{
	  if(num_cand_track_id[num_cand_track-1] != dcData[i].trackID){
	    num_cand_track_id[num_cand_track]=dcData[i].trackID;
	    num_cand_track++;
	    num_hit_track[num_cand_track-1]=1;

	    //initial para for rungekutta
	    kurama_tr[num_cand_track-1].rKInitPara[0]=dcData[i].pos[0]+CLHEP::RandGauss::shoot(0.,1.);
	    kurama_tr[num_cand_track-1].rKInitPara[1]=dcData[i].pos[1]+CLHEP::RandGauss::shoot(0.,1.);
	    kurama_tr[num_cand_track-1].rKInitPara[2]=dcData[i].mom[0]/sqrt(dcData[i].mom[0]*dcData[i].mom[0]+dcData[i].mom[1]*dcData[i].mom[1]+dcData[i].mom[2]*dcData[i].mom[2])*CLHEP::RandGauss::shoot(1.,0.02);
	    kurama_tr[num_cand_track-1].rKInitPara[3]=dcData[i].mom[1]/sqrt(dcData[i].mom[0]*dcData[i].mom[0]+dcData[i].mom[1]*dcData[i].mom[1]+dcData[i].mom[2]*dcData[i].mom[2])*CLHEP::RandGauss::shoot(1.,0.02);
	    kurama_tr[num_cand_track-1].rKInitPara[4]=dcData[i].qqSH/sqrt(dcData[i].mom[0]*dcData[i].mom[0]+dcData[i].mom[1]*dcData[i].mom[1]+dcData[i].mom[2]*dcData[i].mom[2])*CLHEP::RandGauss::shoot(1.,0.0);

	  }else{
	    num_hit_track[num_cand_track-1]++;
	  }

	}
	kurama_tr[num_cand_track-1].numHits=num_hit_track[num_cand_track-1];
	kurama_tr[num_cand_track-1].x[num_hit_track[num_cand_track-1]-1][0]=dcData[i].pos[0]+CLHEP::RandGauss::shoot(0,0.2);
	kurama_tr[num_cand_track-1].x[num_hit_track[num_cand_track-1]-1][1]=dcData[i].pos[1]+CLHEP::RandGauss::shoot(0,0.2);
	kurama_tr[num_cand_track-1].x[num_hit_track[num_cand_track-1]-1][2]=dcData[i].pos[2];

	kurama_tr[num_cand_track-1].uv[num_hit_track[num_cand_track-1]-1][0]=dcData[i].mom[0]/dcData[i].mom[2]+CLHEP::RandGauss::shoot(0.,0.02);
	kurama_tr[num_cand_track-1].uv[num_hit_track[num_cand_track-1]-1][1]=dcData[i].mom[1]/dcData[i].mom[2]+CLHEP::RandGauss::shoot(0.,0.02);



	for(int j=0;j<22;j++){
	  if(dcData[i].detectorID==dc_id[j]){
	    kurama_tr[num_cand_track-1].res[num_hit_track[num_cand_track-1]-1]=dc_res[j];
	    kurama_tr[num_cand_track-1].idwi[num_hit_track[num_cand_track-1]-1]=dc_idwi[j];
	    //	  std::cout<<"plane number:"<<j<<std::endl;
	    //	  std::cout<<"res:"<<kurama_tr[num_cand_track-1].res[num_hit_track[num_cand_track-1]-1]<<std::endl;
	  }
	}
      }

      //    std::cout<<num_cand_track<<std::endl;

      for( G4int kk=0; kk<num_cand_track; kk++){
	//      std::cout<<"number of hits in a track:"<<kurama_tr[kk].numHits<<std::endl;
	//      if(kurama_tr[kk].numHits>9){
	if(kurama_tr[kk].numHits>13){
	  //	RungeKuttaTracker rungekuttatrack(0, kurama_tr+kk);// 1 is C matrix usage.
	}
      }

    }
    //anaRoot.FillTree();
  }//trigger parts
  anaRoot.FillTree();

  return 0;
}

/////shhwang TPC
void TPCAnaManager::SetCounterData( G4int ntrk,G4double time, G4ThreeVector pos,
				    G4ThreeVector mom,
				    G4int track, G4int particle,
				    G4int iLay,  G4int iRow, G4double beta,
				    G4double edep, G4int parentid,
				    G4double /* tlength */, G4double slength )
{
  G4int hitnum = HitNum;
  G4bool flag=true;
  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }

  //  G4ThreeVector tar_pos(0.,0.*mm,-150.*mm);

  G4ThreeVector tar_pos(0.,0.,target_pos_z);
  G4ThreeVector sh_pos(0.,0.,0.);
  sh_pos=pos-tar_pos;

  G4double sh_r = sh_pos.r();
  G4double sh_theta = sh_pos.theta();
  G4double sh_phi = sh_pos.phi();

  G4double sh_x = sh_r*sin(sh_theta)*cos(sh_phi);
  G4double sh_y = sh_r*sin(sh_theta)*sin(sh_phi);
  G4double sh_z = sh_r*cos(sh_theta);

  counterData[hitnum].particleID = particle;
  ////shhwang check, check a multiplicity of layers
  for(G4int i=0;i<hitnum;i++){
    if((counterData[i].iLay == iLay &&
  	counterData[i].ntrk == ntrk)){
      flag = false;
    }
  }
  flag=true;
  if(flag == true){
    counterData[hitnum].ntrk = ntrk;
    counterData[hitnum].time = time;
    counterData[hitnum].beta = beta;
    counterData[hitnum].dedx = edep/slength;
    counterData[hitnum].edep = edep;
    counterData[hitnum].slength = slength;

    //////shhwang position smearing////


    G4double sh_alpha =  atan2(sh_x,sh_z);
    G4double sh_rho =  sqrt(pow(sh_z,2)+pow(sh_x,2));
    //    G4double sh_dalpha = 0.300*mm/sh_rho; // rho * theta = arc --> sigma=300 um
    //    G4double sh_smear_alpha = CLHEP::RandGauss::shoot(sh_alpha, sh_dalpha);
    G4double sh_sigmaY = 0.500*mm; //--> smearing : 400 um

    G4double ang_sh=atan2(sh_pos.getY(),sh_pos.getX());
    // G4double ang_check=0;
    if(ang_sh>acos(-1.)){
      ang_sh=ang_sh-2*acos(-1.);
      // ang_check=-1;
    }else{
      // ang_check=1.;
    }
    // G4double arc_sh=pad_in[iLay]*ang_sh;
    // G4int ith_pad_in=arc_sh/pad_in_width;
    // G4int ith_pad_out=arc_sh/pad_out_width;


    // G4double delta_x=0.;

    // if( env_pad_config ==1 ){
    //   if(iLay<pad_in_num){   // const G4int pad_in_num=10;
    // 	delta_x=arc_sh-ith_pad_in*pad_in_width;
    // 	if(delta_x<0){
    // 	  delta_x=delta_x+pad_in_width;
    // 	}
    //   }else if(iLay>=pad_in_num){
    // 	delta_x=arc_sh-ith_pad_out*pad_out_width;
    // 	if(delta_x<0){
    // 	  delta_x=delta_x+pad_out_width;
    // 	}
    //   }
    // }else if( env_pad_config ==2 ){
    //   if(iLay<pad_in_num){   //

    // 	ith_pad_in=arc_sh/seg_width[iLay];
    // 	delta_x=arc_sh-ith_pad_in*seg_width[iLay];
    // 	if(delta_x<0){
    // 	  delta_x=delta_x+seg_width[iLay];
    // 	}

    //   }else if(iLay>=pad_in_num){
    // 	ith_pad_out=arc_sh/seg_width[iLay];
    // 	delta_x=arc_sh-ith_pad_out*seg_width[iLay];
    // 	if(delta_x<0){
    // 	  delta_x=delta_x+seg_width[iLay];
    // 	}

    //   }
    // }


    ///include saho-san's code

    // G4double compy=0.;
    G4double compx=0.;

    // G4int n_electron=0; //output
    // G4int n_pad=0;//output
    // G4double x_rms=0;//output
    // G4double x_track = delta_x;// pad edge
    // G4double y_track = double(250.+50.+pos.getY());// 50 is office set
    // if(y_track<0 || y_track>600) G4cout<<"Y_track estimation is wrong(y_track:y_pos)"<<y_track<<":"<<pos.getY()<<G4endl;


    // G4double dxdz_track=tan(acos((sh_pos.getX()*mom.getX()+sh_pos.getZ()*mom.getZ())/(sqrt(sh_pos.getX()*sh_pos.getX()+sh_pos.getZ()*sh_pos.getZ())*sqrt(mom.getX()*mom.getX()+mom.getZ()*mom.getZ()))));


    // //pad and particle angle


    // G4double dydz_track=0;

    // if( env_pad_config ==1){
    //   if(iLay<pad_in_num){
    // 	ResHypTPC reshyptpc(pad_in_width, pad_length_in, 0.1,0.18, 0);
    // 	compx = reshyptpc.getXDeviation(n_electron, n_pad, x_rms, x_track, y_track, dxdz_track, dydz_track);
    // 	compy = reshyptpc.getYDeviation(pos.getY());
    //   }else if(iLay>=pad_in_num){
    // 	ResHypTPC reshyptpc(pad_out_width, pad_length_out, 0.1,0.18, 0);
    // 	compx = reshyptpc.getXDeviation(n_electron, n_pad, x_rms, x_track, y_track, dxdz_track, dydz_track);
    // 	compy = reshyptpc.getYDeviation(pos.getY());
    //   }


    // }else if( env_pad_config ==2){
    //   if(iLay<pad_in_num){
    // 	ResHypTPC reshyptpc(double(seg_width[iLay]-0.5), pad_length_in, 0.1,0.18, 0);
    // 	compx = reshyptpc.getXDeviation(n_electron, n_pad, x_rms, x_track, y_track, dxdz_track, dydz_track);
    // 	compy = reshyptpc.getYDeviation(pos.getY());
    //   }else if(iLay>=pad_in_num){
    //ResHypTPC reshyptpc(double(seg_width[iLay]-0.5), pad_length_out, 0.1,0.18, 0);
    // 	compx = reshyptpc.getXDeviation(n_electron, n_pad, x_rms, x_track, y_track, dxdz_track, dydz_track);
    // 	compy = reshyptpc.getYDeviation(pos.getY());
    //   }
    // }
    compx = GetTransverseRes(sh_y);
    double s_compx = CLHEP::RandGauss::shoot(0.,compx);

    // std::cout<<"compx="<<compx<<", sh_sigmaY"<<sh_sigmaY<<std::endl;
    // getchar();
    //G4double sh_dalpha = compx/sh_rho; // rho * theta = arc --> from sako-san's code
    G4double sh_dalpha = s_compx/sh_rho; // rho * theta = arc --> from sako-san's code
    G4double sh_smear_alpha = sh_alpha+sh_dalpha;
    //    G4cout<<compx<<":"<<sh_dalpha<<G4endl;

    counterData[hitnum].resoX = compx;

    counterData[hitnum].pos[ThreeVector::Z] = sh_rho*cos(sh_smear_alpha)+tar_pos.getZ();
    counterData[hitnum].pos[ThreeVector::X] = sh_rho*sin(sh_smear_alpha);
    counterData[hitnum].pos[ThreeVector::Y] = CLHEP::RandGauss::shoot(sh_y,sh_sigmaY);

    counterData[hitnum].pos0[ThreeVector::X] = pos.getX();
    counterData[hitnum].pos0[ThreeVector::Y] = pos.getY();
    counterData[hitnum].pos0[ThreeVector::Z] = pos.getZ();

    counterData[hitnum].mom[ThreeVector::X] = mom.getX();
    counterData[hitnum].mom[ThreeVector::Y] = mom.getY();
    counterData[hitnum].mom[ThreeVector::Z] = mom.getZ();

    counterData[hitnum].trackID = track;
    counterData[hitnum].particleID = particle;
    counterData[hitnum].iLay = iLay;
    G4int iPad=0.;

    if( env_pad_config ==2 ){
    //    G4bool pass_check=false;
      G4bool pass_check=true;
      G4double cur_angle= (acos(-1.)-atan2(sh_x,sh_z))*180./acos(-1.);
    //    G4cout<<"--------------------"<<G4endl;
    //    G4cout<<"currrent angle:"<<cur_angle<<G4endl;
    //    G4cout<<"layer angle:"<<angle[iLay]<<G4endl;
    //    G4cout<<"seg angle:"<<seg_angle[iLay]<<G4endl;
      if(iLay<pad_in_num){
	G4double check_num_pads=(cur_angle)/seg_angle[iLay];
      //      G4cout<<check_num_pads<<G4endl;
	iPad=int(check_num_pads);
      }else if(iLay>=pad_in_num){
	G4double check_num_pads=(cur_angle-angle[iLay])/seg_angle[iLay];
      //      G4cout<<check_num_pads<<G4endl;
	iPad=int(check_num_pads);
      }
      if(iPad>numpads[iLay]){
	G4cout<<"this code has a error(iPad:numpads)-->"<<iPad<<":"<<numpads[iLay]<<G4endl;
      }
      if(pass_check){
	//    G4cout<<"iLay:"<<iLay<<G4endl;
	//        G4cout<<"iLay:"<<iLay<<G4endl;
	//        G4cout<<"iPad:"<<iPad<<G4endl;
	counterData[hitnum].iPad = iPad;
      }else{
	G4cout<<"wrong:"<<iLay<<G4endl;
      }
    }

    counterData[hitnum].iRow = iRow;
    counterData[hitnum].parentID = parentid;
    HitNum++;

    if(particle==321)
      HitNum_K++;

    if(particle==2212)
      HitNum_p++;
  }

  return;
}

//_____________________________________________________________________________
void
TPCAnaManager::SetTPCData( G4int tpctr2, G4int tpcpid2, G4int tpcparentid2,
			   G4int tpcparentid_pid2, G4double tpcpx2,
			   G4double tpcpy2, G4double tpcpz2,
			   G4double /* tpcpp2 */,
			   G4int tpcqq2, G4double tpcpm2, G4double tpcde2,
			   G4double tpclen2, G4int tpclay2,
			   G4double vtxpxtpc2, G4double vtxpytpc2,
			   G4double vtxpztpc2,
			   G4double vtxxtpc2, G4double vtxytpc2,
			   G4double vtxztpc2, G4double vtxene2 )
{
  G4int hitnum = tpctr2;

  // G4double theta=acos(tpcpz2/tpcpp2);

  tpcData[hitnum].tpctr = tpctr2;
  tpcData[hitnum].tpcpid = tpcpid2;
  tpcData[hitnum].tpcparentid = tpcparentid2;
  tpcData[hitnum].tpcparentid_pid = tpcparentid_pid2;

  //// w/o smearing
  tpcData[hitnum].tpcpx = tpcpx2;
  tpcData[hitnum].tpcpy = tpcpy2;
  tpcData[hitnum].tpcpz = tpcpz2;
  //// with smearing
  //    tpcData[hitnum].tpcpx = px;
  //    tpcData[hitnum].tpcpz = pz;
  //    tpcData[hitnum].tpcpy = py;


  //kine E = sqrt(p^2+m^2)-m
  //p=sqrt((E+m)^2-m^2)
  G4double totalmom=sqrt(pow(vtxene2+tpcpm2,2)-pow(tpcpm2,2));
  tpcData[hitnum].tpcvtxpx = totalmom*vtxpxtpc2;
  tpcData[hitnum].tpcvtxpy = totalmom*vtxpytpc2;
  tpcData[hitnum].tpcvtxpz = totalmom*vtxpztpc2;

  tpcData[hitnum].tpcvtxx = vtxxtpc2;
  tpcData[hitnum].tpcvtxy = vtxytpc2;
  tpcData[hitnum].tpcvtxz = vtxztpc2;

  //// with smearing
  //  tpcData[hitnum].tpcpp = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

  //// w/o smearing
  tpcData[hitnum].tpcpp = sqrt(pow(tpcpx2,2)+pow(tpcpy2,2)+pow(tpcpz2,2));
  //  tpcData[hitnum].tpcppfit = sqrt(pow(,2)+pow(tpcpy2,2));

  tpcData[hitnum].tpcqq = tpcqq2;
  tpcData[hitnum].tpcpm = tpcpm2;
  tpcData[hitnum].tpclen = tpclen2;
  tpcData[hitnum].tpcdedx = tpcde2;
  tpcData[hitnum].tpclay = tpclay2;
  tpctrNum++;
  return;
}


void TPCAnaManager::SetPrimaryInfo(G4double mm_d, G4double mm_p, G4double theta, G4double theta_scat, G4double theta_CM)
{
  primaryInfo.mm_d = mm_d;
  primaryInfo.mm_p = mm_p;
  primaryInfo.theta = theta;
  primaryInfo.theta_scat = theta_scat;
  primaryInfo.theta_CM = theta_CM;
}


void TPCAnaManager::SetNumberOfPrimaryParticle(G4int num)
{
  primaryParticle.NumOfParticle = num;
}


void TPCAnaManager::SetPrimaryParticle(G4double px, G4double py, G4double pz)
{
  G4int id = 0;
  primaryParticle.p0[id][0] = px;
  primaryParticle.p0[id][1] = py;
  primaryParticle.p0[id][2] = pz;
  primaryParticle.p0[id][3] = sqrt(pow(px,2.0)+pow(py,2.0)+pow(pz,2.0));
  primaryParticle.p0[id][4] = 0.0;
}


void TPCAnaManager::SetGeneratorID(G4int gen)
{
  primaryBeam.gen = gen;
}

void TPCAnaManager::SetModeID(G4int mode)
{
  primaryBeam.mode = mode;
}


void TPCAnaManager::SetPrimaryParticle(G4int id, G4double px, G4double py, G4double pz, G4double mass)
{
  if(id>primaryParticle.NumOfParticle-1){
    G4cout << "SetPrimaryParticle: Invalid Primary particle ID" << G4endl;
    exit(0);
  } else {
    primaryParticle.p0[id][0] = px;
    primaryParticle.p0[id][1] = py;
    primaryParticle.p0[id][2] = pz;
    primaryParticle.p0[id][3] = sqrt(pow(px,2.0)+pow(py,2.0)+pow(pz,2.0)+pow(mass,2.0));
    primaryParticle.p0[id][4] = mass;
  }
}

void TPCAnaManager::SetPrimaryParticle(G4int id, G4double px, G4double py, G4double pz, G4double mass, G4int pid)
{
  if(id>primaryParticle.NumOfParticle-1){
    G4cout << "SetPrimaryParticle: Invalid Primary particle ID" << G4endl;
    exit(0);
  } else {
    primaryParticle.p0[id][0] = px;
    primaryParticle.p0[id][1] = py;
    primaryParticle.p0[id][2] = pz;
    primaryParticle.p0[id][3] = sqrt(pow(px,2.0)+pow(py,2.0)+pow(pz,2.0)+pow(mass,2.0));
    primaryParticle.p0[id][4] = mass;
    primaryParticle.pid0[id] = pid;
  }
}

void TPCAnaManager::SetPrimaryVertex(G4int id, G4double x, G4double y, G4double z)
{
  if(id>primaryParticle.NumOfParticle-1){
    G4cout << "SetPrimaryVertex: Invalid Primary particle ID" << G4endl;
    exit(0);
  } else {
    //    printf("vertex: %f %f %f\n",x,y,z);
    primaryParticle.x0[id][0] = x;
    primaryParticle.x0[id][1] = y;
    primaryParticle.x0[id][2] = z;
  }
}

//_____________________________________________________________________________
void
TPCAnaManager::SetPrimaryBeam( const G4ThreeVector& p )
{
  SetPrimaryBeam( p.x(), p.y(), p.z() );
}

//_____________________________________________________________________________
void
TPCAnaManager::SetPrimaryBeam( G4double px, G4double py, G4double pz )
{
  primaryBeam.pg[0] = px;
  primaryBeam.pg[1] = py;
  primaryBeam.pg[2] = pz;
  primaryBeam.pg[3] = std::sqrt( px*px + py*py + pz*pz );
}

void TPCAnaManager::SetScintData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
				 G4int track, G4int particle, G4int detector, G4double mass,G4int qq,G4int parentid,G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength)
{

  G4int hitnum = HitNumScint;

  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }

  scintData[hitnum].time = time;
  scintData[hitnum].pos[ThreeVector::X] = pos.getX();
  scintData[hitnum].pos[ThreeVector::Y] = pos.getY();
  scintData[hitnum].pos[ThreeVector::Z] = pos.getZ();
  scintData[hitnum].mom[ThreeVector::X] = mom.getX();
  scintData[hitnum].mom[ThreeVector::Y] = mom.getY();
  scintData[hitnum].mom[ThreeVector::Z] = mom.getZ();
  scintData[hitnum].trackID = track;
  scintData[hitnum].massSH = mass;
  scintData[hitnum].qqSH = qq;
  scintData[hitnum].particleID = particle;
  scintData[hitnum].detectorID = detector;
  scintData[hitnum].parentID = parentid;
  scintData[hitnum].length = tlength;

  //kine E = sqrt(p^2+m^2)-m
  //p=sqrt((E+m)^2-m^2)
  //G4ThreeVecor vtxpos, G4ThreeVecor vtxmom, G4double vtxene
  G4double totalmom=sqrt(pow(vtxene+mass,2)-pow(mass,2));
  scintData[hitnum].scintvtxpx = totalmom*double(vtxmom.getX());
  scintData[hitnum].scintvtxpy = totalmom*double(vtxmom.getY());
  scintData[hitnum].scintvtxpz = totalmom*double(vtxmom.getZ());
  scintData[hitnum].scintvtxx = double(vtxpos.getX());
  scintData[hitnum].scintvtxy = double(vtxpos.getY());
  scintData[hitnum].scintvtxz = double(vtxpos.getZ());

  //  G4cout<<"particle:parentid :"<<particle<<":"<<parentid<<G4endl;
  HitNumScint++;
  return;
}


void TPCAnaManager::SetACData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
				 G4int track, G4int particle, G4int detector, G4double mass,G4int qq,G4int parentid,G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength)
{

  G4int hitnum = HitNumAC;

  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }
  acData[hitnum].time = time;
  acData[hitnum].pos[ThreeVector::X] = pos.getX();
  acData[hitnum].pos[ThreeVector::Y] = pos.getY();
  acData[hitnum].pos[ThreeVector::Z] = pos.getZ();
  acData[hitnum].mom[ThreeVector::X] = mom.getX();
  acData[hitnum].mom[ThreeVector::Y] = mom.getY();
  acData[hitnum].mom[ThreeVector::Z] = mom.getZ();
  acData[hitnum].trackID = track;
  acData[hitnum].massSH = mass;
  acData[hitnum].qqSH = qq;
  acData[hitnum].particleID = particle;
  acData[hitnum].detectorID = detector;
  acData[hitnum].parentID = parentid;
  acData[hitnum].length = tlength;

  //kine E = sqrt(p^2+m^2)-m
  //p=sqrt((E+m)^2-m^2)
  //G4ThreeVecor vtxpos, G4ThreeVecor vtxmom, G4double vtxene
  G4double totalmom=sqrt(pow(vtxene+mass,2)-pow(mass,2));
  acData[hitnum].acvtxpx = totalmom*double(vtxmom.getX());
  acData[hitnum].acvtxpy = totalmom*double(vtxmom.getY());
  acData[hitnum].acvtxpz = totalmom*double(vtxmom.getZ());
  acData[hitnum].acvtxx = double(vtxpos.getX());
  acData[hitnum].acvtxy = double(vtxpos.getY());
  acData[hitnum].acvtxz = double(vtxpos.getZ());

  //  G4cout<<"particle:parentid :"<<particle<<":"<<parentid<<G4endl;
  HitNumAC++;

  if(particle==321)//kaon
    HitNumAC_K++;

  if(particle==2212)//proton
    HitNumAC_p++;
  return;
}




void TPCAnaManager::SetNBARData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
				 G4int track, G4int particle, G4int detector, G4double mass,G4int qq,G4int parentid,G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength)
{

  G4int hitnum = HitNumNBAR;

  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }
  nbarData[hitnum].time = time;
  nbarData[hitnum].pos[ThreeVector::X] = pos.getX();
  nbarData[hitnum].pos[ThreeVector::Y] = pos.getY();
  nbarData[hitnum].pos[ThreeVector::Z] = pos.getZ();
  nbarData[hitnum].mom[ThreeVector::X] = mom.getX();
  nbarData[hitnum].mom[ThreeVector::Y] = mom.getY();
  nbarData[hitnum].mom[ThreeVector::Z] = mom.getZ();
  nbarData[hitnum].trackID = track;
  nbarData[hitnum].massSH = mass;
  nbarData[hitnum].qqSH = qq;
  nbarData[hitnum].particleID = particle;
  nbarData[hitnum].detectorID = detector;
  nbarData[hitnum].parentID = parentid;
  nbarData[hitnum].length = tlength;

  //kine E = sqrt(p^2+m^2)-m
  //p=sqrt((E+m)^2-m^2)
  //G4ThreeVecor vtxpos, G4ThreeVecor vtxmom, G4double vtxene
  G4double totalmom=sqrt(pow(vtxene+mass,2)-pow(mass,2));
  nbarData[hitnum].nbarvtxpx = totalmom*double(vtxmom.getX());
  nbarData[hitnum].nbarvtxpy = totalmom*double(vtxmom.getY());
  nbarData[hitnum].nbarvtxpz = totalmom*double(vtxmom.getZ());
  nbarData[hitnum].nbarvtxx = double(vtxpos.getX());
  nbarData[hitnum].nbarvtxy = double(vtxpos.getY());
  nbarData[hitnum].nbarvtxz = double(vtxpos.getZ());

  //  G4cout<<"particle:parentid :"<<particle<<":"<<parentid<<G4endl;
  HitNumNBAR++;

  if(particle==321)//kaon
    HitNumNBAR_K++;

  if(particle==2212)//proton
    HitNumNBAR_p++;

  return;
}




void TPCAnaManager::SetDCData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
			      G4int track, G4int particle, G4int detector, G4double mass,G4int qq,
			      G4int parentid,G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength)
{

  G4int hitnum = HitNumDC;

  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }
  dcData[hitnum].time = time;
  dcData[hitnum].pos[ThreeVector::X] = pos.getX();
  dcData[hitnum].pos[ThreeVector::Y] = pos.getY();
  dcData[hitnum].pos[ThreeVector::Z] = pos.getZ();
  dcData[hitnum].mom[ThreeVector::X] = mom.getX();
  dcData[hitnum].mom[ThreeVector::Y] = mom.getY();
  dcData[hitnum].mom[ThreeVector::Z] = mom.getZ();
  dcData[hitnum].trackID = track;
  //  std::cout<<"DC"<<std::endl;
  //  std::cout<<track<<std::endl;
  dcData[hitnum].massSH = mass;
  dcData[hitnum].qqSH = qq;
  dcData[hitnum].particleID = particle;
  dcData[hitnum].detectorID = detector;
  dcData[hitnum].parentID = parentid;
  dcData[hitnum].length = tlength;

  //kine E = sqrt(p^2+m^2)-m
  //p=sqrt((E+m)^2-m^2)
  //G4ThreeVecor vtxpos, G4ThreeVecor vtxmom, G4double vtxene
  G4double totalmom=sqrt(pow(vtxene+mass,2)-pow(mass,2));
  dcData[hitnum].vtxpx = totalmom*double(vtxmom.getX());
  dcData[hitnum].vtxpy = totalmom*double(vtxmom.getY());
  dcData[hitnum].vtxpz = totalmom*double(vtxmom.getZ());
  dcData[hitnum].vtxx = double(vtxpos.getX());
  dcData[hitnum].vtxy = double(vtxpos.getY());
  dcData[hitnum].vtxz = double(vtxpos.getZ());

  //  G4cout<<"particle:parentid :"<<particle<<":"<<parentid<<G4endl;
  HitNumDC++;

  if(particle==321)//kaon
    HitNumDC_K++;

  if(particle==2212)//proton
    HitNumDC_p++;

  return;
}


void TPCAnaManager::SetSCHData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
				 G4int track, G4int particle, G4int detector, G4double mass,G4int qq,
			      G4int parentid,G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength)
{

  G4int hitnum = HitNumSCH;

  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }
  schData[hitnum].time = time;
  schData[hitnum].pos[ThreeVector::X] = pos.getX();
  schData[hitnum].pos[ThreeVector::Y] = pos.getY();
  schData[hitnum].pos[ThreeVector::Z] = pos.getZ();
  schData[hitnum].mom[ThreeVector::X] = mom.getX();
  schData[hitnum].mom[ThreeVector::Y] = mom.getY();
  schData[hitnum].mom[ThreeVector::Z] = mom.getZ();
  schData[hitnum].trackID = track;
  //  std::cout<<"SCH"<<std::endl;
  //  std::cout<<track<<std::endl;
  schData[hitnum].massSH = mass;
  schData[hitnum].qqSH = qq;
  schData[hitnum].particleID = particle;
  schData[hitnum].detectorID = detector;
  schData[hitnum].parentID = parentid;
  schData[hitnum].length = tlength;

  //kine E = sqrt(p^2+m^2)-m
  //p=sqrt((E+m)^2-m^2)
  //G4ThreeVecor vtxpos, G4ThreeVecor vtxmom, G4double vtxene
  G4double totalmom=sqrt(pow(vtxene+mass,2)-pow(mass,2));
  schData[hitnum].vtxpx = totalmom*double(vtxmom.getX());
  schData[hitnum].vtxpy = totalmom*double(vtxmom.getY());
  schData[hitnum].vtxpz = totalmom*double(vtxmom.getZ());
  schData[hitnum].vtxx = double(vtxpos.getX());
  schData[hitnum].vtxy = double(vtxpos.getY());
  schData[hitnum].vtxz = double(vtxpos.getZ());

  //  G4cout<<"particle:parentid :"<<particle<<":"<<parentid<<G4endl;
  HitNumSCH++;

  if(particle==321)//kaon
    HitNumSCH_K++;

  if(particle==2212)//proton
    HitNumSCH_p++;

  return;
}


void TPCAnaManager::SetFTOFData(G4double time, G4ThreeVector pos, G4ThreeVector mom,
				 G4int track, G4int particle, G4int detector, G4double mass,
				G4int qq,G4int parentid,G4ThreeVector vtxpos, G4ThreeVector vtxmom, G4double vtxene, G4double tlength)
{

  G4int hitnum = HitNumFTOF;

  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }
  ftofData[hitnum].time = time;
  ftofData[hitnum].pos[ThreeVector::X] = pos.getX();
  ftofData[hitnum].pos[ThreeVector::Y] = pos.getY();
  ftofData[hitnum].pos[ThreeVector::Z] = pos.getZ();
  ftofData[hitnum].mom[ThreeVector::X] = mom.getX();
  ftofData[hitnum].mom[ThreeVector::Y] = mom.getY();
  ftofData[hitnum].mom[ThreeVector::Z] = mom.getZ();
  ftofData[hitnum].trackID = track;
  //  std::cout<<"FTOF"<<std::endl;
  //  std::cout<<track<<std::endl;
  ftofData[hitnum].massSH = mass;
  ftofData[hitnum].qqSH = qq;
  ftofData[hitnum].particleID = particle;
  ftofData[hitnum].detectorID = detector;
  ftofData[hitnum].parentID = parentid;
  ftofData[hitnum].length = tlength;

  //kine E = sqrt(p^2+m^2)-m
  //p=sqrt((E+m)^2-m^2)
  //G4ThreeVecor vtxpos, G4ThreeVecor vtxmom, G4double vtxene
  G4double totalmom=sqrt(pow(vtxene+mass,2)-pow(mass,2));
  ftofData[hitnum].vtxpx = totalmom*double(vtxmom.getX());
  ftofData[hitnum].vtxpy = totalmom*double(vtxmom.getY());
  ftofData[hitnum].vtxpz = totalmom*double(vtxmom.getZ());
  ftofData[hitnum].vtxx = double(vtxpos.getX());
  ftofData[hitnum].vtxy = double(vtxpos.getY());
  ftofData[hitnum].vtxz = double(vtxpos.getZ());

  //  G4cout<<"particle:parentid :"<<particle<<":"<<parentid<<G4endl;
  HitNumFTOF++;

  if(particle==321)//kaon
    HitNumFTOF_K++;

  if(particle==2212)//proton
    HitNumFTOF_p++;

  return;
}

//_____________________________________________________________________________
// AnaManager->SetTargetData(tof, xyz, mom, tid, pid,mass,parentid);
void
TPCAnaManager::SetTargetData( G4int /* nhits */,
			      G4ThreeVector xyz, G4ThreeVector /* mom */,
			      G4int track, G4int particle,
			      G4int parentid,G4ThreeVector vtxpos,
			      G4ThreeVector /* vtxmom */, G4double )
{
  G4int hitnum = HitNumTarget;

  //  if (hitnum >= MaxTrack) {
  //    fprintf(stderr, "TPCAnaManager::SetTargetData Too Much multiplicity %d\n",
  //	    hitnum);
  //    return;
  //  }
  targetData[hitnum].targetparticleid = particle;
  targetData[hitnum].targetparentid = parentid;
  targetData[hitnum].targettrackid = track;
  targetData[hitnum].targetpos = xyz;
  targetData[hitnum].targetvtx = vtxpos;
  //  targetData[hitnum].targetvtxmom = vtxmom;
  //  G4cout<<targetData[hitnum].targetpos<<G4endl;
  //  G4cout<<hitnum<<":"<<nhits<<"th hits :"<<particle<<":"<<parentid<<":"<<xyz<<G4endl;
  //  G4cout<<mom<<G4endl;
  //  G4cout<<"kine energy:"<<vtx
  //  scintData[hitnum].scintvtxx = double(vtxpos.getX());
  //  scintData[hitnum].scintvtxy = double(vtxpos.getY());
  //  scintData[hitnum].scintvtxz = double(vtxpos.getZ());
  //G4cout<<"particle:parentid :"<<particle<<":"<<parentid<<G4endl;
  HitNumTarget++;

  if(particle==321)//kaon
    HitNumTarget_K++;

  if(particle==2212)//proton
    HitNumTarget_p++;

  return;
}



/*************************************
 *************************************/
void initTrack(Track* tracks){
  static const std::string funcname = "[InitTrack]";
  int i,j;
  //  G4cout<<"init track"<<G4endl;
  for( i = 0; i < MAX_TRACK; i++){
    tracks[i].nout   =  0;
    tracks[i].ngood  =  0;
    tracks[i].igroup = -1;
    tracks[i].trkQual = -1;
    tracks[i].numHits = 0;
    ////    tracks[i].numSectors = 0;
    tracks[i].numLayers = 0;
    tracks[i].charge = 1000;
    tracks[i].totalLength = 0.0;
    tracks[i].totalLengthTOF = 0.0; /*NTPC TOF*/
    tracks[i].meanAdc = -1.0;
    tracks[i].chi2Pad = -1.0;
    tracks[i].chi2Z       = -1.0;
    tracks[i].chi2Prob = -1.0;
    tracks[i].chi2 = -1.0;
    tracks[i].radius = 0.0;
    tracks[i].center[0] = tracks[i].center[1] = 1000;
    tracks[i].rKNumIter = -1;
    tracks[i].CrossOuter = -1;
    tracks[i].mom[0] = tracks[i].mom[1] =
      tracks[i].mom[2] = 10000.0;
    tracks[i].mom[3] = -1.0;
    tracks[i].resVirtual[0] = tracks[i].resVirtual[1]
      =tracks[i].resVirtual[2] = tracks[i].resVirtual[3] = -100;
    tracks[i].xOnTrack[0] = tracks[i].xOnTrack[1]
      = tracks[i].xOnTrack[2] = 1000.0;
    tracks[i].RKPFinal[0] =
      tracks[i].RKPFinal[1] = tracks[i].RKPFinal[2]  = -1000;
    //    for( j=0 ; j<NUM_SECTOR ; j++ ){
    //      tracks[i].numLayersinSec[j] = 0;
    //    }
    for(j=0; j < MAX_ITERATION; j++){
      tracks[i].rKChi2[j] = 1000.;
    }

    for(j=0; j < NUM_PARA_RK; j++){
      tracks[i].rKInitPara[j] =
        tracks[i].rKFinalPara[j] = -1000.;
    }

    /*  tracks[i].hitPattern = 0;*/
    for( j = 0; j < MAX_HIT_IN_TRACK; j++){
      tracks[i].ibad[j]  = 10;
      tracks[i].zbad[j]  = 10;
      //      tracks[i].hit[j]    = NULL;
      tracks[i].resPad[j] = 1000;
      tracks[i].resPady[j] = 1000; /*NTPC*/
      tracks[i].resZ[j]   = 1000;
      tracks[i].sector[j] = 100;
      tracks[i].lay[j] = 100;
      tracks[i].x[j][0] = 1000.;/*originally 100.*/
      tracks[i].x[j][1] = 1000.;/*originally 100.*/
      tracks[i].x[j][2] = -10000.;
      tracks[i].err[j][0] = tracks[i].err[j][1]
        = tracks[i].err[j][2] = 1000.;
      tracks[i].arcLen[j] = -1.0;
      tracks[i].resPad[j] = tracks[i].resPady[j] = tracks[i].resZ[j] =
        tracks[i].rKresXYZ[j][0] = tracks[i].rKresXYZ[j][1] =
        tracks[i].rKresXYZ[j][2] = tracks[i].rKresXYZ[j][3] =
        tracks[i].initRes[j][0] = tracks[i].initRes[j][1] =
        tracks[i].initRes[j][2] = tracks[i].phi_local[j]
        = -10000.;

    }
  }
}


void initTrack_ku(Track* tracks){
  static const std::string funcname = "[InitTrack]";
  int i,j;
  //  G4cout<<"init track"<<G4endl;
  for( i = 0; i < MAX_TRACK; i++){
    tracks[i].nout   =  0;
    tracks[i].ngood  =  0;
    tracks[i].igroup = -1;
    tracks[i].trkQual = -1;
    tracks[i].numHits = 0;
    ////    tracks[i].numSectors = 0;
    tracks[i].numLayers = 0;
    tracks[i].charge = 1000;
    tracks[i].totalLength = 0.0;
    tracks[i].totalLengthTOF = 0.0; /*NTPC TOF*/
    tracks[i].meanAdc = -1.0;
    tracks[i].chi2Pad = -1.0;
    tracks[i].chi2Z       = -1.0;
    tracks[i].chi2Prob = -1.0;
    tracks[i].chi2 = -1.0;
    tracks[i].radius = 0.0;
    tracks[i].center[0] = tracks[i].center[1] = 1000;
    tracks[i].rKNumIter = -1;
    tracks[i].CrossOuter = -1;
    tracks[i].mom[0] = tracks[i].mom[1] =
      tracks[i].mom[2] = 10000.0;
    tracks[i].mom[3] = -1.0;
    tracks[i].resVirtual[0] = tracks[i].resVirtual[1]
      =tracks[i].resVirtual[2] = tracks[i].resVirtual[3] = -100;
    tracks[i].xOnTrack[0] = tracks[i].xOnTrack[1]
      = tracks[i].xOnTrack[2] = 1000.0;
    tracks[i].RKPFinal[0] =
      tracks[i].RKPFinal[1] = tracks[i].RKPFinal[2]  = -1000;
    //    for( j=0 ; j<NUM_SECTOR ; j++ ){
    //      tracks[i].numLayersinSec[j] = 0;
    //    }
    for(j=0; j < MAX_ITERATION; j++){
      tracks[i].rKChi2[j] = 1000.;
    }

    for(j=0; j < NUM_PARA_RK; j++){
      tracks[i].rKInitPara[j] =
        tracks[i].rKFinalPara[j] = -1000.;
    }

    /*  tracks[i].hitPattern = 0;*/
    for( j = 0; j < MAX_HIT_IN_TRACK; j++){
      tracks[i].ibad[j]  = 10;
      tracks[i].zbad[j]  = 10;
      //      tracks[i].hit[j]    = NULL;
      tracks[i].resPad[j] = 1000;
      tracks[i].resPady[j] = 1000; /*NTPC*/
      tracks[i].resZ[j]   = 1000;
      tracks[i].sector[j] = 100;
      tracks[i].lay[j] = 100;
      tracks[i].x[j][0] = 1000.;/*originally 100.*/
      tracks[i].x[j][1] = 1000.;/*originally 100.*/
      tracks[i].x[j][2] = -10000.;
      tracks[i].res[j] = 0.;
      tracks[i].err[j][0] = tracks[i].err[j][1]
        = tracks[i].err[j][2] = 1000.;
      tracks[i].arcLen[j] = -1.0;
      tracks[i].resPad[j] = tracks[i].resPady[j] = tracks[i].resZ[j] =
        tracks[i].rKresXYZ[j][0] = tracks[i].rKresXYZ[j][1] =
        tracks[i].rKresXYZ[j][2] = tracks[i].rKresXYZ[j][3] =
        tracks[i].initRes[j][0] = tracks[i].initRes[j][1] =
        tracks[i].initRes[j][2] = tracks[i].phi_local[j]
        = -10000.;

    }
  }
}

/*************************************
 *************************************/
/*void setTrack(Track* tracks,int ntrk, double* x, double* y, double* z, double* ede, double* nhit){
  static const std::string funcname = "[SetTrack]";
  int i,j;
  for( i = 0; i < ntrk; i++){
    for( j = 0; j < nhit[i]; j++){
      tracks[i].lay[j] = 100;
      tracks[i].x[j][0] = x[j];
      tracks[i].x[j][1] = y[j];
      tracks[i].x[j][2] = z[j];
    }
  }
}
*/
/*************************************
 *************************************/

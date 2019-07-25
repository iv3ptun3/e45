// -*- C++ -*-

#include "TPCAnaManager.hh"

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TTree.h>

#include "ConfMan.hh"
#include "FuncName.hh"
#include "ResHypTPC.hh"
#include "RungeKuttaTracker.hh"
#include "minuit2.hh"
#include "switch.h"
#include "track.hh"
#include "VHitInfo.hh"

namespace
{
  const ConfMan& gConf = ConfMan::GetInstance();
  TTree* tree;
  Event event;
  std::map<TString, TH1*> hmap;
}

//_____________________________________________________________________________
TPCAnaManager::TPCAnaManager( void )
  : m_htof_hc()
{
  tree = new TTree( "tree", "GEANT4 simulation for HypTPC" );
  event.pb = new TVector3;
  tree->Branch( "evnum", &event.evnum, "evnum/I" );
  tree->Branch( "pb", "TVector3", event.pb );
  tree->Branch("npid",&event.npid,"npid/I");
  tree->Branch("pid0",event.pid0,"pid0[npid]/I");
  tree->Branch("x0",event.x0,"x0[npid][3]/D");
  tree->Branch("p0",event.p0,"p0[npid][5]/D");
  tree->Branch("pt0",event.pt0,"pt0[npid]/D");
  tree->Branch("mass0",event.mass0,"mass0[npid]/D");
  tree->Branch("theta0",event.theta0,"theta0[npid]/D");

  tree->Branch("mm_d",&event.mm_d,"mm_d/D");
  // tree->Branch("mm_p",&event.mm_p,"mm_p/D");
  tree->Branch("theta",&event.theta,"theta/D");
  tree->Branch("theta_scat",&event.theta_scat,"theta_scat/D");
  tree->Branch("theta_CM",&event.theta_CM,"theta_CM/D");

  // tree->Branch("mm",&event.mm,"mm/D");

  //generator mode
  tree->Branch("gen",&event.gen,"gen/I");
  tree->Branch("mode",&event.mode,"mode/I");

  //Num Of Hit for K+

  ///////shhwang tpc hit step

  // tree->Branch("nttpc",&event.nttpc,"nttpc/I");
  // tree->Branch("ntrk",event.ntrk,"ntrk[nttpc]/I");
  // tree->Branch("ititpc",event.ititpc,"ititpc[nttpc]/I");
  // tree->Branch("idtpc",event.idtpc,"idtpc[nttpc]/I");
  // tree->Branch("xtpc",event.xtpc,"xtpc[nttpc]/D");
  // tree->Branch("ytpc",event.ytpc,"ytpc[nttpc]/D");
  // tree->Branch("ztpc",event.ztpc,"ztpc[nttpc]/D");
  // tree->Branch("x0tpc",event.x0tpc,"x0tpc[nttpc]/D");
  // tree->Branch("y0tpc",event.y0tpc,"y0tpc[nttpc]/D");
  // tree->Branch("z0tpc",event.z0tpc,"z0tpc[nttpc]/D");
  // tree->Branch("resoX",event.resoX,"resoX[nttpc]/D");
  // tree->Branch("pxtpc",event.pxtpc,"pxtpc[nttpc]/D");
  // tree->Branch("pytpc",event.pytpc,"pytpc[nttpc]/D");
  // tree->Branch("pztpc",event.pztpc,"pztpc[nttpc]/D");
  // tree->Branch("pptpc",event.pptpc,"pptpc[nttpc]/D");   // total mometum
  // tree->Branch("masstpc",event.masstpc,"masstpc[nttpc]/D");   // mass TPC
  // //  tree->Branch("betatpc",event.betatpc,"betatpc[nttpc]/D");
  // // tree->Branch("edeptpc",event.edeptpc,"edeptpc[nttpc]/D");
  // //tree->Branch("dedxtpc",event.dedxtpc,"dedxtpc[nttpc]/D");
  // //tree->Branch("slengthtpc",event.slengthtpc,"slengthtpc[nttpc]/D");
  // tree->Branch("laytpc",event.laytpc,"laytpc[nttpc]/I");
  // tree->Branch("rowtpc",event.rowtpc,"rowtpc[nttpc]/I");
  // tree->Branch("parentID",event.parentID,"parentID[nttpc]/I");

  //// study on multiplicity
  // tree->Branch("nthlay",event.nthlay,"nthlay[nttpc]/I");
  // tree->Branch("nthpad",event.nthpad,"nthpad[nttpc]/I");
  // tree->Branch("laypad",event.laypad,"laytpadpc[nttpc][40][250]/I");


  //shhwang ntrtpc --> number of trak in tpc
  tree->Branch("ntrtpc",&event.ntrtpc,"ntrtpc/I");
  tree->Branch("trpmtpc",event.trpmtpc,"trpmtpc[ntrtpc]/D");
  tree->Branch("trqqtpc",event.trqqtpc,"trqqtpc[ntrtpc]/I");
  tree->Branch("trpidtpc",event.trpidtpc,"trpidtpc[ntrtpc]/I");
  tree->Branch("trparentidtpc",event.trparentidtpc,"trparentidtpc[ntrtpc]/I");
  //tree->Branch("trparentid_pid_tpc",event.trparentid_pid_tpc,"trparentid_pid_tpc[ntrtpc]/I");

  tree->Branch("trpxtpc",event.trpxtpc,"trpxtpc[ntrtpc]/D");
  tree->Branch("trpytpc",event.trpytpc,"trpytpc[ntrtpc]/D");
  tree->Branch("trpztpc",event.trpztpc,"trpztpc[ntrtpc]/D");
  tree->Branch("trpptpc",event.trpptpc,"trpptpc[ntrtpc]/D");
  tree->Branch("trpttpc",event.trpttpc,"trpttpc[ntrtpc]/D");

  tree->Branch("trpxtpcfit",event.trpxtpcfit,"trpxtpcfit[ntrtpc]/D");
  tree->Branch("trpytpcfit",event.trpytpcfit,"trpytpcfit[ntrtpc]/D");
  tree->Branch("trpztpcfit",event.trpztpcfit,"trpztpcfit[ntrtpc]/D");
  tree->Branch("trpptpcfit",event.trpptpcfit,"trpptpcfit[ntrtpc]/D");
  tree->Branch("trpttpcfit",event.trpttpcfit,"trpttpcfit[ntrtpc]/D");

  tree->Branch("trvtxpxtpc",event.trvtxpxtpc,"trvtxpxtpc[ntrtpc]/D");
  tree->Branch("trvtxpytpc",event.trvtxpytpc,"trvtxpytpc[ntrtpc]/D");
  tree->Branch("trvtxpztpc",event.trvtxpztpc,"trvtxpztpc[ntrtpc]/D");
  tree->Branch("trvtxpptpc",event.trvtxpptpc,"trvtxpptpc[ntrtpc]/D");

  tree->Branch("trvtxxtpc",event.trvtxxtpc,"trvtxxtpc[ntrtpc]/D");
  tree->Branch("trvtxytpc",event.trvtxytpc,"trvtxytpc[ntrtpc]/D");
  tree->Branch("trvtxztpc",event.trvtxztpc,"trvtxztpc[ntrtpc]/D");

  tree->Branch("trvtxxtpcfit",event.trvtxxtpcfit,"trvtxxtpcfit[ntrtpc]/D");
  tree->Branch("trvtxytpcfit",event.trvtxytpcfit,"trvtxytpcfit[ntrtpc]/D");
  tree->Branch("trvtxztpcfit",event.trvtxztpcfit,"trvtxztpcfit[ntrtpc]/D");

  // tree->Branch("trdetpc",event.trdetpc,"trdetpc[ntrtpc]/D");
  // tree->Branch("trlentpc",event.trlentpc,"trlentpc[ntrtpc]/D");
  // tree->Branch("trdedxtpc",event.trdedxtpc,"trdedxtpc[ntrtpc]/D");
  // tree->Branch("trdedxtrtpc",event.trdedxtrtpc,"trdedxtrtpc[ntrtpc]/D");
  tree->Branch("trlaytpc",event.trlaytpc,"trlaytpc[ntrtpc]/I");
  tree->Branch("cir_r",event.cir_r,"cir_r[ntrtpc]/D");
  tree->Branch("cir_x",event.cir_x,"cir_x[ntrtpc]/D");
  tree->Branch("cir_z",event.cir_z,"cir_z[ntrtpc]/D");
  tree->Branch("cir_fit",event.cir_fit,"cir_fit[ntrtpc]/D");
  tree->Branch("vtx_flag",event.vtx_flag,"vtx_flag[ntrtpc]/I");
  tree->Branch("a_fory",event.a_fory,"a_fory[ntrtpc]/D");
  tree->Branch("b_fory",event.b_fory,"b_fory[ntrtpc]/D");

  //
  // Scinti
  //

  tree->Branch("ntsc",&event.ntsc,"ntsc/I");
  tree->Branch("tidsc",event.tidsc,"tidsc[ntsc]/I");
  tree->Branch("pidsc",event.pidsc,"pidsc[ntsc]/I");
  tree->Branch("didsc",event.didsc,"didsc[ntsc]/I");
  tree->Branch("masssc",event.masssc,"masssc[ntsc]/D");
  tree->Branch("qqsc",event.qqsc,"qqsc[ntsc]/I");
  tree->Branch("xsc",event.xsc,"xsc[ntsc]/D");
  tree->Branch("ysc",event.ysc,"ysc[ntsc]/D");
  tree->Branch("zsc",event.zsc,"zsc[ntsc]/D");
  tree->Branch("pxsc",event.pxsc,"pxsc[ntsc]/D");
  tree->Branch("pysc",event.pysc,"pysc[ntsc]/D");
  tree->Branch("pzsc",event.pzsc,"pzsc[ntsc]/D");
  tree->Branch("ppsc",event.ppsc,"ppsc[ntsc]/D");
  tree->Branch("tofsc",event.tofsc,"tofsc[ntsc]/D");
  tree->Branch("scpID",event.scpID,"scpID[ntsc]/I");

  tree->Branch("trvtxpxscint",event.trvtxpxscint,"trvtxpxscint[ntsc]/D");
  tree->Branch("trvtxpyscint",event.trvtxpyscint,"trvtxpyscint[ntsc]/D");
  tree->Branch("trvtxpzscint",event.trvtxpzscint,"trvtxpzscint[ntsc]/D");
  tree->Branch("trvtxppscint",event.trvtxppscint,"trvtxppscint[ntsc]/D");

  tree->Branch("trvtxxscint",event.trvtxxscint,"trvtxxscint[ntsc]/D");
  tree->Branch("trvtxyscint",event.trvtxyscint,"trvtxyscint[ntsc]/D");
  tree->Branch("trvtxzscint",event.trvtxzscint,"trvtxzscint[ntsc]/D");
  tree->Branch("lengthsc",event.lengthsc,"lengthsc[ntsc]/D");

  //target
  // tree->Branch("targethits",&event.targethits,"targethits/I");
  // tree->Branch("targetpid",event.targetpid,"targetpid[targethits]/I");
  // tree->Branch("targetparentid",event.targetparentid,"targetparentid[targethits]/I");
  // tree->Branch("targettid",event.targettid,"targettid[targethits]/I");
  // tree->Branch("targetpos",event.targetpos,"targetpos[targethits][3]/D");
  // tree->Branch("targetvtx",event.targetvtx,"targetvtx[targethits][3]/D");

  // HTOF
  tree->Branch("nhHtof",&event.nhHtof,"nhHtof/I");
  tree->Branch("tidHtof",event.tidHtof,"tidHtof[nhHtof]/I");
  tree->Branch("pidHtof",event.pidHtof,"pidHtof[nhHtof]/I");
  tree->Branch("didHtof",event.didHtof,"didHtof[nhHtof]/I");
  tree->Branch("massHtof",event.massHtof,"massHtof[nhHtof]/D");
  tree->Branch("qqHtof",event.qqHtof,"qqHtof[nhHtof]/I");
  tree->Branch("xHtof",event.xHtof,"xHtof[nhHtof]/D");
  tree->Branch("yHtof",event.yHtof,"yHtof[nhHtof]/D");
  tree->Branch("zHtof",event.zHtof,"zHtof[nhHtof]/D");
  tree->Branch("pxHtof",event.pxHtof,"pxHtof[nhHtof]/D");
  tree->Branch("pyHtof",event.pyHtof,"pyHtof[nhHtof]/D");
  tree->Branch("pzHtof",event.pzHtof,"pzHtof[nhHtof]/D");
  tree->Branch("ppHtof",event.ppHtof,"ppHtof[nhHtof]/D");
  tree->Branch("tofHtof",event.tofHtof,"tofHtof[nhHtof]/D");
  tree->Branch("HtofpID",event.HtofpID,"HtofpID[nhHtof]/I");
  tree->Branch("trvtxpxHtof",event.trvtxpxHtof,"trvtxpxHtof[nhHtof]/D");
  tree->Branch("trvtxpyHtof",event.trvtxpyHtof,"trvtxpyHtof[nhHtof]/D");
  tree->Branch("trvtxpzHtof",event.trvtxpzHtof,"trvtxpzHtof[nhHtof]/D");
  tree->Branch("trvtxppHtof",event.trvtxppHtof,"trvtxppHtof[nhHtof]/D");
  tree->Branch("trvtxxHtof",event.trvtxxHtof,"trvtxxHtof[nhHtof]/D");
  tree->Branch("trvtxyHtof",event.trvtxyHtof,"trvtxyHtof[nhHtof]/D");
  tree->Branch("trvtxzHtof",event.trvtxzHtof,"trvtxzHtof[nhHtof]/D");
  tree->Branch("lengthHtof",event.lengthHtof,"lengthHtof[nhHtof]/D");
}

//_____________________________________________________________________________
TPCAnaManager::~TPCAnaManager( void )
{
}

//_____________________________________________________________________________
void
TPCAnaManager::BeginOfRunAction( int runnum )
{
  event.evnum = 0;

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

  for( auto& p : hmap ){
    delete p.second;
  }
  hmap.clear();
  TString key;
  key = "Time";
  hmap[key] = new TH1D( key, key, 400, 0.0, 10.0 );
  for( G4int i=0; i<G4ThreeVector::SIZE; ++i ){
    key = Form( "Pos%d", i );
    hmap[key] = new TH1D( key, key, 500, -25.0*CLHEP::cm, 25.0*CLHEP::cm );
    key = Form( "Mom%d", i );
    if( i==2 )
      hmap[key] = new TH1D( key, key, 400, 0.0, 2.0 );
    else
      hmap[key] = new TH1D( key, key, 400, -2.0, 2. );
  }
}

//_____________________________________________________________________________
void
TPCAnaManager::EndOfRunAction( void )
{
  tree->Write();
  for( auto& p : hmap ){
    p.second->Write();
  }
}

//_____________________________________________________________________________
void
TPCAnaManager::BeginOfEventAction( void )
{
  m_htof_hc.clear();

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

  event.nhHtof = 0;
  for( G4int i=0; i<MaxHits; ++i ){
    event.xHtof[i] = -9999.;
    event.yHtof[i] = -9999.;
    event.zHtof[i] = -9999.;
    event.pxHtof[i] = -9999.;
    event.pyHtof[i] = -9999.;
    event.pzHtof[i] = -9999.;
    event.tofHtof[i] = -9999.;
    event.HtofpID[i] = -9999;
    event.tidHtof[i] = -1;
    event.pidHtof[i] = -1;
    event.didHtof[i] = -1;
    event.massHtof[i] = -1;
    event.qqHtof[i] = -1;
    event.trvtxppHtof[i] = -9999.;
    event.trvtxpxHtof[i] = -9999.;
    event.trvtxpyHtof[i] = -9999.;
    event.trvtxpzHtof[i] = -9999.;
    event.trvtxxHtof[i] = -9999.;
    event.trvtxyHtof[i] = -9999.;
    event.trvtxzHtof[i] = -9999.;
  }

  event.nttpc = 0;
  event.ntac = 0;
  event.ntnbar = 0;
  event.ntsc = 0;
  event.ntdc = 0;
  event.ntsch = 0;
  event.ntftof = 0;
  event.npid = 0;
  event.ntrtpc = 0;
  event.targethits = 0;

  event.gen = 0;
  event.mode = 0;

  event.HitNum_K=-1;
  event.HitNumAC_K=-1;
  event.HitNumNBAR_K=-1;
  event.HitNumDC_K=-1;
  event.HitNumFTOF_K=-1;
  event.HitNumSCH_K=-1;
  event.HitNumScint_K=-1;
  event.HitNumTarget_K=-1;


  event.HitNum_p=-1;
  event.HitNumAC_p=-1;
  event.HitNumNBAR_p=-1;
  event.HitNumDC_p=-1;
  event.HitNumFTOF_p=-1;
  event.HitNumSCH_p=-1;
  event.HitNumScint_p=-1;
  event.HitNumTarget_p=-1;


  event.mm_d = 0.;
  event.mm_p = 0.;
  event.theta = 0.;
  event.theta_scat = 0.;
  event.theta_CM = 0.;


  /* ntrtpc initialization */

  for( G4int i=0; i<MaxHitsTPC;++i ){
    event.trpidtpc[i]  = -1;
    event.trparentidtpc[i]  = -1;
    event.trparentid_pid_tpc[i]  = -1;


    event.trpptpc[i]  = -9999.9999;
    event.trpttpc[i]  = -9999.9999;
    event.trpxtpc[i]  = -9999.9999;
    event.trpytpc[i]  = -9999.9999;
    event.trpztpc[i]  = -9999.9999;

    event.trvtxpxtpc[i]  = -9999.9999;
    event.trvtxpytpc[i]  = -9999.9999;
    event.trvtxpztpc[i]  = -9999.9999;

    event.trvtxxtpc[i]  = -9999.9999;
    event.trvtxytpc[i]  = -9999.9999;
    event.trvtxztpc[i]  = -9999.9999;

    event.trpttpcfit[i]  = -9999.9999;

    event.trpptpcfit[i]  = -9999.9999;
    event.trpxtpcfit[i]  = -9999.9999;
    event.trpytpcfit[i]  = -9999.9999;
    event.trpztpcfit[i]  = -9999.9999;

    event.trpmtpc[i]  = -9999.9999;
    event.trqqtpc[i]  = -9999;

    event.trdetpc[i]  = -9999.9999;
    event.trlentpc[i]  = -9999.9999;
    event.trdedxtpc[i]  = -9999.9999;
    event.trdedxtrtpc[i]  = -9999.9999;
    event.trlaytpc[i]  = -9999;

    event.cir_r[i]  = -9999.9999;
    event.cir_x[i]  = -9999.9999;
    event.cir_z[i]  = -9999.9999;
    event.cir_fit[i]  = -9999.9999;

    event.vtx_flag[i]  = -1;
    event.a_fory[i]  = -9999.9999;
    event.b_fory[i]  = -9999.9999;
  }

  for(int i = 0; i< MaxPrimaryParticle;i++){
    for(int j=0;j<3;j++){
      event.x0[i][j] = -9999.9;
    }
    for(int j=0;j<5;j++){
      event.p0[i][j] = -9999.9;
    }
    event.theta0[i] = -9999.9;
    event.pt0[i] = -9999.9;
    event.mass0[i] = -9999.9;
    event.pid0[i] = -9999;
  }


  for(int i=0;i<MaxTrack;i++){

    /// initialization pad multiplicity
    event.nthlay[i]=-9999.;
    event.nthpad[i]=-9999.;
    for(int j = 0; j< MaxNthLay;j++){
      for(int k = 0; k< MaxNthPad;k++){
	event.laypad[i][j][k]  = 0.;
      }
    }
    //////////////


    event.xtpc[i] = -9999.9;
    event.ytpc[i] = -9999.9;
    event.ztpc[i] = -9999.9;

    event.x0tpc[i] = -9999.9;
    event.y0tpc[i] = -9999.9;
    event.z0tpc[i] = -9999.9;
    event.resoX[i] = -9999.9;

    event.pxtpc[i] = -9999.9;
    event.pytpc[i] = -9999.9;
    event.pztpc[i] = -9999.9;
    event.pptpc[i] = -9999.9;

    event.masstpc[i] = -9999.9;

    event.betatpc[i] = -9999.9;

    event.edeptpc[i] = -9999.9;

    event.ititpc[i] = -1;
    event.idtpc[i] = -1;
    event.laytpc[i] = -1;
    event.rowtpc[i] = -1;
    event.parentID[i] = -1;
  }


  for(int i=0; i<MaxHits; i++){
    event.xsc[i] = -9999.9;
    event.ysc[i] = -9999.9;
    event.zsc[i] = -9999.9;
    event.pxsc[i] = -9999.9;
    event.pysc[i] = -9999.9;
    event.pzsc[i] = -9999.9;
    event.tofsc[i] = -9999.9;
    event.scpID[i] = -9999;

    event.tidsc[i] = -1;
    event.pidsc[i] = -1;
    event.didsc[i] = -1;
    event.masssc[i] = -1;
    event.qqsc[i] = -1;

    event.trvtxppscint[i]  = -9999.9999;
    event.trvtxpxscint[i]  = -9999.9999;
    event.trvtxpyscint[i]  = -9999.9999;
    event.trvtxpzscint[i]  = -9999.9999;

    event.trvtxxscint[i]  = -9999.9999;
    event.trvtxyscint[i]  = -9999.9999;
    event.trvtxzscint[i]  = -9999.9999;
  }
  ///ac
  for(int i=0; i<MaxHits; i++){
    event.xac[i] = -9999.9;
    event.yac[i] = -9999.9;
    event.zac[i] = -9999.9;
    event.pxac[i] = -9999.9;
    event.pyac[i] = -9999.9;
    event.pzac[i] = -9999.9;
    event.tofac[i] = -9999.9;
    event.acpID[i] = -9999;

    event.tidac[i] = -1;
    event.pidac[i] = -1;
    event.didac[i] = -1;
    event.massac[i] = -1;
    event.qqac[i] = -1;

    event.trvtxppac[i]  = -9999.9999;
    event.trvtxpxac[i]  = -9999.9999;
    event.trvtxpyac[i]  = -9999.9999;
    event.trvtxpzac[i]  = -9999.9999;

    event.trvtxxac[i]  = -9999.9999;
    event.trvtxyac[i]  = -9999.9999;
    event.trvtxzac[i]  = -9999.9999;
  }

  ///dc
  for(int i=0; i<MaxHits; i++){
    event.xdc[i] = -9999.9;
    event.ydc[i] = -9999.9;
    event.zdc[i] = -9999.9;
    event.pxdc[i] = -9999.9;
    event.pydc[i] = -9999.9;
    event.pzdc[i] = -9999.9;
    event.tofdc[i] = -9999.9;
    event.dcpID[i] = -9999;

    event.tiddc[i] = -1;
    event.piddc[i] = -1;
    event.diddc[i] = -1;
    event.massdc[i] = -1;
    event.qqdc[i] = -1;

    event.trvtxppdc[i]  = -9999.9999;
    event.trvtxpxdc[i]  = -9999.9999;
    event.trvtxpydc[i]  = -9999.9999;
    event.trvtxpzdc[i]  = -9999.9999;

    event.trvtxxdc[i]  = -9999.9999;
    event.trvtxydc[i]  = -9999.9999;
    event.trvtxzdc[i]  = -9999.9999;
  }
  // SCH
  for(int i=0; i<MaxHits; i++){
    event.xsch[i] = -9999.9;
    event.ysch[i] = -9999.9;
    event.zsch[i] = -9999.9;
    event.pxsch[i] = -9999.9;
    event.pysch[i] = -9999.9;
    event.pzsch[i] = -9999.9;
    event.tofsch[i] = -9999.9;
    event.schpID[i] = -9999;

    event.tidsch[i] = -1;
    event.pidsch[i] = -1;
    event.didsch[i] = -1;
    event.masssch[i] = -1;
    event.qqsch[i] = -1;

    event.trvtxppsch[i]  = -9999.9999;
    event.trvtxpxsch[i]  = -9999.9999;
    event.trvtxpysch[i]  = -9999.9999;
    event.trvtxpzsch[i]  = -9999.9999;

    event.trvtxxsch[i]  = -9999.9999;
    event.trvtxysch[i]  = -9999.9999;
    event.trvtxzsch[i]  = -9999.9999;
  }

  ///ftof
  for(int i=0; i<MaxHits; i++){
    event.xftof[i] = -9999.9;
    event.yftof[i] = -9999.9;
    event.zftof[i] = -9999.9;
    event.pxftof[i] = -9999.9;
    event.pyftof[i] = -9999.9;
    event.pzftof[i] = -9999.9;
    event.tofftof[i] = -9999.9;
    event.ftofpID[i] = -9999;

    event.tidftof[i] = -1;
    event.pidftof[i] = -1;
    event.didftof[i] = -1;
    event.massftof[i] = -1;
    event.qqftof[i] = -1;

    event.trvtxppftof[i]  = -9999.9999;
    event.trvtxpxftof[i]  = -9999.9999;
    event.trvtxpyftof[i]  = -9999.9999;
    event.trvtxpzftof[i]  = -9999.9999;

    event.trvtxxftof[i]  = -9999.9999;
    event.trvtxyftof[i]  = -9999.9999;
    event.trvtxzftof[i]  = -9999.9999;
  }

  /*
    for(int i=0; i<MaxTrackFDC; i++){
    event.xfdc[i] = -9999.9;
    event.yfdc[i] = -9999.9;
    event.zfdc[i] = -9999.9;
    event.pxfdc[i] = -9999.9;
    event.pyfdc[i] = -9999.9;
    event.pzfdc[i] = -9999.9;
    event.toffdc[i] = -9999.9;

    event.tidfdc[i] = -1;
    event.pidfdc[i] = -1;
    event.didfdc[i] = -1;
    }
  */

  for(int i=0; i<MaxTrack; i++){
    //    event.targethits = -9999;
    event.targetpid[i]=-9999;
    event.targettid[i]=-9999;
    event.targetparentid[i]=-9999;
    event.targetpos[i][0]=-9999.9999;
    event.targetpos[i][1]=-9999.9999;
    event.targetpos[i][2]=-9999.9999;

    event.targetvtx[i][0]=-9999.9999;
    event.targetvtx[i][1]=-9999.9999;
    event.targetvtx[i][2]=-9999.9999;
  }
}

//_____________________________________________________________________________
int
TPCAnaManager::EndOfEventAction( void )
{
  event.evnum++;
  // HTOF
  event.nhHtof = m_htof_hc.size();
  for( G4int i=0, n=m_htof_hc.size(); i<n; ++i ){
    event.tofHtof[i] = m_htof_hc[i]->GetTime();
    event.xHtof[i] = m_htof_hc[i]->GetPosition()[0];
    event.yHtof[i] = m_htof_hc[i]->GetPosition()[1];
    event.zHtof[i] = m_htof_hc[i]->GetPosition()[2];
    event.pidHtof[i] = m_htof_hc[i]->GetParticleID();
  }

  {
    event.gen = primaryBeam.gen;
    event.mode = primaryBeam.mode;
  }

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
    G4double* x0 = &primaryParticle.x0[id][0];
    G4double* p0 = &primaryParticle.p0[id][0];
    G4int pid = primaryParticle.pid0[id];
    for(int i=0;i<3;i++){
      event.x0[id][i] = x0[i];
    }

    for(int i=0;i<4;i++){
      event.p0[id][i] = p0[i];
    }
    event.pid0[id] = pid;

    event.p0[id][4] = sqrt(pow(p0[0],2.0)+pow(p0[1],2.0)+pow(p0[2],2.0));
    event.pt0[id] = sqrt(pow(p0[0],2.0)+pow(p0[1],2.0));
    event.mass0[id] = p0[4];
    event.theta0[id] =
      180.0*atan2(sqrt(pow(p0[0],2.0)+pow(p0[1],2.0)),p0[2])/CLHEP::pi;
    event.npid++;
  }

  //Fill Primary Infomation for E27
  if( env_Experiment_num ==27 || env_Experiment_num ==45 ){
    event.mm_d = primaryInfo.mm_d;
    event.mm_p = primaryInfo.mm_p;
    event.theta = primaryInfo.theta;
    event.theta_scat = primaryInfo.theta_scat;
    event.theta_CM = primaryInfo.theta_CM;
    event.mm = CLHEP::mm;
  }

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
      event.ntrk[event.nttpc] = counterData[i].ntrk;
      hmap["Time"]->Fill( counterData[i].time );
      for( G4int j=0; j<G4ThreeVector::SIZE; ++j ){
	hmap[Form( "Pos%d", j )]->Fill( counterData[i].pos[j]/CLHEP::mm );
	hmap[Form( "Mom%d", j )]->Fill( counterData[i].mom[j]/CLHEP::GeV );
      }
      event.xtpc[event.nttpc] = counterData[i].pos[0]/CLHEP::mm;
      event.ytpc[event.nttpc] = counterData[i].pos[1]/CLHEP::mm;
      event.ztpc[event.nttpc] = counterData[i].pos[2]/CLHEP::mm;
      event.x0tpc[event.nttpc] = counterData[i].pos0[0]/CLHEP::mm;
      event.y0tpc[event.nttpc] = counterData[i].pos0[1]/CLHEP::mm;
      event.z0tpc[event.nttpc] = counterData[i].pos0[2]/CLHEP::mm;
      event.resoX[event.nttpc] = counterData[i].resoX;
      event.pxtpc[event.nttpc] = counterData[i].mom[0]/CLHEP::GeV;
      event.pytpc[event.nttpc] = counterData[i].mom[1]/CLHEP::GeV;
      event.pztpc[event.nttpc] = counterData[i].mom[2]/CLHEP::GeV;
      event.pptpc[event.nttpc] = sqrt(pow(counterData[i].mom[0], 2) +
				      pow(counterData[i].mom[1], 2) +
				      pow(counterData[i].mom[2], 2))/CLHEP::GeV;
      event.ititpc[event.nttpc] = counterData[i].trackID;
      event.idtpc[event.nttpc] = counterData[i].particleID;
      event.laytpc[event.nttpc] = counterData[i].iLay;
      event.rowtpc[event.nttpc] = counterData[i].iRow;
      event.betatpc[event.nttpc] = counterData[i].beta;
      event.edeptpc[event.nttpc] = counterData[i].edep;
      event.dedxtpc[event.nttpc] = counterData[i].dedx;
      event.slengthtpc[event.nttpc] = counterData[i].slength;
      event.nthlay[event.nttpc] = counterData[i].iLay;
      event.nthpad[event.nttpc] = counterData[i].iPad;
      event.laypad[event.nttpc][event.nthlay[event.nttpc]][event.nthpad[event.nttpc]]
	= event.laypad[event.nttpc][event.nthlay[event.nttpc]][event.nthpad[event.nttpc]]+1.;
      event.nttpc += 1;
    }

    //
    // Scinti.
    //

    for( G4int i=0; i<HitNumScint; i++){
      event.tofsc[event.ntsc] = scintData[i].time;
      event.scpID[event.ntsc] = scintData[i].parentID;
      event.xsc[event.ntsc] = scintData[i].pos[0];
      event.ysc[event.ntsc] = scintData[i].pos[1];
      event.zsc[event.ntsc] = scintData[i].pos[2];
      event.pxsc[event.ntsc] = scintData[i].mom[0]/CLHEP::GeV;
      event.pysc[event.ntsc] = scintData[i].mom[1]/CLHEP::GeV;
      event.pzsc[event.ntsc] = scintData[i].mom[2]/CLHEP::GeV;
      event.ppsc[event.ntsc] = sqrt(pow(scintData[i].mom[0]/CLHEP::GeV,2)+
				    pow(scintData[i].mom[1]/CLHEP::GeV,2)+
				    pow(scintData[i].mom[2]/CLHEP::GeV,2));
      event.tidsc[event.ntsc] = scintData[i].trackID;
      event.pidsc[event.ntsc] = scintData[i].particleID;
      event.masssc[event.ntsc] = scintData[i].massSH/CLHEP::GeV;
      event.qqsc[event.ntsc] = scintData[i].qqSH;
      event.didsc[event.ntsc] = scintData[i].detectorID;
      event.trvtxpxscint[event.ntsc] = scintData[i].scintvtxpx/CLHEP::GeV;
      event.trvtxpyscint[event.ntsc] = scintData[i].scintvtxpy/CLHEP::GeV;
      event.trvtxpzscint[event.ntsc] = scintData[i].scintvtxpz/CLHEP::GeV;
      event.trvtxppscint[event.ntsc] = sqrt(pow(scintData[i].scintvtxpx,2)+
					    pow(scintData[i].scintvtxpy,2)+
					    pow(scintData[i].scintvtxpz,2))/CLHEP::GeV;
      event.trvtxxscint[event.ntsc] = scintData[i].scintvtxx;
      event.trvtxyscint[event.ntsc] = scintData[i].scintvtxy;
      event.trvtxzscint[event.ntsc] = scintData[i].scintvtxz;
      //  event.lengthsc[event.ntsc] = tlength+CLHEP::RandGauss::shoot(0.,20.);
      event.lengthsc[event.ntsc] = scintData[i].length;
      event.ntsc++;
    }

    //
    // AC
    //

    for( G4int i=0; i<HitNumAC; i++){
      event.tofac[event.ntac] = acData[i].time;
      event.acpID[event.ntac] = acData[i].parentID;
      event.xac[event.ntac] = acData[i].pos[0];
      event.yac[event.ntac] = acData[i].pos[1];
      event.zac[event.ntac] = acData[i].pos[2];
      event.pxac[event.ntac] = acData[i].mom[0]/CLHEP::GeV;
      event.pyac[event.ntac] = acData[i].mom[1]/CLHEP::GeV;
      event.pzac[event.ntac] = acData[i].mom[2]/CLHEP::GeV;
      event.ppac[event.ntac] = sqrt(pow(acData[i].mom[0]/CLHEP::GeV,2)+
				    pow(acData[i].mom[1]/CLHEP::GeV,2)+
				    pow(acData[i].mom[2]/CLHEP::GeV,2));
      event.tidac[event.ntac] = acData[i].trackID;
      event.pidac[event.ntac] = acData[i].particleID;
      event.massac[event.ntac] = acData[i].massSH/CLHEP::GeV;
      event.qqac[event.ntac] = acData[i].qqSH;
      event.didac[event.ntac] = acData[i].detectorID;
      event.trvtxpxac[event.ntac] = acData[i].acvtxpx/CLHEP::GeV;
      event.trvtxpyac[event.ntac] = acData[i].acvtxpy/CLHEP::GeV;
      event.trvtxpzac[event.ntac] = acData[i].acvtxpz/CLHEP::GeV;
      event.trvtxppac[event.ntac] = sqrt(pow(acData[i].acvtxpx,2)+
					 pow(acData[i].acvtxpy,2)+
					 pow(acData[i].acvtxpz,2))/CLHEP::GeV;
      event.trvtxxac[event.ntac] = acData[i].acvtxx;
      event.trvtxyac[event.ntac] = acData[i].acvtxy;
      event.trvtxzac[event.ntac] = acData[i].acvtxz;
      event.lengthac[event.ntac] = acData[i].length;// +CLHEP::RandGauss::shoot(0.,20.);
      event.ntac++;
    }

    //
    // nbar
    //

    // for( G4int i=0; i<HitNumNBAR; i++){
    //   anaRoot.FillNBARData(nbarData[i].time, nbarData[i].pos,
    // 			   nbarData[i].mom,
    // 			   nbarData[i].trackID, nbarData[i].particleID,
    // 			   nbarData[i].detectorID,nbarData[i].massSH,nbarData[i].qqSH,nbarData[i].parentID,
    // 			   nbarData[i].nbarvtxpx,nbarData[i].nbarvtxpy,nbarData[i].nbarvtxpz,
    // 			   nbarData[i].nbarvtxx,nbarData[i].nbarvtxy,nbarData[i].nbarvtxz,
    // 			   nbarData[i].length
    // 			   );

    // }

    //
    // DC
    //

    // for( G4int i=0; i<HitNumDC; i++){
    //   anaRoot.FillDCData(dcData[i].time, dcData[i].pos,
    // 			 dcData[i].mom,
    // 			 dcData[i].trackID, dcData[i].particleID,
    // 			 dcData[i].detectorID,dcData[i].massSH,dcData[i].qqSH,dcData[i].parentID,
    // 			 dcData[i].vtxpx,dcData[i].vtxpy,dcData[i].vtxpz,
    // 			 dcData[i].vtxx,dcData[i].vtxy,dcData[i].vtxz,
    // 			 dcData[i].length
    // 			 );

    // }

    //
    // SCH
    //

    // for( G4int i=0; i<HitNumSCH; i++){
    //   anaRoot.FillSCHData(schData[i].time, schData[i].pos,
    // 			 schData[i].mom,
    // 			 schData[i].trackID, schData[i].particleID,
    // 			 schData[i].detectorID,schData[i].massSH,schData[i].qqSH,schData[i].parentID,
    // 			 schData[i].vtxpx,schData[i].vtxpy,schData[i].vtxpz,
    // 			 schData[i].vtxx,schData[i].vtxy,schData[i].vtxz,
    // 			 schData[i].length
    // 			 );

    // }

    //
    // ftof
    //

    //    std::cout<<"tof:"<<HitNumFTOF<<std::endl;
    //    std::cout<<"dc:"<<HitNumDC<<std::endl;
    //    std::cout<<"SCH:"<<HitNumSCH<<std::endl;

    // for( G4int i=0; i<HitNumFTOF; i++){
    //   anaRoot.FillFTOFData(ftofData[i].time, ftofData[i].pos,
    // 			   ftofData[i].mom,
    // 			   ftofData[i].trackID, ftofData[i].particleID,
    // 			   ftofData[i].detectorID,ftofData[i].massSH,ftofData[i].qqSH,ftofData[i].parentID,
    // 			   ftofData[i].vtxpx,ftofData[i].vtxpy,ftofData[i].vtxpz,
    // 			   ftofData[i].vtxx,ftofData[i].vtxy,ftofData[i].vtxz,
    // 			   ftofData[i].length
    // 			   );

    // }

    //
    // Target.
    //

    // for( G4int i=0; i<HitNumTarget; i++){
    //   anaRoot.FillTargetData(i,
    // 			     targetData[i].targetparticleid,
    // 			     targetData[i].targetparentid,
    // 			     targetData[i].targettrackid,
    // 			     targetData[i].targetpos,
    // 			     targetData[i].targetvtx
    // 			     );
    // }

    //
    // TPC
    //
    // for( G4int i=0; i<tpctrNum; i++){
    //   //      G4cout<<"abs"<<abs(env_helm_field)<<G4endl;
    //   //      G4cout<<"fabs"<<fabs(env_helm_field)<<G4endl;
    //   anaRoot.FillTPCData(tpcData[i].tpcpx,
    // 			  tpcData[i].tpcpy,tpcData[i].tpcpz,
    // 			  tpcData[i].tpcpp,
    // 			  tpcData[i].tpcpid, tpcData[i].tpcparentid, tpcData[i].tpcparentid_pid,
    // 			  tpcData[i].tpcqq,
    // 			  tpcData[i].tpcpm,tpcData[i].tpcde,
    // 			  tpcData[i].tpclen,mean[i],trmean[i],
    // 			  tpcData[i].tpclay,
    // 			  tpcData[i].tpcvtxpx,tpcData[i].tpcvtxpy,tpcData[i].tpcvtxpz,
    // 			  tpcData[i].tpcvtxx,tpcData[i].tpcvtxy,tpcData[i].tpcvtxz,
    // 			  vtxxfit[i],vtxyfit[i],vtxzfit[i],
    // 			  vtxpxfit[i],Pz[i],vtxpzfit[i],cir_r[i]*(0.299792458)*fabs(env_helm_field),
    // 			  cir_r[i],cir_x[i],cir_z[i],test[i],
    // 			  vtx_flag[i], a_fory[i], b_fory[i]
    // 			  );
    // }

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
  }//trigger parts

  tree->Fill();

  event.pb->SetMag( 0. );
  return 0;
}

//_____________________________________________________________________________
void
TPCAnaManager::SetCounterData( G4int ntrk,G4double time, G4ThreeVector pos,
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

  //  G4ThreeVector tar_pos(0.,0.*CLHEP::mm,-150.*CLHEP::mm);

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
    //    G4double sh_dalpha = 0.300*CLHEP::mm/sh_rho; // rho * theta = arc --> sigma=300 um
    //    G4double sh_smear_alpha = CLHEP::RandGauss::shoot(sh_alpha, sh_dalpha);
    G4double sh_sigmaY = 0.500*CLHEP::mm; //--> smearing : 400 um

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

    counterData[hitnum].pos[G4ThreeVector::Z] = sh_rho*cos(sh_smear_alpha)+tar_pos.getZ();
    counterData[hitnum].pos[G4ThreeVector::X] = sh_rho*sin(sh_smear_alpha);
    counterData[hitnum].pos[G4ThreeVector::Y] = CLHEP::RandGauss::shoot(sh_y,sh_sigmaY);

    counterData[hitnum].pos0[G4ThreeVector::X] = pos.getX();
    counterData[hitnum].pos0[G4ThreeVector::Y] = pos.getY();
    counterData[hitnum].pos0[G4ThreeVector::Z] = pos.getZ();

    counterData[hitnum].mom[G4ThreeVector::X] = mom.getX();
    counterData[hitnum].mom[G4ThreeVector::Y] = mom.getY();
    counterData[hitnum].mom[G4ThreeVector::Z] = mom.getZ();

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
TPCAnaManager::SetHTOFData( VHitInfo* hit )
{
  if( m_htof_hc.size() > MaxHits )
    G4cerr << FUNC_NAME << " too much nhit" << m_htof_hc.size() << G4endl;
  else
    m_htof_hc.push_back( hit );
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

//_____________________________________________________________________________
void
TPCAnaManager::SetPrimaryInfo( G4double mm_d, G4double mm_p,
			       G4double theta, G4double theta_scat,
			       G4double theta_CM )
{
  primaryInfo.mm_d = mm_d;
  primaryInfo.mm_p = mm_p;
  primaryInfo.theta = theta;
  primaryInfo.theta_scat = theta_scat;
  primaryInfo.theta_CM = theta_CM;
}

//_____________________________________________________________________________
void
TPCAnaManager::SetNumberOfPrimaryParticle( G4int num )
{
  primaryParticle.NumOfParticle = num;
}

//_____________________________________________________________________________
void
TPCAnaManager::SetPrimaryParticle( G4double px, G4double py, G4double pz )
{
  G4int id = 0;
  primaryParticle.p0[id][0] = px;
  primaryParticle.p0[id][1] = py;
  primaryParticle.p0[id][2] = pz;
  primaryParticle.p0[id][3] = sqrt(pow(px,2.0)+pow(py,2.0)+pow(pz,2.0));
  primaryParticle.p0[id][4] = 0.0;
}

//_____________________________________________________________________________
void
TPCAnaManager::SetGeneratorID( G4int gen )
{
  primaryBeam.gen = gen;
}

//_____________________________________________________________________________
void
TPCAnaManager::SetModeID( G4int mode )
{
  primaryBeam.mode = mode;
}

//_____________________________________________________________________________
void
TPCAnaManager::SetPrimaryParticle( G4int id, G4double px, G4double py,
				   G4double pz, G4double mass )
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

//_____________________________________________________________________________
void
TPCAnaManager::SetPrimaryParticle( G4int id, G4double px, G4double py,
				   G4double pz, G4double mass, G4int pid )
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

//_____________________________________________________________________________
void
TPCAnaManager::SetPrimaryVertex( G4int id, G4double x, G4double y, G4double z )
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
  event.pb->SetXYZ( p.x(), p.y(), p.z() );
}

//_____________________________________________________________________________
void
TPCAnaManager::SetPrimaryBeam( G4double px, G4double py, G4double pz )
{
  event.pb->SetXYZ( px, py, pz );
}

//_____________________________________________________________________________
void
TPCAnaManager::SetScintData( G4double time, G4ThreeVector pos,
			     G4ThreeVector mom,
			     G4int track, G4int particle, G4int detector,
			     G4double mass, G4int qq, G4int parentid,
			     G4ThreeVector vtxpos, G4ThreeVector vtxmom,
			     G4double vtxene, G4double tlength )
{
  G4int hitnum = HitNumScint;
  if (hitnum >= MaxTrack) {
    G4cerr << FUNC_NAME << " too much nhit" << hitnum << G4endl;
    return;
  }

  scintData[hitnum].time = time;
  scintData[hitnum].pos[G4ThreeVector::X] = pos.getX();
  scintData[hitnum].pos[G4ThreeVector::Y] = pos.getY();
  scintData[hitnum].pos[G4ThreeVector::Z] = pos.getZ();
  scintData[hitnum].mom[G4ThreeVector::X] = mom.getX();
  scintData[hitnum].mom[G4ThreeVector::Y] = mom.getY();
  scintData[hitnum].mom[G4ThreeVector::Z] = mom.getZ();
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

//_____________________________________________________________________________
void
TPCAnaManager::SetACData( G4double time, G4ThreeVector pos, G4ThreeVector mom,
			  G4int track, G4int particle, G4int detector,
			  G4double mass, G4int qq, G4int parentid,
			  G4ThreeVector vtxpos, G4ThreeVector vtxmom,
			  G4double vtxene, G4double tlength )
{
  G4int hitnum = HitNumAC;
  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }
  acData[hitnum].time = time;
  acData[hitnum].pos[G4ThreeVector::X] = pos.getX();
  acData[hitnum].pos[G4ThreeVector::Y] = pos.getY();
  acData[hitnum].pos[G4ThreeVector::Z] = pos.getZ();
  acData[hitnum].mom[G4ThreeVector::X] = mom.getX();
  acData[hitnum].mom[G4ThreeVector::Y] = mom.getY();
  acData[hitnum].mom[G4ThreeVector::Z] = mom.getZ();
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

//_____________________________________________________________________________
void
TPCAnaManager::SetNBARData( G4double time, G4ThreeVector pos, G4ThreeVector mom,
			    G4int track, G4int particle, G4int detector,
			    G4double mass,G4int qq,G4int parentid,
			    G4ThreeVector vtxpos, G4ThreeVector vtxmom,
			    G4double vtxene, G4double tlength )
{
  G4int hitnum = HitNumNBAR;
  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }
  nbarData[hitnum].time = time;
  nbarData[hitnum].pos[G4ThreeVector::X] = pos.getX();
  nbarData[hitnum].pos[G4ThreeVector::Y] = pos.getY();
  nbarData[hitnum].pos[G4ThreeVector::Z] = pos.getZ();
  nbarData[hitnum].mom[G4ThreeVector::X] = mom.getX();
  nbarData[hitnum].mom[G4ThreeVector::Y] = mom.getY();
  nbarData[hitnum].mom[G4ThreeVector::Z] = mom.getZ();
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

//_____________________________________________________________________________
void
TPCAnaManager::SetDCData( G4double time, G4ThreeVector pos, G4ThreeVector mom,
			  G4int track, G4int particle, G4int detector,
			  G4double mass,G4int qq, G4int parentid,
			  G4ThreeVector vtxpos, G4ThreeVector vtxmom,
			  G4double vtxene, G4double tlength )
{
  G4int hitnum = HitNumDC;
  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }
  dcData[hitnum].time = time;
  dcData[hitnum].pos[G4ThreeVector::X] = pos.getX();
  dcData[hitnum].pos[G4ThreeVector::Y] = pos.getY();
  dcData[hitnum].pos[G4ThreeVector::Z] = pos.getZ();
  dcData[hitnum].mom[G4ThreeVector::X] = mom.getX();
  dcData[hitnum].mom[G4ThreeVector::Y] = mom.getY();
  dcData[hitnum].mom[G4ThreeVector::Z] = mom.getZ();
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

//_____________________________________________________________________________
void
TPCAnaManager::SetSCHData( G4double time, G4ThreeVector pos, G4ThreeVector mom,
			   G4int track, G4int particle, G4int detector,
			   G4double mass,G4int qq, G4int parentid,
			   G4ThreeVector vtxpos, G4ThreeVector vtxmom,
			   G4double vtxene, G4double tlength )
{
  G4int hitnum = HitNumSCH;
  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }
  schData[hitnum].time = time;
  schData[hitnum].pos[G4ThreeVector::X] = pos.getX();
  schData[hitnum].pos[G4ThreeVector::Y] = pos.getY();
  schData[hitnum].pos[G4ThreeVector::Z] = pos.getZ();
  schData[hitnum].mom[G4ThreeVector::X] = mom.getX();
  schData[hitnum].mom[G4ThreeVector::Y] = mom.getY();
  schData[hitnum].mom[G4ThreeVector::Z] = mom.getZ();
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

//_____________________________________________________________________________
void
TPCAnaManager::SetFTOFData( G4double time, G4ThreeVector pos, G4ThreeVector mom,
			    G4int track, G4int particle, G4int detector,
			    G4double mass, G4int qq,G4int parentid,
			    G4ThreeVector vtxpos, G4ThreeVector vtxmom,
			    G4double vtxene, G4double tlength )
{
  G4int hitnum = HitNumFTOF;
  if (hitnum >= MaxTrack) {
    fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	    hitnum);
    return;
  }
  ftofData[hitnum].time = time;
  ftofData[hitnum].pos[G4ThreeVector::X] = pos.getX();
  ftofData[hitnum].pos[G4ThreeVector::Y] = pos.getY();
  ftofData[hitnum].pos[G4ThreeVector::Z] = pos.getZ();
  ftofData[hitnum].mom[G4ThreeVector::X] = mom.getX();
  ftofData[hitnum].mom[G4ThreeVector::Y] = mom.getY();
  ftofData[hitnum].mom[G4ThreeVector::Z] = mom.getZ();
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

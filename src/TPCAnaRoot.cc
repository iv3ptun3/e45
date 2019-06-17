// -*- C++ -*-

#include "TPCAnaRoot.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include "Randomize.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

namespace
{
  using CLHEP::GeV;
  using CLHEP::cm;
  using CLHEP::mm;
}

//_____________________________________________________________________________
TPCAnaRoot::TPCAnaRoot( void )
{
}

//_____________________________________________________________________________
TPCAnaRoot::~TPCAnaRoot( void )
{
}

//_____________________________________________________________________________
void
TPCAnaRoot::BeginOfRunAction( int runnum )
{
  //  char filename[30];


  char *filename;
  int i;
  char histname1[40],histname2[40];

  //  sprintf(filename, "g4E42_%d.root", runnum);
  //  G4String Filename = getenv("Out_ROOT_File_Name");
  filename = getenv("Out_ROOT_File_Name");
  printf("##### Open root file '%s'... #####\n",filename);
  //  filename= Filename.c_str();
  //  filename= Filename;

  /* output file */
  rootfile = gFile; // new TFile(filename,"RECREATE","Output file of TPC for E42");

  /* Time */
  htime = new TH1F("Time","Time",400, 0.0, 10.0);

  /* Hit position */
  for(i=0;i<3;i++){
    sprintf(histname1, "Pos%d", i);
    sprintf(histname2, "Pos%d", i);
    hpos[i] = new TH1F(histname1,histname2,500, -25.0*cm, 25.0*cm);
  }

  /* Momentum */
  for(i=0;i<3;i++){
    sprintf(histname1, "Mom%d", i);
    sprintf(histname2, "Mom%d", i);
    if(i==2){
      hmom[i] = new TH1F(histname1,histname2,400, 0.0, 2.0);
    } else{
      hmom[i] = new TH1F(histname1,histname2,400, -2.0, 2.);
    }
  }
  tree = new TTree("tree","GEANT4 simulation for hypTPC");

  //
  // TPC Pad
  //

  tree->Branch("ev",&tree1ev.ev,"ev/I");
  tree->Branch("pg",tree1ev.pg,"pg[4]/D");
  tree->Branch("npid",&tree1ev.npid,"npid/I");
  tree->Branch("pid0",tree1ev.pid0,"pid0[npid]/I");
  tree->Branch("x0",tree1ev.x0,"x0[npid][3]/D");
  tree->Branch("p0",tree1ev.p0,"p0[npid][5]/D");
  tree->Branch("pt0",tree1ev.pt0,"pt0[npid]/D");
  tree->Branch("mass0",tree1ev.mass0,"mass0[npid]/D");
  tree->Branch("theta0",tree1ev.theta0,"theta0[npid]/D");

  tree->Branch("mm_d",&tree1ev.mm_d,"mm_d/D");
  // tree->Branch("mm_p",&tree1ev.mm_p,"mm_p/D");
  tree->Branch("theta",&tree1ev.theta,"theta/D");
  tree->Branch("theta_scat",&tree1ev.theta_scat,"theta_scat/D");
  tree->Branch("theta_CM",&tree1ev.theta_CM,"theta_CM/D");

  // tree->Branch("mm",&tree1ev.mm,"mm/D");



  //generator mode
  tree->Branch("gen",&tree1ev.gen,"gen/I");
  tree->Branch("mode",&tree1ev.mode,"mode/I");



  //Num Of Hit for K+



  ///////shhwang tpc hit step

  // tree->Branch("nttpc",&tree1ev.nttpc,"nttpc/I");
  // tree->Branch("ntrk",tree1ev.ntrk,"ntrk[nttpc]/I");
  // tree->Branch("ititpc",tree1ev.ititpc,"ititpc[nttpc]/I");
  // tree->Branch("idtpc",tree1ev.idtpc,"idtpc[nttpc]/I");
  // tree->Branch("xtpc",tree1ev.xtpc,"xtpc[nttpc]/D");
  // tree->Branch("ytpc",tree1ev.ytpc,"ytpc[nttpc]/D");
  // tree->Branch("ztpc",tree1ev.ztpc,"ztpc[nttpc]/D");
  // tree->Branch("x0tpc",tree1ev.x0tpc,"x0tpc[nttpc]/D");
  // tree->Branch("y0tpc",tree1ev.y0tpc,"y0tpc[nttpc]/D");
  // tree->Branch("z0tpc",tree1ev.z0tpc,"z0tpc[nttpc]/D");
  // tree->Branch("resoX",tree1ev.resoX,"resoX[nttpc]/D");
  // tree->Branch("pxtpc",tree1ev.pxtpc,"pxtpc[nttpc]/D");
  // tree->Branch("pytpc",tree1ev.pytpc,"pytpc[nttpc]/D");
  // tree->Branch("pztpc",tree1ev.pztpc,"pztpc[nttpc]/D");
  // tree->Branch("pptpc",tree1ev.pptpc,"pptpc[nttpc]/D");   // total mometum
  // tree->Branch("masstpc",tree1ev.masstpc,"masstpc[nttpc]/D");   // mass TPC
  // //  tree->Branch("betatpc",tree1ev.betatpc,"betatpc[nttpc]/D");
  // // tree->Branch("edeptpc",tree1ev.edeptpc,"edeptpc[nttpc]/D");
  // //tree->Branch("dedxtpc",tree1ev.dedxtpc,"dedxtpc[nttpc]/D");
  // //tree->Branch("slengthtpc",tree1ev.slengthtpc,"slengthtpc[nttpc]/D");
  // tree->Branch("laytpc",tree1ev.laytpc,"laytpc[nttpc]/I");
  // tree->Branch("rowtpc",tree1ev.rowtpc,"rowtpc[nttpc]/I");
  // tree->Branch("parentID",tree1ev.parentID,"parentID[nttpc]/I");

  //// study on multiplicity
  // tree->Branch("nthlay",tree1ev.nthlay,"nthlay[nttpc]/I");
  // tree->Branch("nthpad",tree1ev.nthpad,"nthpad[nttpc]/I");
  // tree->Branch("laypad",tree1ev.laypad,"laytpadpc[nttpc][40][250]/I");


  //shhwang ntrtpc --> number of trak in tpc
  tree->Branch("ntrtpc",&tree1ev.ntrtpc,"ntrtpc/I");
  tree->Branch("trpmtpc",tree1ev.trpmtpc,"trpmtpc[ntrtpc]/D");
  tree->Branch("trqqtpc",tree1ev.trqqtpc,"trqqtpc[ntrtpc]/I");
  tree->Branch("trpidtpc",tree1ev.trpidtpc,"trpidtpc[ntrtpc]/I");
  tree->Branch("trparentidtpc",tree1ev.trparentidtpc,"trparentidtpc[ntrtpc]/I");
  //tree->Branch("trparentid_pid_tpc",tree1ev.trparentid_pid_tpc,"trparentid_pid_tpc[ntrtpc]/I");

  tree->Branch("trpxtpc",tree1ev.trpxtpc,"trpxtpc[ntrtpc]/D");
  tree->Branch("trpytpc",tree1ev.trpytpc,"trpytpc[ntrtpc]/D");
  tree->Branch("trpztpc",tree1ev.trpztpc,"trpztpc[ntrtpc]/D");
  tree->Branch("trpptpc",tree1ev.trpptpc,"trpptpc[ntrtpc]/D");
  tree->Branch("trpttpc",tree1ev.trpttpc,"trpttpc[ntrtpc]/D");

  tree->Branch("trpxtpcfit",tree1ev.trpxtpcfit,"trpxtpcfit[ntrtpc]/D");
  tree->Branch("trpytpcfit",tree1ev.trpytpcfit,"trpytpcfit[ntrtpc]/D");
  tree->Branch("trpztpcfit",tree1ev.trpztpcfit,"trpztpcfit[ntrtpc]/D");
  tree->Branch("trpptpcfit",tree1ev.trpptpcfit,"trpptpcfit[ntrtpc]/D");
  tree->Branch("trpttpcfit",tree1ev.trpttpcfit,"trpttpcfit[ntrtpc]/D");

  tree->Branch("trvtxpxtpc",tree1ev.trvtxpxtpc,"trvtxpxtpc[ntrtpc]/D");
  tree->Branch("trvtxpytpc",tree1ev.trvtxpytpc,"trvtxpytpc[ntrtpc]/D");
  tree->Branch("trvtxpztpc",tree1ev.trvtxpztpc,"trvtxpztpc[ntrtpc]/D");
  tree->Branch("trvtxpptpc",tree1ev.trvtxpptpc,"trvtxpptpc[ntrtpc]/D");

  tree->Branch("trvtxxtpc",tree1ev.trvtxxtpc,"trvtxxtpc[ntrtpc]/D");
  tree->Branch("trvtxytpc",tree1ev.trvtxytpc,"trvtxytpc[ntrtpc]/D");
  tree->Branch("trvtxztpc",tree1ev.trvtxztpc,"trvtxztpc[ntrtpc]/D");

  tree->Branch("trvtxxtpcfit",tree1ev.trvtxxtpcfit,"trvtxxtpcfit[ntrtpc]/D");
  tree->Branch("trvtxytpcfit",tree1ev.trvtxytpcfit,"trvtxytpcfit[ntrtpc]/D");
  tree->Branch("trvtxztpcfit",tree1ev.trvtxztpcfit,"trvtxztpcfit[ntrtpc]/D");

  // tree->Branch("trdetpc",tree1ev.trdetpc,"trdetpc[ntrtpc]/D");
  // tree->Branch("trlentpc",tree1ev.trlentpc,"trlentpc[ntrtpc]/D");
  // tree->Branch("trdedxtpc",tree1ev.trdedxtpc,"trdedxtpc[ntrtpc]/D");
  // tree->Branch("trdedxtrtpc",tree1ev.trdedxtrtpc,"trdedxtrtpc[ntrtpc]/D");
  tree->Branch("trlaytpc",tree1ev.trlaytpc,"trlaytpc[ntrtpc]/I");
  tree->Branch("cir_r",tree1ev.cir_r,"cir_r[ntrtpc]/D");
  tree->Branch("cir_x",tree1ev.cir_x,"cir_x[ntrtpc]/D");
  tree->Branch("cir_z",tree1ev.cir_z,"cir_z[ntrtpc]/D");
  tree->Branch("cir_fit",tree1ev.cir_fit,"cir_fit[ntrtpc]/D");
  tree->Branch("vtx_flag",tree1ev.vtx_flag,"vtx_flag[ntrtpc]/I");
  tree->Branch("a_fory",tree1ev.a_fory,"a_fory[ntrtpc]/D");
  tree->Branch("b_fory",tree1ev.b_fory,"b_fory[ntrtpc]/D");

  //
  // Scinti
  //

  tree->Branch("ntsc",&tree1ev.ntsc,"ntsc/I");
  tree->Branch("tidsc",tree1ev.tidsc,"tidsc[ntsc]/I");
  tree->Branch("pidsc",tree1ev.pidsc,"pidsc[ntsc]/I");
  tree->Branch("didsc",tree1ev.didsc,"didsc[ntsc]/I");
  tree->Branch("masssc",tree1ev.masssc,"masssc[ntsc]/D");
  tree->Branch("qqsc",tree1ev.qqsc,"qqsc[ntsc]/I");
  tree->Branch("xsc",tree1ev.xsc,"xsc[ntsc]/D");
  tree->Branch("ysc",tree1ev.ysc,"ysc[ntsc]/D");
  tree->Branch("zsc",tree1ev.zsc,"zsc[ntsc]/D");
  tree->Branch("pxsc",tree1ev.pxsc,"pxsc[ntsc]/D");
  tree->Branch("pysc",tree1ev.pysc,"pysc[ntsc]/D");
  tree->Branch("pzsc",tree1ev.pzsc,"pzsc[ntsc]/D");
  tree->Branch("ppsc",tree1ev.ppsc,"ppsc[ntsc]/D");
  tree->Branch("tofsc",tree1ev.tofsc,"tofsc[ntsc]/D");
  tree->Branch("scpID",tree1ev.scpID,"scpID[ntsc]/I");

  tree->Branch("trvtxpxscint",tree1ev.trvtxpxscint,"trvtxpxscint[ntsc]/D");
  tree->Branch("trvtxpyscint",tree1ev.trvtxpyscint,"trvtxpyscint[ntsc]/D");
  tree->Branch("trvtxpzscint",tree1ev.trvtxpzscint,"trvtxpzscint[ntsc]/D");
  tree->Branch("trvtxppscint",tree1ev.trvtxppscint,"trvtxppscint[ntsc]/D");

  tree->Branch("trvtxxscint",tree1ev.trvtxxscint,"trvtxxscint[ntsc]/D");
  tree->Branch("trvtxyscint",tree1ev.trvtxyscint,"trvtxyscint[ntsc]/D");
  tree->Branch("trvtxzscint",tree1ev.trvtxzscint,"trvtxzscint[ntsc]/D");
  tree->Branch("lengthsc",tree1ev.lengthsc,"lengthsc[ntsc]/D");

  //target
  // tree->Branch("targethits",&tree1ev.targethits,"targethits/I");
  // tree->Branch("targetpid",tree1ev.targetpid,"targetpid[targethits]/I");
  // tree->Branch("targetparentid",tree1ev.targetparentid,"targetparentid[targethits]/I");
  // tree->Branch("targettid",tree1ev.targettid,"targettid[targethits]/I");
  // tree->Branch("targetpos",tree1ev.targetpos,"targetpos[targethits][3]/D");
  // tree->Branch("targetvtx",tree1ev.targetvtx,"targetvtx[targethits][3]/D");

  tree1ev.ev = 0;


}

void TPCAnaRoot::EndOfRunAction()
{
  rootfile->Write();
  rootfile->Close();
  tree1ev.ev = 0;
}


void TPCAnaRoot::BeginOfTrackingAction()
{
  G4cout<< "test tracking action"<<G4endl;
}
void TPCAnaRoot::BeginOfEventAction()
{

  tree1ev.ev++;
  tree1ev.nttpc = 0;
  tree1ev.ntac = 0;
  tree1ev.ntnbar = 0;
  tree1ev.ntsc = 0;
  tree1ev.ntdc = 0;
  tree1ev.ntch = 0;
  tree1ev.ntftof = 0;
  tree1ev.npid = 0;
  tree1ev.ntrtpc = 0;
  tree1ev.targethits = 0;

  tree1ev.gen = 0;
  tree1ev.mode = 0;

  tree1ev.HitNum_K=-1;
  tree1ev.HitNumAC_K=-1;
  tree1ev.HitNumNBAR_K=-1;
  tree1ev.HitNumDC_K=-1;
  tree1ev.HitNumFTOF_K=-1;
  tree1ev.HitNumCH_K=-1;
  tree1ev.HitNumScint_K=-1;
  tree1ev.HitNumTarget_K=-1;


  tree1ev.HitNum_p=-1;
  tree1ev.HitNumAC_p=-1;
  tree1ev.HitNumNBAR_p=-1;
  tree1ev.HitNumDC_p=-1;
  tree1ev.HitNumFTOF_p=-1;
  tree1ev.HitNumCH_p=-1;
  tree1ev.HitNumScint_p=-1;
  tree1ev.HitNumTarget_p=-1;


  tree1ev.mm_d = 0.;
  tree1ev.mm_p = 0.;
  tree1ev.theta = 0.;
  tree1ev.theta_scat = 0.;
  tree1ev.theta_CM = 0.;


  /* ntrtpc initialization */

  for(int i = 0; i< MaxTrackTPC;i++){
    tree1ev.trpidtpc[i]  = -1;
    tree1ev.trparentidtpc[i]  = -1;
    tree1ev.trparentid_pid_tpc[i]  = -1;


    tree1ev.trpptpc[i]  = -9999.9999;
    tree1ev.trpttpc[i]  = -9999.9999;
    tree1ev.trpxtpc[i]  = -9999.9999;
    tree1ev.trpytpc[i]  = -9999.9999;
    tree1ev.trpztpc[i]  = -9999.9999;

    tree1ev.trvtxpxtpc[i]  = -9999.9999;
    tree1ev.trvtxpytpc[i]  = -9999.9999;
    tree1ev.trvtxpztpc[i]  = -9999.9999;

    tree1ev.trvtxxtpc[i]  = -9999.9999;
    tree1ev.trvtxytpc[i]  = -9999.9999;
    tree1ev.trvtxztpc[i]  = -9999.9999;

    tree1ev.trpttpcfit[i]  = -9999.9999;

    tree1ev.trpptpcfit[i]  = -9999.9999;
    tree1ev.trpxtpcfit[i]  = -9999.9999;
    tree1ev.trpytpcfit[i]  = -9999.9999;
    tree1ev.trpztpcfit[i]  = -9999.9999;

    tree1ev.trpmtpc[i]  = -9999.9999;
    tree1ev.trqqtpc[i]  = -9999.9999;

    tree1ev.trdetpc[i]  = -9999.9999;
    tree1ev.trlentpc[i]  = -9999.9999;
    tree1ev.trdedxtpc[i]  = -9999.9999;
    tree1ev.trdedxtrtpc[i]  = -9999.9999;
    tree1ev.trlaytpc[i]  = -9999.9999;

    tree1ev.cir_r[i]  = -9999.9999;
    tree1ev.cir_x[i]  = -9999.9999;
    tree1ev.cir_z[i]  = -9999.9999;
    tree1ev.cir_fit[i]  = -9999.9999;

    tree1ev.vtx_flag[i]  = -1;
    tree1ev.a_fory[i]  = -9999.9999;
    tree1ev.b_fory[i]  = -9999.9999;
  }

  for(int i = 0; i< 4;i++){
    tree1ev.pg[i]  = -9999.9;
  }

  for(int i = 0; i< MaxPrimaryParticle;i++){
    for(int j=0;j<3;j++){
      tree1ev.x0[i][j] = -9999.9;
    }
    for(int j=0;j<5;j++){
      tree1ev.p0[i][j] = -9999.9;
    }
    tree1ev.theta0[i] = -9999.9;
    tree1ev.pt0[i] = -9999.9;
    tree1ev.mass0[i] = -9999.9;
    tree1ev.pid0[i] = -9999;
  }


  for(int i=0;i<MaxTrack;i++){

  /// initialization pad multiplicity
  tree1ev.nthlay[i]=-9999.;
  tree1ev.nthpad[i]=-9999.;
  for(int j = 0; j< MaxNthLay;j++){
    for(int k = 0; k< MaxNthPad;k++){
      tree1ev.laypad[i][j][k]  = 0.;
    }
  }
  //////////////


    tree1ev.xtpc[i] = -9999.9;
    tree1ev.ytpc[i] = -9999.9;
    tree1ev.ztpc[i] = -9999.9;

    tree1ev.x0tpc[i] = -9999.9;
    tree1ev.y0tpc[i] = -9999.9;
    tree1ev.z0tpc[i] = -9999.9;
    tree1ev.resoX[i] = -9999.9;

    tree1ev.pxtpc[i] = -9999.9;
    tree1ev.pytpc[i] = -9999.9;
    tree1ev.pztpc[i] = -9999.9;
    tree1ev.pptpc[i] = -9999.9;

    tree1ev.masstpc[i] = -9999.9;

    tree1ev.betatpc[i] = -9999.9;

    tree1ev.edeptpc[i] = -9999.9;

    tree1ev.ititpc[i] = -1;
    tree1ev.idtpc[i] = -1;
    tree1ev.laytpc[i] = -1;
    tree1ev.rowtpc[i] = -1;
    tree1ev.parentID[i] = -1;
  }

  for(int i=0; i<MaxTrig; i++){
    tree1ev.xsc[i] = -9999.9;
    tree1ev.ysc[i] = -9999.9;
    tree1ev.zsc[i] = -9999.9;
    tree1ev.pxsc[i] = -9999.9;
    tree1ev.pysc[i] = -9999.9;
    tree1ev.pzsc[i] = -9999.9;
    tree1ev.tofsc[i] = -9999.9;
    tree1ev.scpID[i] = -9999;

    tree1ev.tidsc[i] = -1;
    tree1ev.pidsc[i] = -1;
    tree1ev.didsc[i] = -1;
    tree1ev.masssc[i] = -1;
    tree1ev.qqsc[i] = -1;

    tree1ev.trvtxppscint[i]  = -9999.9999;
    tree1ev.trvtxpxscint[i]  = -9999.9999;
    tree1ev.trvtxpyscint[i]  = -9999.9999;
    tree1ev.trvtxpzscint[i]  = -9999.9999;

    tree1ev.trvtxxscint[i]  = -9999.9999;
    tree1ev.trvtxyscint[i]  = -9999.9999;
    tree1ev.trvtxzscint[i]  = -9999.9999;
  }
  ///ac
  for(int i=0; i<MaxTrig; i++){
    tree1ev.xac[i] = -9999.9;
    tree1ev.yac[i] = -9999.9;
    tree1ev.zac[i] = -9999.9;
    tree1ev.pxac[i] = -9999.9;
    tree1ev.pyac[i] = -9999.9;
    tree1ev.pzac[i] = -9999.9;
    tree1ev.tofac[i] = -9999.9;
    tree1ev.acpID[i] = -9999;

    tree1ev.tidac[i] = -1;
    tree1ev.pidac[i] = -1;
    tree1ev.didac[i] = -1;
    tree1ev.massac[i] = -1;
    tree1ev.qqac[i] = -1;

    tree1ev.trvtxppac[i]  = -9999.9999;
    tree1ev.trvtxpxac[i]  = -9999.9999;
    tree1ev.trvtxpyac[i]  = -9999.9999;
    tree1ev.trvtxpzac[i]  = -9999.9999;

    tree1ev.trvtxxac[i]  = -9999.9999;
    tree1ev.trvtxyac[i]  = -9999.9999;
    tree1ev.trvtxzac[i]  = -9999.9999;
  }

  ///dc
  for(int i=0; i<MaxTrig; i++){
    tree1ev.xdc[i] = -9999.9;
    tree1ev.ydc[i] = -9999.9;
    tree1ev.zdc[i] = -9999.9;
    tree1ev.pxdc[i] = -9999.9;
    tree1ev.pydc[i] = -9999.9;
    tree1ev.pzdc[i] = -9999.9;
    tree1ev.tofdc[i] = -9999.9;
    tree1ev.dcpID[i] = -9999;

    tree1ev.tiddc[i] = -1;
    tree1ev.piddc[i] = -1;
    tree1ev.diddc[i] = -1;
    tree1ev.massdc[i] = -1;
    tree1ev.qqdc[i] = -1;

    tree1ev.trvtxppdc[i]  = -9999.9999;
    tree1ev.trvtxpxdc[i]  = -9999.9999;
    tree1ev.trvtxpydc[i]  = -9999.9999;
    tree1ev.trvtxpzdc[i]  = -9999.9999;

    tree1ev.trvtxxdc[i]  = -9999.9999;
    tree1ev.trvtxydc[i]  = -9999.9999;
    tree1ev.trvtxzdc[i]  = -9999.9999;
  }


  ///ch
  for(int i=0; i<MaxTrig; i++){
    tree1ev.xch[i] = -9999.9;
    tree1ev.ych[i] = -9999.9;
    tree1ev.zch[i] = -9999.9;
    tree1ev.pxch[i] = -9999.9;
    tree1ev.pych[i] = -9999.9;
    tree1ev.pzch[i] = -9999.9;
    tree1ev.tofch[i] = -9999.9;
    tree1ev.chpID[i] = -9999;

    tree1ev.tidch[i] = -1;
    tree1ev.pidch[i] = -1;
    tree1ev.didch[i] = -1;
    tree1ev.massch[i] = -1;
    tree1ev.qqch[i] = -1;

    tree1ev.trvtxppch[i]  = -9999.9999;
    tree1ev.trvtxpxch[i]  = -9999.9999;
    tree1ev.trvtxpych[i]  = -9999.9999;
    tree1ev.trvtxpzch[i]  = -9999.9999;

    tree1ev.trvtxxch[i]  = -9999.9999;
    tree1ev.trvtxych[i]  = -9999.9999;
    tree1ev.trvtxzch[i]  = -9999.9999;
  }

  ///ftof
  for(int i=0; i<MaxTrig; i++){
    tree1ev.xftof[i] = -9999.9;
    tree1ev.yftof[i] = -9999.9;
    tree1ev.zftof[i] = -9999.9;
    tree1ev.pxftof[i] = -9999.9;
    tree1ev.pyftof[i] = -9999.9;
    tree1ev.pzftof[i] = -9999.9;
    tree1ev.tofftof[i] = -9999.9;
    tree1ev.ftofpID[i] = -9999;

    tree1ev.tidftof[i] = -1;
    tree1ev.pidftof[i] = -1;
    tree1ev.didftof[i] = -1;
    tree1ev.massftof[i] = -1;
    tree1ev.qqftof[i] = -1;

    tree1ev.trvtxppftof[i]  = -9999.9999;
    tree1ev.trvtxpxftof[i]  = -9999.9999;
    tree1ev.trvtxpyftof[i]  = -9999.9999;
    tree1ev.trvtxpzftof[i]  = -9999.9999;

    tree1ev.trvtxxftof[i]  = -9999.9999;
    tree1ev.trvtxyftof[i]  = -9999.9999;
    tree1ev.trvtxzftof[i]  = -9999.9999;
  }

  /*
  for(int i=0; i<MaxTrackFDC; i++){
    tree1ev.xfdc[i] = -9999.9;
    tree1ev.yfdc[i] = -9999.9;
    tree1ev.zfdc[i] = -9999.9;
    tree1ev.pxfdc[i] = -9999.9;
    tree1ev.pyfdc[i] = -9999.9;
    tree1ev.pzfdc[i] = -9999.9;
    tree1ev.toffdc[i] = -9999.9;

    tree1ev.tidfdc[i] = -1;
    tree1ev.pidfdc[i] = -1;
    tree1ev.didfdc[i] = -1;
  }
  */

  for(int i=0; i<MaxTrack; i++){
    //    tree1ev.targethits = -9999;
    tree1ev.targetpid[i]=-9999;
    tree1ev.targettid[i]=-9999;
    tree1ev.targetparentid[i]=-9999;
    tree1ev.targetpos[i][0]=-9999.9999;
    tree1ev.targetpos[i][1]=-9999.9999;
    tree1ev.targetpos[i][2]=-9999.9999;

    tree1ev.targetvtx[i][0]=-9999.9999;
    tree1ev.targetvtx[i][1]=-9999.9999;
    tree1ev.targetvtx[i][2]=-9999.9999;
  }
}

//void TPCAnaRoot::FillScintData(G4double time, G4double *pos, G4double *mom,
//			  G4int tid, G4int pid, G4int did)
//{
//  tree1ev.tofsc[tree1ev.ntsc] = time;
//  tree1ev.xsc[tree1ev.ntsc] = pos[0];
//  tree1ev.ysc[tree1ev.ntsc] = pos[1];
//  tree1ev.zsc[tree1ev.ntsc] = pos[2];


void TPCAnaRoot::FillGenMode(G4int gen_in,G4int mode_in)
{
  //  tree1ev.slengthtpc[tree1ev.nttpc] = slength;
  tree1ev.gen = gen_in;
  tree1ev.mode = mode_in;
}


void TPCAnaRoot::FillLayerPad(G4int nlay,G4int npad)
{
  //  tree1ev.slengthtpc[tree1ev.nttpc] = slength;

  tree1ev.nthlay[tree1ev.nttpc] = nlay;
  tree1ev.nthpad[tree1ev.nttpc] = npad;
  //  G4cout<<nlay<<":"<<npad<<G4endl;
  tree1ev.laypad[tree1ev.nttpc][tree1ev.nthlay[tree1ev.nttpc]][tree1ev.nthpad[tree1ev.nttpc]] = tree1ev.laypad[tree1ev.nttpc][tree1ev.nthlay[tree1ev.nttpc]][tree1ev.nthpad[tree1ev.nttpc]]+1.;
  //  G4cout<<tree1ev.laypad[tree1ev.nttpc][tree1ev.nthlay[tree1ev.nttpc]][tree1ev.nthpad[tree1ev.nttpc]]<<G4endl;
}



void TPCAnaRoot::FillTPCData(G4double tpcpx,G4double tpcpy,
			     G4double tpcpz,G4double tpcpp,
			     G4int tpcpid, G4int tpcparentid, G4int tpcparentid_pid,
			     G4int tpcqq,
			     G4double tpcpm, G4double tpcde, G4double tpclen,
			     G4double tpcdedx,G4double tpcdedxtr, G4int tpclay,
			     G4double tpcvtxpx,G4double tpcvtxpy,G4double tpcvtxpz,
			     G4double tpcvtxx,G4double tpcvtxy,G4double tpcvtxz,
			     G4double tpcvtxxfit,G4double tpcvtxyfit,G4double tpcvtxzfit,
			     G4double tpcpxfit,G4double tpcpyfit,G4double tpcpzfit,G4double tpcptfit,
			     G4double cir_r, G4double cir_x, G4double cir_z, G4double cir_fit,
			     G4int vtx_flag, G4double a_fory, G4double b_fory)//--> should be changed

{
  tree1ev.trpidtpc[tree1ev.ntrtpc] = tpcpid;
  tree1ev.trparentidtpc[tree1ev.ntrtpc] = tpcparentid;
  tree1ev.trparentid_pid_tpc[tree1ev.ntrtpc] = tpcparentid_pid;

  tree1ev.trpxtpc[tree1ev.ntrtpc] = tpcpx/GeV;
  tree1ev.trpytpc[tree1ev.ntrtpc] = tpcpy/GeV;
  tree1ev.trpztpc[tree1ev.ntrtpc] = tpcpz/GeV;
  tree1ev.trpptpc[tree1ev.ntrtpc] = sqrt(pow(tpcpx,2)+pow(tpcpy,2)+pow(tpcpz,2))/GeV;
  tree1ev.trpttpc[tree1ev.ntrtpc] = sqrt(pow(tpcpx,2)+pow(tpcpz,2))/GeV;

  tree1ev.trpxtpcfit[tree1ev.ntrtpc] = tpcpxfit/GeV;
  tree1ev.trpytpcfit[tree1ev.ntrtpc] = tpcpyfit/GeV;
  tree1ev.trpztpcfit[tree1ev.ntrtpc] = tpcpzfit/GeV;
  //  tree1ev.trpptpcfit[tree1ev.ntrtpc] = sqrt(pow(tpcpxfit,2)+pow(tpcpyfit,2)+pow(tpcpzfit,2))/GeV;
  tree1ev.trpptpcfit[tree1ev.ntrtpc] = sqrt(pow(tpcpyfit,2)+pow(tpcptfit,2))/GeV;
  tree1ev.trpttpcfit[tree1ev.ntrtpc] = tpcptfit/GeV;

  tree1ev.trvtxpxtpc[tree1ev.ntrtpc] = tpcvtxpx/GeV;
  tree1ev.trvtxpytpc[tree1ev.ntrtpc] = tpcvtxpy/GeV;
  tree1ev.trvtxpztpc[tree1ev.ntrtpc] = tpcvtxpz/GeV;
  tree1ev.trvtxpptpc[tree1ev.ntrtpc] = sqrt(pow(tpcvtxpx,2)+pow(tpcvtxpy,2)+pow(tpcvtxpz,2))/GeV;

  //  G4cout<<tpcpxfit/GeV<<" "<<tpcvtxpx/GeV<<" "<<tpcvtxpx/GeV<<G4endl;

  tree1ev.trvtxxtpc[tree1ev.ntrtpc] = tpcvtxx;
  tree1ev.trvtxytpc[tree1ev.ntrtpc] = tpcvtxy;
  tree1ev.trvtxztpc[tree1ev.ntrtpc] = tpcvtxz;

  tree1ev.trvtxxtpcfit[tree1ev.ntrtpc] = tpcvtxxfit;
  tree1ev.trvtxytpcfit[tree1ev.ntrtpc] = tpcvtxyfit;
  tree1ev.trvtxztpcfit[tree1ev.ntrtpc] = tpcvtxzfit;

  //  G4cout<<"anaroot"<<G4endl;
  //  G4cout<<tpcvtxz<<G4endl;

  tree1ev.trqqtpc[tree1ev.ntrtpc] = tpcqq;
  tree1ev.trpmtpc[tree1ev.ntrtpc] = tpcpm/GeV;

  tree1ev.trdetpc[tree1ev.ntrtpc] = tpcde;
  tree1ev.trlentpc[tree1ev.ntrtpc] = tpclen;
  tree1ev.trdedxtpc[tree1ev.ntrtpc] = tpcdedx;

  tree1ev.trdedxtrtpc[tree1ev.ntrtpc] = tpcdedxtr;
  //  G4cout<<"dedx:"<<tpclen<<G4endl;

  //  G4cout <<"check fill :"<<tpcdedx<< G4endl;
  tree1ev.trlaytpc[tree1ev.ntrtpc] = tpclay;
  tree1ev.cir_r[tree1ev.ntrtpc] = cir_r;
  tree1ev.cir_x[tree1ev.ntrtpc] = cir_x;
  tree1ev.cir_z[tree1ev.ntrtpc] = cir_z;
  tree1ev.cir_fit[tree1ev.ntrtpc] = cir_fit;

  tree1ev.vtx_flag[tree1ev.ntrtpc] = vtx_flag;
  tree1ev.a_fory[tree1ev.ntrtpc] = a_fory;
  tree1ev.b_fory[tree1ev.ntrtpc] = b_fory;

  tree1ev.ntrtpc++;

  //  G4cout <<"check fill :"<<tpctr<< G4endl;
  //  G4cout <<"check fill11 :"<<tree1ev.ntrtpc<< G4endl;
  //  G4cout <<"check fill11111 :"<<tree1ev.trpxtpc[tree1ev.ntrtpc-1]<< G4endl;
}



void TPCAnaRoot::FillMom(G4double *mom)
{
  G4int i;
  for(i=0;i<3;i++){
    hmom[i]->Fill(mom[i]/GeV);
  }

  tree1ev.pxtpc[tree1ev.nttpc] = mom[0]/GeV;
  tree1ev.pytpc[tree1ev.nttpc] = mom[1]/GeV;
  tree1ev.pztpc[tree1ev.nttpc] = mom[2]/GeV;
  tree1ev.pptpc[tree1ev.nttpc] = sqrt(pow(2,mom[0])+pow(2,mom[1])+pow(2,mom[2]))/GeV;
}

void TPCAnaRoot::FillPos(G4double *pos)
{
  G4int i;
  for(i=0;i<3;i++){
    hpos[i]->Fill(pos[i]/mm);
  }

  tree1ev.xtpc[tree1ev.nttpc] = pos[0]/mm;
  tree1ev.ytpc[tree1ev.nttpc] = pos[1]/mm;
  tree1ev.ztpc[tree1ev.nttpc] = pos[2]/mm;

}

void TPCAnaRoot::FillNtrk(G4int ntrk)
{
  tree1ev.ntrk[tree1ev.nttpc] = ntrk;
}


void TPCAnaRoot::FillPos0(G4double *pos, G4double resoX)
{
  tree1ev.x0tpc[tree1ev.nttpc] = pos[0]/mm;
  tree1ev.y0tpc[tree1ev.nttpc] = pos[1]/mm;
  tree1ev.z0tpc[tree1ev.nttpc] = pos[2]/mm;
  tree1ev.resoX[tree1ev.nttpc] = resoX;
}

void TPCAnaRoot::FillTime(G4double time)
{
  //  printf("time = %f\n",time);
  htime->Fill(time);

}

void TPCAnaRoot::FillBeta(G4double beta)
{
  tree1ev.betatpc[tree1ev.nttpc] = beta;
}

void TPCAnaRoot::FillEdep(G4double edep)
{
  tree1ev.edeptpc[tree1ev.nttpc] = edep;
}

void TPCAnaRoot::FilldEdx(G4double dedx)
{
  tree1ev.dedxtpc[tree1ev.nttpc] = dedx;
}

void TPCAnaRoot::FillsLength(G4double slength)
{
  tree1ev.slengthtpc[tree1ev.nttpc] = slength;
}

void TPCAnaRoot::FillTree()
{
  tree->Fill();
}

void TPCAnaRoot::FillPrimaryParticle(int id, double* x0, double* p0, int pid)
{

  for(int i=0;i<3;i++){
    tree1ev.x0[id][i] = x0[i];
  }

  for(int i=0;i<4;i++){
    tree1ev.p0[id][i] = p0[i];
  }
  tree1ev.pid0[id] = pid;

  tree1ev.p0[id][4] = sqrt(pow(p0[0],2.0)+pow(p0[1],2.0)+pow(p0[2],2.0));
  tree1ev.pt0[id] = sqrt(pow(p0[0],2.0)+pow(p0[1],2.0));
  tree1ev.mass0[id] = p0[4];
  tree1ev.theta0[id] =
    180.0*atan2(sqrt(pow(p0[0],2.0)+pow(p0[1],2.0)),p0[2])/PI;
  tree1ev.npid++;
  //  G4cout<<"p0 test:"<<p0[id][4]<<G4endl;
}

void TPCAnaRoot::FillBeam(double px, double py, double pz)
{

  tree1ev.pg[0] = px;
  tree1ev.pg[1] = py;
  tree1ev.pg[2] = pz;
  tree1ev.pg[3] = sqrt(pow(px,2.0)+ pow(py,2.0) +  pow(pz,2.0));

}


void TPCAnaRoot::FillPrimaryInfo(double mm_d, double mm_p, double theta,
				 double theta_scat, double theta_CM)
{

  tree1ev.mm_d = mm_d;
  tree1ev.mm_p = mm_p;
  tree1ev.theta = theta;
  tree1ev.theta_scat = theta_scat;
  tree1ev.theta_CM = theta_CM;

  tree1ev.mm = mm;
}

void TPCAnaRoot::FillNumOfK(int HitNum_K, int HitNumAC_K, int HitNumNBAR_K,
			    int HitNumDC_K, int HitNumCH_K, int HitNumFTOF_K,
			    int HitNumScint_K, int HitNumTarget_K)
{
  tree1ev.HitNum_K = HitNum_K;
  tree1ev.HitNumAC_K = HitNumAC_K;
  tree1ev.HitNumNBAR_K = HitNumNBAR_K;
  tree1ev.HitNumDC_K = HitNumDC_K;
  tree1ev.HitNumCH_K = HitNumCH_K;
  tree1ev.HitNumFTOF_K = HitNumFTOF_K;
  tree1ev.HitNumScint_K = HitNumScint_K;
  tree1ev.HitNumTarget_K = HitNumTarget_K;
}



void TPCAnaRoot::FillNumOfp(int HitNum_p, int HitNumAC_p, int HitNumNBAR_p,
			    int HitNumDC_p, int HitNumCH_p, int HitNumFTOF_p,
			    int HitNumScint_p, int HitNumTarget_p)
{
  tree1ev.HitNum_p = HitNum_p;
  tree1ev.HitNumAC_p = HitNumAC_p;
  tree1ev.HitNumNBAR_p = HitNumNBAR_p;
  tree1ev.HitNumDC_p = HitNumDC_p;
  tree1ev.HitNumCH_p = HitNumCH_p;
  tree1ev.HitNumFTOF_p = HitNumFTOF_p;
  tree1ev.HitNumScint_p = HitNumScint_p;
  tree1ev.HitNumTarget_p = HitNumTarget_p;
}


void TPCAnaRoot::FillScintData(G4double time, G4double *pos, G4double *mom,
			       G4int tid, G4int pid, G4int did,G4double mass, G4int qq,G4int parentID,
			       G4double scintvtxpx,G4double scintvtxpy,G4double scintvtxpz,
			       G4double scintvtxx,G4double scintvtxy,G4double scintvtxz, G4double tlength
)
{

  //  tree1ev.tofsc[tree1ev.ntsc] = time+CLHEP::RandGauss::shoot(0.,0.130);
  tree1ev.tofsc[tree1ev.ntsc] = time;
  tree1ev.scpID[tree1ev.ntsc] = parentID;
  //  G4cout<<"particle:parentid :"<<pid<<":"<<parentID<<G4endl;
  tree1ev.xsc[tree1ev.ntsc] = pos[0];
  tree1ev.ysc[tree1ev.ntsc] = pos[1];
  tree1ev.zsc[tree1ev.ntsc] = pos[2];

  tree1ev.pxsc[tree1ev.ntsc] = mom[0]/GeV;
  tree1ev.pysc[tree1ev.ntsc] = mom[1]/GeV;
  tree1ev.pzsc[tree1ev.ntsc] = mom[2]/GeV;
  tree1ev.ppsc[tree1ev.ntsc] = sqrt(pow(mom[0]/GeV,2)+pow(mom[1]/GeV,2)+pow(mom[2]/GeV,2));

  tree1ev.tidsc[tree1ev.ntsc] = tid;
  tree1ev.pidsc[tree1ev.ntsc] = pid;
  tree1ev.masssc[tree1ev.ntsc] = mass/GeV;
  tree1ev.qqsc[tree1ev.ntsc] = qq;
  tree1ev.didsc[tree1ev.ntsc] = did;

  tree1ev.trvtxpxscint[tree1ev.ntsc] = scintvtxpx/GeV;
  tree1ev.trvtxpyscint[tree1ev.ntsc] = scintvtxpy/GeV;
  tree1ev.trvtxpzscint[tree1ev.ntsc] = scintvtxpz/GeV;
  tree1ev.trvtxppscint[tree1ev.ntsc] = sqrt(pow(scintvtxpx,2)+pow(scintvtxpy,2)+pow(scintvtxpz,2))/GeV;

  tree1ev.trvtxxscint[tree1ev.ntsc] = scintvtxx;
  tree1ev.trvtxyscint[tree1ev.ntsc] = scintvtxy;
  tree1ev.trvtxzscint[tree1ev.ntsc] = scintvtxz;

  //  tree1ev.lengthsc[tree1ev.ntsc] = tlength+CLHEP::RandGauss::shoot(0.,20.);
  tree1ev.lengthsc[tree1ev.ntsc] = tlength;
  tree1ev.ntsc++;

}


void TPCAnaRoot::FillACData(G4double time, G4double *pos, G4double *mom,
			       G4int tid, G4int pid, G4int did,G4double mass, G4int qq,G4int parentID,
			       G4double acvtxpx,G4double acvtxpy,G4double acvtxpz,
			       G4double acvtxx,G4double acvtxy,G4double acvtxz, G4double tlength
)
{

  tree1ev.tofac[tree1ev.ntac] = time;
  tree1ev.acpID[tree1ev.ntac] = parentID;
  //  G4cout<<"particle:parentid :"<<pid<<":"<<parentID<<G4endl;
  tree1ev.xac[tree1ev.ntac] = pos[0];
  tree1ev.yac[tree1ev.ntac] = pos[1];
  tree1ev.zac[tree1ev.ntac] = pos[2];

  tree1ev.pxac[tree1ev.ntac] = mom[0]/GeV;
  tree1ev.pyac[tree1ev.ntac] = mom[1]/GeV;
  tree1ev.pzac[tree1ev.ntac] = mom[2]/GeV;
  tree1ev.ppac[tree1ev.ntac] = sqrt(pow(mom[0]/GeV,2)+pow(mom[1]/GeV,2)+pow(mom[2]/GeV,2));

  tree1ev.tidac[tree1ev.ntac] = tid;
  tree1ev.pidac[tree1ev.ntac] = pid;
  tree1ev.massac[tree1ev.ntac] = mass/GeV;
  tree1ev.qqac[tree1ev.ntac] = qq;
  tree1ev.didac[tree1ev.ntac] = did;

  tree1ev.trvtxpxac[tree1ev.ntac] = acvtxpx/GeV;
  tree1ev.trvtxpyac[tree1ev.ntac] = acvtxpy/GeV;
  tree1ev.trvtxpzac[tree1ev.ntac] = acvtxpz/GeV;
  tree1ev.trvtxppac[tree1ev.ntac] = sqrt(pow(acvtxpx,2)+pow(acvtxpy,2)+pow(acvtxpz,2))/GeV;

  tree1ev.trvtxxac[tree1ev.ntac] = acvtxx;
  tree1ev.trvtxyac[tree1ev.ntac] = acvtxy;
  tree1ev.trvtxzac[tree1ev.ntac] = acvtxz;

  tree1ev.lengthac[tree1ev.ntac] = tlength+CLHEP::RandGauss::shoot(0.,20.);

  tree1ev.ntac++;
}





void TPCAnaRoot::FillNBARData(G4double time, G4double *pos, G4double *mom,
			       G4int tid, G4int pid, G4int did,G4double mass, G4int qq,G4int parentID,
			       G4double nbarvtxpx,G4double nbarvtxpy,G4double nbarvtxpz,
			       G4double nbarvtxx,G4double nbarvtxy,G4double nbarvtxz, G4double tlength
)
{

  tree1ev.tofnbar[tree1ev.ntnbar] = time;
  tree1ev.nbarpID[tree1ev.ntnbar] = parentID;
  //  G4cout<<"particle:parentid :"<<pid<<":"<<parentID<<G4endl;
  tree1ev.xnbar[tree1ev.ntnbar] = pos[0];
  tree1ev.ynbar[tree1ev.ntnbar] = pos[1];
  tree1ev.znbar[tree1ev.ntnbar] = pos[2];

  tree1ev.pxnbar[tree1ev.ntnbar] = mom[0]/GeV;
  tree1ev.pynbar[tree1ev.ntnbar] = mom[1]/GeV;
  tree1ev.pznbar[tree1ev.ntnbar] = mom[2]/GeV;
  tree1ev.ppnbar[tree1ev.ntnbar] = sqrt(pow(mom[0]/GeV,2)+pow(mom[1]/GeV,2)+pow(mom[2]/GeV,2));

  tree1ev.tidnbar[tree1ev.ntnbar] = tid;
  tree1ev.pidnbar[tree1ev.ntnbar] = pid;
  tree1ev.massnbar[tree1ev.ntnbar] = mass/GeV;
  tree1ev.qqnbar[tree1ev.ntnbar] = qq;
  tree1ev.didnbar[tree1ev.ntnbar] = did;

  tree1ev.trvtxpxnbar[tree1ev.ntnbar] = nbarvtxpx/GeV;
  tree1ev.trvtxpynbar[tree1ev.ntnbar] = nbarvtxpy/GeV;
  tree1ev.trvtxpznbar[tree1ev.ntnbar] = nbarvtxpz/GeV;
  tree1ev.trvtxppnbar[tree1ev.ntnbar] = sqrt(pow(nbarvtxpx,2)+pow(nbarvtxpy,2)+pow(nbarvtxpz,2))/GeV;

  tree1ev.trvtxxnbar[tree1ev.ntnbar] = nbarvtxx;
  tree1ev.trvtxynbar[tree1ev.ntnbar] = nbarvtxy;
  tree1ev.trvtxznbar[tree1ev.ntnbar] = nbarvtxz;

  tree1ev.lengthnbar[tree1ev.ntnbar] = tlength+CLHEP::RandGauss::shoot(0.,20.);

  tree1ev.ntnbar++;
}




void TPCAnaRoot::FillDCData(G4double time, G4double *pos, G4double *mom,
			       G4int tid, G4int pid, G4int did,G4double mass, G4int qq,G4int parentID,
			       G4double vtxpx,G4double vtxpy,G4double vtxpz,
			       G4double vtxx,G4double vtxy,G4double vtxz, G4double tlength
)
{

  tree1ev.tofdc[tree1ev.ntdc] = time;
  tree1ev.dcpID[tree1ev.ntdc] = parentID;
  //  G4cout<<"particle:parentid :"<<pid<<":"<<parentID<<G4endl;
  tree1ev.xdc[tree1ev.ntdc] = pos[0];
  tree1ev.ydc[tree1ev.ntdc] = pos[1];
  tree1ev.zdc[tree1ev.ntdc] = pos[2];

  tree1ev.pxdc[tree1ev.ntdc] = mom[0]/GeV;
  tree1ev.pydc[tree1ev.ntdc] = mom[1]/GeV;
  tree1ev.pzdc[tree1ev.ntdc] = mom[2]/GeV;
  tree1ev.ppdc[tree1ev.ntdc] = sqrt(pow(mom[0]/GeV,2)+pow(mom[1]/GeV,2)+pow(mom[2]/GeV,2));

  tree1ev.tiddc[tree1ev.ntdc] = tid;
  tree1ev.piddc[tree1ev.ntdc] = pid;
  tree1ev.massdc[tree1ev.ntdc] = mass/GeV;
  tree1ev.qqdc[tree1ev.ntdc] = qq;
  tree1ev.diddc[tree1ev.ntdc] = did;

  tree1ev.trvtxpxdc[tree1ev.ntdc] = vtxpx/GeV;
  tree1ev.trvtxpydc[tree1ev.ntdc] = vtxpy/GeV;
  tree1ev.trvtxpzdc[tree1ev.ntdc] = vtxpz/GeV;
  tree1ev.trvtxppdc[tree1ev.ntdc] = sqrt(pow(vtxpx,2)+pow(vtxpy,2)+pow(vtxpz,2))/GeV;

  tree1ev.trvtxxdc[tree1ev.ntdc] = vtxx;
  tree1ev.trvtxydc[tree1ev.ntdc] = vtxy;
  tree1ev.trvtxzdc[tree1ev.ntdc] = vtxz;

  tree1ev.lengthdc[tree1ev.ntdc] = tlength+CLHEP::RandGauss::shoot(0.,20.);

  tree1ev.ntdc++;
}


void TPCAnaRoot::FillCHData(G4double time, G4double *pos, G4double *mom,
			       G4int tid, G4int pid, G4int did,G4double mass, G4int qq,G4int parentID,
			       G4double vtxpx,G4double vtxpy,G4double vtxpz,
			       G4double vtxx,G4double vtxy,G4double vtxz, G4double tlength
)
{

  tree1ev.tofch[tree1ev.ntch] = time;
  tree1ev.chpID[tree1ev.ntch] = parentID;
  //  G4cout<<"particle:parentid :"<<pid<<":"<<parentID<<G4endl;
  tree1ev.xch[tree1ev.ntch] = pos[0];
  tree1ev.ych[tree1ev.ntch] = pos[1];
  tree1ev.zch[tree1ev.ntch] = pos[2];

  tree1ev.pxch[tree1ev.ntch] = mom[0]/GeV;
  tree1ev.pych[tree1ev.ntch] = mom[1]/GeV;
  tree1ev.pzch[tree1ev.ntch] = mom[2]/GeV;
  tree1ev.ppch[tree1ev.ntch] = sqrt(pow(mom[0]/GeV,2)+pow(mom[1]/GeV,2)+pow(mom[2]/GeV,2));

  tree1ev.tidch[tree1ev.ntch] = tid;
  tree1ev.pidch[tree1ev.ntch] = pid;
  tree1ev.massch[tree1ev.ntch] = mass/GeV;
  tree1ev.qqch[tree1ev.ntch] = qq;
  tree1ev.didch[tree1ev.ntch] = did;

  tree1ev.trvtxpxch[tree1ev.ntch] = vtxpx/GeV;
  tree1ev.trvtxpych[tree1ev.ntch] = vtxpy/GeV;
  tree1ev.trvtxpzch[tree1ev.ntch] = vtxpz/GeV;
  tree1ev.trvtxppch[tree1ev.ntch] = sqrt(pow(vtxpx,2)+pow(vtxpy,2)+pow(vtxpz,2))/GeV;

  tree1ev.trvtxxch[tree1ev.ntch] = vtxx;
  tree1ev.trvtxych[tree1ev.ntch] = vtxy;
  tree1ev.trvtxzch[tree1ev.ntch] = vtxz;

  tree1ev.lengthch[tree1ev.ntch] = tlength+CLHEP::RandGauss::shoot(0.,20.);

  tree1ev.ntch++;
}



void TPCAnaRoot::FillFTOFData(G4double time, G4double *pos, G4double *mom,
			       G4int tid, G4int pid, G4int did,G4double mass, G4int qq,G4int parentID,
			       G4double vtxpx,G4double vtxpy,G4double vtxpz,
			       G4double vtxx,G4double vtxy,G4double vtxz, G4double tlength
)
{

  tree1ev.tofftof[tree1ev.ntftof] = time;
  tree1ev.ftofpID[tree1ev.ntftof] = parentID;
  //  G4cout<<"particle:parentid :"<<pid<<":"<<parentID<<G4endl;
  tree1ev.xftof[tree1ev.ntftof] = pos[0];
  tree1ev.yftof[tree1ev.ntftof] = pos[1];
  tree1ev.zftof[tree1ev.ntftof] = pos[2];

  tree1ev.pxftof[tree1ev.ntftof] = mom[0]/GeV;
  tree1ev.pyftof[tree1ev.ntftof] = mom[1]/GeV;
  tree1ev.pzftof[tree1ev.ntftof] = mom[2]/GeV;
  tree1ev.ppftof[tree1ev.ntftof] = sqrt(pow(mom[0]/GeV,2)+pow(mom[1]/GeV,2)+pow(mom[2]/GeV,2));

  tree1ev.tidftof[tree1ev.ntftof] = tid;
  tree1ev.pidftof[tree1ev.ntftof] = pid;
  tree1ev.massftof[tree1ev.ntftof] = mass/GeV;
  tree1ev.qqftof[tree1ev.ntftof] = qq;
  tree1ev.didftof[tree1ev.ntftof] = did;

  tree1ev.trvtxpxftof[tree1ev.ntftof] = vtxpx/GeV;
  tree1ev.trvtxpyftof[tree1ev.ntftof] = vtxpy/GeV;
  tree1ev.trvtxpzftof[tree1ev.ntftof] = vtxpz/GeV;
  tree1ev.trvtxppftof[tree1ev.ntftof] = sqrt(pow(vtxpx,2)+pow(vtxpy,2)+pow(vtxpz,2))/GeV;

  tree1ev.trvtxxftof[tree1ev.ntftof] = vtxx;
  tree1ev.trvtxyftof[tree1ev.ntftof] = vtxy;
  tree1ev.trvtxzftof[tree1ev.ntftof] = vtxz;
  tree1ev.lengthftof[tree1ev.ntftof] = tlength;

  tree1ev.ntftof++;
  //  std::cout<<  tree1ev.ntftof<<std::endl;
}



/*void TPCAnaRoot::FillFDCData(G4double time, G4double *pos, G4double *mom,
			  G4int tid, G4int pid, G4int did)
{
  tree1ev.toffdc[tree1ev.ntfdc] = time;

  tree1ev.xfdc[tree1ev.ntfdc] = pos[0];
  tree1ev.yfdc[tree1ev.ntfdc] = pos[1];
  tree1ev.zfdc[tree1ev.ntfdc] = pos[2];

  tree1ev.pxfdc[tree1ev.ntfdc] = mom[0]/GeV;
  tree1ev.pyfdc[tree1ev.ntfdc] = mom[1]/GeV;
  tree1ev.pzfdc[tree1ev.ntfdc] = mom[2]/GeV;

  tree1ev.tidfdc[tree1ev.ntfdc] = tid;
  tree1ev.pidfdc[tree1ev.ntfdc] = pid;
  tree1ev.didfdc[tree1ev.ntfdc] = did;

  tree1ev.ntfdc++;

}
*/

void TPCAnaRoot::FillTargetData(G4int nhit, G4int particleid, G4int parentid, G4int trackid, G4ThreeVector pos, G4ThreeVector vtx)
{
  //  tree1ev.targethits[tree1ev.ntfdc] = nhit;
  tree1ev.targetpid[tree1ev.targethits] = particleid;
  tree1ev.targetparentid[tree1ev.targethits] = parentid;
  tree1ev.targettid[tree1ev.targethits] = trackid;
  tree1ev.targetpos[tree1ev.targethits][0] = pos.getX();
  tree1ev.targetpos[tree1ev.targethits][1] = pos.getY();
  tree1ev.targetpos[tree1ev.targethits][2] = pos.getZ();
  tree1ev.targetvtx[tree1ev.targethits][0] = vtx.getX();
  tree1ev.targetvtx[tree1ev.targethits][1] = vtx.getY();
  tree1ev.targetvtx[tree1ev.targethits][2] = vtx.getZ();

  tree1ev.targethits++;
  //  G4cout<<tree1ev.targethits<<G4endl;

}

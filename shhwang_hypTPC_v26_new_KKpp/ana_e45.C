//void ana_e45(Int_t beam, char generator[10]){
void ana_e45(){
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Fri Oct 11 08:34:35 2013 by ROOT version5.34/05)
//   from TTree tree/GEANT4 simulation for hypTPC
//   found on file: ./grun_e45_elastic.root
//////////////////////////////////////////////////////////

  char name1[100];
  char name2[100];
  char name3[100];
  char name4[100];
    char generator[10]="mode4";
    Int_t beam=8;
  //  Int_t beam=20;
//Reset ROOT and connect tree file
   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./grun_123456.root");
   if (!f) {
     //           f = new TFile("./grun_e45_elastic.root");
     //     f = new TFile("./grun1.root");
     //     f = new TFile(Form("./raw_root/grun_e45_elastic_%02d_%s.root",beam,generator));     
     f = new TFile("grun1.root");     
     //     f = new TFile("./grun_e45_elastic_20_mode1.root");     
     //     f = new TFile("./grun_e45_elastic_pip.root");
     //      f = new TFile("../shhwang_hypTPC_v18.1/root_40/grun_gen24_1035.root");
   }
   //    f->GetObject("tree",tree);
  TTree *tree = (TTree*)gDirectory->Get("tree");

  cout<<"ana root file:"<<Form("./raw_root/grun_e45_elastic_%02d_%s.root",beam,generator)<<endl;
//Declaration of leaves types
   Int_t           ev;
   Double_t        pg[4];
   Int_t           npid;
   Double_t        x0[10][3];
   Double_t        p0[10][5];
   Double_t        pt0[10];
   Double_t        mass0[10];
   Double_t        theta0[10];
   Int_t           gen;
   Int_t           mode;
   Int_t           nttpc;
   Int_t           ntrk[500];
   Int_t           ititpc[500];
   Int_t           idtpc[500];
   Double_t        xtpc[500];
   Double_t        ytpc[500];
   Double_t        ztpc[500];
   Double_t        x0tpc[500];
   Double_t        y0tpc[500];
   Double_t        z0tpc[500];
   Double_t        resoX[500];
   Double_t        pxtpc[500];
   Double_t        pytpc[500];
   Double_t        pztpc[500];
   Double_t        pptpc[500];
   Double_t        masstpc[500];
   Double_t        betatpc[500];
   Double_t        edeptpc[500];
   Double_t        dedxtpc[500];
   Double_t        slengthtpc[500];
   Int_t           laytpc[500];
   Int_t           rowtpc[500];
   Int_t           parentID[500];
   Int_t           nthlay[500];
   Int_t           nthpad[500];
   Int_t           laypad[500][40][250];
   Int_t           ntrtpc;
   Double_t        trpmtpc[10];
   Int_t           trqqtpc[10];
   Double_t        trpxtpc[10];
   Double_t        trpytpc[10];
   Double_t        trpztpc[10];
   Double_t        trpptpc[10];
   Double_t        trpttpc[10];
   Double_t        trpxtpcfit[10];
   Double_t        trpytpcfit[10];
   Double_t        trpztpcfit[10];
   Double_t        trpptpcfit[10];
   Double_t        trpttpcfit[10];
   Double_t        trvtxpxtpc[10];
   Double_t        trvtxpytpc[10];
   Double_t        trvtxpztpc[10];
   Double_t        trvtxpptpc[10];
   Double_t        trvtxxtpc[10];
   Double_t        trvtxytpc[10];
   Double_t        trvtxztpc[10];
   Double_t        trvtxxtpcfit[10];
   Double_t        trvtxytpcfit[10];
   Double_t        trvtxztpcfit[10];
   Double_t        trdetpc[10];
   Double_t        trlentpc[10];
   Double_t        trdedxtpc[10];
   Double_t        trdedxtrtpc[10];
   Int_t           trlaytpc[10];
   Double_t        cir_r[10];
   Double_t        cir_x[10];
   Double_t        cir_z[10];
   Double_t        cir_fit[10];
   Int_t           ntsc;
   Int_t           tidsc[40];
   Int_t           pidsc[40];
   Int_t           didsc[40];
   Double_t        masssc[40];
   Int_t           qqsc[40];
   Double_t        xsc[40];
   Double_t        ysc[40];
   Double_t        zsc[40];
   Double_t        pxsc[40];
   Double_t        pysc[40];
   Double_t        pzsc[40];
   Double_t        ppsc[40];
   Double_t        tofsc[40];
   Int_t           scpID[40];
   Double_t        trvtxpxscint[40];
   Double_t        trvtxpyscint[40];
   Double_t        trvtxpzscint[40];
   Double_t        trvtxppscint[40];
   Double_t        trvtxxscint[40];
   Double_t        trvtxyscint[40];
   Double_t        trvtxzscint[40];
   Double_t        lengthsc[40];
   Int_t           ntac;
   Int_t           tidac[0];
   Int_t           pidac[0];
   Int_t           didac[0];
   Double_t        massac[0];
   Int_t           qqac[0];
   Double_t        xac[0];
   Double_t        yac[0];
   Double_t        zac[0];
   Double_t        pxac[0];
   Double_t        pyac[0];
   Double_t        pzac[0];
   Double_t        ppac[0];
   Double_t        tofac[0];
   Int_t           acpID[0];
   Double_t        trvtxpxac[0];
   Double_t        trvtxpyac[0];
   Double_t        trvtxpzac[0];
   Double_t        trvtxppac[0];
   Double_t        trvtxxac[0];
   Double_t        trvtxyac[0];
   Double_t        trvtxzac[0];
   Double_t        lengthac[0];
   Int_t           ntnbar;
   Int_t           tidnbar[100];
   Int_t           pidnbar[100];
   Int_t           didnbar[100];
   Double_t        massnbar[100];
   Int_t           qqnbar[100];
   Double_t        xnbar[100];
   Double_t        ynbar[100];
   Double_t        znbar[100];
   Double_t        pxnbar[100];
   Double_t        pynbar[100];
   Double_t        pznbar[100];
   Double_t        ppnbar[100];
   Double_t        tofnbar[100];
   Int_t           nbarpID[100];
   Double_t        trvtxpxnbar[100];
   Double_t        trvtxpynbar[100];
   Double_t        trvtxpznbar[100];
   Double_t        trvtxppnbar[100];
   Double_t        trvtxxnbar[100];
   Double_t        trvtxynbar[100];
   Double_t        trvtxznbar[100];
   Double_t        lengthnbar[100];
   Int_t           ntch;
   Int_t           tidch[0];
   Int_t           pidch[0];
   Int_t           didch[0];
   Double_t        massch[0];
   Int_t           qqch[0];
   Double_t        xch[0];
   Double_t        ych[0];
   Double_t        zch[0];
   Double_t        pxch[0];
   Double_t        pych[0];
   Double_t        pzch[0];
   Double_t        ppch[0];
   Double_t        tofch[0];
   Int_t           chpID[0];
   Double_t        trvtxpxch[0];
   Double_t        trvtxpych[0];
   Double_t        trvtxpzch[0];
   Double_t        trvtxppch[0];
   Double_t        trvtxxch[0];
   Double_t        trvtxych[0];
   Double_t        trvtxzch[0];
   Double_t        lengthch[0];
   Int_t           ntftof;
   Int_t           tidftof[0];
   Int_t           pidftof[0];
   Int_t           didftof[0];
   Double_t        massftof[0];
   Int_t           qqftof[0];
   Double_t        xftof[0];
   Double_t        yftof[0];
   Double_t        zftof[0];
   Double_t        pxftof[0];
   Double_t        pyftof[0];
   Double_t        pzftof[0];
   Double_t        ppftof[0];
   Double_t        tofftof[0];
   Int_t           ftofpID[0];
   Double_t        trvtxpxftof[0];
   Double_t        trvtxpyftof[0];
   Double_t        trvtxpzftof[0];
   Double_t        trvtxppftof[0];
   Double_t        trvtxxftof[0];
   Double_t        trvtxyftof[0];
   Double_t        trvtxzftof[0];
   Double_t        lengthftof[0];
   Int_t           ntdc;
   Int_t           tiddc[0];
   Int_t           piddc[0];
   Int_t           diddc[0];
   Double_t        massdc[0];
   Int_t           qqdc[0];
   Double_t        xdc[0];
   Double_t        ydc[0];
   Double_t        zdc[0];
   Double_t        pxdc[0];
   Double_t        pydc[0];
   Double_t        pzdc[0];
   Double_t        ppdc[0];
   Double_t        tofdc[0];
   Int_t           dcpID[0];
   Double_t        trvtxpxdc[0];
   Double_t        trvtxpydc[0];
   Double_t        trvtxpzdc[0];
   Double_t        trvtxppdc[0];
   Double_t        trvtxxdc[0];
   Double_t        trvtxydc[0];
   Double_t        trvtxzdc[0];
   Double_t        lengthdc[0];
   Int_t           targethits;
   Int_t           targetpid[10];
   Int_t           targetparentid[10];
   Int_t           targettid[10];
   Double_t        targetpos[10][3];
   Double_t        targetvtx[10][3];

   // Set branch addresses.
   tree->SetBranchAddress("ev",&ev);
   tree->SetBranchAddress("pg",pg);
   tree->SetBranchAddress("npid",&npid);
   tree->SetBranchAddress("x0",x0);
   tree->SetBranchAddress("p0",p0);
   tree->SetBranchAddress("pt0",pt0);
   tree->SetBranchAddress("mass0",mass0);
   tree->SetBranchAddress("theta0",theta0);
   tree->SetBranchAddress("gen",&gen);
   tree->SetBranchAddress("mode",&mode);
   tree->SetBranchAddress("nttpc",&nttpc);
   tree->SetBranchAddress("ntrk",ntrk);
   tree->SetBranchAddress("ititpc",ititpc);
   tree->SetBranchAddress("idtpc",idtpc);
   tree->SetBranchAddress("xtpc",xtpc);
   tree->SetBranchAddress("ytpc",ytpc);
   tree->SetBranchAddress("ztpc",ztpc);
   tree->SetBranchAddress("x0tpc",x0tpc);
   tree->SetBranchAddress("y0tpc",y0tpc);
   tree->SetBranchAddress("z0tpc",z0tpc);
   tree->SetBranchAddress("resoX",resoX);
   tree->SetBranchAddress("pxtpc",pxtpc);
   tree->SetBranchAddress("pytpc",pytpc);
   tree->SetBranchAddress("pztpc",pztpc);
   tree->SetBranchAddress("pptpc",pptpc);
   tree->SetBranchAddress("masstpc",masstpc);
   tree->SetBranchAddress("betatpc",betatpc);
   tree->SetBranchAddress("edeptpc",edeptpc);
   tree->SetBranchAddress("dedxtpc",dedxtpc);
   tree->SetBranchAddress("slengthtpc",slengthtpc);
   tree->SetBranchAddress("laytpc",laytpc);
   tree->SetBranchAddress("rowtpc",rowtpc);
   tree->SetBranchAddress("parentID",parentID);
   tree->SetBranchAddress("nthlay",nthlay);
   tree->SetBranchAddress("nthpad",nthpad);
   tree->SetBranchAddress("laypad",laypad);
   tree->SetBranchAddress("ntrtpc",&ntrtpc);
   tree->SetBranchAddress("trpmtpc",trpmtpc);
   tree->SetBranchAddress("trqqtpc",trqqtpc);
   tree->SetBranchAddress("trpxtpc",trpxtpc);
   tree->SetBranchAddress("trpytpc",trpytpc);
   tree->SetBranchAddress("trpztpc",trpztpc);
   tree->SetBranchAddress("trpptpc",trpptpc);
   tree->SetBranchAddress("trpttpc",trpttpc);
   tree->SetBranchAddress("trpxtpcfit",trpxtpcfit);
   tree->SetBranchAddress("trpytpcfit",trpytpcfit);
   tree->SetBranchAddress("trpztpcfit",trpztpcfit);
   tree->SetBranchAddress("trpptpcfit",trpptpcfit);
   tree->SetBranchAddress("trpttpcfit",trpttpcfit);
   tree->SetBranchAddress("trvtxpxtpc",trvtxpxtpc);
   tree->SetBranchAddress("trvtxpytpc",trvtxpytpc);
   tree->SetBranchAddress("trvtxpztpc",trvtxpztpc);
   tree->SetBranchAddress("trvtxpptpc",trvtxpptpc);
   tree->SetBranchAddress("trvtxxtpc",trvtxxtpc);
   tree->SetBranchAddress("trvtxytpc",trvtxytpc);
   tree->SetBranchAddress("trvtxztpc",trvtxztpc);
   tree->SetBranchAddress("trvtxxtpcfit",trvtxxtpcfit);
   tree->SetBranchAddress("trvtxytpcfit",trvtxytpcfit);
   tree->SetBranchAddress("trvtxztpcfit",trvtxztpcfit);
   tree->SetBranchAddress("trdetpc",trdetpc);
   tree->SetBranchAddress("trlentpc",trlentpc);
   tree->SetBranchAddress("trdedxtpc",trdedxtpc);
   tree->SetBranchAddress("trdedxtrtpc",trdedxtrtpc);
   tree->SetBranchAddress("trlaytpc",trlaytpc);
   tree->SetBranchAddress("cir_r",cir_r);
   tree->SetBranchAddress("cir_x",cir_x);
   tree->SetBranchAddress("cir_z",cir_z);
   tree->SetBranchAddress("cir_fit",cir_fit);
   tree->SetBranchAddress("ntsc",&ntsc);
   tree->SetBranchAddress("tidsc",tidsc);
   tree->SetBranchAddress("pidsc",pidsc);
   tree->SetBranchAddress("didsc",didsc);
   tree->SetBranchAddress("masssc",masssc);
   tree->SetBranchAddress("qqsc",qqsc);
   tree->SetBranchAddress("xsc",xsc);
   tree->SetBranchAddress("ysc",ysc);
   tree->SetBranchAddress("zsc",zsc);
   tree->SetBranchAddress("pxsc",pxsc);
   tree->SetBranchAddress("pysc",pysc);
   tree->SetBranchAddress("pzsc",pzsc);
   tree->SetBranchAddress("ppsc",ppsc);
   tree->SetBranchAddress("tofsc",tofsc);
   tree->SetBranchAddress("scpID",scpID);
   tree->SetBranchAddress("trvtxpxscint",trvtxpxscint);
   tree->SetBranchAddress("trvtxpyscint",trvtxpyscint);
   tree->SetBranchAddress("trvtxpzscint",trvtxpzscint);
   tree->SetBranchAddress("trvtxppscint",trvtxppscint);
   tree->SetBranchAddress("trvtxxscint",trvtxxscint);
   tree->SetBranchAddress("trvtxyscint",trvtxyscint);
   tree->SetBranchAddress("trvtxzscint",trvtxzscint);
   tree->SetBranchAddress("lengthsc",lengthsc);
   tree->SetBranchAddress("ntac",&ntac);
   tree->SetBranchAddress("tidac",&tidac);
   tree->SetBranchAddress("pidac",&pidac);
   tree->SetBranchAddress("didac",&didac);
   tree->SetBranchAddress("massac",&massac);
   tree->SetBranchAddress("qqac",&qqac);
   tree->SetBranchAddress("xac",&xac);
   tree->SetBranchAddress("yac",&yac);
   tree->SetBranchAddress("zac",&zac);
   tree->SetBranchAddress("pxac",&pxac);
   tree->SetBranchAddress("pyac",&pyac);
   tree->SetBranchAddress("pzac",&pzac);
   tree->SetBranchAddress("ppac",&ppac);
   tree->SetBranchAddress("tofac",&tofac);
   tree->SetBranchAddress("acpID",&acpID);
   tree->SetBranchAddress("trvtxpxac",&trvtxpxac);
   tree->SetBranchAddress("trvtxpyac",&trvtxpyac);
   tree->SetBranchAddress("trvtxpzac",&trvtxpzac);
   tree->SetBranchAddress("trvtxppac",&trvtxppac);
   tree->SetBranchAddress("trvtxxac",&trvtxxac);
   tree->SetBranchAddress("trvtxyac",&trvtxyac);
   tree->SetBranchAddress("trvtxzac",&trvtxzac);
   tree->SetBranchAddress("lengthac",&lengthac);
   tree->SetBranchAddress("ntnbar",&ntnbar);
   tree->SetBranchAddress("tidnbar",&tidnbar);
   tree->SetBranchAddress("pidnbar",&pidnbar);
   tree->SetBranchAddress("didnbar",&didnbar);
   tree->SetBranchAddress("massnbar",&massnbar);
   tree->SetBranchAddress("qqnbar",&qqnbar);
   tree->SetBranchAddress("xnbar",&xnbar);
   tree->SetBranchAddress("ynbar",&ynbar);
   tree->SetBranchAddress("znbar",&znbar);
   tree->SetBranchAddress("pxnbar",&pxnbar);
   tree->SetBranchAddress("pynbar",&pynbar);
   tree->SetBranchAddress("pznbar",&pznbar);
   tree->SetBranchAddress("ppnbar",&ppnbar);
   tree->SetBranchAddress("tofnbar",&tofnbar);
   tree->SetBranchAddress("nbarpID",&nbarpID);
   tree->SetBranchAddress("trvtxpxnbar",&trvtxpxnbar);
   tree->SetBranchAddress("trvtxpynbar",&trvtxpynbar);
   tree->SetBranchAddress("trvtxpznbar",&trvtxpznbar);
   tree->SetBranchAddress("trvtxppnbar",&trvtxppnbar);
   tree->SetBranchAddress("trvtxxnbar",&trvtxxnbar);
   tree->SetBranchAddress("trvtxynbar",&trvtxynbar);
   tree->SetBranchAddress("trvtxznbar",&trvtxznbar);
   tree->SetBranchAddress("lengthnbar",&lengthnbar);
   tree->SetBranchAddress("ntch",&ntch);
   tree->SetBranchAddress("tidch",&tidch);
   tree->SetBranchAddress("pidch",&pidch);
   tree->SetBranchAddress("didch",&didch);
   tree->SetBranchAddress("massch",&massch);
   tree->SetBranchAddress("qqch",&qqch);
   tree->SetBranchAddress("xch",&xch);
   tree->SetBranchAddress("ych",&ych);
   tree->SetBranchAddress("zch",&zch);
   tree->SetBranchAddress("pxch",&pxch);
   tree->SetBranchAddress("pych",&pych);
   tree->SetBranchAddress("pzch",&pzch);
   tree->SetBranchAddress("ppch",&ppch);
   tree->SetBranchAddress("tofch",&tofch);
   tree->SetBranchAddress("chpID",&chpID);
   tree->SetBranchAddress("trvtxpxch",&trvtxpxch);
   tree->SetBranchAddress("trvtxpych",&trvtxpych);
   tree->SetBranchAddress("trvtxpzch",&trvtxpzch);
   tree->SetBranchAddress("trvtxppch",&trvtxppch);
   tree->SetBranchAddress("trvtxxch",&trvtxxch);
   tree->SetBranchAddress("trvtxych",&trvtxych);
   tree->SetBranchAddress("trvtxzch",&trvtxzch);
   tree->SetBranchAddress("lengthch",&lengthch);
   tree->SetBranchAddress("ntftof",&ntftof);
   tree->SetBranchAddress("tidftof",&tidftof);
   tree->SetBranchAddress("pidftof",&pidftof);
   tree->SetBranchAddress("didftof",&didftof);
   tree->SetBranchAddress("massftof",&massftof);
   tree->SetBranchAddress("qqftof",&qqftof);
   tree->SetBranchAddress("xftof",&xftof);
   tree->SetBranchAddress("yftof",&yftof);
   tree->SetBranchAddress("zftof",&zftof);
   tree->SetBranchAddress("pxftof",&pxftof);
   tree->SetBranchAddress("pyftof",&pyftof);
   tree->SetBranchAddress("pzftof",&pzftof);
   tree->SetBranchAddress("ppftof",&ppftof);
   tree->SetBranchAddress("tofftof",&tofftof);
   tree->SetBranchAddress("ftofpID",&ftofpID);
   tree->SetBranchAddress("trvtxpxftof",&trvtxpxftof);
   tree->SetBranchAddress("trvtxpyftof",&trvtxpyftof);
   tree->SetBranchAddress("trvtxpzftof",&trvtxpzftof);
   tree->SetBranchAddress("trvtxppftof",&trvtxppftof);
   tree->SetBranchAddress("trvtxxftof",&trvtxxftof);
   tree->SetBranchAddress("trvtxyftof",&trvtxyftof);
   tree->SetBranchAddress("trvtxzftof",&trvtxzftof);
   tree->SetBranchAddress("lengthftof",&lengthftof);
   tree->SetBranchAddress("ntdc",&ntdc);
   tree->SetBranchAddress("tiddc",&tiddc);
   tree->SetBranchAddress("piddc",&piddc);
   tree->SetBranchAddress("diddc",&diddc);
   tree->SetBranchAddress("massdc",&massdc);
   tree->SetBranchAddress("qqdc",&qqdc);
   tree->SetBranchAddress("xdc",&xdc);
   tree->SetBranchAddress("ydc",&ydc);
   tree->SetBranchAddress("zdc",&zdc);
   tree->SetBranchAddress("pxdc",&pxdc);
   tree->SetBranchAddress("pydc",&pydc);
   tree->SetBranchAddress("pzdc",&pzdc);
   tree->SetBranchAddress("ppdc",&ppdc);
   tree->SetBranchAddress("tofdc",&tofdc);
   tree->SetBranchAddress("dcpID",&dcpID);
   tree->SetBranchAddress("trvtxpxdc",&trvtxpxdc);
   tree->SetBranchAddress("trvtxpydc",&trvtxpydc);
   tree->SetBranchAddress("trvtxpzdc",&trvtxpzdc);
   tree->SetBranchAddress("trvtxppdc",&trvtxppdc);
   tree->SetBranchAddress("trvtxxdc",&trvtxxdc);
   tree->SetBranchAddress("trvtxydc",&trvtxydc);
   tree->SetBranchAddress("trvtxzdc",&trvtxzdc);
   tree->SetBranchAddress("lengthdc",&lengthdc);
   tree->SetBranchAddress("targethits",&targethits);
   tree->SetBranchAddress("targetpid",&targetpid);
   tree->SetBranchAddress("targetparentid",&targetparentid);
   tree->SetBranchAddress("targettid",&targettid);
   tree->SetBranchAddress("targetpos",&targetpos);
   tree->SetBranchAddress("targetvtx",&targetvtx);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// tree->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   TFile *test=new TFile(Form("./root/out_e45_%02d_%s.root",beam,generator),"recreate");


   TH2F *phi2d = new TH2F("phi2d","",100,-2.,2., 100,-2.,2.);
   TH2F *id2d = new TH2F("id2d","",32,0.,32., 32,0.,32.);
   TH2F *id2d_all = new TH2F("id2d_all","",32,0.,32., 32,0.,32.);


   TH2F *phi2d_nbar = new TH2F("phi2d_nbar","",100,-2.,2., 100,-2.,2.);
   TH2F *id2d_nbar = new TH2F("id2d_nbar","",32,0.,32., 32,0.,32.);
   TH2F *id2d_all_nbar = new TH2F("id2d_all_nbar","",32,0.,32., 32,0.,32.);


   Long64_t nentries = tree->GetEntries();

   Long64_t nbytes = 0;
   for (Long64_t i=0; i<nentries;i++) {
      nbytes += tree->GetEntry(i);
      Int_t check_pi=0;
      Int_t check_pro=0;
      Int_t check_histo=0;
      Double_t phi_pro,phi_pi;
      Int_t id_pro,id_pi;
      Double_t ppsc_tmp;
      if(ntsc>1.){
      for(Int_t j=0;j<ntsc;j++){
	if( ((masssc[j]>0.11 && masssc[j]<0.2) || (masssc[j]>0.9 && masssc[j]<0.99))
	    && check_pi==0. && check_pro==0. && qqsc[j]!=0){
	  phi_pi=atan2(xsc[j],zsc[j]);
	  id_pi=didsc[j];
	  ppsc_tmp=trvtxppscint[j];
	  check_pi=1.;
	}else if(((masssc[j]>0.11 && masssc[j]<0.2) || (masssc[j]>0.9 && masssc[j]<0.99))
	    && check_pi==1. && check_pro==0. && qqsc[j]!=0 && ppsc_tmp!=trvtxppscint[j]){
	  phi_pro=atan2(xsc[j],zsc[j]);
	  id_pro=didsc[j];
	  check_pro=1.;
	}

	if(check_pi==1 && check_pro==1 && check_histo==0){
	  if(id_pi>=id_pro){
	    id2d_all->Fill(id_pro,id_pi);	    
	  }else{
	    id2d_all->Fill(id_pi,id_pro);
	  }

	  if(gRandom->Rndm()<0.5){
	    phi2d->Fill(phi_pi,phi_pro);
	    id2d->Fill(id_pi,id_pro);
	  }else if{
	    phi2d->Fill(phi_pro,phi_pi);
	    id2d->Fill(id_pro,id_pi);
	  }
	  //	  cout<<"phi_pi:phi_pro-->"<<phi_pi<<", "<<phi_pro<<endl;
	  check_histo=1.;
	}

	
      }
      }


      ///////////////nbar

      Int_t check_pi_nbar=0;
      Int_t check_pro_nbar=0;
      Int_t check_histo_nbar=0;
      Double_t phi_pro_nbar,phi_pi_nbar;
      Int_t id_pro_nbar,id_pi_nbar;
      Double_t ppnbar_tmp;



      if(ntnbar>2 && check_histo==1){
	for(Int_t j=0;j<ntnbar;j++){


	  if( ( (massnbar[j]>0.11 && massnbar[j]<0.2 && qqnbar[j]!=0) || (massnbar[j]>0.9 && massnbar[j]<0.99 && qqnbar[j]==0)) && check_pi_nbar==0  && check_pro_nbar==0 ){


	    phi_pi_nbar=atan2(xnbar[j],znbar[j]);
	    id_pi_nbar=didnbar[j];
	    ppnbar_tmp=trvtxppnbar[j];
	    check_pi_nbar=1.;

	  }else if(((massnbar[j]>0.11 && massnbar[j]<0.2 && qqnbar[j]!=0) || (massnbar[j]>0.9 && massnbar[j]<0.99 && qqnbar[j]==0))
	    && check_pi_nbar==1. && check_pro_nbar==0.  && ppnbar_tmp!=trvtxppnbar[j]){
	  phi_pro_nbar=atan2(xnbar[j],znbar[j]);
	  id_pro_nbar=didnbar[j];
	  check_pro_nbar=1.;
	}


	if(check_pi_nbar==1 && check_pro_nbar==1 && check_histo_nbar==0){
	  if(id_pi_nbar>=id_pro_nbar){
	    id2d_all_nbar->Fill(id_pro_nbar,id_pi_nbar);	    
	  }else{
	    id2d_all_nbar->Fill(id_pi_nbar,id_pro_nbar);
	  }

	  if(gRandom->Rndm()<0.5){
	    phi2d_nbar->Fill(phi_pi_nbar,phi_pro_nbar);
	    id2d_nbar->Fill(id_pi_nbar,id_pro_nbar);
	  }else if{
	    phi2d_nbar->Fill(phi_pro_nbar,phi_pi_nbar);
	    id2d_nbar->Fill(id_pro_nbar,id_pi_nbar);
	  }
	  //	  cout<<"phi_pi:phi_pro-->"<<phi_pi<<", "<<phi_pro<<endl;
	  check_histo_nbar=1.;
	}

	
      }
      }



   }

   ildStyle->SetOptStat(0);
   new TCanvas("c1","",500.,500.);
   phi2d->Draw("colz");
   phi2d->GetXaxis()->SetTitle("#phi of 1st hit on x-z plane");
   phi2d->GetYaxis()->SetTitle("#phi of 2nd hit on x-z plane");
   c1->Print(Form("./ps/phi2d_%02d_%s.eps",beam,generator));

   new TCanvas("c2","",500.,500.);
   id2d->Draw("colz");
   id2d->GetXaxis()->SetTitle("ID of 1st hit for TPC-Hodo");
   id2d->GetYaxis()->SetTitle("ID of 2nd hit for TPC-Hodo");
   c2->Print(Form("./ps/id2d_%02d_%s.eps",beam,generator));

   test->Write();
   test->Close();
   f->Close();
}

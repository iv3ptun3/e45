// -*- C++ -*-

#include "TPCAnaManager.hh"
#include "MatrixReader.hh"

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "ConfMan.hh"
#include "TPCParamMan.hh"
#include "TPCSteppingAction.hh"
#include "FuncName.hh"
#include "ResHypTPC.hh"
#include "RungeKuttaTracker.hh"
#include "switch.h"
#include "track.hh"
#include "VHitInfo.hh"
#include "padHelper.hh"
#include "TString.h"
namespace
{
	Event event;
  const ConfMan& gConf = ConfMan::GetInstance();
  const TPCParamMan& gTPC    = TPCParamMan::GetInstance();
	auto& gTrackBuffer = TPCTrackBuffer::GetInstance();
  //TTree* tree;
  TTree* TPC_g;
  const auto& ResParamInnerLayerHSOn = gTPC.TPCResolutionParams(true, false); //B=1 T, Inner layers
  const auto& ResParamOuterLayerHSOn = gTPC.TPCResolutionParams(true, true); //B=1 T, Outer layers
  const auto& ResParamInnerLayerHSOff = gTPC.TPCResolutionParams(false, false); //B=0, Inner layers
  const auto& ResParamOuterLayerHSOff = gTPC.TPCResolutionParams(false, true); //B=0, Outer layers
 	const auto& DiscardData = gConf.Get<G4bool> ("DiscardData");
// 	const ng& Matrix_2D  = gConf.Get<G4String> ("MTX2D");
}

//_____________________________________________________________________________
TPCAnaManager::TPCAnaManager( void )
{
	const TString Matrix_2D = "param/Matrix/mtx2d1/mtx2d1_e42_Kaon_20210602";
	MatrixReader::mat2d = Matrix_2D;
	MatrixReader::ReadMatrix();
//	if(Matrix_2D)
  TPC_g = new TTree( "TPC_g", "GEANT4 simulation for HypTPC" );
  event.pb = new TVector3;
  TPC_g->Branch( "evnum", &event.evnum, "evnum/I" );
  TPC_g->Branch( "pb", "TVector3", event.pb );
  TPC_g->Branch( "nhPrm", &event.nhPrm, "nhPrm/I" );
  TPC_g->Branch( "data_runnum", &event.data_runnum, "data_runnum/I" );
  TPC_g->Branch( "data_evnum", &event.data_evnum, "data_evnum/I" );
  TPC_g->Branch( "NumberOfTracks",&event.NumberOfTracks, "NumberOfTracks/I" );
  TPC_g->Branch( "PIDOfTrack",event.PIDOfTrack, "PIDOfTrack[1000]/I" );
  TPC_g->Branch( "ParentIDOfTrack",event.ParentIDOfTrack, "ParentIDOfTrack[1000]/I" );
  TPC_g->Branch( "VertexOfTrack_x",event.VertexOfTrack_x, "VertexOfTrack_x[1000]/D" );
  TPC_g->Branch( "VertexOfTrack_y",event.VertexOfTrack_y, "VertexOfTrack_y[1000]/D" );
  TPC_g->Branch( "VertexOfTrack_z",event.VertexOfTrack_z, "VertexOfTrack_z[1000]/D" );
  TPC_g->Branch( "MomentumOfTrack",event.MomentumOfTrack, "MomentumOfTrack[1000]/D" );
  TPC_g->Branch( "MomentumOfTrack_x",event.MomentumOfTrack_x, "MomentumOfTrack_x[1000]/D" );
  TPC_g->Branch( "MomentumOfTrack_y",event.MomentumOfTrack_y, "MomentumOfTrack_y[1000]/D" );
  TPC_g->Branch( "MomentumOfTrack_z",event.MomentumOfTrack_z, "MomentumOfTrack_z[1000]/D" );
  TPC_g->Branch( "trigpat", event.trigpat, "trigpat[32]/I" );
  TPC_g->Branch( "pidPrm", event.pidPrm, "pidPrm[nhPrm]/I" );
  TPC_g->Branch( "xPrm", event.xPrm, "xPrm[nhPrm]/D" );
  TPC_g->Branch( "yPrm", event.yPrm, "yPrm[nhPrm]/D" );
  TPC_g->Branch( "zPrm", event.zPrm, "zPrm[nhPrm]/D" );
  TPC_g->Branch( "pxPrm", event.pxPrm, "pxPrm[nhPrm]/D" );
  TPC_g->Branch( "pyPrm", event.pyPrm, "pyPrm[nhPrm]/D" );
  TPC_g->Branch( "pzPrm", event.pzPrm, "pzPrm[nhPrm]/D" );
  TPC_g->Branch( "ppPrm", event.ppPrm, "ppPrm[nhPrm]/D" );
  TPC_g->Branch( "mPrm", event.mPrm, "mPrm[nhPrm]/D" );
  TPC_g->Branch( "thetaPrm", event.thetaPrm, "thetaPrm[nhPrm]/D" );
  TPC_g->Branch( "phiPrm", event.phiPrm, "phiPrm[nhPrm]/D" );

  TPC_g->Branch("mm_d",&event.mm_d,"mm_d/D");
  // TPC_g->Branch("mm_p",&event.mm_p,"mm_p/D");
  TPC_g->Branch("theta",&event.theta,"theta/D");
  TPC_g->Branch("theta_scat",&event.theta_scat,"theta_scat/D");
  TPC_g->Branch("theta_CM",&event.theta_CM,"theta_CM/D");

  // TPC_g->Branch("mm",&event.mm,"mm/D");

  //generator mode
  TPC_g->Branch("gen",&event.gen,"gen/I");
  TPC_g->Branch("mode",&event.mode,"mode/I");
  TPC_g->Branch("inc",&event.inc,"inc/I");
  // BH2
  TPC_g->Branch( "nhBh2", &event.nhBh2, "nhBh2/I" );
  TPC_g->Branch( "tidBh2", event.tidBh2, "tidBh2[nhBh2]/I" );
  TPC_g->Branch( "pidBh2", event.pidBh2, "pidBh2[nhBh2]/I" );
  TPC_g->Branch( "didBh2", event.didBh2, "didBh2[nhBh2]/I" );
  TPC_g->Branch( "prtBh2", event.prtBh2, "prtBh2[nhBh2]/I" );
  TPC_g->Branch( "qBh2", event.qBh2, "qBh2[nhBh2]/I" );
  TPC_g->Branch( "massBh2", event.massBh2, "massBh2[nhBh2]/D" );
  TPC_g->Branch( "xBh2", event.xBh2, "xBh2[nhBh2]/D" );
  TPC_g->Branch( "yBh2", event.yBh2, "yBh2[nhBh2]/D" );
  TPC_g->Branch( "zBh2", event.zBh2, "zBh2[nhBh2]/D" );
  TPC_g->Branch( "pxBh2", event.pxBh2, "pxBh2[nhBh2]/D" );
  TPC_g->Branch( "pyBh2", event.pyBh2, "pyBh2[nhBh2]/D" );
  TPC_g->Branch( "pzBh2", event.pzBh2, "pzBh2[nhBh2]/D" );
  TPC_g->Branch( "ppBh2", event.ppBh2, "ppBh2[nhBh2]/D" );
  TPC_g->Branch( "deBh2", event.deBh2, "deBh2[nhBh2]/D" );
  TPC_g->Branch( "tBh2", event.tBh2, "tBh2[nhBh2]/D" );
  TPC_g->Branch( "vtpxBh2", event.vtpxBh2, "vtpxBh2[nhBh2]/D" );
  TPC_g->Branch( "vtpyBh2", event.vtpyBh2, "vtpyBh2[nhBh2]/D" );
  TPC_g->Branch( "vtpzBh2", event.vtpzBh2, "vtpzBh2[nhBh2]/D" );
  TPC_g->Branch( "vtppBh2", event.vtppBh2, "vtppBh2[nhBh2]/D" );
  TPC_g->Branch( "vtxBh2", event.vtxBh2, "vtxBh2[nhBh2]/D" );
  TPC_g->Branch( "vtyBh2", event.vtyBh2, "vtyBh2[nhBh2]/D" );
  TPC_g->Branch( "vtzBh2", event.vtzBh2, "vtzBh2[nhBh2]/D" );
  TPC_g->Branch( "lengthBh2", event.lengthBh2, "lengthBh2[nhBh2]/D" );

  ///////shhwang tpc hit step

  // comment out for trigger study
  TPC_g->Branch("nhittpc",&event.nhittpc,"nhittpc/I");
  TPC_g->Branch("ntrk",event.ntrk,"ntrk[nhittpc]/I");
  TPC_g->Branch("ititpc",event.ititpc,"ititpc[nhittpc]/I");
  TPC_g->Branch("idtpc",event.idtpc,"idtpc[nhittpc]/I");
  TPC_g->Branch("xtpc",event.xtpc,"xtpc[nhittpc]/D");//after smeared by resolution
  TPC_g->Branch("ytpc",event.ytpc,"ytpc[nhittpc]/D");//after smeared by resolution
  TPC_g->Branch("ztpc",event.ztpc,"ztpc[nhittpc]/D");//after smeared by resolution
  TPC_g->Branch("x0tpc",event.x0tpc,"x0tpc[nhittpc]/D");
  TPC_g->Branch("y0tpc",event.y0tpc,"y0tpc[nhittpc]/D");
  TPC_g->Branch("z0tpc",event.z0tpc,"z0tpc[nhittpc]/D");
  TPC_g->Branch("resoX",event.resoX,"resoX[nhittpc]/D");
  TPC_g->Branch("resxtpc",event.resxtpc,"resxtpc[nhittpc]/D");
  TPC_g->Branch("resytpc",event.resytpc,"resytpc[nhittpc]/D");
  TPC_g->Branch("resztpc",event.resztpc,"resztpc[nhittpc]/D");
  TPC_g->Branch("pxtpc",event.pxtpc,"pxtpc[nhittpc]/D");
  TPC_g->Branch("pytpc",event.pytpc,"pytpc[nhittpc]/D");
  TPC_g->Branch("pztpc",event.pztpc,"pztpc[nhittpc]/D");
  TPC_g->Branch("pptpc",event.pptpc,"pptpc[nhittpc]/D");   // total mometum TPC_g->Branch("masstpc",event.masstpc,"masstpc[nhittpc]/D");   // mass TPC
  TPC_g->Branch("timetpc",event.timetpc,"timetpc[nhittpc]/D");
  TPC_g->Branch("betatpc",event.betatpc,"betatpc[nhittpc]/D");
  TPC_g->Branch("edeptpc",event.edeptpc,"edeptpc[nhittpc]/D");
  TPC_g->Branch("dedxtpc",event.dedxtpc,"dedxtpc[nhittpc]/D");
  TPC_g->Branch("slengthtpc",event.slengthtpc,"slengthtpc[nhittpc]/D");
  TPC_g->Branch("tlengthtpc",event.tlengthtpc,"tlengthtpc[nhittpc]/D");
  TPC_g->Branch("iPadtpc",event.iPadtpc,"iPadtpc[nhittpc]/I");
  TPC_g->Branch("laytpc",event.laytpc,"laytpc[nhittpc]/I");
  TPC_g->Branch("rowtpc",event.rowtpc,"rowtpc[nhittpc]/I");
  TPC_g->Branch("parentID",event.parentID,"parentID[nhittpc]/I");
  TPC_g->Branch("xtpc_pad",event.xtpc_pad,"xtpc_pad[nhittpc]/D");//pad center position
  TPC_g->Branch("ytpc_pad",event.ytpc_pad,"ytpc_pad[nhittpc]/D");//pad center position (dummy = ytpc)
  TPC_g->Branch("ztpc_pad",event.ztpc_pad,"ztpc_pad[nhittpc]/D");//pad center position
  TPC_g->Branch("dxtpc_pad",event.dxtpc_pad,"dxtpc_pad[nhittpc]/D");//x0tpc - xtpc_pad
  TPC_g->Branch("dytpc_pad",event.dytpc_pad,"dytpc_pad[nhittpc]/D");//y0tpc - ytpc_pad (dummy = 0)
  TPC_g->Branch("dztpc_pad",event.dztpc_pad,"dztpc_pad[nhittpc]/D");//z0tpc - ztpc_pad


	//Spectrometer for Databased simulation//
  TPC_g->Branch( "ntK18",      &event.ntK18);
  TPC_g->Branch( "xvpHS",    &event.xvpHS);
  TPC_g->Branch( "yvpHS",    &event.yvpHS);
  TPC_g->Branch( "zvpHS",    &event.zvpHS);
  TPC_g->Branch( "xtgtHS",    &event.xtgtHS);
  TPC_g->Branch( "ytgtHS",    &event.ytgtHS);
  TPC_g->Branch( "ztgtHS",    &event.ztgtHS);
  TPC_g->Branch( "p_3rd",       &event.p_3rd);
  TPC_g->Branch( "xoutK18",    &event.xoutK18);
  TPC_g->Branch( "youtK18",    &event.youtK18);
  TPC_g->Branch( "uoutK18",    &event.uoutK18);
  TPC_g->Branch( "voutK18",    &event.voutK18);
  TPC_g->Branch( "layerK18",    &event.layerK18);
  TPC_g->Branch( "wireK18",    &event.wireK18);
  TPC_g->Branch( "localhitposK18",    &event.localhitposK18);
	TPC_g->Branch( "ntKurama",      &event.ntKurama);
  TPC_g->Branch( "xout",    &event.xout);
  TPC_g->Branch( "yout",    &event.yout);
  TPC_g->Branch( "zout",    &event.zout);
  TPC_g->Branch( "pxout",    &event.pxout);
  TPC_g->Branch( "pyout",    &event.pyout);
  TPC_g->Branch( "pzout",    &event.pzout);
  TPC_g->Branch( "xvpKurama",     &event.xvpKurama);
  TPC_g->Branch( "yvpKurama",     &event.yvpKurama);
  TPC_g->Branch( "zvpKurama",     &event.zvpKurama);
  TPC_g->Branch( "xtgtKurama",     &event.xtgtKurama);
  TPC_g->Branch( "ytgtKurama",     &event.ytgtKurama);
  TPC_g->Branch( "layer",    &event.layer);
  TPC_g->Branch( "wire",    &event.wire);
  TPC_g->Branch( "localhitpos",    &event.localhitpos);


  //// Study on multiplicity
  // TPC_g->Branch("nthlay",event.nthlay,"nthlay[nhittpc]/I");
  // TPC_g->Branch("nthpad",event.nthpad,"nthpad[nhittpc]/I");
  // TPC_g->Branch("laypad",event.laypad,"laytpadpc[nhittpc][40][250]/I");


  //shhwang ntrtpc --> number of trak in tpc
  // TPC_g->Branch("ntrtpc",&event.ntrtpc,"ntrtpc/I");
  // TPC_g->Branch("trpmtpc",event.trpmtpc,"trpmtpc[ntrtpc]/D");
  // TPC_g->Branch("trqqtpc",event.trqqtpc,"trqqtpc[ntrtpc]/I");
  // TPC_g->Branch("trpidtpc",event.trpidtpc,"trpidtpc[ntrtpc]/I");
  // TPC_g->Branch("trparentidtpc",event.trparentidtpc,"trparentidtpc[ntrtpc]/I");
  // //TPC_g->Branch("trparentid_pid_tpc",event.trparentid_pid_tpc,"trparentid_pid_tpc[ntrtpc]/I");

  // TPC_g->Branch("trpxtpc",event.trpxtpc,"trpxtpc[ntrtpc]/D");
  // TPC_g->Branch("trpytpc",event.trpytpc,"trpytpc[ntrtpc]/D");
  // TPC_g->Branch("trpztpc",event.trpztpc,"trpztpc[ntrtpc]/D");
  // TPC_g->Branch("trpptpc",event.trpptpc,"trpptpc[ntrtpc]/D");
  // TPC_g->Branch("trpttpc",event.trpttpc,"trpttpc[ntrtpc]/D");

  // TPC_g->Branch("trpxtpcfit",event.trpxtpcfit,"trpxtpcfit[ntrtpc]/D");
  // TPC_g->Branch("trpytpcfit",event.trpytpcfit,"trpytpcfit[ntrtpc]/D");
  // TPC_g->Branch("trpztpcfit",event.trpztpcfit,"trpztpcfit[ntrtpc]/D");
  // TPC_g->Branch("trpptpcfit",event.trpptpcfit,"trpptpcfit[ntrtpc]/D");
  // TPC_g->Branch("trpttpcfit",event.trpttpcfit,"trpttpcfit[ntrtpc]/D");

  // TPC_g->Branch("vtpxtpc",event.vtpxtpc,"vtpxtpc[ntrtpc]/D");
  // TPC_g->Branch("vtpytpc",event.vtpytpc,"vtpytpc[ntrtpc]/D");
  // TPC_g->Branch("vtpztpc",event.vtpztpc,"vtpztpc[ntrtpc]/D");
  // TPC_g->Branch("vtpptpc",event.vtpptpc,"vtpptpc[ntrtpc]/D");

  // TPC_g->Branch("vtxtpc",event.vtxtpc,"vtxtpc[ntrtpc]/D");
  // TPC_g->Branch("vtytpc",event.vtytpc,"vtytpc[ntrtpc]/D");
  // TPC_g->Branch("vtztpc",event.vtztpc,"vtztpc[ntrtpc]/D");

  // TPC_g->Branch("vtxtpcfit",event.vtxtpcfit,"vtxtpcfit[ntrtpc]/D");
  // TPC_g->Branch("vtytpcfit",event.vtytpcfit,"vtytpcfit[ntrtpc]/D");
  // TPC_g->Branch("vtztpcfit",event.vtztpcfit,"vtztpcfit[ntrtpc]/D");

  // TPC_g->Branch("trdetpc",event.trdetpc,"trdetpc[ntrtpc]/D");
  // TPC_g->Branch("trlentpc",event.trlentpc,"trlentpc[ntrtpc]/D");
  // TPC_g->Branch("trdedxtpc",event.trdedxtpc,"trdedxtpc[ntrtpc]/D");
  // TPC_g->Branch("trdedxtrtpc",event.trdedxtrtpc,"trdedxtrtpc[ntrtpc]/D");
  // TPC_g->Branch("trlaytpc",event.trlaytpc,"trlaytpc[ntrtpc]/I");
  // TPC_g->Branch("cir_r",event.cir_r,"cir_r[ntrtpc]/D");
  // TPC_g->Branch("cir_x",event.cir_x,"cir_x[ntrtpc]/D");
  // TPC_g->Branch("cir_z",event.cir_z,"cir_z[ntrtpc]/D");
  // TPC_g->Branch("cir_fit",event.cir_fit,"cir_fit[ntrtpc]/D");
  // TPC_g->Branch("vtx_flag",event.vtx_flag,"vtx_flag[ntrtpc]/I");
  // TPC_g->Branch("a_fory",event.a_fory,"a_fory[ntrtpc]/D");
  // TPC_g->Branch("b_fory",event.b_fory,"b_fory[ntrtpc]/D");

  // TARGET
  TPC_g->Branch( "nhTgt", &event.nhTgt, "nhTgt/I" );
  TPC_g->Branch( "tidTgt", event.tidTgt, "tidTgt[nhTgt]/I" );
  TPC_g->Branch( "pidTgt", event.pidTgt, "pidTgt[nhTgt]/I" );
  TPC_g->Branch( "prtTgt", event.prtTgt, "prtTgt[nhTgt]/I" );
  TPC_g->Branch( "xTgt", event.xTgt, "xTgt[nhTgt]/D" );
  TPC_g->Branch( "yTgt", event.yTgt, "yTgt[nhTgt]/D" );
  TPC_g->Branch( "zTgt", event.zTgt, "zTgt[nhTgt]/D" );
  TPC_g->Branch( "uTgt", event.uTgt, "uTgt[nhTgt]/D" );
  TPC_g->Branch( "vTgt", event.vTgt, "vTgt[nhTgt]/D" );
  TPC_g->Branch( "vtxTgt", event.vtxTgt, "vtxTgt[nhTgt]/D" );
  TPC_g->Branch( "vtyTgt", event.vtyTgt, "vtyTgt[nhTgt]/D" );
  TPC_g->Branch( "vtzTgt", event.vtzTgt, "vtzTgt[nhTgt]/D" );
  // HTOF
  TPC_g->Branch( "nhHtof", &event.nhHtof, "nhHtof/I" );
  TPC_g->Branch( "tidHtof", event.tidHtof, "tidHtof[nhHtof]/I" );
  TPC_g->Branch( "pidHtof", event.pidHtof, "pidHtof[nhHtof]/I" );
  TPC_g->Branch( "didHtof", event.didHtof, "didHtof[nhHtof]/I" );
  TPC_g->Branch( "prtHtof", event.prtHtof, "prtHtof[nhHtof]/I" );
  TPC_g->Branch( "qHtof", event.qHtof, "qHtof[nhHtof]/I" );
  TPC_g->Branch( "massHtof", event.massHtof, "massHtof[nhHtof]/D" );
  TPC_g->Branch( "xHtof", event.xHtof, "xHtof[nhHtof]/D" );
  TPC_g->Branch( "yHtof", event.yHtof, "yHtof[nhHtof]/D" );
  TPC_g->Branch( "zHtof", event.zHtof, "zHtof[nhHtof]/D" );
  TPC_g->Branch( "pxHtof", event.pxHtof, "pxHtof[nhHtof]/D" );
  TPC_g->Branch( "pyHtof", event.pyHtof, "pyHtof[nhHtof]/D" );
  TPC_g->Branch( "pzHtof", event.pzHtof, "pzHtof[nhHtof]/D" );
  TPC_g->Branch( "ppHtof", event.ppHtof, "ppHtof[nhHtof]/D" );
  TPC_g->Branch( "deHtof", event.deHtof, "deHtof[nhHtof]/D" );
  TPC_g->Branch( "tHtof", event.tHtof, "tHtof[nhHtof]/D" );
  TPC_g->Branch( "vtpxHtof", event.vtpxHtof, "vtpxHtof[nhHtof]/D" );
  TPC_g->Branch( "vtpyHtof", event.vtpyHtof, "vtpyHtof[nhHtof]/D" );
  TPC_g->Branch( "vtpzHtof", event.vtpzHtof, "vtpzHtof[nhHtof]/D" );
  TPC_g->Branch( "vtppHtof", event.vtppHtof, "vtppHtof[nhHtof]/D" );
  TPC_g->Branch( "vtxHtof", event.vtxHtof, "vtxHtof[nhHtof]/D" );
  TPC_g->Branch( "vtyHtof", event.vtyHtof, "vtyHtof[nhHtof]/D" );
  TPC_g->Branch( "vtzHtof", event.vtzHtof, "vtzHtof[nhHtof]/D" );
  TPC_g->Branch( "lengthHtof", event.lengthHtof, "lengthHtof[nhHtof]/D" );
  // SDC
  TPC_g->Branch( "nhSdc", &event.nhSdc, "nhSdc/I" );
  TPC_g->Branch( "tidSdc", event.tidSdc, "tidSdc[nhSdc]/I" );
  TPC_g->Branch( "pidSdc", event.pidSdc, "pidSdc[nhSdc]/I" );
  TPC_g->Branch( "didSdc", event.didSdc, "didSdc[nhSdc]/I" );
  TPC_g->Branch( "prtSdc", event.prtSdc, "prtSdc[nhSdc]/I" );
  TPC_g->Branch( "qSdc", event.qSdc, "qSdc[nhSdc]/I" );
  TPC_g->Branch( "massSdc", event.massSdc, "massSdc[nhSdc]/D" );
  TPC_g->Branch( "xSdc", event.xSdc, "xSdc[nhSdc]/D" );
  TPC_g->Branch( "ySdc", event.ySdc, "ySdc[nhSdc]/D" );
  TPC_g->Branch( "zSdc", event.zSdc, "zSdc[nhSdc]/D" );
  TPC_g->Branch( "pxSdc", event.pxSdc, "pxSdc[nhSdc]/D" );
  TPC_g->Branch( "pySdc", event.pySdc, "pySdc[nhSdc]/D" );
  TPC_g->Branch( "pzSdc", event.pzSdc, "pzSdc[nhSdc]/D" );
  TPC_g->Branch( "ppSdc", event.ppSdc, "ppSdc[nhSdc]/D" );
  TPC_g->Branch( "deSdc", event.deSdc, "deSdc[nhSdc]/D" );
  TPC_g->Branch( "tSdc", event.tSdc, "tSdc[nhSdc]/D" );
  TPC_g->Branch( "vtpxSdc", event.vtpxSdc, "vtpxSdc[nhSdc]/D" );
  TPC_g->Branch( "vtpySdc", event.vtpySdc, "vtpySdc[nhSdc]/D" );
  TPC_g->Branch( "vtpzSdc", event.vtpzSdc, "vtpzSdc[nhSdc]/D" );
  TPC_g->Branch( "vtppSdc", event.vtppSdc, "vtppSdc[nhSdc]/D" );
  TPC_g->Branch( "vtxSdc", event.vtxSdc, "vtxSdc[nhSdc]/D" );
  TPC_g->Branch( "vtySdc", event.vtySdc, "vtySdc[nhSdc]/D" );
  TPC_g->Branch( "vtzSdc", event.vtzSdc, "vtzSdc[nhSdc]/D" );
  TPC_g->Branch( "lengthSdc", event.lengthSdc, "lengthSdc[nhSdc]/D" );
  // SCH
  TPC_g->Branch( "nhSch", &event.nhSch, "nhSch/I" );
  TPC_g->Branch( "tidSch", event.tidSch, "tidSch[nhSch]/I" );
  TPC_g->Branch( "pidSch", event.pidSch, "pidSch[nhSch]/I" );
  TPC_g->Branch( "didSch", event.didSch, "didSch[nhSch]/I" );
  TPC_g->Branch( "prtSch", event.prtSch, "prtSch[nhSch]/I" );
  TPC_g->Branch( "qSch", event.qSch, "qSch[nhSch]/I" );
  TPC_g->Branch( "massSch", event.massSch, "massSch[nhSch]/D" );
  TPC_g->Branch( "xSch", event.xSch, "xSch[nhSch]/D" );
  TPC_g->Branch( "ySch", event.ySch, "ySch[nhSch]/D" );
  TPC_g->Branch( "zSch", event.zSch, "zSch[nhSch]/D" );
  TPC_g->Branch( "pxSch", event.pxSch, "pxSch[nhSch]/D" );
  TPC_g->Branch( "pySch", event.pySch, "pySch[nhSch]/D" );
  TPC_g->Branch( "pzSch", event.pzSch, "pzSch[nhSch]/D" );
  TPC_g->Branch( "ppSch", event.ppSch, "ppSch[nhSch]/D" );
  TPC_g->Branch( "deSch", event.deSch, "deSch[nhSch]/D" );
  TPC_g->Branch( "tSch", event.tSch, "tSch[nhSch]/D" );
  TPC_g->Branch( "vtpxSch", event.vtpxSch, "vtpxSch[nhSch]/D" );
  TPC_g->Branch( "vtpySch", event.vtpySch, "vtpySch[nhSch]/D" );
  TPC_g->Branch( "vtpzSch", event.vtpzSch, "vtpzSch[nhSch]/D" );
  TPC_g->Branch( "vtppSch", event.vtppSch, "vtppSch[nhSch]/D" );
  TPC_g->Branch( "vtxSch", event.vtxSch, "vtxSch[nhSch]/D" );
  TPC_g->Branch( "vtySch", event.vtySch, "vtySch[nhSch]/D" );
  TPC_g->Branch( "vtzSch", event.vtzSch, "vtzSch[nhSch]/D" );
  TPC_g->Branch( "lengthSch", event.lengthSch, "lengthSch[nhSch]/D" );
  // FTOF
  TPC_g->Branch( "nhFtof", &event.nhFtof, "nhFtof/I" );
  TPC_g->Branch( "tidFtof", event.tidFtof, "tidFtof[nhFtof]/I" );
  TPC_g->Branch( "pidFtof", event.pidFtof, "pidFtof[nhFtof]/I" );
  TPC_g->Branch( "didFtof", event.didFtof, "didFtof[nhFtof]/I" );
  TPC_g->Branch( "prtFtof", event.prtFtof, "prtFtof[nhFtof]/I" );
  TPC_g->Branch( "qFtof", event.qFtof, "qFtof[nhFtof]/I" );
  TPC_g->Branch( "massFtof", event.massFtof, "massFtof[nhFtof]/D" );
  TPC_g->Branch( "xFtof", event.xFtof, "xFtof[nhFtof]/D" );
  TPC_g->Branch( "yFtof", event.yFtof, "yFtof[nhFtof]/D" );
  TPC_g->Branch( "zFtof", event.zFtof, "zFtof[nhFtof]/D" );
  TPC_g->Branch( "pxFtof", event.pxFtof, "pxFtof[nhFtof]/D" );
  TPC_g->Branch( "pyFtof", event.pyFtof, "pyFtof[nhFtof]/D" );
  TPC_g->Branch( "pzFtof", event.pzFtof, "pzFtof[nhFtof]/D" );
  TPC_g->Branch( "ppFtof", event.ppFtof, "ppFtof[nhFtof]/D" );
  TPC_g->Branch( "deFtof", event.deFtof, "deFtof[nhFtof]/D" );
  TPC_g->Branch( "tFtof", event.tFtof, "tFtof[nhFtof]/D" );
  TPC_g->Branch( "vtpxFtof", event.vtpxFtof, "vtpxFtof[nhFtof]/D" );
  TPC_g->Branch( "vtpyFtof", event.vtpyFtof, "vtpyFtof[nhFtof]/D" );
  TPC_g->Branch( "vtpzFtof", event.vtpzFtof, "vtpzFtof[nhFtof]/D" );
  TPC_g->Branch( "vtppFtof", event.vtppFtof, "vtppFtof[nhFtof]/D" );
  TPC_g->Branch( "vtxFtof", event.vtxFtof, "vtxFtof[nhFtof]/D" );
  TPC_g->Branch( "vtyFtof", event.vtyFtof, "vtyFtof[nhFtof]/D" );
  TPC_g->Branch( "vtzFtof", event.vtzFtof, "vtzFtof[nhFtof]/D" );
  TPC_g->Branch( "lengthFtof", event.lengthFtof, "lengthFtof[nhFtof]/D" );
  // LAC
  TPC_g->Branch( "nhLac", &event.nhLac, "nhLac/I" );
  TPC_g->Branch( "tidLac", event.tidLac, "tidLac[nhLac]/I" );
  TPC_g->Branch( "pidLac", event.pidLac, "pidLac[nhLac]/I" );
  TPC_g->Branch( "didLac", event.didLac, "didLac[nhLac]/I" );
  TPC_g->Branch( "prtLac", event.prtLac, "prtLac[nhLac]/I" );
  TPC_g->Branch( "qLac", event.qLac, "qLac[nhLac]/I" );
  TPC_g->Branch( "massLac", event.massLac, "massLac[nhLac]/D" );
  TPC_g->Branch( "xLac", event.xLac, "xLac[nhLac]/D" );
  TPC_g->Branch( "yLac", event.yLac, "yLac[nhLac]/D" );
  TPC_g->Branch( "zLac", event.zLac, "zLac[nhLac]/D" );
  TPC_g->Branch( "pxLac", event.pxLac, "pxLac[nhLac]/D" );
  TPC_g->Branch( "pyLac", event.pyLac, "pyLac[nhLac]/D" );
  TPC_g->Branch( "pzLac", event.pzLac, "pzLac[nhLac]/D" );
  TPC_g->Branch( "ppLac", event.ppLac, "ppLac[nhLac]/D" );
  TPC_g->Branch( "deLac", event.deLac, "deLac[nhLac]/D" );
  TPC_g->Branch( "tLac", event.tLac, "tLac[nhLac]/D" );
  TPC_g->Branch( "vtpxLac", event.vtpxLac, "vtpxLac[nhLac]/D" );
  TPC_g->Branch( "vtpyLac", event.vtpyLac, "vtpyLac[nhLac]/D" );
  TPC_g->Branch( "vtpzLac", event.vtpzLac, "vtpzLac[nhLac]/D" );
  TPC_g->Branch( "vtppLac", event.vtppLac, "vtppLac[nhLac]/D" );
  TPC_g->Branch( "vtxLac", event.vtxLac, "vtxLac[nhLac]/D" );
  TPC_g->Branch( "vtyLac", event.vtyLac, "vtyLac[nhLac]/D" );
  TPC_g->Branch( "vtzLac", event.vtzLac, "vtzLac[nhLac]/D" );
  TPC_g->Branch( "lengthLac", event.lengthLac, "lengthLac[nhLac]/D" );
  // WC
  TPC_g->Branch( "nhWc", &event.nhWc, "nhWc/I" );
  TPC_g->Branch( "tidWc", event.tidWc, "tidWc[nhWc]/I" );
  TPC_g->Branch( "pidWc", event.pidWc, "pidWc[nhWc]/I" );
  TPC_g->Branch( "didWc", event.didWc, "didWc[nhWc]/I" );
  TPC_g->Branch( "prtWc", event.prtWc, "prtWc[nhWc]/I" );
  TPC_g->Branch( "qWc", event.qWc, "qWc[nhWc]/I" );
  TPC_g->Branch( "massWc", event.massWc, "massWc[nhWc]/D" );
  TPC_g->Branch( "xWc", event.xWc, "xWc[nhWc]/D" );
  TPC_g->Branch( "yWc", event.yWc, "yWc[nhWc]/D" );
  TPC_g->Branch( "zWc", event.zWc, "zWc[nhWc]/D" );
  TPC_g->Branch( "pxWc", event.pxWc, "pxWc[nhWc]/D" );
  TPC_g->Branch( "pyWc", event.pyWc, "pyWc[nhWc]/D" );
  TPC_g->Branch( "pzWc", event.pzWc, "pzWc[nhWc]/D" );
  TPC_g->Branch( "ppWc", event.ppWc, "ppWc[nhWc]/D" );
  TPC_g->Branch( "deWc", event.deWc, "deWc[nhWc]/D" );
  TPC_g->Branch( "tWc", event.tWc, "tWc[nhWc]/D" );
  TPC_g->Branch( "vtpxWc", event.vtpxWc, "vtpxWc[nhWc]/D" );
  TPC_g->Branch( "vtpyWc", event.vtpyWc, "vtpyWc[nhWc]/D" );
  TPC_g->Branch( "vtpzWc", event.vtpzWc, "vtpzWc[nhWc]/D" );
  TPC_g->Branch( "vtppWc", event.vtppWc, "vtppWc[nhWc]/D" );
  TPC_g->Branch( "vtxWc", event.vtxWc, "vtxWc[nhWc]/D" );
  TPC_g->Branch( "vtyWc", event.vtyWc, "vtyWc[nhWc]/D" );
  TPC_g->Branch( "vtzWc", event.vtzWc, "vtzWc[nhWc]/D" );
  TPC_g->Branch( "lengthWc", event.lengthWc, "lengthWc[nhWc]/D" );
	//BVH 
  TPC_g->Branch( "nhBvh", &event.nhBvh, "nhBvh/I" );
  TPC_g->Branch( "tidBvh", event.tidBvh, "tidBvh[nhBvh]/I" );
  TPC_g->Branch( "pidBvh", event.pidBvh, "pidBvh[nhBvh]/I" );
  TPC_g->Branch( "didBvh", event.didBvh, "didBvh[nhBvh]/I" );
  TPC_g->Branch( "prtBvh", event.prtBvh, "prtBvh[nhBvh]/I" );
  TPC_g->Branch( "qBvh", event.qBvh, "qBvh[nhBvh]/I" );
  TPC_g->Branch( "massBvh", event.massBvh, "massBvh[nhBvh]/D" );
  TPC_g->Branch( "xBvh", event.xBvh, "xBvh[nhBvh]/D" );
  TPC_g->Branch( "yBvh", event.yBvh, "yBvh[nhBvh]/D" );
  TPC_g->Branch( "zBvh", event.zBvh, "zBvh[nhBvh]/D" );
  TPC_g->Branch( "pxBvh", event.pxBvh, "pxBvh[nhBvh]/D" );
  TPC_g->Branch( "pyBvh", event.pyBvh, "pyBvh[nhBvh]/D" );
  TPC_g->Branch( "pzBvh", event.pzBvh, "pzBvh[nhBvh]/D" );
  TPC_g->Branch( "ppBvh", event.ppBvh, "ppBvh[nhBvh]/D" );
  TPC_g->Branch( "deBvh", event.deBvh, "deBvh[nhBvh]/D" );
  TPC_g->Branch( "tBvh", event.tBvh, "tBvh[nhBvh]/D" );
  TPC_g->Branch( "vtpxBvh", event.vtpxBvh, "vtpxBvh[nhBvh]/D" );
  TPC_g->Branch( "vtpyBvh", event.vtpyBvh, "vtpyBvh[nhBvh]/D" );
  TPC_g->Branch( "vtpzBvh", event.vtpzBvh, "vtpzBvh[nhBvh]/D" );
  TPC_g->Branch( "vtppBvh", event.vtppBvh, "vtppBvh[nhBvh]/D" );
  TPC_g->Branch( "vtxBvh", event.vtxBvh, "vtxBvh[nhBvh]/D" );
  TPC_g->Branch( "vtyBvh", event.vtyBvh, "vtyBvh[nhBvh]/D" );
  TPC_g->Branch( "vtzBvh", event.vtzBvh, "vtzBvh[nhBvh]/D" );
  TPC_g->Branch( "lengthBvh", event.lengthBvh, "lengthBvh[nhBvh]/D" );
	// VP
  TPC_g->Branch( "nhVp", &event.nhVp, "nhVp/I" );
  TPC_g->Branch( "tidVp", event.tidVp, "tidVp[nhVp]/I" );
  TPC_g->Branch( "pidVp", event.pidVp, "pidVp[nhVp]/I" );
  TPC_g->Branch( "didVp", event.didVp, "didVp[nhVp]/I" );
  TPC_g->Branch( "prtVp", event.prtVp, "prtVp[nhVp]/I" );
  TPC_g->Branch( "qVp", event.qVp, "qVp[nhVp]/I" );
  TPC_g->Branch( "massVp", event.massVp, "massVp[nhVp]/D" );
  TPC_g->Branch( "xVp", event.xVp, "xVp[nhVp]/D" );
  TPC_g->Branch( "yVp", event.yVp, "yVp[nhVp]/D" );
  TPC_g->Branch( "zVp", event.zVp, "zVp[nhVp]/D" );
  TPC_g->Branch( "pxVp", event.pxVp, "pxVp[nhVp]/D" );
  TPC_g->Branch( "pyVp", event.pyVp, "pyVp[nhVp]/D" );
  TPC_g->Branch( "pzVp", event.pzVp, "pzVp[nhVp]/D" );
  TPC_g->Branch( "ppVp", event.ppVp, "ppVp[nhVp]/D" );
  TPC_g->Branch( "deVp", event.deVp, "deVp[nhVp]/D" );
  TPC_g->Branch( "tVp", event.tVp, "tVp[nhVp]/D" );
  TPC_g->Branch( "vtpxVp", event.vtpxVp, "vtpxVp[nhVp]/D" );
  TPC_g->Branch( "vtpyVp", event.vtpyVp, "vtpyVp[nhVp]/D" );
  TPC_g->Branch( "vtpzVp", event.vtpzVp, "vtpzVp[nhVp]/D" );
  TPC_g->Branch( "vtppVp", event.vtppVp, "vtppVp[nhVp]/D" );
  TPC_g->Branch( "vtxVp", event.vtxVp, "vtxVp[nhVp]/D" );
  TPC_g->Branch( "vtyVp", event.vtyVp, "vtyVp[nhVp]/D" );
  TPC_g->Branch( "vtzVp", event.vtzVp, "vtzVp[nhVp]/D" );
  TPC_g->Branch( "lengthVp", event.lengthVp, "lengthVp[nhVp]/D" );

	TPC_g->Branch( "SpinXi", &event.SpinXi,"SpinXi/D");
	TPC_g->Branch( "SpinXi_x", &event.SpinXi_x,"SpinXi_x/D");
	TPC_g->Branch( "SpinXi_y", &event.SpinXi_y,"SpinXi_y/D");
	TPC_g->Branch( "SpinXi_z", &event.SpinXi_z,"SpinXi_x/D");
	TPC_g->Branch( "MomXi", &event.MomXi,"MomXi/D");
	TPC_g->Branch( "MomXi_x", &event.MomXi_x,"MomXi_x/D");
	TPC_g->Branch( "MomXi_y", &event.MomXi_y,"MomXi_y/D");
	TPC_g->Branch( "MomXi_z", &event.MomXi_z,"MomXi_z/D");
	TPC_g->Branch( "SpinLd", &event.SpinLd,"SpinLd/D");
	TPC_g->Branch( "SpinLd_x", &event.SpinLd_x,"SpinLd_x/D");
	TPC_g->Branch( "SpinLd_y", &event.SpinLd_y,"SpinLd_y/D");
	TPC_g->Branch( "SpinLd_z", &event.SpinLd_z,"SpinLd_x/D");
	TPC_g->Branch( "MomLd", &event.MomLd,"MomLd/D");
	TPC_g->Branch( "MomLd_x", &event.MomLd_x,"MomLd_x/D");
	TPC_g->Branch( "MomLd_y", &event.MomLd_y,"MomLd_y/D");
	TPC_g->Branch( "MomLd_z", &event.MomLd_z,"MomLd_z/D");
	TPC_g->Branch( "CM_E", &event.CM_E,"CM_E/D");
	TPC_g->Branch( "CM_x", &event.CM_x,"CM_x/D");
	TPC_g->Branch( "CM_y", &event.CM_y,"CM_y/D");
	TPC_g->Branch( "CM_z", &event.CM_z,"CM_z/D");
	TPC_g->Branch( "ThXi_CM", &event.ThXi_CM,"ThXi_CM/D");
	TPC_g->Branch( "MomP", &event.MomP,"MomP/D");
	TPC_g->Branch( "MomP_x", &event.MomP_x,"MomP_x/D");
	TPC_g->Branch( "MomP_y", &event.MomP_y,"MomP_y/D");
	TPC_g->Branch( "MomP_z", &event.MomP_z,"MomP_z/D");
	
	TPC_g->Branch( "ThLd_CM", &event.ThLd_CM,"ThLd_CM/D");
}

//_____________________________________________________________________________
TPCAnaManager::~TPCAnaManager( void )
{
}

//_____________________________________________________________________________
void
TPCAnaManager::BeginOfRunAction( G4int /* runnum */ )
{
  event.evnum = 0;

  G4double target_pos_z=-143.;
  truncated_mean_cut = gConf.Get<G4double>("TruncatedMeanCut");
  m_experiment = gConf.Get<G4int>("Experiment");
  //out side less 100 mm. 10+5*x < 100 mm is pad_in_num
  pad_length_in = gConf.Get<G4double>("PadLengthIn");
  pad_length_out = gConf.Get<G4double>("PadLengthOut");
  pad_gap = gConf.Get<G4double>("PadGap");

  ////pad configure
  m_pad_config = gConf.Get<G4int>("PadConfigure");
  pad_in_num = gConf.Get<G4int>("PadNumIn");
  pad_out_num = gConf.Get<G4int>("PadNumOut");
  pad_in_width = gConf.Get<G4double>("PadWidthOut");
  pad_out_width = gConf.Get<G4double>("PadWidthOut");

  m_on_off_helm = gConf.Get<G4int>("ShsFieldMap");

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

  if( m_pad_config ==1 ){
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


  }else if( m_pad_config ==2 ){
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
  for( auto& p : hmap2d ){
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
    key = Form( "Fermi%d", i );
    hmap[key] = new TH1D( key, key, 500, -1.0*CLHEP::GeV, 1.0*CLHEP::GeV );
    hmap[key]->GetXaxis()->SetTitle( "[MeV/c]" );
  }
	key ="BeamGenThetaP";
	hmap2d[key] = new TH2D(key,key,300,0,30,320,0.4,2.);
	key ="BeamGenCosTP";
	hmap2d[key] = new TH2D(key,key,100,0.85,1,80,0.4,2.);
	key ="BeamGenCosTPhi";
	hmap2d[key] = new TH2D(key,key,100,0.85,1,100,-3.15,3.15);
	key ="BeamGenPhiP";
	hmap2d[key] = new TH2D(key,key,100,-3.15,3.15,80,0.4,2.);
	key ="BeamGenXThetaP";
	hmap2d[key] = new TH2D(key,key,100,-30,30,80,0.4,2.);
	key ="BeamGenYThetaP";
	hmap2d[key] = new TH2D(key,key,100,-30,30,80,0.4,2.);


	key ="Sdc1Hitpat";
	hmap2d[key] = new TH2D(key,key,1000,-1000,1000,1000,-500,500);
	key ="Sdc2Hitpat";
	hmap2d[key] = new TH2D(key,key,1000,-1000,1000,1000,-500,500);
	key ="Sdc3Hitpat";
	hmap2d[key] = new TH2D(key,key,1000,-1000,1000,1000,-500,500);
	key ="Sdc4Hitpat";
	hmap2d[key] = new TH2D(key,key,1000,-1000,1000,1000,-500,500);


	key ="VP1Hitpat";
	hmap2d[key] = new TH2D(key,key,1000,-1000,1000,1000,-500,500);
	key ="VP2Hitpat";
	hmap2d[key] = new TH2D(key,key,1000,-1000,1000,1000,-500,500);
	key ="VP3Hitpat";
	hmap2d[key] = new TH2D(key,key,1000,-1000,1000,1000,-500,500);
	key ="VP4Hitpat";
	hmap2d[key] = new TH2D(key,key,1000,-1000,1000,1000,-500,500);
	key ="VP5Hitpat";
	hmap2d[key] = new TH2D(key,key,1000,-1000,1000,1000,-500,500);

	key ="BeamGenZP_th_0_5";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);
	key ="BeamGenZP_th_5_10";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);
	key ="BeamGenZP_th_10_15";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);
	key ="BeamGenZP_th_15_20";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);
	key ="BeamGenZP_th_20_25";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);
	key ="BeamGenZP_th_25_30";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);

	
	
	key ="BeamAcptThetaP";
	hmap2d[key] = new TH2D(key,key,300,0,30,320,0.4,2.);
	key ="BeamAcptCosTP";
	hmap2d[key] =hmap2d[key] = new TH2D(key,key,100,0.85,1,80,0.4,2.);
	key ="BeamAcptCosTPhi";
	hmap2d[key] = new TH2D(key,key,100,0.85,1,100,-3.15,3.15);
	key ="BeamAcptPhiP";
	hmap2d[key] = new TH2D(key,key,100,-3.15,3.15,80,0.4,2.);
	key ="BeamAcptXThetaP";
	hmap2d[key] = new TH2D(key,key,100,-30,30,80,0.4,2.);
	key ="BeamAcptYThetaP";
	hmap2d[key] = new TH2D(key,key,100,-30,30,80,0.4,2.);




	key ="BeamAcptZP_th_0_5";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);
	key ="BeamAcptZP_th_5_10";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);
	key ="BeamAcptZP_th_10_15";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);
	key ="BeamAcptZP_th_15_20";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);
	key ="BeamAcptZP_th_20_25";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);
	key ="BeamAcptZP_th_25_30";
	hmap2d[key] = new TH2D(key,key,100,-153,-133,80,0.4,2);


}

//_____________________________________________________________________________
void
TPCAnaManager::EndOfRunAction( void )
{
		TPC_g->Write();
	for( auto& p : hmap ){
    p.second->Write();
  }
	for( auto& p : hmap2d ){
    p.second->Write();
  }
}

//_____________________________________________________________________________
void
TPCAnaManager::BeginOfEventAction( void )
{
  HitNum=0;
  tpctrNum=0;

  //for K+
  HitNum_K=0;
  //  tpctrNum_K=0;

  //for proton
  HitNum_p=0;
  //  tpctrNum_K=0;

  event.nhBh2 = 0;
  event.nhTgt = 0;
  event.nhHtof = 0;
  event.nhSdc = 0;
  event.nhSch = 0;
  event.nhFtof = 0;
  event.nhLac = 0;
  event.nhWc = 0;
  event.nhBvh = 0;
  event.nhVp = 0;
  for( G4int i=0; i<MaxHits; ++i ){
    // BH2
    event.tidBh2[i] = -9999;
    event.pidBh2[i] = -9999;
    event.didBh2[i] = -9999;
    event.prtBh2[i] = -9999;
    event.qBh2[i] = -9999;
    event.massBh2[i] = -9999.;
    event.xBh2[i] = -9999.;
    event.yBh2[i] = -9999.;
    event.zBh2[i] = -9999.;
    event.pxBh2[i] = -9999.;
    event.pyBh2[i] = -9999.;
    event.pzBh2[i] = -9999.;
    event.tBh2[i] = -9999.;
    event.vtppBh2[i] = -9999.;
    event.vtpxBh2[i] = -9999.;
    event.vtpyBh2[i] = -9999.;
    event.vtpzBh2[i] = -9999.;
    event.vtxBh2[i] = -9999.;
    event.vtyBh2[i] = -9999.;
    event.vtzBh2[i] = -9999.;
    // TARGET
    event.pidTgt[i] = -9999;
    event.tidTgt[i] = -9999;
    event.prtTgt[i] = -9999;
    event.xTgt[i] = -9999.;
    event.yTgt[i] = -9999.;
    event.zTgt[i] = -9999.;
    event.uTgt[i] = -9999.;
    event.vTgt[i] = -9999.;
    event.vtxTgt[i] = -9999.;
    event.vtyTgt[i] = -9999.;
    event.vtzTgt[i] = -9999.;
    // HTOF
    event.tidHtof[i] = -9999;
    event.pidHtof[i] = -9999;
    event.didHtof[i] = -9999;
    event.prtHtof[i] = -9999;
    event.qHtof[i] = -9999;
    event.massHtof[i] = -9999.;
    event.xHtof[i] = -9999.;
    event.yHtof[i] = -9999.;
    event.zHtof[i] = -9999.;
    event.pxHtof[i] = -9999.;
    event.pyHtof[i] = -9999.;
    event.pzHtof[i] = -9999.;
    event.ppHtof[i] = -9999.;
    event.tHtof[i] = -9999.;
    event.vtppHtof[i] = -9999.;
    event.vtpxHtof[i] = -9999.;
    event.vtpyHtof[i] = -9999.;
    event.vtpzHtof[i] = -9999.;
    event.vtxHtof[i] = -9999.;
    event.vtyHtof[i] = -9999.;
    event.vtzHtof[i] = -9999.;
    event.lengthHtof[i] = -9999.;
    // SDC
    event.tidSdc[i] = -9999;
    event.pidSdc[i] = -9999;
    event.didSdc[i] = -9999;
    event.prtSdc[i] = -9999;
    event.qSdc[i] = -9999;
    event.massSdc[i] = -9999.;
    event.xSdc[i] = -9999.;
    event.ySdc[i] = -9999.;
    event.zSdc[i] = -9999.;
    event.pxSdc[i] = -9999.;
    event.pySdc[i] = -9999.;
    event.pzSdc[i] = -9999.;
    event.tSdc[i] = -9999.;
    event.vtppSdc[i] = -9999.;
    event.vtpxSdc[i] = -9999.;
    event.vtpySdc[i] = -9999.;
    event.vtpzSdc[i] = -9999.;
    event.vtxSdc[i] = -9999.;
    event.vtySdc[i] = -9999.;
    event.vtzSdc[i] = -9999.;
    // SCH
    event.tidSch[i] = -9999;
    event.pidSch[i] = -9999;
    event.didSch[i] = -9999;
    event.prtSch[i] = -9999;
    event.qSch[i] = -9999;
    event.massSch[i] = -9999.;
    event.xSch[i] = -9999.;
    event.ySch[i] = -9999.;
    event.zSch[i] = -9999.;
    event.pxSch[i] = -9999.;
    event.pySch[i] = -9999.;
    event.pzSch[i] = -9999.;
    event.tSch[i] = -9999.;
    event.vtppSch[i] = -9999.;
    event.vtpxSch[i] = -9999.;
    event.vtpySch[i] = -9999.;
    event.vtpzSch[i] = -9999.;
    event.vtxSch[i] = -9999.;
    event.vtySch[i] = -9999.;
    event.vtzSch[i] = -9999.;
    // FTOF
    event.tidFtof[i] = -9999;
    event.pidFtof[i] = -9999;
    event.didFtof[i] = -9999;
    event.prtFtof[i] = -9999;
    event.qFtof[i] = -9999;
    event.massFtof[i] = -9999.;
    event.xFtof[i] = -9999.;
    event.yFtof[i] = -9999.;
    event.zFtof[i] = -9999.;
    event.pxFtof[i] = -9999.;
    event.pyFtof[i] = -9999.;
    event.pzFtof[i] = -9999.;
    event.tFtof[i] = -9999.;
    event.vtppFtof[i] = -9999.;
    event.vtpxFtof[i] = -9999.;
    event.vtpyFtof[i] = -9999.;
    event.vtpzFtof[i] = -9999.;
    event.vtxFtof[i] = -9999.;
    event.vtyFtof[i] = -9999.;
    event.vtzFtof[i] = -9999.;
    // LAC
    event.tidLac[i] = -9999;
    event.pidLac[i] = -9999;
    event.didLac[i] = -9999;
    event.prtLac[i] = -9999;
    event.qLac[i] = -9999;
    event.massLac[i] = -9999.;
    event.xLac[i] = -9999.;
    event.yLac[i] = -9999.;
    event.zLac[i] = -9999.;
    event.pxLac[i] = -9999.;
    event.pyLac[i] = -9999.;
    event.pzLac[i] = -9999.;
    event.tLac[i] = -9999.;
    event.vtppLac[i] = -9999.;
    event.vtpxLac[i] = -9999.;
    event.vtpyLac[i] = -9999.;
    event.vtpzLac[i] = -9999.;
    event.vtxLac[i] = -9999.;
    event.vtyLac[i] = -9999.;
    event.vtzLac[i] = -9999.;
    // WC
    event.tidWc[i] = -9999;
    event.pidWc[i] = -9999;
    event.didWc[i] = -9999;
    event.prtWc[i] = -9999;
    event.qWc[i] = -9999;
    event.massWc[i] = -9999.;
    event.xWc[i] = -9999.;
    event.yWc[i] = -9999.;
    event.zWc[i] = -9999.;
    event.pxWc[i] = -9999.;
    event.pyWc[i] = -9999.;
    event.pzWc[i] = -9999.;
    event.tWc[i] = -9999.;
    event.vtppWc[i] = -9999.;
    event.vtpxWc[i] = -9999.;
    event.vtpyWc[i] = -9999.;
    event.vtpzWc[i] = -9999.;
    event.vtxWc[i] = -9999.;
    event.vtyWc[i] = -9999.;
    event.vtzWc[i] = -9999.;
    //Bvh
		event.tidBvh[i] = -9999;
    event.pidBvh[i] = -9999;
    event.didBvh[i] = -9999;
    event.prtBvh[i] = -9999;
    event.qBvh[i] = -9999;
    event.massBvh[i] = -9999.;
    event.xBvh[i] = -9999.;
    event.yBvh[i] = -9999.;
    event.zBvh[i] = -9999.;
    event.pxBvh[i] = -9999.;
    event.pyBvh[i] = -9999.;
    event.pzBvh[i] = -9999.;
    event.tBvh[i] = -9999.;
    event.vtppBvh[i] = -9999.;
    event.vtpxBvh[i] = -9999.;
    event.vtpyBvh[i] = -9999.;
    event.vtpzBvh[i] = -9999.;
    event.vtxBvh[i] = -9999.;
    event.vtyBvh[i] = -9999.;
    event.vtzBvh[i] = -9999.;
    // VP
    event.tidVp[i] = -9999;
    event.pidVp[i] = -9999;
    event.didVp[i] = -9999;
    event.prtVp[i] = -9999;
    event.qVp[i] = -9999;
    event.massVp[i] = -9999.;
    event.xVp[i] = -9999.;
    event.yVp[i] = -9999.;
    event.zVp[i] = -9999.;
    event.pxVp[i] = -9999.;
    event.pyVp[i] = -9999.;
    event.pzVp[i] = -9999.;
    event.tVp[i] = -9999.;
    event.vtppVp[i] = -9999.;
    event.vtpxVp[i] = -9999.;
    event.vtpyVp[i] = -9999.;
    event.vtpzVp[i] = -9999.;
    event.vtxVp[i] = -9999.;
    event.vtyVp[i] = -9999.;
    event.vtzVp[i] = -9999.;
  }

  event.nhittpc = 0;
  event.ntrtpc = 0;

  event.HitNum_K=-1;

  event.HitNum_p=-1;


  event.mm_d = 0.;
  event.mm_p = 0.;
  event.theta = 0.;
  event.theta_scat = 0.;
  event.theta_CM = 0.;


  /* ntrtpc initialization */

  for( G4int i=0; i<MaxHitsTPC;++i ){
    event.trpidtpc[i]  = -1;


    event.trpptpc[i]  = -9999.9999;
    event.trpttpc[i]  = -9999.9999;
    event.trpxtpc[i]  = -9999.9999;
    event.trpytpc[i]  = -9999.9999;
    event.trpztpc[i]  = -9999.9999;

    event.vtpxtpc[i]  = -9999.9999;
    event.vtpytpc[i]  = -9999.9999;
    event.vtpztpc[i]  = -9999.9999;

    event.vtxtpc[i]  = -9999.9999;
    event.vtytpc[i]  = -9999.9999;
    event.vtztpc[i]  = -9999.9999;

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

    event.xtpc_pad[i] = -9999.9;
    event.ytpc_pad[i] = -9999.9;
    event.ztpc_pad[i] = -9999.9;
    event.dxtpc_pad[i] = -9999.9;
    event.dytpc_pad[i] = -9999.9;
    event.dztpc_pad[i] = -9999.9;



    event.x0tpc[i] = -9999.9;
    event.y0tpc[i] = -9999.9;
    event.z0tpc[i] = -9999.9;
    event.resoX[i] = -9999.9;
    event.resxtpc[i] = -9999.9;
    event.resytpc[i] = -9999.9;
    event.resztpc[i] = -9999.9;

    event.pxtpc[i] = -9999.9;
    event.pytpc[i] = -9999.9;
    event.pztpc[i] = -9999.9;
    event.pptpc[i] = -9999.9;

    event.masstpc[i] = -9999.9;

    event.timetpc[i] = -9999.9;
    event.betatpc[i] = -9999.9;

    event.edeptpc[i] = -9999.9;

    event.ititpc[i] = -1;
    event.idtpc[i] = -1;
    event.iPadtpc[i] = -1;
    event.laytpc[i] = -1;
    event.rowtpc[i] = -1;
    event.parentID[i] = -1;
  }
	event.SpinXi=0;
	event.SpinXi_x=0;
	event.SpinXi_y=0;
	event.SpinXi_z=0;
	event.MomXi=0;
	event.MomXi_x=0;
	event.MomXi_y=0;
	event.MomXi_z=0;
	event.SpinLd=0;
	event.SpinLd_x=0;
	event.SpinLd_y=0;
	event.SpinLd_z=0;
	event.MomLd=0;
	event.MomLd_x=0;
	event.MomLd_y=0;
	event.MomLd_z=0;
	event.MomP=0;
	event.MomP_x=0;
	event.MomP_y=0;
	event.MomP_z=0;
	event.CM_E=0;
	event.CM_x=0;
	event.CM_y=0;
	event.CM_z=0;
	event.ThXi_CM=0;
	event.ThLd_CM=0;
}

//_____________________________________________________________________________
int
TPCAnaManager::EndOfEventAction( void )
{
  event.evnum++;
	auto Mat2D = MatrixReader::Mat2D;

    if(tpctrNum>9){
      G4cout<<"Error--> over the number of tracks in the TPC:"<<tpctrNum<<G4endl;
		}
  //Fill Primary Infomation for E27
  if( m_experiment == 27 || m_experiment == 45 ){
    event.mm_d = primaryInfo.mm_d;
    event.mm_p = primaryInfo.mm_p;
    event.theta = primaryInfo.theta;
    event.theta_scat = primaryInfo.theta_scat;
    event.theta_CM = primaryInfo.theta_CM;
    event.mm = CLHEP::mm;
  }
	for(int it=0;it<1000;++it){
		event.NumberOfTracks = gTrackBuffer.GetNumberOfTracks();
		event.PIDOfTrack[it] = gTrackBuffer.GetPIDOfTrack()[it];
		event.ParentIDOfTrack[it] = gTrackBuffer.GetParentIDOfTrack()[it];
		event.VertexOfTrack_x[it] = gTrackBuffer.GetVertexOfTrack_x()[it];
		event.VertexOfTrack_y[it] = gTrackBuffer.GetVertexOfTrack_y()[it];
		event.VertexOfTrack_z[it] = gTrackBuffer.GetVertexOfTrack_z()[it];
		event.MomentumOfTrack[it] = gTrackBuffer.GetMomentumOfTrack()[it];
		event.MomentumOfTrack_x[it] = gTrackBuffer.GetMomentumOfTrack_x()[it];
		event.MomentumOfTrack_y[it] = gTrackBuffer.GetMomentumOfTrack_y()[it];
		event.MomentumOfTrack_z[it] = gTrackBuffer.GetMomentumOfTrack_z()[it];
	}
  auto SXi = gTrackBuffer.GetPolarity(0);
  auto PXi = gTrackBuffer.GetMomentum(0);
	event.SpinXi = SXi.mag();
	event.SpinXi_x = SXi.x();
	event.SpinXi_y = SXi.y();
	event.SpinXi_z = SXi.z();
	event.MomXi = PXi.mag();
	event.MomXi_x = PXi.x();
	event.MomXi_y = PXi.y();
	event.MomXi_z = PXi.z();
  
	auto PLd = gTrackBuffer.GetVertexMomentum(1);
	event.MomLd = PLd.mag();
	event.MomLd_x = PLd.x();
	event.MomLd_y = PLd.y();
	event.MomLd_z = PLd.z();

	auto XiLV = gTrackBuffer.GetLV(0);
	auto XiU = XiLV.boostVector();
	auto CMLV = gTrackBuffer.GetCMLV();
	event.CM_E = CMLV.t();
	event.CM_x = CMLV.x();
	event.CM_y = CMLV.y();
	event.CM_z = CMLV.z();
	auto LdVLV = gTrackBuffer.GetVertexLV(1);
	auto PLd_Lab = LdVLV.vect();
	double ThXi_Lab = acos(PLd_Lab*SXi*(1./SXi.mag()/PLd_Lab.mag()));
	LdVLV.boost(-XiU);
	auto PLd_CM = LdVLV.vect();
	
	double ThXi_CM = acos(PLd_CM*SXi*(1./SXi.mag()/PLd_CM.mag()));
	event.ThXi_CM = ThXi_CM;



	auto LdLV = gTrackBuffer.GetLV(1);
  auto SLd = gTrackBuffer.GetPolarity(1);
	auto LdU = LdLV.boostVector();

	auto PVLV = gTrackBuffer.GetVertexLV(2);
	auto PP = PVLV.vect();
	event.SpinLd = SLd.mag();
	event.SpinLd_x = SLd.x();
	event.SpinLd_y = SLd.y();
	event.SpinLd_z = SLd.z();
	event.MomP = PP.mag();
	event.MomP_x = PP.x();
	event.MomP_y = PP.y();
	event.MomP_z = PP.z();
	PVLV.boost(-LdU);
	auto PP_CM = PVLV.vect();
	double ThLd_CM = acos(PP_CM*SLd*(1./SLd.mag()/PP_CM.mag()));
	event.ThLd_CM = ThLd_CM;

	if( HitNum > 0){
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



    G4int sh_paID[MAX_TRACK] = {};
    G4int sh_paPID[MAX_TRACK] = {};
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
	// test[kk]=circleFit(x[kk],z[kk],y[kk],c[kk],&cx[kk],&cz[kk],&rad[kk],&Pz[kk],
	// 		   &a_fory[kk], &b_fory[kk], &theta0_fory[kk]);
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
      G4double cx2 = event.xPrm[0];
      G4double cz2 = event.zPrm[0];
      G4double theta12=atan2(cz2-cz1, cx2-cx1);
      G4double ca1=a_fory[i];
      G4double cb1=b_fory[i];
      G4double ct01=theta0_fory[i];


      G4double cent_dist=sqrt(pow(cx1-cx2,2)+pow(cz1-cz2,2));

      vtxxfit[i]=cos(theta12)*cent_dist+cx1;
      vtxzfit[i]=sin(theta12)*cent_dist+cz1;
      vtxyfit[i]=-1.*tpcData[i].tpcqq*ca1*rho1*(theta12-ct01)+cb1;

      mom_theta[i]=atan2(vtxzfit[i]-cz[i],vtxxfit[i]-cx[i])-acos(-1.)/2;

      // vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
      // vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);

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


	      // vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      // vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      // vtxpxfit[j]=cos(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
	      // vtxpzfit[j]=sin(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
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

	      // vtxpxfit[i]=cos(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      // vtxpzfit[i]=sin(mom_theta[i])*(cir_r[i])*(-0.299792458)*(env_helm_field)*(tpcData[i].tpcqq);
	      // vtxpxfit[j]=cos(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
	      // vtxpzfit[j]=sin(mom_theta[j])*(cir_r[j])*(-0.299792458)*(env_helm_field)*(tpcData[j].tpcqq);
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
    if( HitNum > MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit (TPC) " << HitNum << G4endl;
    } else {

      for( G4int i=0; i<HitNum; i++){
	event.ntrk[event.nhittpc] = counterData[i].ntrk;
	hmap["Time"]->Fill( counterData[i].time );
	for( G4int j=0; j<G4ThreeVector::SIZE; ++j ){
	  hmap[Form( "Pos%d", j )]->Fill( counterData[i].pos[j]/CLHEP::mm );
	  hmap[Form( "Mom%d", j )]->Fill( counterData[i].mom[j]/CLHEP::GeV );
	}
	event.xtpc[event.nhittpc] = counterData[i].pos[0]/CLHEP::mm;
	event.ytpc[event.nhittpc] = counterData[i].pos[1]/CLHEP::mm;
	event.ztpc[event.nhittpc] = counterData[i].pos[2]/CLHEP::mm;

	event.x0tpc[event.nhittpc] = counterData[i].pos0[0]/CLHEP::mm;
	event.y0tpc[event.nhittpc] = counterData[i].pos0[1]/CLHEP::mm;
	event.z0tpc[event.nhittpc] = counterData[i].pos0[2]/CLHEP::mm;

	event.resoX[event.nhittpc] = counterData[i].resoX;
	event.resxtpc[event.nhittpc] = counterData[i].resxtpc;
	event.resytpc[event.nhittpc] = counterData[i].resytpc;
	event.resztpc[event.nhittpc] = counterData[i].resztpc;
	event.pxtpc[event.nhittpc] = counterData[i].mom[0]/CLHEP::GeV;
	event.pytpc[event.nhittpc] = counterData[i].mom[1]/CLHEP::GeV;
	event.pztpc[event.nhittpc] = counterData[i].mom[2]/CLHEP::GeV;
	event.pptpc[event.nhittpc] = sqrt(pow(counterData[i].mom[0], 2) +
					pow(counterData[i].mom[1], 2) +
					pow(counterData[i].mom[2], 2))/CLHEP::GeV;
	event.ititpc[event.nhittpc] = counterData[i].trackID;
	event.idtpc[event.nhittpc] = counterData[i].particleID;
	event.laytpc[event.nhittpc] = counterData[i].iLay;

	event.rowtpc[event.nhittpc] = counterData[i].iRow;
	event.iPadtpc[event.nhittpc] = padHelper::getPadID(event.laytpc[event.nhittpc], event.rowtpc[event.nhittpc]);
	TVector3 Point = padHelper::getPoint(event.iPadtpc[event.nhittpc]);
	event.xtpc_pad[event.nhittpc] = Point.x();
	event.ytpc_pad[event.nhittpc] = event.ytpc[event.nhittpc];
	event.ztpc_pad[event.nhittpc] = Point.z();

	event.dxtpc_pad[event.nhittpc] = event.x0tpc[event.nhittpc] - event.xtpc_pad[event.nhittpc];
	event.dytpc_pad[event.nhittpc] = event.y0tpc[event.nhittpc] - event.ytpc_pad[event.nhittpc];
	event.dztpc_pad[event.nhittpc] = event.z0tpc[event.nhittpc] - event.ztpc_pad[event.nhittpc];


	event.timetpc[event.nhittpc] = counterData[i].time/CLHEP::ns;
	event.betatpc[event.nhittpc] = counterData[i].beta;
	event.edeptpc[event.nhittpc] = counterData[i].edep/(CLHEP::MeV/CLHEP::mm);
	event.dedxtpc[event.nhittpc] = counterData[i].dedx;
	event.slengthtpc[event.nhittpc] = counterData[i].slength/CLHEP::mm;
	event.tlengthtpc[event.nhittpc] = counterData[i].tlength/CLHEP::mm;
	event.nthlay[event.nhittpc] = counterData[i].iLay;
	event.nthpad[event.nhittpc] = counterData[i].iPad;
	event.laypad[event.nhittpc][event.nthlay[event.nhittpc]][event.nthpad[event.nhittpc]]
	  = event.laypad[event.nhittpc][event.nthlay[event.nhittpc]][event.nthpad[event.nhittpc]]+1.;
	event.parentID[event.nhittpc] = counterData[i].parentID;
	event.nhittpc += 1;
      }
    }
    //
    // TPC
    //
    for( G4int i=0; i<tpctrNum; i++){
      //      G4cout<<"abs"<<abs(env_helm_field)<<G4endl;
      //      G4cout<<"fabs"<<fabs(env_helm_field)<<G4endl;
      // anaRoot.FillTPCData(tpcData[i].tpcpx,
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
    }
  }//trigger parts
	if(DiscardData){
		event.Clear();
	}
	else{
  	TPC_g->Fill();
		event.Clear();
	}
  event.pb->SetXYZ( 0., 0., 0. );
  event.nhPrm = 0;
  for( Int_t i=0; i<MaxPrimaryParticle; ++i ){
    event.pidPrm[i] = -9999;
    event.xPrm[i] = -9999.;
    event.yPrm[i] = -9999.;
    event.zPrm[i] = -9999.;
    event.pxPrm[i] = -9999.;
    event.pyPrm[i] = -9999.;
    event.pzPrm[i] = -9999.;
    event.ppPrm[i] = -9999.;
    event.mPrm[i] = -9999.;
    event.thetaPrm[i] = -9999.;
    event.phiPrm[i] = -9999.;
  }
	bool Trig = false;
	bool SDC = false;
	int nhSdc1=0;
	int nhSdc2=0;
	int nhSdc3=0;
	int nhSdc4=0;
	double xSdc1=-9999;
	double ySdc1=-9999;
	double xSdc2=-9999;
	double ySdc2=-9999;
	double xSdc3=-9999;
	double ySdc3=-9999;
	double xSdc4=-9999;
	double ySdc4=-9999;
	for(int ih=0;ih<event.nhSdc;++ih){
		if(event.tidSdc[ih]!=1) continue;
		double zsdc = event.zSdc[ih];
		if(zsdc<0) continue;
		if(zsdc<1000){
			nhSdc1++;
			if(zsdc<787){
				xSdc1 = event.xSdc[ih];
				ySdc1 = event.ySdc[ih];
			}
		}
		else if(zsdc<1500){
			nhSdc2++;
			if(zsdc < 1240){
				xSdc2 = event.xSdc[ih];
				ySdc2 = event.ySdc[ih];
			}
		}
		else if(zsdc<2700){
			nhSdc3++;
			if(zsdc>2644){
				xSdc3 = event.xSdc[ih];
				ySdc3 = event.ySdc[ih];
			}
		}
		else if(zsdc<3000){
			nhSdc4++;
			if(zsdc>2905){
				xSdc4 = event.xSdc[ih];
				ySdc4 = event.ySdc[ih];
			}
		}
	}
	int nhSdcIn= nhSdc1+nhSdc2;
	int nhSdcOut= nhSdc3+nhSdc4;
	if(nhSdcIn>=8 and nhSdcOut>=6)SDC = true;
	int ToFHit = -1;
	int ToFHit2 = -1;
	int SchHit = -1;
	for(int ih=0;ih<event.nhFtof;++ih){
		if(event.tidFtof[ih] !=1) continue;
		ToFHit = event.didFtof[ih] ;
	}
	for(int ih=0;ih<event.nhBvh;++ih){
		if(event.tidBvh[ih] !=1) continue;
		int ibbb= event.didBvh[ih];
		if(ibbb < 2) ToFHit2 = 24;
		else if(ibbb < 5) ToFHit2 = 25;
		else if(ibbb < 9) ToFHit2 = 26;
		else if(ibbb < 14) ToFHit2 = 27;
	}
	for(int ih=0;ih<event.nhSch;++ih){
		if(event.tidSch[ih] !=1) continue;
		SchHit = event.didSch[ih];
	}
	/*	for(int itof=0;itof<24;++itof){
		for(int isch=0;isch<64;++isch){
		}
	}
	*/
	if(ToFHit > 0 and SchHit > 0){
		if(Mat2D[ToFHit][SchHit])Trig = true;	
	}
	if(ToFHit2 > 0 and SchHit > 0){
		if(Mat2D[ToFHit2][SchHit])Trig = true;	
	}

	double xVP1=-9999;
	double yVP1=-9999;
	double xVP2=-9999;
	double yVP2=-9999;
	double xVP3=-9999;
	double yVP3=-9999;
	double xVP4=-9999;
	double yVP4=-9999;
	double xVP5=-9999;
	double yVP5=-9999;
	for(int ih=0;ih<event.nhVp;++ih){
		if(event.tidVp[ih]!=1)continue;
		TVector3 KuramaPos(50,0,1719.5);
		TVector3 VPPos(event.xVp[ih],event.yVp[ih],event.zVp[ih]);
		VPPos-=KuramaPos;
		if(abs(VPPos.Z()+400)<1){
			xVP1=VPPos.X();
			yVP1=VPPos.y();
		}
		if(abs(VPPos.Z()+200)<1){
			xVP2=VPPos.X();
			yVP2=VPPos.y();
		}
		if(abs(VPPos.Z()+0)<1){
			xVP3=VPPos.X();
			yVP3=VPPos.y();
		}
		if(abs(VPPos.Z()-200)<1){
			xVP4=VPPos.X();
			yVP4=VPPos.y();
		}
		if(abs(VPPos.Z()-400)<1){
			xVP5=VPPos.X();
			yVP5=VPPos.y();
		}
	}


	
	double pk = event.MomentumOfTrack[1]*0.001;
	double pkx = event.MomentumOfTrack_x[1]*0.001;
	double pky = event.MomentumOfTrack_y[1]*0.001;
	double pkz = event.MomentumOfTrack_z[1]*0.001;
	G4ThreeVector pkk(pkx,pky,pkz);
	double pkth = cos(pkk.theta());
	double pkph = pkk.phi();
	double vtx_z = event.VertexOfTrack_z[1];	
	double pkangle = acos(pkth)*180./acos(-1); 
	
	TString key = "BeamGenThetaP"; 
	auto H0 = hmap2d[key];
	key = "BeamGenCosTP"; 
	auto H1 = hmap2d[key];
	key ="BeamGenCosTPhi";
	auto H2 = hmap2d[key];
	key ="BeamGenPhiP";
	auto H3 = hmap2d[key];

	key ="BeamGenZP_th_0_5";
	auto H4 = hmap2d[key];
	key ="BeamGenZP_th_5_10";
	auto H5 = hmap2d[key];
	key ="BeamGenZP_th_10_15";
	auto H6 = hmap2d[key];
	key ="BeamGenZP_th_15_20";
	auto H7 = hmap2d[key];
	key ="BeamGenZP_th_20_25";
	auto H8 = hmap2d[key];
	key ="BeamGenZP_th_25_30";
	auto H9 = hmap2d[key];


	key ="BeamGenXThetaP";
	auto H10 = hmap2d[key];
	key ="BeamGenYThetaP";
	auto H11 = hmap2d[key];



	key ="Sdc1Hitpat";
	auto Sdc1Hit = hmap2d[key];
	key ="Sdc2Hitpat";
	auto Sdc2Hit = hmap2d[key];
	key ="Sdc3Hitpat";
	auto Sdc3Hit = hmap2d[key];
	key ="Sdc3Hitpat";
	auto Sdc4Hit = hmap2d[key];
	
	key ="VP1Hitpat";
	auto VP1Hit = hmap2d[key];
	key ="VP2Hitpat";
	auto VP2Hit = hmap2d[key];
	key ="VP3Hitpat";
	auto VP3Hit = hmap2d[key];
	key ="VP4Hitpat";
	auto VP4Hit = hmap2d[key];
	key ="VP5Hitpat";
	auto VP5Hit = hmap2d[key];


	key ="BeamAcptThetaP";
	auto HA0 = hmap2d[key];	
	key ="BeamAcptCosTP";
	auto HA1 = hmap2d[key];	
	key ="BeamAcptCosTPhi";
	auto HA2 = hmap2d[key];	
	key ="BeamAcptPhiP";
	auto HA3 = hmap2d[key];	
	
	key ="BeamAcptZP_th_0_5";
	auto HA4 = hmap2d[key];
	key ="BeamAcptZP_th_5_10";
	auto HA5 = hmap2d[key];
	key ="BeamAcptZP_th_10_15";
	auto HA6 = hmap2d[key];
	key ="BeamAcptZP_th_15_20";
	auto HA7 = hmap2d[key];
	key ="BeamAcptZP_th_20_25";
	auto HA8 = hmap2d[key];
	key ="BeamAcptZP_th_25_30";
	auto HA9 = hmap2d[key];
	
	key ="BeamAcptXThetaP";
	auto HA10 = hmap2d[key];
	key ="BeamAcptYThetaP";
	auto HA11 = hmap2d[key];

	H0->Fill(pkangle,pk);
	H1->Fill(pkth,pk);
	H2->Fill(pkth,pkph);
	H3->Fill(pkph,pk);
	H10->Fill(180*asin(pkx/pk)/acos(-1),pk);
	H11->Fill(180*asin(pky/pk)/acos(-1),pk);
{
	if(pkangle < 5){
		H4->Fill(vtx_z,pk);
	}
	else if(pkangle < 10){
		H5->Fill(vtx_z,pk);
	}
	else if(pkangle<15){
		H6->Fill(vtx_z,pk);
	}
	else if(pkangle < 20){
		H7->Fill(vtx_z,pk);
	}
	else if(pkangle < 25){
		H8->Fill(vtx_z,pk);
	}
	else{
		H9->Fill(vtx_z,pk);
	};
}

	if(Trig and SDC){
		HA0->Fill(pkangle,pk);
		HA1->Fill(pkth,pk);
		HA2->Fill(pkth,pkph);
		HA3->Fill(pkph,pk);
		if(pkangle < 5){
			HA4->Fill(vtx_z,pk);
		}
		else if(pkangle < 10){
			HA5->Fill(vtx_z,pk);
		}
		else if(pkangle<15){
			HA6->Fill(vtx_z,pk);
		}
		else if(pkangle < 20){
			HA7->Fill(vtx_z,pk);
		}
		else if(pkangle < 25){
			HA8->Fill(vtx_z,pk);
		}
		else{
			HA9->Fill(vtx_z,pk);
		}
		HA10->Fill(180*asin(pkx/pk)/acos(-1),pk);
		HA11->Fill(180*asin(pky/pk)/acos(-1),pk);
		Sdc1Hit->Fill(xSdc1,ySdc1);	
		Sdc2Hit->Fill(xSdc2,ySdc2);	
		Sdc3Hit->Fill(xSdc3,ySdc3);	
		Sdc4Hit->Fill(xSdc4,ySdc4);	
		VP1Hit->Fill(xVP1,yVP1);	
		VP2Hit->Fill(xVP2,yVP2);	
		VP3Hit->Fill(xVP3,yVP3);	
		VP4Hit->Fill(xVP4,yVP4);	
		VP5Hit->Fill(xVP5,yVP5);	
	}
  return 0;
}

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetBH2Data( const VHitInfo* hit )
  {
    if( event.nhBh2 >= MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit " << event.nhBh2 << G4endl;
    } else {
      Int_t i = event.nhBh2;
      event.tidBh2[i] = hit->GetTrackID();
      event.pidBh2[i] = hit->GetParticleID();
      event.didBh2[i] = hit->GetDetectorID();
      event.prtBh2[i] = hit->GetParentID();
      event.xBh2[i] = hit->GetPosition().x();
      event.yBh2[i] = hit->GetPosition().y();
      event.zBh2[i] = hit->GetPosition().z();
      event.pxBh2[i] = hit->GetMomentum().x();
      event.pyBh2[i] = hit->GetMomentum().y();
      event.pzBh2[i] = hit->GetMomentum().z();
      event.ppBh2[i] = hit->GetMomentum().mag();
      event.deBh2[i] = hit->GetEnergyDeposit();
      event.tBh2[i] = hit->GetTime();
      event.nhBh2++;
    }
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetCounterData( G4int ntrk, G4double time, G4ThreeVector pos,
				   G4ThreeVector mom,
				   G4int track, G4int particle,
				   G4int iLay,  G4int iRow, G4double beta,
				   G4double edep, G4int parentid,
				   G4double tlength, G4double slength )
  {
    G4int hitnum = HitNum;
    G4bool flag=true;
    if (hitnum > MaxTrack) {
      fprintf(stderr, "TPCAnaManager::SetCounterData Too Much multiplicity %d\n",
	      hitnum);
      return;
    }
		
    //  G4ThreeVector tar_pos(0.,0.*CLHEP::mm,-150.*CLHEP::mm);

    G4ThreeVector tar_pos(0.,0., -143);
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
      counterData[hitnum].tlength = tlength;



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
			std::vector<double>ResPar;
			if(iLay < 10){
				ResPar = ResParamInnerLayerHSOn;
			}
			else{
				ResPar = ResParamOuterLayerHSOn;
			}
			double par_t[6]={
				ResPar[0],ResPar[1],ResPar[2],ResPar[3],ResPar[4],ResPar[5]};
			double par_y[4] = {
				ResPar[6],ResPar[1],ResPar[7],ResPar[8]};
      // G4double compy=0.;
      G4double compx=0.;
			auto SmearingVector =
				GetSmearingVector(sh_pos,mom,par_y,par_t);

			auto ResVector =
				GetResVector(sh_pos,mom,par_y,par_t);
			compx = ResVector.mag();
  /*
			compx = GetTransverseRes(sh_y);
      double s_compx = CLHEP::RandGauss::shoot(0.,compx);

      // std::cout<<"compx="<<compx<<", sh_sigmaY"<<sh_sigmaY<<std::endl;
      // getchar();
      //G4double sh_dalpha = compx/sh_rho; // rho * theta = arc --> from sako-san's code
      //      G4double sh_dalpha = s_compx/sh_rho; // rho * theta = arc --> from sako-san's code
      G4double sh_dalpha = atan2(s_compx, sh_rho); // rho * theta = arc --> from sako-san's code
      G4double sh_smear_alpha = sh_alpha+sh_dalpha;
      //    G4cout<<compx<<":"<<sh_dalpha<<G4endl;
      G4double randx = sh_rho*(sin(sh_smear_alpha)-sin(sh_alpha));
      G4double randz = sh_rho*(cos(sh_smear_alpha)-cos(sh_alpha));
      //G4double s0 = 0.204;// mm HIMAC result
      G4double s0 = 0.199;// mm HIMAC result
      G4double randx_com = CLHEP::RandGauss::shoot(0.,s0);
      G4double randz_com = CLHEP::RandGauss::shoot(0.,s0);
      G4double smear_x = sqrt(randx*randx + randx_com*randx_com )*(randx_com/fabs(randx_com));
      G4double smear_z = sqrt(randz*randz + randz_com*randz_com )*(randz_com/fabs(randz_com));
*/
      counterData[hitnum].resoX = compx;
      counterData[hitnum].resxtpc = ResVector.x();
      counterData[hitnum].resytpc = ResVector.y();
      counterData[hitnum].resztpc = ResVector.z();

     // counterData[hitnum].pos[G4ThreeVector::Z] = sh_rho*cos(sh_alpha)+smear_z+tar_pos.getZ();
      // counterData[hitnum].pos[G4ThreeVector::X] = sh_rho*sin(sh_alpha)+smear_x;
 
			counterData[hitnum].pos[G4ThreeVector::X] = SmearingVector.x()+pos.x();
			counterData[hitnum].pos[G4ThreeVector::Y] = SmearingVector.y()+pos.y();
			counterData[hitnum].pos[G4ThreeVector::Z] = SmearingVector.z()+pos.z();

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

      if( m_pad_config == 2 ){
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
    TPCAnaManager::SetFermiMomentum( const G4ThreeVector& p )
  {
    for( G4int i=0; i<G4ThreeVector::SIZE; ++i ){
      TString key = Form( "Fermi%d", i );
      hmap[key]->Fill( p[i]*CLHEP::GeV );
    }
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetFTOFData( const VHitInfo* hit )
  {
    if( event.nhFtof >= MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit " << event.nhFtof << G4endl;
    } else {
      Int_t i = event.nhFtof;
      event.tidFtof[i] = hit->GetTrackID();
      event.pidFtof[i] = hit->GetParticleID();
      event.didFtof[i] = hit->GetDetectorID();
      event.prtFtof[i] = hit->GetParentID();
      event.xFtof[i] = hit->GetPosition().x();
      event.yFtof[i] = hit->GetPosition().y();
      event.zFtof[i] = hit->GetPosition().z();
      event.pxFtof[i] = hit->GetMomentum().x();
      event.pyFtof[i] = hit->GetMomentum().y();
      event.pzFtof[i] = hit->GetMomentum().z();
      event.ppFtof[i] = hit->GetMomentum().mag();
      event.deFtof[i] = hit->GetEnergyDeposit();
      event.tFtof[i] = hit->GetTime();
      event.nhFtof++;
    }
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetHTOFData( const VHitInfo* hit )
  {
    if( event.nhHtof >= MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit " << event.nhHtof << G4endl;
    } else {
      Int_t i = event.nhHtof;
      event.tidHtof[i] = hit->GetTrackID();
      event.pidHtof[i] = hit->GetParticleID();
      event.didHtof[i] = hit->GetDetectorID();
      event.prtHtof[i] = hit->GetParentID();
      event.qHtof[i] = hit->GetCharge();
      event.massHtof[i] = hit->GetMass();
      event.xHtof[i] = hit->GetPosition().x();
      event.yHtof[i] = hit->GetPosition().y();
      event.zHtof[i] = hit->GetPosition().z();
      event.pxHtof[i] = hit->GetMomentum().x();
      event.pyHtof[i] = hit->GetMomentum().y();
      event.pzHtof[i] = hit->GetMomentum().z();
      event.ppHtof[i] = hit->GetMomentum().mag();
      event.deHtof[i] = hit->GetEnergyDeposit();
      event.tHtof[i] = hit->GetTime();
      event.vtppHtof[i] = hit->GetVertexMomentum().mag();
      event.vtpxHtof[i] = hit->GetVertexMomentum().x();
      event.vtpyHtof[i] = hit->GetVertexMomentum().y();
      event.vtpzHtof[i] = hit->GetVertexMomentum().z();
      event.vtxHtof[i] = hit->GetVertexPosition().x();
      event.vtyHtof[i] = hit->GetVertexPosition().y();
      event.vtzHtof[i] = hit->GetVertexPosition().z();
      event.lengthHtof[i] = hit->GetTrackLength();
      event.nhHtof++;
    }
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetLACData( const VHitInfo* hit )
  {
    if( event.nhLac >= MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit " << event.nhLac << G4endl;
    } else {
      Int_t i = event.nhLac;
      event.tidLac[i] = hit->GetTrackID();
      event.pidLac[i] = hit->GetParticleID();
      event.didLac[i] = hit->GetDetectorID();
      event.prtLac[i] = hit->GetParentID();
      event.xLac[i] = hit->GetPosition().x();
      event.yLac[i] = hit->GetPosition().y();
      event.zLac[i] = hit->GetPosition().z();
      event.pxLac[i] = hit->GetMomentum().x();
      event.pyLac[i] = hit->GetMomentum().y();
      event.pzLac[i] = hit->GetMomentum().z();
      event.ppLac[i] = hit->GetMomentum().mag();
      event.deLac[i] = hit->GetEnergyDeposit();
      event.tLac[i] = hit->GetTime();
      event.nhLac++;
    }
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetSCHData( const VHitInfo* hit )
  {
    if( event.nhSch >= MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit " << event.nhSch << G4endl;
    } else {
      Int_t i = event.nhSch;
      event.tidSch[i] = hit->GetTrackID();
      event.pidSch[i] = hit->GetParticleID();
      event.didSch[i] = hit->GetDetectorID();
      event.prtSch[i] = hit->GetParentID();
      event.xSch[i] = hit->GetPosition().x();
      event.ySch[i] = hit->GetPosition().y();
      event.zSch[i] = hit->GetPosition().z();
      event.pxSch[i] = hit->GetMomentum().x();
      event.pySch[i] = hit->GetMomentum().y();
      event.pzSch[i] = hit->GetMomentum().z();
      event.ppSch[i] = hit->GetMomentum().mag();
      event.deSch[i] = hit->GetEnergyDeposit();
      event.tSch[i] = hit->GetTime();
      event.nhSch++;
    }
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetSDCData( const VHitInfo* hit )
  {
    if( event.nhSdc >= MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit " << event.nhSdc << G4endl;
    } else {
      Int_t i = event.nhSdc;
      event.tidSdc[i] = hit->GetTrackID();
      event.pidSdc[i] = hit->GetParticleID();
      event.didSdc[i] = hit->GetDetectorID();
      event.prtSdc[i] = hit->GetParentID();
      event.xSdc[i] = hit->GetPosition().x();
      event.ySdc[i] = hit->GetPosition().y();
      event.zSdc[i] = hit->GetPosition().z();
      event.pxSdc[i] = hit->GetMomentum().x();
      event.pySdc[i] = hit->GetMomentum().y();
      event.pzSdc[i] = hit->GetMomentum().z();
      event.ppSdc[i] = hit->GetMomentum().mag();
      event.deSdc[i] = hit->GetEnergyDeposit();
      event.tSdc[i] = hit->GetTime();
      event.nhSdc++;
    }
  }

  void
    TPCAnaManager::SetBVHData( const VHitInfo* hit )
  {
    if( event.nhBvh >= MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit " << event.nhBvh << G4endl;
    } else {
      Int_t i = event.nhBvh;
      event.tidBvh[i] = hit->GetTrackID();
      event.pidBvh[i] = hit->GetParticleID();
      event.didBvh[i] = hit->GetDetectorID();
      event.prtBvh[i] = hit->GetParentID();
      event.xBvh[i] = hit->GetPosition().x();
      event.yBvh[i] = hit->GetPosition().y();
      event.zBvh[i] = hit->GetPosition().z();
      event.pxBvh[i] = hit->GetMomentum().x();
      event.pyBvh[i] = hit->GetMomentum().y();
      event.pzBvh[i] = hit->GetMomentum().z();
      event.ppBvh[i] = hit->GetMomentum().mag();
      event.deBvh[i] = hit->GetEnergyDeposit();
      event.tBvh[i] = hit->GetTime();
      event.nhBvh++;
    }
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetVPData( const VHitInfo* hit )
  {
    if( event.nhVp >= MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit " << event.nhVp << G4endl;
    } else {
      Int_t i = event.nhVp;
      event.tidVp[i] = hit->GetTrackID();
      event.pidVp[i] = hit->GetParticleID();
      event.didVp[i] = hit->GetDetectorID();
      event.prtVp[i] = hit->GetParentID();
      event.xVp[i] = hit->GetPosition().x();
      event.yVp[i] = hit->GetPosition().y();
      event.zVp[i] = hit->GetPosition().z();
      event.pxVp[i] = hit->GetMomentum().x();
      event.pyVp[i] = hit->GetMomentum().y();
      event.pzVp[i] = hit->GetMomentum().z();
      event.ppVp[i] = hit->GetMomentum().mag();
      event.deVp[i] = hit->GetEnergyDeposit();
      event.tVp[i] = hit->GetTime();
      event.nhVp++;
    }
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetWCData( const VHitInfo* hit )
  {
    if( event.nhWc >= MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit " << event.nhWc << G4endl;
    } else {
      Int_t i = event.nhWc;
      event.tidWc[i] = hit->GetTrackID();
      event.pidWc[i] = hit->GetParticleID();
      event.didWc[i] = hit->GetDetectorID();
      event.prtWc[i] = hit->GetParentID();
      event.xWc[i] = hit->GetPosition().x();
      event.yWc[i] = hit->GetPosition().y();
      event.zWc[i] = hit->GetPosition().z();
      event.pxWc[i] = hit->GetMomentum().x();
      event.pyWc[i] = hit->GetMomentum().y();
      event.pzWc[i] = hit->GetMomentum().z();
      event.ppWc[i] = hit->GetMomentum().mag();
      event.deWc[i] = hit->GetEnergyDeposit();
      event.tWc[i] = hit->GetTime();
      event.nhWc++;
    }
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
    TPCAnaManager::SetNumberOfPrimaryParticle( G4int n )
  {
    event.nhPrm = n;
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetPrimaryParticle( G4double px, G4double py, G4double pz )
  {
    Int_t id = 0;
    event.pxPrm[id] = px;
    event.pyPrm[id] = py;
    event.pzPrm[id] = pz;
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetGeneratorID( G4int gen )
  {
    event.gen = gen;
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetModeID( G4int mode )
  {
    event.mode = mode;
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetIncID( G4int inc )
  {
    event.inc = inc;
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetPrimaryParticle( G4int id, const G4ThreeVector& p,
				       G4double mass, G4int pid )
  {
    if( id >= event.nhPrm ){
      G4cerr << FUNC_NAME << " Invalid Primary particle ID" << G4endl;
    } else {
      event.pxPrm[id] = p.x();
      event.pyPrm[id] = p.y();
      event.pzPrm[id] = p.z();
      event.ppPrm[id] = p.mag();
      event.mPrm[id] = mass;
      event.thetaPrm[id] = p.theta();
      event.phiPrm[id] = p.phi();
      event.pidPrm[id] = pid;
    }
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetPrimaryParticle( G4int id, G4double px, G4double py,
				       G4double pz, G4double mass, G4int pid )
  {
    SetPrimaryParticle( id, G4ThreeVector( px, py, pz ), mass, pid );
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetPrimaryVertex( G4int id, G4double x, G4double y, G4double z )
  {
    if( id >= event.nhPrm ){
      G4cerr << FUNC_NAME << " Invalid Primary particle ID" << G4endl;
    } else {
      event.xPrm[id] = x;
      event.yPrm[id] = y;
      event.zPrm[id] = z;
    }
  }

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetPrimaryVertex( G4int id, const G4ThreeVector& x )
  {
    if( id >= event.nhPrm ){
      G4cerr << FUNC_NAME << " Invalid Primary particle ID" << G4endl;
    } else {
      event.xPrm[id] = x.x();
      event.yPrm[id] = x.y();
      event.zPrm[id] = x.z();
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
	void
		TPCAnaManager::SetRealBeamData(G4int runnum, G4int evnum, G4int* trigpat){
			event.data_runnum = runnum;
			event.data_evnum = evnum;
			for(int i=0;i<32;++i){
				event.trigpat[i] = trigpat[i];
			}
		};

  //_____________________________________________________________________________
  void
    TPCAnaManager::SetTargetData( const VHitInfo* hit )
  {
    if( event.nhTgt >= MaxHits ){
      G4cerr << FUNC_NAME << " too much nhit " << event.nhTgt << G4endl;
    } else {
      Int_t i = event.nhTgt;
      event.tidTgt[i] = hit->GetTrackID();
      event.pidTgt[i] = hit->GetParticleID();
      event.prtTgt[i] = hit->GetParentID();
      event.xTgt[i] = hit->GetPosition().x();
      event.yTgt[i] = hit->GetPosition().y();
      event.zTgt[i] = hit->GetPosition().z();
      auto pTgt = hit->GetMomentum();
			event.uTgt[i] = pTgt.x()/pTgt.z();
			event.vTgt[i] = pTgt.y()/pTgt.z();
			event.vtxTgt[i] = hit->GetVertexPosition().x();
      event.vtyTgt[i] = hit->GetVertexPosition().y();
      event.vtzTgt[i] = hit->GetVertexPosition().z();
      event.nhTgt++;
    }
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

	void TPCAnaManager::SetMMVertex(MMVertex* vert){
		event.ntK18 = vert->ntK18;
		event.xvpHS = vert->xvpHS;
		event.yvpHS = vert->yvpHS;
		event.zvpHS = vert->zvpHS;
		event.xtgtHS = vert->xtgtHS;
		event.ytgtHS = vert->ytgtHS;
		event.ztgtHS = vert->ztgtHS;
		event.xoutK18 = vert->xoutK18;
		event.youtK18 = vert->youtK18;
		event.uoutK18 = vert->uoutK18;
		event.voutK18 = vert->voutK18;
		event.p_3rd = vert->p_3rd;
		event.layerK18 = vert->layerK18;
		event.wireK18 = vert->wireK18;
		event.localhitposK18 = vert->localhitposK18;
		event.ntKurama = vert->ntKurama;
		event.xvpKurama = vert->xvpKurama;
		event.yvpKurama = vert->yvpKurama;
		event.zvpKurama = vert->zvpKurama;
		event.xtgtKurama = vert->xtgtKurama;
		event.ytgtKurama = vert->ytgtKurama;
		event.xout = vert->xout;
		event.yout = vert->yout;
		event.zout = vert->zout;
		event.pxout = vert->pxout;
		event.pyout = vert->pyout;
		event.pzout = vert->pzout;
		event.layer = vert->layer;
		event.wire = vert->wire;
		event.localhitpos = vert->localhitpos;
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

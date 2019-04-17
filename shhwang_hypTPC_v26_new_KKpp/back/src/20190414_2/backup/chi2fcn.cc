#include "chi2fcn.hh"

double chi2fcn::chi2fcn::operator()(int np, double *g, double*u, int iflag){

  //rungekutta integration has 5 parameters
  //p, q, u=(x0,y0), u'=(x0',y0'), z
  double xi=u[0], yi=u[1], ui=u[2], vi=u[3], delta=u[4];
  double xo, yo, uo, vo;
  //  int nin=??//hit in--> before magnet
  //  int nout=//hit out--> after magnet
  double chi=0.0;
  int nh=0;
  int nin=10;
  for(int i=0;i<nin;++i){

    /*
    int lnum=hit->GetLayer();
    double z=geomMan.GetLocalZ( lnum );
    double dd=geomMan.GetResolution( lnum );
    double pos=hit->GetLocalHitPos();
    double aa=hit->GetHit()->GetTiltAngle()*Deg2Rad;
    double cpx=xo+uo*z, cpy=yo+vo*z;
    double c=cpx*cos(aa)+cpy*sin(aa);
    chi += ((pos-c)*(pos-c)/(dd*dd)); ++nh;
    if( flag==3 ) hit->SetCalLPos(c);
    */
  }
  if(nh>5) chi /=double(nh-5);
  return chi;
}



#ifndef __E42_CHAMBER_H__
#define __E42_CHAMBER_H__

#include "chamber_const.hh"
#include "common.hh"

typedef struct {
  double x0,y0,z0;
  double shiftx,shifty,shiftz;
  double sinth,costh;
  double sinphi1,cosphi1;
  double sinphi2,cosphi2;
} Position_Par;

class E42_chamber {
 private:
  Position_Par BPCpos[BPCMAX];
  Position_Par BDCpos[BDCMAX];
  Position_Par PCpos[PCMAX];
  Position_Par DCpos[DCMAX];

  double BPCreso[BPCPlaneMax];
  double BDCreso[BDCPlaneMax];
  double PCreso[PCPlaneMax];
  double DCreso[DCPlaneMax];
  
  int WireMax[4];

  double PCPlaneEdge1[PCPlaneMax][XYZ];
  double PCPlaneEdge2[PCPlaneMax][XYZ];
  double PCPlaneCenter[PCPlaneMax][XYZ];
  
  double DCPlaneEdge1[DCPlaneMax][XYZ];
  double DCPlaneEdge2[DCPlaneMax][XYZ];
  double DCPlaneCenter[DCPlaneMax][XYZ];

 public:
  E42_chamber(void);
  // Read Position Parameter
  void ReadParPos(void);
  void DumpParPos(void);
  // Read Resolution Parameter
  int ReadParReso(void);
  void DumpParReso(void);
  double get_BPCreso(int plane);
  double get_PCreso(int plane);
  double get_DCreso(int plane);
  // BPC
  int defbpc(int id,double *orgn,double *shift, double *angle);
  int calbpx(double *hitpos,int plane,double wire, int *stat);
  int calbpy(double *hitpos,int plane,double wire, int *stat);
  // BDC
  int defbdc(int id,double *orgn,double *shift, double *angle);
  // PC
  int defpc(int id,double *orgn,double *shift, double *angle);
  int calpc(double *Pos, int iplane, double wire);
  int calWirePosPC(double *intpos, int dir, double cent, double z, double span,
		   int maxwire, double wire);
  void sgpc(int chamid, double *pos, double *intpos);
  void PCPlaneEdgedef(void);
  void DumpPCPlaneEdge(void);
  double get_PCPlaneEdge1(int plane, int coord);
  double get_PCPlaneEdge2(int plane, int coord);
  double get_PCPlaneCenter(int plane, int coord);
  void get_PCPlaneCenter(int plane, double *pos);
  // DC
  int defdc(int id,double *orgn,double *shift, double *angle);
  //DC1
  int caldc1(double *hitpos, int plane, int wire, int *stat);
  void caldc1u(double *hitpos, double *intpos, int plane, int wire, int *istat);
  void sudc1(double *intpos, int plane, int wire, int *stat);
  void sxdc1(double *intpos, int plane, int wire, int *stat);
  void sydc1(double *intpos, int plane, int wire, int *stat);
  void sgdc1(double *hitpos, double *intpos);
  //DC2
  int caldc2(double *hitpos, int plane, int wire, int *stat);
  void sxdc2(double *intpos, int plane, int wire, int *stat);
  void sydc2(double *intpos, int plane, int wire, int *stat);
  void sgdc2(double *hitpos, double *intpos);
  //DC3
  int caldc3(double *hitpos, int plane, int wire, int *stat);
  void sxdc3(double *intpos, int plane, int wire, int *stat);
  void sydc3(double *intpos, int plane, int wire, int *stat);
  void sgdc3(double *hitpos, double *intpos);

  void DCPlaneEdgedef(void);
  void DumpDCPlaneEdge(void);
  void PlaneEdgeDC1(void);
  void PlaneEdgeDC2(void);
  void PlaneEdgeDC3(void);
  double get_DCPlaneEdge1(int plane, int coord);
  double get_DCPlaneEdge2(int plane, int coord);
  double get_DCPlaneCenter(int plane, int coord);
  void get_DCPlaneCenter(int plane, double *pos);

  int getPlaneID(char *planenam, int *plane);
  void get_PlaneName(char *chamnam, int plane, char *planename);
};
  
#endif

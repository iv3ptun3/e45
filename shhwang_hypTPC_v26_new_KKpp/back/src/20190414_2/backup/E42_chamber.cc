#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "E42_chamber.hh"

E42_chamber::E42_chamber(void)
{
  int i,j;

  printf("ChamRESOinit::initialising Chamber Resolution  ...");
  for (i = 0; i < BPCPlaneMax; i++)
    BPCreso[i] = 9999.0;
  for (i = 0; i < BDCPlaneMax; i++)
    BDCreso[i] = 9999.0;
  for (i = 0; i < PCPlaneMax; i++)
    PCreso[i] = 9999.0;
  for (i = 0; i < DCPlaneMax; i++)
    DCreso[i] = 9999.0;

  for (i=0; i<PCPlaneMax; i++){
    for (j=0; j<3; j++) {
      PCPlaneEdge1[i][j] = 0.0;
      PCPlaneEdge2[i][j] = 0.0;
      PCPlaneCenter[i][j] = 0.0;
    }
  }
  for (i=0; i<DCPlaneMax; i++){
    for (j=0; j<3; j++) {
      DCPlaneEdge1[i][j] = 0.0;
      DCPlaneEdge2[i][j] = 0.0;
      DCPlaneCenter[i][j] = 0.0;
    }
  }

  WireMax[VHx] = 32;
  WireMax[VHy] = 18;
  WireMax[CHx] = 24;
  WireMax[FTOFx] = 24;

  ReadParPos();
  ReadParReso();

  return; 
}

/////////Read Position Parameter//////////
void E42_chamber::ReadParPos(void)
{
  int i, id, tokennum;
  double orgn[3],shift[3],angle[3];
  char dummy[256],cham[5],name[10];
  char fnam[30];
  FILE *fchampos;

  sprintf(fnam, "champosE42.d");
  if ((fchampos = fopen(fnam, "r")) == NULL) {
    fprintf(stderr, "ChamPOSinit::Cannot open champos file: %s\n", fnam);
    exit(-1);
  }

  printf("ChamPOSinit::initialising Chamber Position ...");
  for (i=0;;) {
    printf(".");
    if (fgets(dummy,sizeof(dummy),fchampos)==NULL)
      break;
    if (dummy[0] == '#' || dummy[0] == '\n')
      continue;
    /* orginal position */
    tokennum = sscanf(dummy, "%s %d %s %lf %lf %lf",
		      cham, &id, name, &orgn[0], &orgn[1], &orgn[2]);
    printf("orgn=(%f,%f,%f)\n", orgn[0], orgn[1], orgn[2]);
    if (tokennum != 6) {
      fprintf(stderr, "ChamPOSinit::invarid data structure\n");
      fprintf(stderr, "%s\n", dummy);
      exit(-1);
    }
    /* position shift */
    fgets(dummy, sizeof(dummy), fchampos);
    sscanf(dummy, "%s %d %s %lf %lf %lf",
	   cham, &id, name, &shift[0], &shift[1], &shift[2]);
    printf("shift=(%f,%f,%f)\n", shift[0], shift[1], shift[2]);
    if (tokennum != 6) {
      fprintf(stderr, "ChamPOSinit::invarid data structure\n");
      fprintf(stderr, "%s\n", dummy);
      exit(-1);
    }
    /* rotation */
    fgets(dummy, sizeof(dummy), fchampos);
    sscanf(dummy, "%s %d %s %lf %lf %lf",
	   cham, &id, name, &angle[0], &angle[1], &angle[2]);
    printf("angle=(%f,%f,%f)\n", angle[0], angle[1], angle[2]);
    if (tokennum != 6) {
      fprintf(stderr, "ChamPOSinit::invalid data structure\n");
      fprintf(stderr, "%s\n", dummy);
      exit(-1);
    }

    if ( strcmp(cham, "BPC") == 0 ) {
      if (id < 1 || id > 5) {
	fprintf(stderr, "ChamPOSinit::no such BPC:id=%d\n", id);
	exit(-1);
      }

      defbpc(id, orgn, shift, angle);
    } else if ( strcmp(cham, "BDC") == 0 ) {
      if (id < 1 || id > 3) {
	fprintf(stderr, "ChamPOSinit::no such BDC:id=%d\n", id);
	exit(-1);
      }
      defbdc(id-1, orgn, shift, angle);
    } else if ((strcmp(cham, "VH") && strcmp(cham, "CH")
		&& strcmp(cham, "FTOF")) == 0) {
      if (id != 1) {
	fprintf(stderr, "ChamPOSinit::no such %s:id=%d\n", cham, id);
	exit(-1);
      }
      if (strcmp(cham, "VH") == 0)
	id = VH;
      if (strcmp(cham, "CH") == 0)
	id = CHPC;
      if (strcmp(cham, "FTOF") == 0)
	id = FTOFPC;

      defpc(id, orgn, shift, angle);
    
    } else if (strcmp(cham, "DC") == 0) {
      defdc(id - 1, orgn, shift, angle);
    } else {
      fprintf(stderr, "ChamPOSinit::no such Chamber type=%s\n", cham);
      exit(-1);
    }
  }
  /* Define Plane Edge */
  PCPlaneEdgedef();
  DumpPCPlaneEdge();
  DCPlaneEdgedef();
  DumpDCPlaneEdge();

  printf("done !\n");
  fclose(fchampos);

  return ;
}

void E42_chamber::DumpParPos(void)
{
  int i;
  printf("E42_chamber::DumpPosPar BPC\n");
  for (i=0; i<BPCMAX; i++) {
    printf("--BPC%d---\n",i+1);
    printf("x0:%f y0:%f z0:%f\n",BPCpos[i].x0,BPCpos[i].y0,BPCpos[i].z0);
    printf("shiftx:%f shifty:%f shiftz:%f\n",
	   BPCpos[i].shiftx,BPCpos[i].shifty,BPCpos[i].shiftz);
    printf("sinth:%f costh:%f\n",BPCpos[i].sinth,BPCpos[i].costh);
    printf("sinphi1:%f cosphi1:%f\n",BPCpos[i].sinphi1,BPCpos[i].cosphi1);
    printf("sinphi2:%f cosphi2:%f\n",BPCpos[i].sinphi2,BPCpos[i].cosphi2);
  }
  for (i=0; i<BDCMAX; i++) {
    printf("--BDC%d---\n",i+1);
    printf("x0:%f y0:%f z0:%f\n",BDCpos[i].x0,BDCpos[i].y0,BDCpos[i].z0);
    printf("shiftx:%f shifty:%f shiftz:%f\n",
	   BDCpos[i].shiftx,BDCpos[i].shifty,BDCpos[i].shiftz);
    printf("sinth:%f costh:%f\n",BDCpos[i].sinth,BDCpos[i].costh);
    printf("sinphi1:%f cosphi1:%f\n",BDCpos[i].sinphi1,BDCpos[i].cosphi1);
    printf("sinphi2:%f cosphi2:%f\n",BDCpos[i].sinphi2,BDCpos[i].cosphi2);
  }
  for (i=0; i<PCMAX; i++) {
    printf("--PC%d---\n",i+1);
    printf("x0:%f y0:%f z0:%f\n",PCpos[i].x0,PCpos[i].y0,PCpos[i].z0);
    printf("shiftx:%f shifty:%f shiftz:%f\n",
	   PCpos[i].shiftx,PCpos[i].shifty,PCpos[i].shiftz);
    printf("sinth:%f costh:%f\n",PCpos[i].sinth,PCpos[i].costh);
    printf("sinphi1:%f cosphi1:%f\n",PCpos[i].sinphi1,PCpos[i].cosphi1);
    printf("sinphi2:%f cosphi2:%f\n",PCpos[i].sinphi2,PCpos[i].cosphi2);
  }
  for (i=0; i<DCMAX; i++) {
    printf("--DC%d---\n",i+1);
    printf("x0:%f y0:%f z0:%f\n",DCpos[i].x0,DCpos[i].y0,DCpos[i].z0);
    printf("shiftx:%f shifty:%f shiftz:%f\n",
	   DCpos[i].shiftx,DCpos[i].shifty,DCpos[i].shiftz);
    printf("sinth:%f costh:%f\n",DCpos[i].sinth,DCpos[i].costh);
    printf("sinphi1:%f cosphi1:%f\n",DCpos[i].sinphi1,DCpos[i].cosphi1);
    printf("sinphi2:%f cosphi2:%f\n",DCpos[i].sinphi2,DCpos[i].cosphi2);
  }
  return;
}


////////E42_chamberReso////////////

int E42_chamber::ReadParReso(void)
{
  int i, tokennum;
  FILE *fchamreso;
  char fnam[30], cham[10], dummy[256];
  double reso;

  sprintf(fnam, "chamreso.d");
  if ((fchamreso = fopen(fnam, "r")) == NULL) {
    fprintf(stderr, "ChamRESOinit::Cannot open chamreso file: %s\n", fnam);
    exit(-1);
  }

  while (fgets(dummy, sizeof(dummy), fchamreso) != NULL) {
    printf(".");
    if (dummy[0] == '#' || dummy[0] == '\n')
      continue;
    tokennum = sscanf(dummy, "%s %lf", cham, &reso);
    if (tokennum != 2) {
      fprintf(stderr, "ChamRESOinit::invarid data structure\n");
      fprintf(stderr, "%s\n", dummy);
      exit(-1);
    }
    if (strcmp(cham, "BPC1X") == 0)
      BPCreso[BPC1X] = reso;
    else if (strcmp(cham, "BPC1Y") == 0)
      BPCreso[BPC1Y] = reso;
    else if (strcmp(cham, "BPC2X") == 0)
      BPCreso[BPC2X] = reso;
    else if (strcmp(cham, "BPC3X") == 0)
      BPCreso[BPC3X] = reso;
    else if (strcmp(cham, "BPC3Y") == 0)
      BPCreso[BPC3Y] = reso;
    else if (strcmp(cham, "BPC4X") == 0)
      BPCreso[BPC4X] = reso;
    else if (strcmp(cham, "BPC4Y") == 0)
      BPCreso[BPC4Y] = reso;
    else if (strcmp(cham, "BPC5X") == 0)
      BPCreso[BPC5X] = reso;
    else if (strcmp(cham, "BPC5Y") == 0)
      BPCreso[BPC5Y] = reso;
    else if (strcmp(cham, "BDC1X") == 0)
      BDCreso[BDC1X] = reso;
    else if (strcmp(cham, "BDC1Xp") == 0)
      BDCreso[BDC1Xp] = reso;
    else if (strcmp(cham, "BDC2X") == 0)
      BDCreso[BDC2X] = reso;
    else if (strcmp(cham, "BDC2Xp") == 0)
      BDCreso[BDC2Xp] = reso;
    else if (strcmp(cham, "BDC3X") == 0)
      BDCreso[BDC3X] = reso;
    else if (strcmp(cham, "BDC3Xp") == 0)
      BDCreso[BDC3Xp] = reso;
    else if (strcmp(cham, "VHx") == 0)
      PCreso[VHx] = reso;
    else if (strcmp(cham, "VHy") == 0)
      PCreso[VHy] = reso;
    else if (strcmp(cham, "CHx") == 0)
      PCreso[CHx] = reso;
    else if (strcmp(cham, "DC1X") == 0)
      DCreso[DC1X] = reso;
    else if (strcmp(cham, "DC1Xp") == 0)
      DCreso[DC1Xp] = reso;
    else if (strcmp(cham, "DC1Y") == 0)
      DCreso[DC1Y] = reso;
    else if (strcmp(cham, "DC1U") == 0)
      DCreso[DC1U] = reso;
    else if (strcmp(cham, "DC2X") == 0)
      DCreso[DC2X] = reso;
    else if (strcmp(cham, "DC2Xp") == 0)
      DCreso[DC2Xp] = reso;
    else if (strcmp(cham, "DC2Y") == 0)
      DCreso[DC2Y] = reso;
    else if (strcmp(cham, "DC2Yp") == 0)
      DCreso[DC2Yp] = reso;
    else if (strcmp(cham, "DC3X") == 0)
      DCreso[DC3X] = reso;
    else if (strcmp(cham, "DC3Xp") == 0)
      DCreso[DC3Xp] = reso;
    else if (strcmp(cham, "DC3Y") == 0)
      DCreso[DC3Y] = reso;
    else if (strcmp(cham, "DC3Yp") == 0)
      DCreso[DC3Yp] = reso;
  }
  printf("\n");

  return 0;
}

void E42_chamber::DumpParReso(void)
{
  int plane;
  char planename[30];

  for (plane=0; plane<BPCPlaneMax; plane++ ) {
    get_PlaneName("BPC", plane, planename);
    printf("BPCreso[%s]:%f\n",planename,BPCreso[plane]);
  }
  for (plane=0; plane<PCPlaneMax; plane++ ) {
    get_PlaneName("PC", plane, planename);
    printf("PCreso[%s]:%f\n",planename,PCreso[plane]);
  }
  for (plane=0; plane<DCPlaneMax; plane++ ) {
    get_PlaneName("DC", plane, planename);
    printf("DCreso[%s]:%f\n",planename,DCreso[plane]);
  }
  return;
}



double E42_chamber::get_DCreso(int plane)
{
  return DCreso[plane];
}


double E42_chamber::get_PCreso(int plane)
{
  return PCreso[plane];
}


double E42_chamber::get_BPCreso(int plane)
{
  return BPCreso[plane];
}

///////////////// For Beam Propotional Chamber ////////////////////////////

int E42_chamber::defbpc(int id,double *orgn,double *shift, double *angle)
{
  int chamID;   /* BPC1 ---> 0,,,,BPC5 ---->4 */

  chamID = id-1;

  if (chamID >= 0 && chamID < BPCMAX) {
    BPCpos[chamID].x0 = orgn[0];
    BPCpos[chamID].y0 = orgn[1];
    BPCpos[chamID].z0 = orgn[2];
    
    BPCpos[chamID].shiftx = shift[0];
    BPCpos[chamID].shifty = shift[1];
    BPCpos[chamID].shiftz = shift[2];
    
    BPCpos[chamID].sinth = sin(angle[0]*3.141592654/180.0);
    BPCpos[chamID].costh = cos(angle[0]*3.141592654/180.0);
    
    BPCpos[chamID].sinphi1 = sin(angle[1]*3.141592654/180.0);
    BPCpos[chamID].cosphi1 = cos(angle[1]*3.141592654/180.0);
    
    BPCpos[chamID].sinphi2 = sin(angle[2]*3.141592654/180.0);
    BPCpos[chamID].cosphi2 = cos(angle[2]*3.141592654/180.0);
 
    return 0;
  } else {
    fprintf(stderr,"E42_chamber::defpc chamID is abnormal:%d\n",chamID);
    exit(1);
  }
}


int E42_chamber::calbpx(double *hitpos,int plane,double wire, int *stat)
{
  int minw[5]={1, 1, 17, 17, 33};
  int maxw[5]={128, 128, 128, 128, 128};
  double center[5]={72.5, 72.5, 72.5, 72.5, 72.5};
  double x,xx,zz;

  if (plane < 0 || plane >= 5) {
    *stat = -1;
    fprintf(stderr,"E42_chamber::calbpx invalid plane:%d \n",plane);
    exit(1);
  }
  if (wire < minw[plane] || wire > maxw[plane]) {
    *stat = -1;
    fprintf(stderr,"E42_chamber::calbpx invalid wire:%f",wire);
    exit(1);
  }
  if((plane == 2) || (plane == 3) || (plane == 4)) {
    x = center[plane] - wire;
  } else {
    x = wire - center[plane];
  }

  xx = x * BPCpos[plane].costh + BPCpos[plane] .shiftx;
  zz = -x * BPCpos[plane].sinth + BPCpos[plane] .shiftz;

  hitpos[XCOORD] = xx + BPCpos[plane].x0;
  hitpos[YCOORD] = 0.0;            // dummy
  hitpos[ZCOORD] = zz + BPCpos[plane].z0;

  *stat = 1;
  return 0;
}


int E42_chamber::calbpy(double *hitpos,int plane,double wire, int *stat)
{
  int minw[5]={1, 1, 1, 1, 1};
  int maxw[5]={64, 64, 64, 64, 64};
  double center[5]={32.5, 32.5, 32.5, 32.5, 32.5};
  double y,yy,zz;
  double intposz;
  double x,xx;

  if (plane < 0 || plane >= 5) {
    *stat = -1;
    fprintf(stderr,"E42_chamber::calbpy invalid plane:%d \n",plane);
    exit(1);
  }
  if (wire < minw[plane] || wire > maxw[plane]) {
    *stat = -1;
    fprintf(stderr,"E42_chamber::calbpy invalid wire:%f",wire);
    exit(1);
  }

  if (plane == 2) {
    intposz = 16.0;
  } else {
    intposz = -16.0;
  }

  if (plane == 3 || plane == 4) {
    y = wire - center[plane];
  } else {
    y = center[plane] - wire;
  }

  x = 0.0;
  xx = x * BPCpos[plane].costh + BPCpos[plane] .shiftx;
  yy = y * BPCpos[plane].cosphi2 + BPCpos[plane] .shifty;
  zz = -y * BPCpos[plane].sinphi2 + intposz + BPCpos[plane] .shiftz;

  hitpos[XCOORD] = xx + BPCpos[plane].x0;            // dummy
  hitpos[YCOORD] = yy + BPCpos[plane].y0;
  hitpos[ZCOORD] = zz + BPCpos[plane].z0;

  *stat = 1;
  return 0;
}


///////////////// For Beam Drift Chamber ////////////////////////////

int E42_chamber::defbdc(int id,double *orgn,double *shift, double *angle)
{
  int chamID;   /* BDC1 ---> 0,,,,BDC3 ---->2 */
               
  chamID = id;  /* これは*idが0から渡される */  

  if (chamID >= 0 && chamID < BDCMAX) {
    BDCpos[chamID].x0 = orgn[0];
    BDCpos[chamID].y0 = orgn[1];
    BDCpos[chamID].z0 = orgn[2];
    
    BDCpos[chamID].shiftx = shift[0];
    BDCpos[chamID].shifty = shift[1];
    BDCpos[chamID].shiftz = shift[2];
    
    BDCpos[chamID].sinth = sin(angle[0]*3.141592654/180.0);
    BDCpos[chamID].costh = cos(angle[0]*3.141592654/180.0);
    
    BDCpos[chamID].sinphi1 = sin(angle[1]*3.141592654/180.0);
    BDCpos[chamID].cosphi1 = cos(angle[1]*3.141592654/180.0);
    
    BDCpos[chamID].sinphi2 = sin(angle[2]*3.141592654/180.0);
    BDCpos[chamID].cosphi2 = cos(angle[2]*3.141592654/180.0);
 
    return 0;
  } else {
    fprintf(stderr,"E42_chamber::defbdc chamID is abnormal:%d\n",chamID);
    exit(1);
  }
}



///////////////// For Propotional Chamber ////////////////////////////

int E42_chamber::defpc(int id,double *orgn,double *shift, double *angle)
{
  int chamID;  /* これは*idが0から渡される */

  chamID = id;
  
  if (chamID >= 0 && chamID < PCMAX) {
    PCpos[chamID].x0 = orgn[0];
    PCpos[chamID].y0 = orgn[1];
    PCpos[chamID].z0 = orgn[2];
    
    PCpos[chamID].shiftx = shift[0];
    PCpos[chamID].shifty = shift[1];
    PCpos[chamID].shiftz = shift[2];
    
    PCpos[chamID].sinth = sin(angle[0]*3.141592654/180.0);
    PCpos[chamID].costh = cos(angle[0]*3.141592654/180.0);
    
    PCpos[chamID].sinphi1 = sin(angle[1]*3.141592654/180.0);
    PCpos[chamID].cosphi1 = cos(angle[1]*3.141592654/180.0);
    
    PCpos[chamID].sinphi2 = sin(angle[2]*3.141592654/180.0);
    PCpos[chamID].cosphi2 = cos(angle[2]*3.141592654/180.0);
    
    return 0;
  } else {
    fprintf(stderr,"E42_chamber::defpc chamID is abnormal:%d\n",chamID);
    exit(1);
  }
}

int E42_chamber::calpc(double *Pos, int iplane, double wire)
{
  /********************************************************/
  /*  PC wire or cluster center position caliicuration    */
  /*      Pos[3]    :  calicurated position               */
  /*      intpos[3] :  internal position                  */
  /*      wire     :   wire # (from 1)                    */
  /********************************************************/
  int istat;
  double intpos[3];
  int chamid=-1;

  switch (iplane) {
  case VHx:
    istat = calWirePosPC(intpos, XCOORD, 16.5, 20.0, 4.4, 32, wire);
    if ((int)wire%2 == 0)
      intpos[ZCOORD] -= 1.0;
    else
      intpos[ZCOORD] += 1.0;
    chamid = VH;
    break;
  case VHy:
    istat = calWirePosPC(intpos, YCOORD, 9.5, 14.0, 6.0, 18, wire);
    if ((int)wire%2 == 0)
      intpos[ZCOORD] -= 1.0;
    else
      intpos[ZCOORD] += 1.0;
    chamid = VH;
    break;
  case CHx:
    istat = calWirePosPC(intpos, XCOORD, 12.5, 0.0, 17.5, 24, wire);
    if ((int)wire%2 == 0)
      intpos[ZCOORD] -= 1.0;
    else
      intpos[ZCOORD] += 1.0;
    chamid = CHPC;
    break;
  case FTOFx:
    istat = calWirePosPC(intpos, XCOORD, 12.5, -121.5, 75.0, 24, wire);
    //istat = calWirePosPC(intpos, XCOORD, 12.5, 0.0, 75.0, 24, wire);
    /*
    if (((int) wire) % 2 == 0)
      intpos[ZCOORD] -= 35.0;
    */
    if (((int) wire) % 2 == 0)
      intpos[ZCOORD] -= 17.5;
    else
      intpos[ZCOORD] += 17.5;
    chamid = FTOFPC;
    break;
  default:
    fprintf(stderr, "calpc::Invarid plane id;iplane=%d\n", iplane);
    istat = -1;
    break;
  }

  if (istat >= 0)
    sgpc(chamid, Pos, intpos);
  else
    printf("PC Invalid wire id chamid=%d, wire=%lf\n", chamid, wire);

  return (istat);

}

int E42_chamber::calWirePosPC(double *intpos, int dir, double cent,
			      double z, double span,int maxwire, double wire)
{
  int istat;

  if (wire >= 1.0 && wire <= (double)maxwire) {
    intpos[XCOORD] = intpos[YCOORD] = 0;
    intpos[dir] = (cent - wire) * span;
    intpos[ZCOORD] = z;
    istat = 1;
  } else
    istat = -1;

  return (istat);
}

void E42_chamber::sgpc(int chamid, double *pos, double *intpos)
{
  double intpos1[3];

  if (0 <= chamid && chamid < PCMAX) {
    intpos1[XCOORD] = intpos[XCOORD] + PCpos[chamid].shiftx;
    intpos1[YCOORD] = intpos[YCOORD] * PCpos[chamid].cosphi2
      + intpos[ZCOORD] * PCpos[chamid].sinphi2 + PCpos[chamid].shifty;
    intpos1[ZCOORD] = -intpos[YCOORD] * PCpos[chamid].sinphi2
      + intpos[ZCOORD] * PCpos[chamid].cosphi2 + PCpos[chamid].shiftz;
    pos[XCOORD] =
      intpos1[XCOORD] * PCpos[chamid].costh
      + intpos1[ZCOORD] * PCpos[chamid].sinth 
      + PCpos[chamid].x0;
    pos[ZCOORD] =
      - intpos1[XCOORD] * PCpos[chamid].sinth
      + intpos1[ZCOORD] * PCpos[chamid].costh
      + PCpos[chamid].z0;
    
    pos[YCOORD] = intpos1[YCOORD] + PCpos[chamid].y0;
    
    return;
  } else {
    fprintf(stderr,"E42_chamber::sgpc chamID is abnormal:%d\n",chamid);
    exit(1);
  }
}

void E42_chamber::PCPlaneEdgedef(void)
{
  int plane,chamid=-1;
  double wire;

  double intpos[XYZ];
  double hitpos[XYZ];

  for (plane=0; plane<PCPlaneMax; plane++) {
    wire = 1.0;
    calpc(hitpos, plane, wire);
    PCPlaneEdge1[plane][XCOORD] = hitpos[XCOORD];
    PCPlaneEdge1[plane][YCOORD] = hitpos[YCOORD];
    PCPlaneEdge1[plane][ZCOORD] = hitpos[ZCOORD];
    
    switch (plane) {
    case VHx:
      intpos[XCOORD] = 0.0;
      intpos[YCOORD] = 0.0;
      intpos[ZCOORD] = 20.0;
      chamid = VH;
      break;
    case VHy:
      intpos[XCOORD] = 0.0;
      intpos[YCOORD] = 0.0;
      intpos[ZCOORD] = 14.0;
      chamid = VH;
      break;
    case CHx:
      intpos[XCOORD] = 0.0;
      intpos[YCOORD] = 0.0;
      intpos[ZCOORD] = 0.0;
      chamid = CHPC;
      break;
    case FTOFx:
      intpos[XCOORD] = 0.0;
      intpos[YCOORD] = 0.0;
      intpos[ZCOORD] = -121.5;
      chamid = FTOFPC;
      break;
    default:
      fprintf(stderr, "PCPlaneEdgedef::Invarid plane id;plane=%d\n", plane);
      break;
    }
    sgpc(chamid, hitpos, intpos);
    PCPlaneCenter[plane][XCOORD] = hitpos[XCOORD];
    PCPlaneCenter[plane][YCOORD] = hitpos[YCOORD];
    PCPlaneCenter[plane][ZCOORD] = hitpos[ZCOORD];

    wire = (double) WireMax[plane];
    calpc(hitpos, plane, wire);
    PCPlaneEdge2[plane][XCOORD] = hitpos[XCOORD];
    PCPlaneEdge2[plane][YCOORD] = hitpos[YCOORD];
    PCPlaneEdge2[plane][ZCOORD] = hitpos[ZCOORD];
  }
  return;
}

void E42_chamber::DumpPCPlaneEdge(void)
{
  int i;
  char planename[30];

  printf("---DumpPCPlaneEdge---\n");
  for (i=0; i<PCPlaneMax; i++) {
    get_PlaneName("PC",i,planename);
    printf("PCPlaneEdge1[%s][XYZ]=(%f, %f, %f)\n", planename,
	   PCPlaneEdge1[i][0], PCPlaneEdge1[i][1], PCPlaneEdge1[i][2]);
    printf("PCPlaneEdge2[%s][XYZ]=(%f, %f, %f)\n", planename,
	   PCPlaneEdge2[i][0],PCPlaneEdge2[i][1], PCPlaneEdge2[i][2]);
    printf("PCPlaneCenter[%s][XYZ]=(%f, %f, %f)\n", planename,
	   PCPlaneCenter[i][0],PCPlaneCenter[i][1], PCPlaneCenter[i][2]);
  }
  
  return;
}

double E42_chamber::get_PCPlaneEdge1(int plane, int coord)
{
  return PCPlaneEdge1[plane][coord];
}

double E42_chamber::get_PCPlaneEdge2(int plane, int coord)
{
  return PCPlaneEdge2[plane][coord];
}

double E42_chamber::get_PCPlaneCenter(int plane, int coord)
{
  return PCPlaneCenter[plane][coord];
}

void E42_chamber::get_PCPlaneCenter(int plane, double *pos)
{
  pos[XCOORD] = PCPlaneCenter[plane][XCOORD];
  pos[YCOORD] = PCPlaneCenter[plane][YCOORD];
  pos[ZCOORD] = PCPlaneCenter[plane][ZCOORD];
  return;
}


///////////////// For Drift Chamber ////////////////////////////
 
int E42_chamber::defdc(int id,double *orgn,double *shift, double *angle)
{
  int chamID;  /* これも*idが0から渡される */

  chamID = id;
  
  if (chamID >= 0 && chamID < DCMAX) {
    DCpos[chamID].x0 = orgn[0];
    DCpos[chamID].y0 = orgn[1];
    DCpos[chamID].z0 = orgn[2];
    
    DCpos[chamID].shiftx = shift[0];
    DCpos[chamID].shifty = shift[1];
    DCpos[chamID].shiftz = shift[2];
    
    DCpos[chamID].sinth = sin(angle[0]*3.141592654/180.0);
    DCpos[chamID].costh = cos(angle[0]*3.141592654/180.0);
    
    DCpos[chamID].sinphi1 = sin(angle[1]*3.141592654/180.0);
    DCpos[chamID].cosphi1 = cos(angle[1]*3.141592654/180.0);
    
    DCpos[chamID].sinphi2 = sin(angle[2]*3.141592654/180.0);
    DCpos[chamID].cosphi2 = cos(angle[2]*3.141592654/180.0);
    
    return 0;
  } else {
    fprintf(stderr,"E42_chamber::defdc chamID is abnormal:%d\n",chamID);
    exit(1);
  }
}

int E42_chamber::caldc1(double *hitpos, int plane, int wire, int *stat)
{
  double intpos[3];
  
  if (plane == 1 || plane == 2) {
    sxdc1(intpos, plane, wire, stat);
    if (*stat != 1) return -1;
  } else if (plane == 3) {
    sydc1(intpos, plane, wire, stat);
    if (*stat != 1) return -1;
  } else if (plane == 4) {
    printf("DC1- U plane not available caldc1!\n");
  } else {
    printf("E42_chamber::caldc1 invalid plane:%d\n",plane);
    exit(1);
  }
  
  sgdc1(hitpos, intpos);

  return 0;
}  

void E42_chamber::caldc1u(double *hitpos, double *intpos, int plane, int wire, int *istat)
{
  if (plane == 4) {
    sudc1(intpos, plane, wire, istat);
    if (*istat != 1)
      return;
  } else {
    *istat = -999;
    return;
  }
  
  sgdc1(hitpos, intpos);

  return;
}

void E42_chamber::sxdc1(double *intpos, int plane, int wire, int *stat)
{
  /*****************************************************/
  /* calculation of the position in internal geometry  */
  /* with some parameters for DC1-X,X' plane           */
  /*****************************************************/

  /*	------------------------------------------     
	internal geometry for DC1  
       
	INTPOS(1) == X direction in the DC1-plane
	INTPOS(2) == Y direction == no usage (dummy)
	INTPOS(3) == Z direction in the DC1-plane

	  !!! caution  * the unit is now "mm" 

       Z-direction (in the DC1 plane)......
	  1) 0 point is X-plane
	  2) X'-plane will be -4.8 mm
	  3) Y -plane will be X' plane -19.6 mm 
	  4) U -plane will be y  plane -4.8  mm 

 ===================================================================
       X-direction (in the DC1 plane -- magnet plane )
         1) 0 point is X-plane 24.0th (24.5th for X') wire position 

 ===================================================================
 -------------------------------------------  */

  const int MXWIRE = 48;
  const double DX = 10.0, ZXX = -4.8;
  const double SHIFTXP = 0.0,SHIFTX = -0.32;

  intpos[1] = 0.0;
  if(plane == 1) {
    if(wire>=1 && wire <= MXWIRE){
      intpos[0] = -(double(wire)-24.0)*DX + SHIFTX;
      intpos[2] = 0.0;
      *stat  = 1;
    } else {
      *stat = -1;
      return;
    }
  } else if(plane == 2) {
    if(wire >= 1 && wire<=MXWIRE) {
      intpos[0] = (double(wire)-24.5)*DX + SHIFTXP;
      intpos[2] = ZXX;
      *stat = 1;
    } else {
      *stat = -1;
      return;
    }
  } else {
    *stat = -999;
    return;
  }
}

void E42_chamber::sydc1(double *intpos, int plane, int wire, int *stat)
{
  /******************************************************/
  /** Calculation of the Position in internal geometry **/
  /** with some parameters for DC1-Y plane             **/
  /******************************************************/

  /*	------------------------------------------
	internal geometry for DC1 - Y Y'
       
	INTPOS(1) == X direction (dummy)
	INTPOS(2) == Y direction in the DC1 plane
	INTPOS(3) == Z direction in the DC1 plane

	  !!! caution  * the unit is now "mm" 

       Z-direction (in the DC1 plane)......
	  1) 0 point is X-plane
	  2) X'-plane will be -4.8 mm
	  3) Y -plane will be X' plane -19.6 mm 
	  4) U -plane will be y  plane -4.8  mm 

       Y-direction 
         1) 0 point is Y-plane 16.5th wire position 
         
	 ------------------------------------------- */
  const double DY  = 10.0;
  const double YX  = -4.8 -19.6;
  //	     YX is the length from X-plane to Y-plane in Z-direction
  const int MXWIRE = 32;
  const double SHIFTY = -0.4;
  //        parameter (SHIFTY = -0.58)

  intpos[0] = 0.0;
  if(plane == 3) {
    if(wire >=1 && wire <= MXWIRE) {
      intpos[1] = (double(wire)-16.5)*DY + SHIFTY;
      intpos[2] = YX;
      *stat = 1;
    } else {
      *stat = -1;
      return;
    }
  } else {
    *stat = -999;
    return;
  }        
  return;
}
      

void E42_chamber::sudc1(double *intpos, int plane, int wire, int *stat)
{
  /*  
	** Calculation of the Position in internal geometry 
	** with some parameters for DC1-U plane
	***************************************************
	
	------------------------------------------
	  !!! caution  * the unit is now "mm" 
       Z-direction (in the DC1 plane)......
	  1) 0 point is X-plane
	  2) X'-plane will be -4.8 mm
	  3) Y -plane will be X' plane -19.6 mm 
	  4) U -plane will be y  plane -4.8  mm 

       Y-direction 
         1) 0 point is Y-plane 1st wire position 
         
       -------------------------------------------
  */

  double  a,b,X0;

  const double DU = 10.0;
  const double tan60 = 1.732050;
  const double sin30 = 0.50;
  const double cos30 = 0.866025;

  const double XU = -4.8 -19.6 -4.8;
  const int MXWIRE = 48;

  if(plane == 4) {
    if(wire >= 1 && wire <= MXWIRE){
      X0=((double)wire-24.5)*DU;
      a=tan60;
      b=X0*(sin30+tan60*cos30);
      intpos[1]=-a*intpos[0]+b-0.84;
      intpos[2]=XU;
      *stat = 1;
    } else {
      *stat = -1;
      return;
    } 
  } else {
    *stat = -999;
    return;
  }

  return;
}


void E42_chamber::sgdc1(double *hitpos, double *intpos)
{
  /**************************************************************/
  /** calcualation of position in LAB-system from internal pos **/
  /**************************************************************/

  double dumpos[3];

  dumpos[0] = intpos[0] + DCpos[DC1].shiftx;
  dumpos[1] = intpos[1]*DCpos[DC1].cosphi2+intpos[2]*DCpos[DC1].sinphi2+DCpos[DC1].shifty;
  dumpos[2] = -intpos[1]*DCpos[DC1].sinphi2+intpos[2]*DCpos[DC1].cosphi2+DCpos[DC1].shiftz;

  hitpos[0] = dumpos[0]*DCpos[DC1].costh+dumpos[2]*DCpos[DC1].sinth
    + DCpos[DC1].x0;
  hitpos[2] = -dumpos[0]*DCpos[DC1].sinth+dumpos[2]*DCpos[DC1].costh
    + DCpos[DC1].z0;
  hitpos[1] = dumpos[1] + DCpos[DC1].y0;

  return;
}


///////////////////// For DC2 //////////////////////////////////
int E42_chamber::caldc2(double *hitpos, int plane, int wire, int *stat)
{
  double intpos[3];
  
  if (plane == 1 || plane == 2) {
    sxdc2(intpos, plane, wire, stat);
    sgdc2(hitpos, intpos);
  } else if (plane == 3 || plane == 4) {
    sydc2(intpos, plane, wire, stat);
    sgdc2(hitpos, intpos);
  } else {
    printf("E42_chamber::caldc2 invalid plane:%d\n",plane);
    *stat = -1;
    exit(1);
  }
  return 0;
}  

void E42_chamber::sxdc2(double *intpos, int plane, int wire, int *stat)
{
  /*****************************************************/
  /* calculation of the position in internal geometry  */
  /* with some parameters for DC2-X,X' plane           */
  /*****************************************************/
  /*	------------------------------------------
	internal geometry for DC2  
       
	INTPOS(1) == X direction in the DC2-plane
	INTPOS(2) == Y direction == no usage (dummy)
	INTPOS(3) == Z direction in the DC2-plane

	  !!! caution  * the unit is now "mm" 

       Z-direction (in the DC2 plane)......
	  1) 0 point is X-plane
	  2) X'-plane will be about   7.5 mm
	  3) Y -plane will be about + 25 mm 
	  4) Y'-plane will be about   y + 7.5 mm 
ccc     I have to have design of DC2  2/26 92' Satoru Yamashita   cccc
ccc		But I think
        X' is 7.8mm
        Y  is 25.2mm
        Y' is 33mm
ccc    96/1/12   Y.Kondo

       X-direction (in the DC2 plane -- magnet plane )
         1) 0 point is X-plane 64.75th (64.25th for X') wire position 
         
	 ------------------------------------------- */

  const int MXWIRE = 128;
  const double DX = 9.0, ZXX = 7.8;

  intpos[1] = 0.0;
  if(plane == 1) {
    if(wire>=1 && wire <= MXWIRE){
      intpos[0] = (double(wire)-64.75)*DX;
      intpos[2] = 0.0;
      *stat  = 1;
    } else {
      *stat = -1;
      return;
    }
  } else if(plane == 2) {
    if(wire >= 1 && wire<=MXWIRE) {
      intpos[0] = (double(wire)-64.25)*DX;
      intpos[2] = ZXX;
      *stat = 1;
    } else {
      *stat = -1;
      return;
    }
  } else {
    *stat = -999;
    return;
  }
}

void E42_chamber::sydc2(double *intpos, int plane, int wire, int *stat)
{
  /******************************************************/
  /** Calculation of the Position in internal geometry **/
  /** with some parameters for DC2-Y plane             **/
  /******************************************************/
  /*	------------------------------------------
	internal geometry for DC2 - Y Y'
       
	INTPOS(1) == X direction (dummy)
	INTPOS(2) == Y direction in the DC2 plane
	INTPOS(3) == Z direction in the DC2 plane

	  !!! caution  * the unit is now "mm" 

       Z-direction (in the DC2 plane)......
	  1) 0 point is X-plane
	  2) X'-plane will be about   7.5 mm
	  3) Y -plane will be about + 25 mm 
	  4) Y'-plane will be about   y + 7.5 mm 
ccc     I have to have design of DC2  2/26 92' Satoru Yamashita   cccc
ccc		But I think
        X' is 7.8mm
	 Y  is 25.2mm
        Y' is 33mm
ccc    96/1/12   Y.Kondo

       Y-direction 
         1) 0 point is Y-plane 48.75th (48.25th for Y') wire position 
         
	 -------------------------------------------*/

  const double DY  = 9.0;
  const double YX  = 25.2;
  //	     YX is the length from X-plane to Y-plane in Z-direction
  const double ZYY  = 7.8;
  const int MXWIRE = 96;

  intpos[0] = 0.0;
  if(plane == 3) {
    if(wire >=1 && wire <= MXWIRE) {
      intpos[1] = (double(wire)-48.75)*DY;
      intpos[2] = YX;
      *stat = 1;
    } else {
      *stat = -1;
      return;
    }
  } else if(plane == 4) {
    if(wire >=1 && wire <= MXWIRE) {
      intpos[1] = (double(wire)-48.25)*DY;
      intpos[2] = YX + ZYY;
      *stat = 1;
    } else {
      *stat = -1;
      return;
    }
  } else {
    *stat = -999;
    return;
  }        
  return;
}
      
void E42_chamber::sgdc2(double *hitpos, double *intpos)
{
  /**************************************************************/
  /** calcualation of position in LAB-system from internal pos **/ 
  /** for DC2   Ver 2.0  Satoru Yamashita   92' 2/26           **/ 
  /** modified by Satoshi Mihara for E251  93' 6/23            **/ 
  /** using measured values and shifted along z-direction (+200.0mm) **/
  /** rewrittened by Y.Kondo for E289                          **/ 
  /**************************************************************/
  /** All geometry parameter is measured value                 **/
  /**************************************************************/

  double dumpos[3];

  dumpos[0] = intpos[0] + DCpos[DC2].shiftx;
  dumpos[1] = intpos[1]*DCpos[DC2].cosphi2+intpos[2]*DCpos[DC2].sinphi2+DCpos[DC2].shifty;
  dumpos[2] = -intpos[1]*DCpos[DC2].sinphi2+intpos[2]*DCpos[DC2].cosphi2+DCpos[DC2].shiftz;

  hitpos[0] = dumpos[0]*DCpos[DC2].costh+dumpos[2]*DCpos[DC2].sinth
    + DCpos[DC2].x0;
  hitpos[2] = -dumpos[0]*DCpos[DC2].sinth+dumpos[2]*DCpos[DC2].costh
    + DCpos[DC2].z0;
  hitpos[1] = dumpos[1] + DCpos[DC2].y0;

  return;
}

///////////////////// For DC3 //////////////////////////////////
int E42_chamber::caldc3(double *hitpos, int plane, int wire, int *stat)
{
  double intpos[3];
  
  if (plane == 1 || plane == 2) {
    sxdc3(intpos, plane, wire, stat);
    sgdc3(hitpos, intpos);
  } else if (plane == 3 || plane == 4) {
    sydc3(intpos, plane, wire, stat);
    sgdc3(hitpos, intpos);
  } else {
    printf("E42_chamber::caldc3 invalid plane:%d\n",plane);
    *stat = -1;
    exit(1);
  }
  return 0;
}  

void E42_chamber::sxdc3(double *intpos, int plane, int wire, int *stat)
{
  /*****************************************************/
  /* calculation of the position in internal geometry  */
  /* with some parameters for DC3-X,X' plane           */
  /*****************************************************/
  /*	------------------------------------------
	internal geometry for DC3  
       
	INTPOS(1) == X direction in the DC3-plane
	INTPOS(2) == Y direction == no usage (dummy)
	INTPOS(3) == Z direction in the DC3-plane

	  !!! caution  * the unit is now "mm" 

       Z-direction (in the DC3 plane)......
	  1) 0 point is X-plane
	  2) X'-plane will be 7+9+9+7 mm
	  3) Y -plane will be X'-plane + 7+9+9+6 mm (caution 6mm)
	  4) Y'-plane will be Y -plane + 6+9+9+6 mm  

       X-direction (in the DC3 plane -- similer to magnet plane )
         1) 0 point is X-plane 16.75 th wire position (center of X,X') 
 
	 ------------------------------------------- */
  const int MXWIRE = 32;
  const double DX = 56.0, ZXX = 7.0 + 9.0 + 9.0 + 7.0;

  intpos[1] = 0.0;
  if(plane == 1) {
    if(wire>=1 && wire <= MXWIRE){
      intpos[0] = (double(wire)-16.75)*DX;
      intpos[2] = 0.0;
      *stat  = 1;
    } else {
      *stat = -1;
      return;
    }
  } else if(plane == 2) {
    if(wire >= 1 && wire<=MXWIRE) {
      intpos[0] = (double(wire)-16.25)*DX;
      intpos[2] = ZXX;
      *stat = 1;
    } else {
      *stat = -1;
      return;
    }
  } else {
    *stat = -999;
    return;
  }
}

void E42_chamber::sydc3(double *intpos, int plane, int wire, int *stat)
{
  /******************************************************/
  /** Calculation of the Position in internal geometry **/
  /** with some parameters for DC3-Y plane             **/
  /******************************************************/
  /*	------------------------------------------
	internal geometry for DC3 - Y Y'
       
	INTPOS(1) == X direction (dummy)
	INTPOS(2) == Y direction in the DC3 plane
	INTPOS(3) == Z direction in the DC3 plane

	  !!! caution  * the unit is now "mm" 

       Z-direction (in the DC3 plane)......
	  1) 0 point is X-plane
	  2) X'-plane will be 7+9+9+7 mm
	  3) Y -plane will be X'-plane + 7+9+9+6 mm (caution 6mm)
	  4) Y'-plane will be Y -plane + 6+9+9+6 mm  

       Y-direction 
         1) 0 point is Y-plane 8.25th (8.75th for Y') wire position 
         
	 -------------------------------------------*/
  const double DY  = 60.0;
  const double YX  = 7.0+9.0+9.0+7.0+7.0+9.0+9.0+6.0;
  //	     YX is the length from X-plane to Y-plane in Z-direction
  const double ZYY  = 6.0 + 9.0 + 9.0 + 6.0;
  const int MXWIRE = 16;

  intpos[0] = 0.0;
  if(plane == 3) {
    if(wire >=1 && wire <= MXWIRE) {
      intpos[1] = (double(wire)-8.25)*DY;
      intpos[2] = YX;
      *stat = 1;
    } else {
      *stat = -1;
      return;
    }
  } else if(plane == 4) {
    if(wire >=1 && wire <= MXWIRE) {
      intpos[1] = (double(wire)-8.75)*DY;
      intpos[2] = YX + ZYY;
      *stat = 1;
    } else {
      *stat = -1;
      return;
    }
  } else {
    *stat = -999;
    return;
  }        
  return;
}
      
void E42_chamber::sgdc3(double *hitpos, double *intpos)
{
  /******************************************************************/
  /** calcualation of position in LAB-system from internal pos     **/
  /** for DC3  Ver 2.0   Satoru Yamashita                          **/
  /*------  1993 2/26 Satoru Yamashita new param by measured value **/
  /*------  1993 3/9  Satoru Yamashita shift parameter defined     **/
  /*  1993 10/4 Satoshi Mihara changed the orignz value for E251 (+200mm) **/
  /*  1996 1/12 Yasuhiro Kondo rewritten for E289                  **/
  /******************************************************************/
  /*	all geometry parameters is measured value .                **/
  /******************************************************************/

  double dumpos[3];

  dumpos[0] = intpos[0] + DCpos[DC3].shiftx;
  dumpos[1] = intpos[1]*DCpos[DC3].cosphi2+intpos[2]*DCpos[DC3].sinphi2+DCpos[DC3].shifty;
  dumpos[2] = -intpos[1]*DCpos[DC3].sinphi2+intpos[2]*DCpos[DC3].cosphi2+DCpos[DC3].shiftz;

  hitpos[0] = dumpos[0]*DCpos[DC3].costh+dumpos[2]*DCpos[DC3].sinth
    + DCpos[DC3].x0;
  hitpos[2] = -dumpos[0]*DCpos[DC3].sinth+dumpos[2]*DCpos[DC3].costh
    + DCpos[DC3].z0;
  hitpos[1] = dumpos[1] + DCpos[DC3].y0;

  return;
}

void E42_chamber::DCPlaneEdgedef(void)
{
/**************************************************************/
/*   --------|---------                                       */
/*  |        |         |        1 : PlaneEdge1                */
/* ----------C-----------       2 : PlaneEdge2                */
/*  |        |         |        C : PlaneCenter               */
/*   ------------------                                       */
/*  ^                  ^                                      */
/*  1                  2                                      */
/**************************************************************/

  PlaneEdgeDC1();
  PlaneEdgeDC2();
  PlaneEdgeDC3();

  return;
  
}

void E42_chamber::DumpDCPlaneEdge(void)
{
  int i;
  char planename[30];

  printf("---DumpDCPlaneEdge---\n");
  for (i=0; i<DCPlaneMax; i++) {
    get_PlaneName("DC",i,planename);
    printf("DCPlaneEdge1[%s][XYZ]=(%f, %f, %f)\n", planename,
	   DCPlaneEdge1[i][0],DCPlaneEdge1[i][1], DCPlaneEdge1[i][2]);
    printf("DCPlaneEdge2[%s][XYZ]=(%f, %f, %f)\n", planename,
	   DCPlaneEdge2[i][0], DCPlaneEdge2[i][1], DCPlaneEdge2[i][2]);
    printf("DCPlaneCenter[%s][XYZ]=(%f, %f, %f)\n", planename,
	   DCPlaneCenter[i][0],DCPlaneCenter[i][1], DCPlaneCenter[i][2]);
  }
  
  return;
}



void E42_chamber::PlaneEdgeDC1(void)
{
  int plane;
  int wire;
  int stat;

  double hitpos[3];
  double intpos[3];

  /* DC1X */
  plane = 1;
  wire = 1;
  caldc1(hitpos, plane, wire, &stat);
  DCPlaneEdge1[DC1X][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC1X][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0.0;
  intpos[YCOORD] = 0.0;
  intpos[ZCOORD] = 0.0;
  sgdc1(hitpos, intpos);
  DCPlaneCenter[DC1X][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC1X][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC1X][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 48;
  caldc1(hitpos, plane, wire, &stat);
  DCPlaneEdge2[DC1X][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC1X][ZCOORD] = (double) hitpos[ZCOORD];

  /* DC1Xp */
  plane = 2;
  wire = 48;
  caldc1(hitpos, plane, wire, &stat);
  DCPlaneEdge1[DC1Xp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC1Xp][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0.0;
  intpos[YCOORD] = 0.0;
  intpos[ZCOORD] = -4.8;
  sgdc1(hitpos, intpos);
  DCPlaneCenter[DC1Xp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC1Xp][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC1Xp][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 1;
  caldc1(hitpos, plane, wire, &stat);
  DCPlaneEdge2[DC1Xp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC1Xp][ZCOORD] = (double) hitpos[ZCOORD];

  /* DC1Y */
  plane = 1;
  wire = 1;
  sxdc1(intpos, plane, wire, &stat);
  intpos[ZCOORD] = -24.4;
  sgdc1(hitpos, intpos);
  DCPlaneEdge1[DC1Y][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC1Y][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0.0;
  intpos[YCOORD] = 0.0;
  intpos[ZCOORD] = -24.4;
  sgdc1(hitpos, intpos);
  DCPlaneCenter[DC1Y][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC1Y][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC1Y][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 48;
  sxdc1(intpos, plane, wire, &stat);
  intpos[ZCOORD] = -24.4;
  sgdc1(hitpos, intpos);
  DCPlaneEdge2[DC1Y][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC1Y][ZCOORD] = (double) hitpos[ZCOORD];

  /* DC1U */
  plane = 1;
  wire = 1;
  sxdc1(intpos, plane, wire, &stat);
  intpos[ZCOORD] = -29.2;
  sgdc1(hitpos, intpos);
  DCPlaneEdge1[DC1U][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC1U][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0.0;
  intpos[YCOORD] = 0.0;
  intpos[ZCOORD] = -29.2;
  sgdc1(hitpos, intpos);
  DCPlaneCenter[DC1U][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC1U][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC1U][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 48;
  sxdc1(intpos, plane, wire, &stat);
  intpos[ZCOORD] = -29.2;
  sgdc1(hitpos, intpos);
  DCPlaneEdge2[DC1U][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC1U][ZCOORD] = (double) hitpos[ZCOORD];

  /*printf("DC1X DCPlaneEdge1(x,z)=(%f,%f)\n",
     DCPlaneEdge1[DC1X][XCOORD],DCPlaneEdge1[DC1X][ZCOORD]);
     printf("DC1X DCPlaneEdge2(x,z)=(%f,%f)\n",
     DCPlaneEdge2[DC1X][XCOORD],DCPlaneEdge2[DC1X][ZCOORD]);
     printf("DC1Xp DCPlaneEdge1(x,z)=(%f,%f)\n",
     DCPlaneEdge1[DC1Xp][XCOORD],DCPlaneEdge1[DC1Xp][ZCOORD]);
     printf("DC1Xp DCPlaneEdge2(x,z)=(%f,%f)\n",
     DCPlaneEdge2[DC1Xp][XCOORD],DCPlaneEdge2[DC1Xp][ZCOORD]);
     printf("DC1Y DCPlaneEdge1(x,z)=(%f,%f)\n",
     DCPlaneEdge1[DC1Y][XCOORD],DCPlaneEdge1[DC1Y][ZCOORD]);
     printf("DC1Y DCPlaneEdge2(x,z)=(%f,%f)\n",
     DCPlaneEdge2[DC1Y][XCOORD],DCPlaneEdge2[DC1Y][ZCOORD]);
     printf("DC1Y DCPlaneCenter(x,y,z)=(%f,%f,%f)\n",
     DCPlaneCenter[DC1Y][XCOORD],DCPlaneCenter[DC1Y][YCOORD],DCPlaneCenter[DC1Y][ZCOORD]);
     printf("DC1U DCPlaneEdge1(x,z)=(%f,%f)\n",
     DCPlaneEdge1[DC1U][XCOORD],DCPlaneEdge1[DC1U][ZCOORD]);
     printf("DC1U DCPlaneEdge2(x,z)=(%f,%f)\n",
     DCPlaneEdge2[DC1U][XCOORD],DCPlaneEdge2[DC1U][ZCOORD]);
     printf("DC1U DCPlaneCenter(x,y,z)=(%f,%f,%f)\n",
     DCPlaneCenter[DC1U][XCOORD],DCPlaneCenter[DC1U][YCOORD],DCPlaneCenter[DC1U][ZCOORD]); */

  return;
}

void E42_chamber::PlaneEdgeDC2(void)
{
  int plane;
  int wire;
  int stat;

  double hitpos[3];
  double intpos[3];

  /* DC2X */
  plane = 1;
  wire = 1;
  caldc2(hitpos, plane, wire, &stat);
  DCPlaneEdge1[DC2X][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC2X][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0.0;
  intpos[YCOORD] = 0.0;
  intpos[ZCOORD] = 0.0;
  sgdc2(hitpos, intpos);
  DCPlaneCenter[DC2X][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC2X][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC2X][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 128;
  caldc2(hitpos, plane, wire, &stat);
  DCPlaneEdge2[DC2X][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC2X][ZCOORD] = (double) hitpos[ZCOORD];

  /* DC2Xp */
  plane = 2;
  wire = 1;
  caldc2(hitpos, plane, wire, &stat);
  DCPlaneEdge1[DC2Xp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC2Xp][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0;
  intpos[YCOORD] = 0;
  intpos[ZCOORD] = 7.8;
  sgdc2(hitpos, intpos);
  DCPlaneCenter[DC2Xp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC2Xp][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC2Xp][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 128;
  caldc2(hitpos, plane, wire, &stat);
  DCPlaneEdge2[DC2Xp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC2Xp][ZCOORD] = (double) hitpos[ZCOORD];

  /* DC2Y */
  plane = 1;
  wire = 1;
  sxdc2(intpos, plane, wire, &stat);
  intpos[ZCOORD] = 25.2;
  sgdc2(hitpos, intpos);
  DCPlaneEdge1[DC2Y][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC2Y][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0;
  intpos[YCOORD] = 0;
  intpos[ZCOORD] = 25.2;
  sgdc2(hitpos, intpos);
  DCPlaneCenter[DC2Y][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC2Y][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC2Y][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 128;
  sxdc2(intpos, plane, wire, &stat);
  intpos[ZCOORD] = 25.2;
  sgdc2(hitpos, intpos);
  DCPlaneEdge2[DC2Y][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC2Y][ZCOORD] = (double) hitpos[ZCOORD];

  /* DC2Yp */
  plane = 1;
  wire = 1;
  sxdc2(intpos, plane, wire, &stat);
  intpos[ZCOORD] = 33.0;
  sgdc2(hitpos, intpos);
  DCPlaneEdge1[DC2Yp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC2Yp][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0;
  intpos[YCOORD] = 0;
  intpos[ZCOORD] = 33.0;
  sgdc2(hitpos, intpos);
  DCPlaneCenter[DC2Yp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC2Yp][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC2Yp][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 128;
  sxdc2(intpos, plane, wire, &stat);
  intpos[ZCOORD] = 33.0;
  sgdc2(hitpos, intpos);
  DCPlaneEdge2[DC2Yp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC2Yp][ZCOORD] = (double) hitpos[ZCOORD];

  /*printf("DC2Y DCPlaneEdge1(x,z)=(%f,%f)\n",
     DCPlaneEdge1[DC2Y][XCOORD],DCPlaneEdge1[DC2Y][ZCOORD]);
     printf("DC2Y DCPlaneEdge2(x,z)=(%f,%f)\n",
     DCPlaneEdge2[DC2Y][XCOORD],DCPlaneEdge2[DC2Y][ZCOORD]);
     printf("DC2Y DCPlaneCenter(x,y,z)=(%f,%f,%f)\n",
     DCPlaneCenter[DC2Y][XCOORD],DCPlaneCenter[DC2Y][YCOORD],DCPlaneCenter[DC2Y][ZCOORD]);
     printf("DC2Yp DCPlaneEdge1(x,z)=(%f,%f)\n",
     DCPlaneEdge1[DC2Yp][XCOORD],DCPlaneEdge1[DC2Yp][ZCOORD]);
     printf("DC2Yp DCPlaneEdge2(x,z)=(%f,%f)\n",
     DCPlaneEdge2[DC2Yp][XCOORD],DCPlaneEdge2[DC2Yp][ZCOORD]);
     printf("DC2Yp DCPlaneCenter(x,y,z)=(%f,%f,%f)\n",
     DCPlaneCenter[DC2Yp][XCOORD],DCPlaneCenter[DC2Yp][YCOORD],DCPlaneCenter[DC2Yp][ZCOORD]); */

  return;
}

void E42_chamber::PlaneEdgeDC3(void)
{
  int plane;
  int wire;
  int stat;

  double hitpos[3];
  double intpos[3];

  /* DC3X */
  plane = 1;
  wire = 1;
  caldc3(hitpos, plane, wire, &stat);
  DCPlaneEdge1[DC3X][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC3X][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0;
  intpos[YCOORD] = 0;
  intpos[ZCOORD] = 0;
  sgdc3(hitpos, intpos);
  DCPlaneCenter[DC3X][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC3X][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC3X][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 32;
  caldc3(hitpos, plane, wire, &stat);
  DCPlaneEdge2[DC3X][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC3X][ZCOORD] = (double) hitpos[ZCOORD];

  /* DC3Xp */
  plane = 2;
  wire = 1;
  caldc3(hitpos, plane, wire, &stat);
  DCPlaneEdge1[DC3Xp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC3Xp][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0;
  intpos[YCOORD] = 0;
  intpos[ZCOORD] = 32.0;
  sgdc3(hitpos, intpos);
  DCPlaneCenter[DC3Xp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC3Xp][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC3Xp][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 32;
  caldc3(hitpos, plane, wire, &stat);
  DCPlaneEdge2[DC3Xp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC3Xp][ZCOORD] = (double) hitpos[ZCOORD];

  /* DC3Y */
  plane = 1;
  wire = 1;
  sxdc3(intpos, plane, wire, &stat);
  intpos[ZCOORD] = 63.0;
  sgdc3(hitpos, intpos);
  DCPlaneEdge1[DC3Y][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC3Y][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0;
  intpos[YCOORD] = 0;
  intpos[ZCOORD] = 63.0;
  sgdc3(hitpos, intpos);
  DCPlaneCenter[DC3Y][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC3Y][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC3Y][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 32;
  sxdc3(intpos, plane, wire, &stat);
  intpos[ZCOORD] = 63.0;
  sgdc3(hitpos, intpos);
  DCPlaneEdge2[DC3Y][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC3Y][ZCOORD] = (double) hitpos[ZCOORD];

  /* DC3Yp */
  plane = 1;
  wire = 1;
  sxdc3(intpos, plane, wire, &stat);
  intpos[ZCOORD] = 93.0;
  sgdc3(hitpos, intpos);
  DCPlaneEdge1[DC3Yp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge1[DC3Yp][ZCOORD] = (double) hitpos[ZCOORD];
  intpos[XCOORD] = 0;
  intpos[YCOORD] = 0;
  intpos[ZCOORD] = 93.0;
  sgdc3(hitpos, intpos);
  DCPlaneCenter[DC3Yp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneCenter[DC3Yp][YCOORD] = (double) hitpos[YCOORD];
  DCPlaneCenter[DC3Yp][ZCOORD] = (double) hitpos[ZCOORD];
  wire = 32;
  sxdc3(intpos, plane, wire, &stat);
  intpos[ZCOORD] = 93.0;
  sgdc3(hitpos, intpos);
  DCPlaneEdge2[DC3Yp][XCOORD] = (double) hitpos[XCOORD];
  DCPlaneEdge2[DC3Yp][ZCOORD] = (double) hitpos[ZCOORD];

  /*printf("DC3Y DCPlaneEdge1(x,z)=(%f,%f)\n",
     DCPlaneEdge1[DC3Y][XCOORD],DCPlaneEdge1[DC3Y][ZCOORD]);
     printf("DC3Y DCPlaneEdge2(x,z)=(%f,%f)\n",
     DCPlaneEdge2[DC3Y][XCOORD],DCPlaneEdge2[DC3Y][ZCOORD]);
     printf("DC3Y DCPlaneCenter(x,y,z)=(%f,%f,%f)\n",
     DCPlaneCenter[DC3Y][XCOORD],DCPlaneCenter[DC3Y][YCOORD],DCPlaneCenter[DC3Y][ZCOORD]);
     printf("DC3Yp DCPlaneEdge1(x,z)=(%f,%f)\n",
     DCPlaneEdge1[DC3Yp][XCOORD],DCPlaneEdge1[DC3Yp][ZCOORD]);
     printf("DC3Yp DCPlaneEdge2(x,z)=(%f,%f)\n",
     DCPlaneEdge2[DC3Yp][XCOORD],DCPlaneEdge2[DC3Yp][ZCOORD]);
     printf("DC3Yp DCPlaneCenter(x,y,z)=(%f,%f,%f)\n",
     DCPlaneCenter[DC3Yp][XCOORD],DCPlaneCenter[DC3Yp][YCOORD],DCPlaneCenter[DC3Yp][ZCOORD]); */

  return;
}

double E42_chamber::get_DCPlaneEdge1(int plane, int coord)
{
  return DCPlaneEdge1[plane][coord];
}

double E42_chamber::get_DCPlaneEdge2(int plane, int coord)
{
  return DCPlaneEdge2[plane][coord];
}

double E42_chamber::get_DCPlaneCenter(int plane, int coord)
{
  return DCPlaneCenter[plane][coord];
}

void E42_chamber::get_DCPlaneCenter(int plane, double *pos)
{
  pos[XCOORD] = DCPlaneCenter[plane][XCOORD];
  pos[YCOORD] = DCPlaneCenter[plane][YCOORD];
  pos[ZCOORD] = DCPlaneCenter[plane][ZCOORD];
  return;
}

int E42_chamber::getPlaneID(char *planenam, int *plane)
{
/*******************************************************/
/* Convert the Plane Name to Plane ID                  */
/*                        of datastructure BDC or DC   */
/*      planenam : Plane Name(for ex. DC1X)            */
/*      plane    : Plane ID                            */
/*      chamid   : chamber id which identify BDC or DC */
/*                   0:BDC 1:DC -1:Error               */
/*******************************************************/

  int chamid;			/* BDC:0 DC:1 */

  chamid = 0;
  if (strcmp(planenam, "BDC1X") == 0) {
    *plane = BDC1X;
    return (chamid);
  }
  if (strcmp(planenam, "BDC1Xp") == 0) {
    *plane = BDC1Xp;
    return (chamid);
  }
  if (strcmp(planenam, "BDC2X") == 0) {
    *plane = BDC2X;
    return (chamid);
  }
  if (strcmp(planenam, "BDC2Xp") == 0) {
    *plane = BDC2Xp;
    return (chamid);
  }
  if (strcmp(planenam, "BDC3X") == 0) {
    *plane = BDC3X;
    return (chamid);
  }
  if (strcmp(planenam, "BDC3Xp") == 0) {
    *plane = BDC3Xp;
    return (chamid);
  }

  chamid = 1;
  if (strcmp(planenam, "DC1X") == 0) {
    *plane = DC1X;
    return (chamid);
  }
  if (strcmp(planenam, "DC1Y") == 0) {
    *plane = DC1Y;
    return (chamid);
  }
  if (strcmp(planenam, "DC1Xp") == 0) {
    *plane = DC1Xp;
    return (chamid);
  }
  if (strcmp(planenam, "DC1U") == 0) {
    *plane = DC1U;
    return (chamid);
  }

  if (strcmp(planenam, "DC2X") == 0) {
    *plane = DC2X;
    return (chamid);
  }
  if (strcmp(planenam, "DC2Xp") == 0) {
    *plane = DC2Xp;
    return (chamid);
  }
  if (strcmp(planenam, "DC2Y") == 0) {
    *plane = DC2Y;
    return (chamid);
  }
  if (strcmp(planenam, "DC2Yp") == 0) {
    *plane = DC2Yp;
    return (chamid);
  }

  if (strcmp(planenam, "DC3X") == 0) {
    *plane = DC3X;
    return (chamid);
  }
  if (strcmp(planenam, "DC3Xp") == 0) {
    *plane = DC3Xp;
    return (chamid);
  }
  if (strcmp(planenam, "DC3Y") == 0) {
    *plane = DC3Y;
    return (chamid);
  }
  if (strcmp(planenam, "DC3Yp") == 0) {
    *plane = DC3Yp;
    return (chamid);
  }

  fprintf(stderr, "getPlaneID::Invarid plane %s\n", planenam);
  fflush(stderr);
  return (-1);
}


void E42_chamber::get_PlaneName(char *chamnam, int plane, char *planename)
{
  
  if (strcmp(chamnam,"BPC")==0) {
    switch (plane) {
    case BPC1X:
      sprintf(planename,"BPC1X");
      break;
    case BPC1Y:
      sprintf(planename,"BPC1Y");
      break;
    case BPC2X:
      sprintf(planename,"BPC2X");
      break;
    case BPC3X:
      sprintf(planename,"BPC3X");
      break;
    case BPC3Y:
      sprintf(planename,"BPC3Y");
      break;
    case BPC4X:
      sprintf(planename,"BPC4X");
      break;
    case BPC4Y:
      sprintf(planename,"BPC4Y");
      break;
    case BPC5X:
      sprintf(planename,"BPC5X");
      break;
    case BPC5Y:
      sprintf(planename,"BPC5Y");
      break;
    default:
      fprintf(stderr,"BPC do not have such plane%d\n",plane);
      break;
    }
  } else if (strcmp(chamnam,"PC") == 0) {
    switch(plane) {
    case VHx:
      sprintf(planename,"VHx");
      break;
    case VHy:
      sprintf(planename,"VHy");
      break;
    case CHx:
      sprintf(planename,"CHx");
      break;
    case FTOFx:
      sprintf(planename,"FTOFx");
      break;
    default:
      fprintf(stderr,"PC do not have such plane%d\n",plane);
    }
  } else if (strcmp(chamnam,"DC") == 0) {
    switch (plane) {
    case DC1X:
      sprintf(planename,"DC1X");
      break;
    case DC1Y:
      sprintf(planename,"DC1Y");
      break;
    case DC1Xp:
      sprintf(planename,"DC1Xp");
      break;
    case DC1U:
      sprintf(planename,"DC1U");
      break;
    case DC2X:
      sprintf(planename,"DC2X");
      break;
    case DC2Xp:
      sprintf(planename,"DC2Xp");
      break;
    case DC2Y:
      sprintf(planename,"DC2Y");
      break;
    case DC2Yp:
      sprintf(planename,"DC2Yp");
      break;
    case DC3X:
      sprintf(planename,"DC3X");
      break;
    case DC3Xp:
      sprintf(planename,"DC3Xp");
      break;
    case DC3Y:
      sprintf(planename,"DC3Y");
      break;
    case DC3Yp:
      sprintf(planename,"DC3Yp");
      break;
    default:
      fprintf(stderr,"DC do not have such plane%d\n",plane);
    }
  } else if (strcmp(chamnam,"BDC") == 0) {
    printf("now BDC is not analyzed\n");
  } else {
    fprintf(stderr,"%s does not exist\n",chamnam);
  }
  
  return;
}

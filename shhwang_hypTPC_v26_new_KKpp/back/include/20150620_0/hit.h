#ifndef __HIT_H_
#define __HIT_H_
//#include <cluster.h>
//#include <tpcParameter.h>
//#include <switch.h>

typedef struct Hit{
  int peakPad;           /* padId, which has largest pulse height  */
  int peakTdc;           /* tdc of peaking time on peakPad         */
  float peakTdcfit;      /* tdc of peaking time on peakPad  fitting*/
  float secondpeakTdc;
  float secondpeakTdcfit;
  int sector; /*NTPC*/
  int layer;
  int numPads;            /* # of pads in this hit                  */
  int quality;
  int ithPulse;           /* average of ithPulse in Pulse structure */
  int numPulse;
  Pulse* pulse[MAX_NUM_PULSE_IN_HIT];

  /*  PulseInWire wirepulse[NUM_WIRE_CLUST];*/
  /*  using pointer version -- memory leak problem
      int wireIdinClust[NUM_WIRE_CLUST]; */
  PulseInWire* wirepulse[NUM_WIRE_CLUST]; /*NTPC wire*/
  /*test wire line fit*/
  PulseInWire* wirepulseOnTime[NUM_WIRE_CLUST][MAX_ONTIME_PULSE];
  
  int numOnTimePulse[NUM_WIRE_CLUST]; /*NTPC wire(test)*/
  float wireweight[NUM_WIRE_CLUST]; /*NTPC wire(test)*/
  double dXdZdrift[2];    /* tangent vector on the drift line at this hit
                            dR/dZ, dPhi/dZ */
  double xlocal[3];      /* cartesian coordinate of this hit in local sector*/
  double x[3];           /* cartesian coordinate of this hit in global system*/
  double err[3];          /* error in (x,y,z)                    */
  float adc[MAX_NUM_PAD_IN_HIT];
                         /* pulse heights on each pad              */
                         /*     in decreasing order                */
  float maxAdc;          /* maximum adc value in this hit          */
  float maxAdcfit;       /* maximum adc value in this hit   fitting*/
  float meanTdc;
  float driftTime;
  float sumAdcPeakPad;   /* adc summed at the PeakPad              */
  float sumAdc;          /* adc value summed up in this hit        */
  float gaussHei;        /* height of gaussian for center of gravity  */
  //  float res;             /* resolution, shhwang 2013. 6. 10.                 */
  float prw;             /* pad response width (mm)                */
  float dx[2];           /* distance between center of gravity     */
                         /* and pad center (x and Y direction NTPC)*/
  float rmsPad;          /* rms of dist. of Adc_i*(xPeakPad-x_i) */
  float rmsZ;            /* rms of dist. of Adc_i*(zPeakTdc-z_i) */
  float tandip; /*obtained from wire TDC diff   for dip angle correction*/
  float tandippad; /*obtained from wire TDC diff   for dip angle correction*/
  float tandipRK;   /* z/xy-plane in the local coordianate for each sector*/
  float tandipyzRK; /* z/y */
  float tanphiRK;   /* x/y */
  int pidsim; /*g3 pid*/
  int tidsim; /*g3 track id*/
}Hit;

#define PILEUP             -4
#define GOOD_HIT            1  
#define TOO_SMALL_PULSE    -2
#define MULTI_PEAK_IN_HIT  -3  

//int anaCluster( Cluster*** clusters, Hit hits[][NUM_LAY][MAX_HIT_IN_LAY],
//                Hit selectedHits[][NUM_LAY][MAX_SELECTED_HITS], Wire* wires,TpcParam*, Switch*);

/*** NTPC test ***/
int nrawhit; /* hit w/o selection */
/*****************/

#endif

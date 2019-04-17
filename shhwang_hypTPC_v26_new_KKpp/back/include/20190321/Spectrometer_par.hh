/**********************************************/
/*  copy of Miwa-san's code                   */
/*  Will update soon(2012. 7. 29.)            */
/*                                            */
/**********************************************/

#ifndef __E42_SPECTROMETER_CONST_H__
#define __E42_SPECTROMETER_CONST_H__

#define MULTIMAXCHAM 30
#define ClustMax   15

#define BPCMAX 5
#define BDCMAX 3
#define PCMAX 3
#define DCMAX 3

#define VH 0
#define CHPC 1
#define FTOFPC 2

#define DC1 0
#define DC2 1
#define DC3 2

#define BPCPlaneMax 9
#define BPCwireMax  144
#define BPCYwireMax 64
#define BPCwireSpace 1.0        /* mm */
#define BPC1X 0
#define BPC1Y 1
#define BPC2X 2
#define BPC3X 3
#define BPC3Y 4
#define BPC4X 5
#define BPC4Y 6
#define BPC5X 7
#define BPC5Y 8

#define PCPlaneMax 4
#define PCwireMax 144
#define CounterMax 5
#define PCMultiMax 10
#define VHx 0
#define VHy 1
#define CHx 2
#define FTOFx 3
#define YHy 4

#define VHXMAX 32
#define VHYMAX 18
#define CHMAX 64 ///e42
//#define CHMAX 48 ///e07
#define FTOFMAX 24
#define YHMAX 6

#define BDCPlaneMax 6
#define BDCwireMax 32
#define BDCwireSpace 5.0        /* mm */
#define BDCMultiMax 5
#define BDC1X  0
#define BDC1Xp 1
#define BDC2X  2
#define BDC2Xp 3
#define BDC3X  4
#define BDC3Xp 5

#define DCPlaneMax 12
#define DCwireMax 128
#define DCMultiMax 10
#define DC1X  0
#define DC1Y  1
#define DC1Xp 2
#define DC1U  3
#define DC2X  4
#define DC2Xp 5
#define DC2Y  6
#define DC2Yp 7
#define DC3X  8
#define DC3Xp 9
#define DC3Y  10
#define DC3Yp 11


#define MaxWireBDC 32

#define MaxWireDC1X  48
#define MaxWireDC1Xp 48
#define MaxWireDC1Y  32
#define MaxWireDC1U  48

#define MaxWireDC2X  128
#define MaxWireDC2Xp 128
#define MaxWireDC2Y  96
#define MaxWireDC2Yp 96

#define MaxWireDC3X  32
#define MaxWireDC3Xp 32
#define MaxWireDC3Y  16
#define MaxWireDC3Yp 16

//void Spectrometer_par(){

//}

#endif

#ifndef PADHELPER_HH
#define PADHELPER_HH
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <TMath.h>
#include <TVector3.h>

namespace padHelper
{

	
	
	
	
	
	
	
	
	
	
	
	//#PadID is defined as 0 origin
  //#OfPad #division #radius padLength
  static const Double_t padParameter[32][6]=
    {{0, 48,    14.75, 48, 0,  9.},
     {1, 48,    24.25, 48, 0,  9.},
     {2, 72,    33.75, 72, 0,  9.},
     {3, 96,    43.25, 96, 0,  9.},
     {4, 120,    52.75,120,0,   9.},
     {5, 144,    62.25,144,0,   9.},
     {6, 168,    71.75,168,0,   9.},
     {7, 192,    81.25,192,0,   9.},
     {8, 216,    90.75,216,0,   9.},
     {9, 240,    100.25,240,0,  9.},
     {10,208,    111.5,241, 0,  12.5},
     {11,218,    124.5,271, 0,  12.5},
     {12,230,    137.5,300, 0,  12.5},
     {13,214,    150.5,330, 0,  12.5},
     {14,212,    163.5,360, 0,  12.5},
     {15,214,    176.5,390, 0,  12.5}, 
     {16,220,    189.5,420, 0,  12.5},
     {17,224,    202.5,449, 0,  12.5},
     {18,232,    215.5,479, 0,  12.5},
     {19,238,    228.5,509, 0,  12.5},
     {20,244,    241.5,539, 0,  12.5},
     {21,232,    254.5,569, 0,  12.5},
     {22,218,    267.5,599, 0,  12.5},
     {23,210,    280.5,628, 0,  12.5},
     {24,206,    293.5,658, 0,  12.5},
     {25,202,    306.5,688, 0,  12.5},
     {26,200,    319.5,718, 0,  12.5},
     {27,196,    332.5,748, 0,  12.5},
     {28,178,    345.5,777, 0,  12.5},
     {29,130,    358.5,807, 0,  12.5},
     {30,108,    371.5,837, 0,  12.5},
     {31,90,     384.5,867, 0, 12.5}};
	static const Int_t padOnCenterFrame[] =
	{

		//Pads on the frame
		965,966,967,968,969,970,971,972,973,974,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1272,1395,1396,1397,1398,1399,1400,1401,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1428,1429,1430,1431,1432,1436,1437,1438,1439,1440,1441,1442,1455,1456,1457,1458,1459,1460,1461,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,1485,1486,1487,1488,1489,1490,1491,1492,1496,1497,1498,1499,1500,1501,1502,1595,1596,1597,1598,1604,1605,1606,1607,1608,1647,1648,1649,1650,1651,1652,1658,1659,1660,1663,1664,1665,1671,1672,1673,1674,1675,1676,1715,1716,1717,1718,1719,1725,1726,1727,1728,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1818,1879,1880,1881,1882,1889,1890,1891,1892,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,1964,2018,2019,2020,2025,2026,2027,2028,2100,2101,2102,2103,2104,2111,2112,2113,2114,2115,2187,2188,2189,2190,2195,2196,2197,2219,2220,2221,2222,2223,2224,2225,2226,2227,2228,2309,2310,2311,2317,2318,2323,2324,2330,2331,2332,2413,2414,2415,2416,2417,2418,2419,2420,2421,2422,2427,2428,2429,2430,2518,2519,2520,2521,2522,2523,2524,2525,2526,2527,2540,2541,2542,2543,2544,2545,2546,2547,2548,2549,2637,2638,2639,2640,2732,2733,2739,2740,2761,2762,2768,2769,2950,2951,2952,2953,2954,2955,2956,2957,2958,2987,2988,2989,2990,2991,2992,2993,2994,2995,3174,3175,3176,3177,3178,3179,3180,3181,3182,3219,3220,3221,3222,3223,3224,3225,3226,3227,3405,3406,3407,3408,3409,3410,3411,3412,3413,3458,3459,3460,3461,3462,3463,3464,3465,3466,3642,3643,3644,3645,3646,3647,3648,3649,3650,3703,3704,3705,3706,3707,3708,3709,3710,3711,3877,3878,3879,3880,3881,3882,3883,3884,3945,3946,3947,3948,3949,3950,3951,3952,4098,4099,4105,4174,4180,4181,4308,4309,4310,4311,4312,4313,4314,4315,4316,4391,4392,4393,4394,4395,4396,4397,4398,4399,4512,4513,4520,4603,4610,4611,4713,4714,4715,4716,4717,4718,4719,4720,4811,4812,4813,4814,4815,4816,4817,4818,4910,4911,4912,4913,4914,4915,4916,4917,5016,5017,5018,5019,5020,5021,5022,5023,5104,5105,5108,5111,5112,5217,5218,5221,5224,5225,5287,5288,5289,5290,5291,5292,5293,5294,5295,5408,5409,5410,5411,5412,5413,5414,5415,5416,5441,5442,5443,5444,5445,5566,5567,5568,5569,5570,

		//Empty Pads
		1017,1018,
		1394,1402,1403,1404,1433,1434,1435,1462,1463,1464,1493,1494,1495,1503,1599,1600,1601,1602,1603,1604,1653,1654,1655,1656,1657,1666,1667,1668,1669,1670,1720,1721,1722,1723,1724,1877,1878,1883,1884,1885,1886,1887,1888,1893,1894,2021,2022,2023,2024,2105,2106,2107,2108,2109,2110,2191,2192,2193,2194,2312,2313,2314,2315,2316,2325,2326,2327,2328,2329,2734,2735,2736,2737,2738,2763,2764,2765,2766,2767,4100,4101,4102,4103,4104,4175,4176,4177,4178,4179,4514,4515,4516,4517,4518,4519,4604,4605,4606,4607,4608,4609,5106,5107,5109,5110,5219,5220,5222,5223
	};

	static const Int_t deadChannel[32][5] =
	{
		//18,179,333,3809, //Dead pads
		//2748, //Empty pads
		//Pads on the gem supporing frame
		//73,135,136,223,333,409,468,627,810,60,111,186,285,556,727,922,1141,1364,1565,1566,1779,1819,1820,2037,2038,2246,2247,2248,2456,2457,2669,2670,2888,3112,3113,3343,3344,3580,3581,3815,4036,1630,1851,2069,2278,2487,2488,2700,2701,2919,3143,3375,3612,3846,4067,4278,4482,4682,4880,5074
		{18, -1, -1, -1, -1}, //0
		{60, 73, -1, -1, -1},
		{111, 135, 136, -1, -1},
		{179, 186, 223, -1, -1},
		{285, 333, -1, -1, -1},
		{409, 468, -1, -1, -1}, //5
		{556, 627, -1, -1, -1},
		{727, 810, -1, -1, -1},
		{922, -1, -1, -1, -1},
		{1141, -1, -1, -1, -1},
		{1364, -1, -1, -1, -1}, //10
		{1565, 1566, 1630, -1, -1},
		{1779, 1819, 1820, 1851, -1},
		{2037, 2038, 2069, -1, -1},
		{2246, 2247, 2248, 2278, -1},
		{2456, 2457, 2487, 2488, -1}, //15
		{2669, 2670, 2700, 2701, 2748},
		{2888, 2919, -1, -1, -1},
		{3112, 3113, 3143, -1, -1},
		{3343, 3344, 3375, -1, -1},
		{3580, 3581, 3612, -1, -1}, //20
		{3809, 3815, 3846, -1, -1},
		{4036, 4067, -1, -1, -1},
		{4278, -1, -1, -1, -1},
		{4481, 4482, -1, -1, -1},
		{4682, -1, -1, -1, -1}, //25
		{4880, -1, -1, -1, -1},
		{5074, -1, -1, -1, -1},
		{-1, -1, -1, -1, -1},
		{-1, -1, -1, -1, -1},
		{-1, -1, -1, -1, -1}, //30
		{-1, -1, -1, -1, -1}

	};

  inline Double_t getDTheta(Int_t layerID){
    return (360./padParameter[layerID][3]);
  }

  inline Double_t getsTheta(Int_t layerID)
  {
    Double_t sTheta = 180.-(360./padParameter[layerID][3])*padParameter[layerID][1]/2.;
    return sTheta;
  }

  inline Double_t getRadius(Int_t layerID)
  {
    return padParameter[layerID][2];
  }

  inline Double_t getLength(Int_t layerID)
  {
    return padParameter[layerID][5];
  }

  inline Int_t getPadID(Int_t layerID, Int_t rowID)
  {
    Int_t padID=0;
    for(int layi = 0 ; layi<layerID; layi++) padID += padParameter[layi][1];
    padID+=rowID;
    return padID;

  }

  inline Int_t getLayerID(Int_t padID)
  {
    //    padID-=1;
    int layer;
    int sum = 0;

    for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
      {
	sum += padParameter[layer][1];
      }
    return layer;
  }

  inline Int_t getRowID(Int_t padID)
  {
    //    padID-=1;
    int layer, row;
    int sum = 0;

    for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
      {
	sum += padParameter[layer][1];
      }
    row = padID - sum;
    return row;
  }
  /*
    Double_t getTheta(Int_t layerID, Int_t rowID)
    {
    Double_t sTheta = 180.-(360./padParameter[layerID][3])*padParameter[layerID][1]/2.;
    Double_t theta = sTheta+(rowID+0.5)*(360.-2*sTheta)/padParameter[layerID][1];
    return theta;
    }
  */


  inline int findPadID(double z, double x) 
  { 
    z += 143;
    double radius = sqrt(x*x + z*z);
    double angle;
    if (z == 0)
      {
	if (x > 0)   angle = 1.5*TMath::Pi();
	else if (x < 0)   angle = 0.5*TMath::Pi();
	else return -1000; // no padID if (0,0)
      }
    else
      {
				if (z > 0) angle = TMath::Pi()+atan(x / z);
				else if( z < 0&&x<0) angle = atan(x / z);
				else angle = 2*TMath::Pi()+ atan(x / z);
//	angle = TMath::Pi() - atan(-x / z);
      }

    int layer, row;
    // find layer_num. 
    for (layer = 0; !(padParameter[layer][2]+padParameter[layer][5]*0.5 >= radius
		      && padParameter[layer][2]-padParameter[layer][5]*0.5 <= radius); layer++)
      {
	if (layer >= 32) return -1000;
	if (layer != 0)
	  {
	    if (padParameter[layer][2] - padParameter[layer][5] * 0.5 >= radius && 
		padParameter[layer - 1][2] + padParameter[layer - 1][5] * 0.5 <= radius) return -layer;
	  }
      }
    
    
    //std::cout<<"padHelper:: layer="<<layer<<", angle="<<angle<<", "<<(getsTheta(layer)*TMath::Pi()/180.)<<std::endl;
    // find row_num
    //  if (angle - (padParameter[layer][4]*TMath::Pi()/180.) < 0) return -1000;
    if (angle - (getsTheta(layer)*TMath::Pi()/180.) < 0) return -2000;

    //    double a, b, c;
    //row = (int)((angle-(padParameter[layer][4]*TMath::Pi()/180.))/(padParameter[layer][3]*TMath::Pi()/180.));
    //    row = (int)((angle-(getsTheta(layer)*TMath::Pi()/180.))/(getDTheta(layer)*TMath::Pi()/180.));
    
    //row = (int)((angle-(getsTheta(layer)*TMath::Pi()/180.))/(getDTheta(layer)*TMath::Pi()/180.))+1;
    row = (int)((angle-(getsTheta(layer)*TMath::Pi()/180.))/(getDTheta(layer)*TMath::Pi()/180.));
    if (row > padParameter[layer][1]) return -1000;

    return getPadID(layer, row);
  }
	inline Bool_t Dead(Int_t padID){

		Bool_t centerframe = std::find(std::begin(padOnCenterFrame), std::end(padOnCenterFrame), padID) != std::end(padOnCenterFrame);
		Int_t layer = getLayerID(padID);
		Bool_t dead = std::find(std::begin(deadChannel[layer]), std::end(deadChannel[layer]), padID) != std::end(deadChannel[layer]);
		if(centerframe||dead) return true;
		else return false;
	}

	//_____________________________________________________________________________
	inline Bool_t Dead(Int_t layer, Int_t row){

		Int_t padID = getPadID(layer, row);
		return Dead(padID);
	}




  inline Double_t getTheta(Int_t padID)
  {
    //    padID-=1;
    int layer, row;
    int sum = 0;

    for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
      {
	sum += padParameter[layer][1];
      }
    row = padID - sum;
    Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
    Double_t theta = sTheta+(row+0.5)*(360.-2*sTheta)/padParameter[layer][1];
    return theta;
  }

  inline Double_t getR(Int_t padID)
  {
    //    padID-=1;
    int layer, row;
    int sum = 0;

    for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
      {
	sum += padParameter[layer][1];
      }
    row = padID - sum;
    Double_t R = padParameter[layer][2];
    return R;
  }

  inline TVector3 getPoint(int padID)
  {
    //padID-=1;
    int layer, row;
    int sum = 0;

    for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
      {
	sum += padParameter[layer][1];
      }
    row = padID - sum;

    TVector3 result;
    if (row > padParameter[layer][1]){ // out of range
      result.SetX(0);
      result.SetY(-1);
      result.SetZ(0);
    }
    else{
      double x, z;
      Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
      x = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
      z = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) - 143.0;
      result.SetX(x);
      result.SetY(0);
      result.SetZ(z);
    }
    return result;
  }

}
#endif

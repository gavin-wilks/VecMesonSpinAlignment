#ifndef StEventPlaneCons_h
#define StEventPlaneCons_h

#include <string>
#include "TString.h"
// #include "StarClassLibrary/SystemOfUnits.h"

namespace recoEP
{
  //--------------------------------------------------
  // used in Event Plane Reconstruction
  const int NumBeamEnergy = 2;
  const std::string mBeamEnergy[NumBeamEnergy] = {"14p5GeV_2019","19p6GeV_2019"};
  const double mEnergyValue[NumBeamEnergy] = {14.5,19.6};
  const int mBeamYear[NumBeamEnergy] = {2019,2019};

  // event cut
  const int Vz_bins = 4;
  const double mVzMaxMap[NumBeamEnergy] = {70.0,70.0}; // 0: 200GeV_2014 | 1: 54GeV_2017 | 2: 27GeV_2018 
  const double mVrMax[NumBeamEnergy] = {2.0,2.0};
  const double mVrMin[NumBeamEnergy] = {0.0,0.0};
  const double mVzVpdDiffMax[NumBeamEnergy] = {10.0,10.0}; // 3.0
  const int mMatchedToFMin[NumBeamEnergy] = {2,2}; // 2

  // track cut
  const double mSigScaleMap[NumBeamEnergy] = {1.0,1.0};
  const double mDcaEPMax[NumBeamEnergy] = {1.0,1.0}; // for event plane reconstruction: 1.0 for BES
  const double mDcaTrMax = 1.0; // for pion, kaon, proton mDcaTrMax = 1.0 for flow
  const int mHitsDedxMin = 5;
  const int mHitsFitTPCMin = 15;
  const int mHitsMaxTPCMin = 0;
  const double mHitsRatioTPCMin = 0.52;
  const double mEtaMax = 1.5;
  const double mEtaGap = 0.05;
  const double mPrimPtMin[NumBeamEnergy] = {0.15,0.15}; // for event plane reconstruction and for pion, kaon, proton: 0.2 for BES
  const double mPrimPtMax = 2.0;
  const double mPrimPtWeight = 2.0;
  const double mPrimMomMax = 10.0; // also use for gMom
  const double mGlobPtMin = 0.1; // for phi, Lambda, K0s
  const int mTrackMin = 2;
  const int mTrackMinFull = 4;

  const int mNumOfRunIndex = 2000;

  // ZDC-SMD Event Plane
  const std::string mEastWest[2] = {"East","West"};
  const std::string mVertHori[2] = {"Vertical","Horizontal"};

  const std::string mVStr[4] = {"neg70","neg30","pos30","pos70"};
  const std::string mOrder = "2nd";
  const int mNumShiftOrder = 5;

  const double pl19[5][5] = {{-22.706,0.863809,-0.00036812,5.97123e-06,-1.27439e-08},   //-145 < vz < -87  cm
                       {-15.5924,0.604207,0.00131806,-2.04754e-06,5.73182e-10},   //-87  < vz < -29  cm 
                       {-13.1158,0.504708,0.00187998,-4.7317e-06,4.82342e-09},    //-29  < vz <  29  cm
                       {-15.9411,0.615061,0.00118242,-1.48902e-06,-2.29371e-10},  // 29  < vz <  87  cm 
                       {-22.468,0.839062,-0.000106213,4.93946e-06,-1.14503e-08}}; // 87  < vz <  145 cm

  const double ph19[5][5] = {{33.7733,1.75938,-0.00285868,8.50344e-06,-1.19215e-08},  //-145 < vz < -87  cm
                     {20.5875,1.67371,-0.00307534,7.93756e-06,-8.46257e-09},  //-87  < vz < -29  cm 
                     {15.1015,1.53929,-0.00269203,6.48876e-06,-6.06074e-09},  //-29  < vz <  29  cm
                     {20.7719,1.67316,-0.00315093,8.35824e-06,-9.14822e-09},  // 27  < vz <  87  cm 
                     {33.4926,1.79373,-0.00319461,9.56613e-06,-1.31049e-08}}; // 87  < vz <  145 cm

  //******** POINTWISE VZ CORRECTION *********//
  const double vz_corr[29] = {
  0.998171,  // (-145,-135)cm
  0.99364,   // (-135,-125)cm
  0.991578,  // (-125,-115)cm
  0.990252,  // (-115,-105)cm
  0.990494,  // (-105,-95)cm
  0.990065,  // (-95,-85)cm
  0.990332,  // (-85,-75)cm
  0.996478,  // (-75,-65)cm
  0.999687,  // (-65,-55)cm
  0.998645,  // (-55,-45)cm
  0.993835, // (-45,-35)cm
  0.996273, // (-35,-25)cm
  0.998307, // (-25,-15)cm
  0.999295, // (-15,-5)cm
  1,        // (-5,5)cm
  1.00056,  // (5,15)cm
  1.00019,  // (15,25)cm
  0.999894, // (25,35)cm
  0.998907, // (35,45)cm
  1.00468,  // (45,55)cm
  1.0055,   // (55,65)cm
  1.00142,  // (65,75)cm
  0.996247, // (75,85)cm
  0.995789, // (85,95)cm
  0.996513, // (95,105)cm
  0.996928, // (105,115)cm
  0.998196, // (115,125)cm
  1.00097,  // (125,135)cm
  1.00699}; // (135,145)cm 
 
  // EPD Event Plane
  const int mEpdEpOrder = 3;
}

#endif

#ifndef StToFMatchCons_h
#define StToFMatchCons_h

#include "TMath.h"
#include <string>

namespace tof
{
  std::string const mPID_ToF[3] = {"K","Pi","P"};
  std::string const mCharge[2] = {"plus","minus"};
  int const mPID_Start = 0;
  int const mPID_Stop  = 3;
  int const pidQA = 0;
  int const chargeQA = 0;

  // int const BinPt  = 160; // for efficiency
  int const BinPt  = 40; // for efficiency
  float const ptMin = 0.2;
  float const ptMax = 8.2;

  // int const BinEta = 20;
  int const BinEta = 10;
  float const etaMin = -1.0;
  float const etaMax = 1.0;
  int const etaQA = 0;

  // int const BinPhi =36 ;
  int const BinPhi = 12;
  float const phiMin = -1.0*TMath::Pi();
  float const phiMax = 1.0*TMath::Pi();
  int const phiQA = 0;

  int const NCentMax = 9; 
  float const weight[NCentMax] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.5};
  int const centQA = 5;

  int const BinCos = 7;
}

#endif

#ifndef StToFMatchCons_h
#define StToFMatchCons_h

#include "TMath.h"
#include <string>

namespace tof
{
  std::string const mPID_ToF[3] = {"K","Pi","P"};
  std::string const mCharge[2] = {"plus","minus"};
  int const mPID_Start = 0;
  int const mPID_Stop  = 2;
  int const pidQA = 0;
  int const chargeQA = 1;

  // int const BinPt  = 160; // for efficiency
  int const BinPt  = 80; // for efficiency
  float const ptMin = 0.2;
  float const ptMax = 8.0;

  int const BinEta = 15;
  // int const BinEta = 1;
  float const etaMin = -1.5;
  float const etaMax = 1.5;
  int const etaQA = 5;

  int const BinPhi = 12;
  // int const BinPhi = 1;
  float const phiMin = -1.0*TMath::Pi();
  float const phiMax = 1.0*TMath::Pi();
  int const phiQA = 4;

  int const NCentMax = 9; 
  float const weight[NCentMax] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.5};
  int const centQA = 5;

  int const BinCos = 7;
}

#endif

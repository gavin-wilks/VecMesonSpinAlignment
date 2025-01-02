#ifndef StToFMatchCons_h
#define StToFMatchCons_h

#include "TMath.h"
#include <string>

namespace tof
{
  std::string const mPID_ToF[3] = {"K","Pi","P"};
  std::string const mCharge[2] = {"plus","minus"};
  int const mPID_Start = 0;
  int const mPID_Stop  = 1;
  int const pidQA = 0;
  int const chargeQA = 1;

  // int const BinPt  = 160; // for efficiency
  int const BinPt  = 40; // for efficiency
  float const ptMin = 0.1;
  float const ptMax = 5.0;

  int const BinEta = 24;
  int const BinPhi = 240;
  // int const BinPhi = 1;
  //float const phiMin = -11./12.*TMath::Pi();
  //float const phiMax =  13./12.*TMath::Pi();
  float const phiMin = -TMath::Pi();
  float const phiMax =  TMath::Pi();
  int const phiQA = 4;

  int const NCentMax = 9; 
  float const weight[NCentMax] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.5};
  int const centQA = 5;

  int const BinCos = 7;
}

#endif

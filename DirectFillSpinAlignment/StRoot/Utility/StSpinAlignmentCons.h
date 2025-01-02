#ifndef StSpinAlignmentCons_h
#define StSpinAlignmentCons_h

#include <string>
#include "TString.h"
#include "TMath.h"

// #include "StarClassLibrary/SystemOfUnits.h"

namespace vmsa
{
  //--------------------------------------------------
  // used in TreeProduction and FillSpinAlginment
  int const NumBeamEnergy = 5;
  // event cut
  float const mVzMaxMap[NumBeamEnergy] = {70.0,70.0,70.0,70.0,70.0}; // 7.7, 9.1, 11.5, 14.6, 19.6 (BES-II)
  float const mVrMax = 2.0;
  float const mVzVpdDiffMax = 10.0;
  int const mMatchedToFMin = 2;

  // track cut
  float const mSigScaleMap[NumBeamEnergy] = {1.0,1.0,1.0,1.0,1.0};
  float const mDcaEPMax[NumBeamEnergy] = {1.0,1.0,1.0,1.0,1.0}; // for event plane reconstruction: 1.0 for BES-II
  float const mDcaTrMax = 1.0; // for pion, kaon, proton mDcaTrMax = 1.0 for flow
  float const mDcaTrMax_pid[3] = {3.0,2.0,2.0}; // cuts for daughter pions and kaons to make mother particles
  float const mDcaTrMax_phi = 3.0;
  float const mDcaTrMax_KStar = 2.0; // for phi meson mDcaTrMax = 2.0 to fill a tree and apply an additional cut
  int const mHitsDedxMin = 5;
  int const mHitsFitTPCMin = 15;
  int const mHitsMaxTPCMin = 0;
  float const mHitsRatioTPCMin = 0.52;
  float const mEtaMax = 1.0;
  float const mPrimPtMin[NumBeamEnergy] = {0.15,0.15,0.15,0.15,0.15}; // for event plane reconstruction and for pion, kaon, proton
  float const mGlobPtMin = 0.1; // for phi, Lambda, K0s
  float const mGlobPtMax = 10.0; // for phi, Lambda, K0s
  float const mPrimPtMax = 2.0;
  float const mPrimPtWeight = 2.0;
  float const mPrimMomMax = 10.0; // also use for gMom
  float const mMass2Min = -10.0;
  // double const MAGFIELDFACTOR = kilogauss;

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
 

  int const mTrackMin = 2;
  int const mTrackMin_Full = 4;
  float const mToFYLocalMax = 1.8;
  float const mToFZLocalMax = 1.8;
  float const mNSigmaElectronMax = 2.5;
  float const mNSigmaPionMax[3] = {0.0,2.0,2.0};
  float const mNSigmaKaonMax[3] = {3.0,0.0,2.0};
  float const mNSigmaProtonMax = 2.5;
  float const mNSigmaMax[3] = {3.0,2.0,2.0};
  float const mMassPion = 0.13957;
  float const mMassKaon = 0.49368;
  float const mMassProton = 0.93827;
  float const mSigKaon = 2.5;
  float const mNSigmaToF = 0.6;

  // used constapt
  // float const mEta_Gap[4] = {0.05,0.10,0.20,0.50};
  float const mEta_Gap = 0.05;
  float const mShiftOrder[5] = {2.0, 4.0, 6.0, 8.0, 10.0};

  int const pt_total = 24; // pt bin
  int const pt_start = 0;
  int const pt_stop  = 24;
  float const ptRawStart[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5};
  float const ptRawStop[pt_total]  = {0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0};
  double const pt_bin[pt_total+1]  = {0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0};

  // mix event
  int const Bin_Centrality = 9;
  int const Bin_VertexZ = 10;
  int const Bin_Phi_Psi = 5;
  int const Bin_Rho_Psi = 5;
  int const Bin_KStar_Psi = 5;
  int const Buffer_depth = 5;
  int const Buffer_depth_KStar = 3; 
  TString const MixEvent[2] = {"SE","ME"};

  TString const vm_tree[3]  = {"PhiMesonEvent","RhoMesonEvent","KStarMesonEvent"};
  TString const vm_branch[3] = {"phi_SpinAlignment_branch","rho_SpinAlignment_branch","KStar_SpinAlignment_branch"};
  int const mList_Delta = 20;
  //--------------------------------------------------

  // used in CalSpinAlginment
  int const pt_rebin = 9; // maximum pt binning
  float const pt_low[NumBeamEnergy][pt_rebin] = {
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2},
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2},
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2},
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2},
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2}
  };
  float const pt_up[NumBeamEnergy][pt_rebin]  = { 
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0},
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0},
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0},
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0},
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0}
  };
  int const pt_rebin_start[NumBeamEnergy][pt_rebin] = {
    {1,3,5, 8,11,14,17,20,24},
    {1,3,5, 8,11,14,17,20,24},
    {1,3,5, 8,11,14,17,20,24},
    {1,3,5, 8,11,14,17,20,24},
    {1,3,5, 8,11,14,17,20,24}
  };
  int const pt_rebin_stop[NumBeamEnergy][pt_rebin]  = {
    {2,4,7,10,13,16,19,23,24},
    {2,4,7,10,13,16,19,23,24},
    {2,4,7,10,13,16,19,23,24},
    {2,4,7,10,13,16,19,23,24},
    {2,4,7,10,13,16,19,23,24}
  };
  int const pt_rebin_first[NumBeamEnergy] = {0,0,0,0,0};
  int const pt_rebin_last[NumBeamEnergy]  = {8,6,6,6,6};
  int const pt_QA[NumBeamEnergy]    = {1,1,2,2,3};
  int const pt_RawQA[NumBeamEnergy]    = {2,4,6,3,10};

  std::string const Centrality[9] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%"}; // Centrality bin
  int const Cent_Total = 10; // centrality: 9 = 20%-60%, 0-8 from RefMultCorr
  int const Cent_start = 0;
  int const Cent_stop  = 10;
  int const cent_low[5] = {2,0,7,4,0}; // 0 = 20%-60%, 1 = 0-80%, 2 = 0-10%, 3 = 10-40%, 4 = 40-80%
  int const cent_up[5]  = {5,8,8,6,3}; // 0 = 20%-60%, 1 = 0-80%, 2 = 0-10%, 3 = 10-40%, 4 = 40-80%
  int const Cent_QA    = 0;

  double const Pi = TMath::Pi();
  double const PiOver10 = Pi/10.0;
  double const PiOver7 = Pi/7.0;
  int const PhiPsi_total = 7;
  double const PhiPsi_low[10] = {0.0,1.0*PiOver7/2.0, 2.0*PiOver7/2.0, 3.0*PiOver7/2.0, 4.0*PiOver7/2.0, 5.0*PiOver7/2.0, 6.0*PiOver7/2.0};
  double const PhiPsi_up[10]  = {1.0*PiOver7/2.0, 2.0*PiOver7/2.0, 3.0*PiOver7/2.0, 4.0*PiOver7/2.0, 5.0*PiOver7/2.0, 6.0*PiOver7/2.0, 7.0*PiOver7/2.0};

  int const CTS_total = 7; // cos(theta*) bin
  int const CTS_start = 0;
  int const CTS_stop  = 7;
  float const CTS_low[7] = {0.0/7.0,1.0/7.0,2.0/7.0,3.0/7.0,4.0/7.0,5.0/7.0,6.0/7.0};
  float const CTS_up[7]  = {1.0/7.0,2.0/7.0,3.0/7.0,4.0/7.0,5.0/7.0,6.0/7.0,7.0/7.0};

  int const CTS_total_KS = 5; // cos(theta*) bin
  int const CTS_start_KS = 0;
  int const CTS_stop_KS  = 5;
  float const CTS_low_KS[5] = {0.0/5.0,1.0/5.0,2.0/5.0,3.0/5.0,4.0/5.0};
  float const CTS_up_KS[5]  = {1.0/5.0,2.0/5.0,3.0/5.0,4.0/5.0,5.0/5.0};

  int const Dca_start = 0;
  int const Dca_stop = 3;
  float const mDcaSys[3] = {2.0,2.5,3.0}; // dca sys. errors
  float const mDcaSysKS[3] = {2.0,1.6,1.8}; // dca sys. errors
  float const mDcaSys_pid[3][3] = {{2.0,2.5,3.0},
                                   {2.0,1.6,1.8},
                                   {2.0,1.6,1.8}};
  
  int const mNHit_start = 0;
  int const mNHit_stop = 3;
  int const mNHitSys[3] = {15,16,18};
  
  int const mNHitSys_pid[3][3] = {{15,15,15},
                                  {15,16,18},
                                  {15,16,18}};

  int const nSigKaon_start = 0;
  int const nSigKaon_stop = 3;
  float const mNSigmaKaonSys[3] = {2.5,2.0,3.0}; // nSigKaon sys. errors
  float const mNSigmaKaonSysKStar[3] = {2.0,1.6,1.8}; // nSigKaon sys. errors

  float const mNSig_start = 0;
  float const mNSig_stop = 3;
  float const mNSigmaSys_pid[3][3] = {{2.5,2.0,3.0},
                                      {2.0,1.6,1.8},
                                      {2.0,1.6,1.8}};

  int const nSigPion_start = 0;
  int const nSigPion_stop = 3;
  float const mNSigmaPionSys[3] = {2.0,1.6,1.8}; // nSigPion sys. errors

  int const EtaGap_total = 4; // EtaGap bin
  int const Eta_start = 0; // EtaGap bin
  int const Eta_stop  = 1;
  int const Eta_QA    = 0;

  int const Charge_total = 2;
  int const Charge_start = 0;
  int const Charge_stop  = 2;

  int const Sys_start = 0; // SysError bin
  int const Sys_stop  = 1;
  int const Sys_QA    = 0;

  int const Norm_start = 0;
  int const Norm_stop  = 3;
  int const Norm_QA    = 0;

  std::string const mInteMethod[2] = {"Count","Inte"};
  int const Method_start = 0;
  int const Method_stop  = 2;
  int const Method_QA    = 0;

  float const nSigVecSys[6] = {2.0,2.5,3.0,1.5,1.0,0.5};
  int const Sig_start = 0;
  int const Sig_stop  = 6;
  int const Sig_QA    = 0;

  // used for systematic errors
  int const FuncParNum[4] = {5,6,6,6};
  int const Func_start = 0;
  int const Func_stop  = 1;
  int const Func_QA    = 0;

  // shared constant
  std::string const mBeamEnergy[NumBeamEnergy] = {"7GeV","9GeV","11GeV","14GeV","19GeV"};
  float const mEnergyValue[NumBeamEnergy] = {7.7,9.1,11.5,14.6,19.6};
  int const mBeamYear[NumBeamEnergy] = {2019,2019,2019,2019,2019};

  std::string const mPID[3]   = {"Phi","Rho","KStar"};
  float const Norm_Start[3][2]  = {{1.04,0.99},{0.41,0.30},{0.41,0.30}}; // normalise to right and left
  float const Norm_Stop[3][2]   = {{1.05,1.00},{0.46,0.31},{0.46,0.31}};
  float const BW_Start[3]     = {0.994,1.0,1.0};
  // float const BW_Start[3]     = {0.99,1.0,1.0}; // for RooFit
  float const BW_Stop[3]      = {1.050,1.0,1.0};
  float const Width[3]        = {0.00426,0.0487,0.0487};
  float const InvMass_low[3]  = {0.98,0.4,0.6};
  float const InvMass_high[3] = {1.05,0.6,1.20};
  float const nSigVec = 2.0;

  float const ptEffMax = 8.0;
  float const ptMin = 0.2;
  float const ptMax = 5.0;
  int const BinPt  = 80; // for efficiency
  // int const BinPt  = 20;
  int const BinEta = 10;
  int const BinY = 30;
  int const BinPhi = 12;
  int const BinCos = 7;

  // used in McPhiResCorr
  double const acceptanceRapidity = 1.0;
  float const rhoDelta = 0.01;
  float const InvMass[3] = {1.01946,0.89594,0.49761}; // 0: phi, 1: K*, 2 K0S
  int const decayChannelsFirst[3] = {656,617,613};
  int const decayChannelsSecond[3] = {666,619,614};
  int const decayMother[3] = {333,313,310};
  int const decayChannels[3] = {656,617,613}; // 0: phi->K+K-, 1: K*->Kpi, 2 K0S->pi+pi-
  float const McEtaBin[20] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,3.5,4.0};

  // used in RcPhiEffCorr
  std::string const mParType[2] = {"Kplus","Kminus"};
  std::string const mYear[2] = {"run11","run10"};
  std::string const mCuts[2] = {"pr","gl"};
  int const NCentMax = 9; 
  float const weight[NCentMax] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.5};

  // plotting
  int const Color[pt_rebin] = {kGray+2,kBlack,kRed,kCyan,kMagenta,kAzure,kViolet,kBlue,kRed};
  int const Style[pt_rebin] = {20,21,22,23,24,25,26,28,29};

  // ZDC-SMD constant
  std::string const mEastWest[2] = {"East","West"};
  std::string const mVertHori[2] = {"Vertical","Horizontal"};
}

#endif

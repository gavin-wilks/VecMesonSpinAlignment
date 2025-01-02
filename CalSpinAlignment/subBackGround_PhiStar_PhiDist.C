#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TProfile2D.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "phi_data_constants_19GeV.h"
//#ifdef MAKECINT
//#pragma link C++ class std::map<std::string,TH1F*>+;
//#endif

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

#ifndef _SaveQA_
#define _SaveQA_  0
#endif

using namespace std;

void subBackGround_PhiStar_PhiDist(int energy = 4, int pid = 0, int year = 0, string date = "20240426", bool random3D = false, int order = 2, string etamode = "eta1_eta1", int deltaonly = 0)
{
  std::string EP[2] = {"","2nd"};
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  //string InPutFile_SE = Form("../data/Yields_Phi_SE_19GeV_20220527.root"); //original eta < 1.0
  string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_deltar_fixedagain.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  //if(energy == 3) InPutFile_SE = Form("../data/Yields_Phi_SE_14GeV_%s.root",etamode.c_str());
  if(order == 1) InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(random3D) InPutFile_SE = Form("../data/3DRandom/Yields_Phi_SE_%s_%s_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_19GeV_20220408.root"); //original eta < 1.0
  string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_deltar_fixedagain.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  //if(energy == 3) InPutFile_ME = Form("../data/Yields_Phi_ME_14GeV_%s.root",etamode.c_str());
  if(order == 1) InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(random3D) InPutFile_ME = Form("../data/3DRandom/Yields_Phi_ME_%s_%s_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  //string InPutFile_RC = "effaccfiles/Phi/19GeV/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240326_Acc_pT0p1_TPC_TOF_phis_20240327/Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0.root";  
  string folder = "CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_phis_20240426_withRho_deltar";
  string folderdelta = "CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_Acc1Rc_pT0p1_TPC_TOF_phis_20240416_delta";
  string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0.root",folder.c_str());  
  TFile *File_RC = TFile::Open(InPutFile_RC.c_str());

  cout << "Loaded All Files" << endl;

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH2FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;
  TH2FMap h_mMass_MC;
  TH2FMap h_mMass_RC;
  TH1FMap h_mMass_SE_1D, h_mMass_ME_1D, h_mMass_SM_1D;
  TH1FMap h_mMass_MC_1D;
  TH1FMap h_mMass_RC_1D;
  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++) // pt bin 
  {
    cout << "i_pt_phi = " << i_pt_phi << endl;
    for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
    {
      cout << "i_phi = " << i_phi << endl;
      string KEY_RC = Form("rcphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_RC[KEY_RC] = (TH2F*)((TH2F*)File_RC->Get(KEY_RC.c_str()))->Clone();
      cout << KEY_RC << endl;

      string KEY_MC = Form("mcphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_MC[KEY_MC] = (TH2F*)((TH2F*)File_RC->Get(KEY_MC.c_str()))->Clone();
      cout << KEY_MC << endl;

      string KEY_SE = Form("phi_phipt_%d_cos2phistarphi_%d_SE",i_pt_phi,i_phi);
      h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();
      cout << KEY_SE << endl;

      string KEY_ME = Form("phi_phipt_%d_cos2phistarphi_%d_ME",i_pt_phi,i_phi);
      h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();
      cout << KEY_ME << endl;

      string KEY_SM = Form("phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
      h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
      h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     
      cout << KEY_SM << endl;
    }

    string KEY_RC = Form("rcphi_phistarphi_phipt_%d",i_pt_phi);
    h_mMass_RC[KEY_RC] = (TH2F*)((TH2F*)File_RC->Get(KEY_RC.c_str()))->Clone();

    string KEY_MC = Form("mcphi_phistarphi_phipt_%d",i_pt_phi);
    h_mMass_MC[KEY_MC] = (TH2F*)((TH2F*)File_RC->Get(KEY_MC.c_str()))->Clone();

    string KEY_SE = Form("phi_phistarphi_phipt_%d_SE",i_pt_phi);
    h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();

    string KEY_ME = Form("phi_phistarphi_phipt_%d_ME",i_pt_phi);
    h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();

    string KEY_SM = Form("phi_phistarphi_phipt_%d",i_pt_phi);
    h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
    h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
    h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     

    KEY_RC = Form("rcphi_phistarmphi_phipt_%d",i_pt_phi);
    h_mMass_RC_1D[KEY_RC] = (TH1F*)((TH1F*)File_RC->Get(KEY_RC.c_str()))->Clone();

    KEY_MC = Form("mcphi_phistarmphi_phipt_%d",i_pt_phi);
    h_mMass_MC_1D[KEY_MC] = (TH1F*)((TH1F*)File_RC->Get(KEY_MC.c_str()))->Clone();

    KEY_SE = Form("phi_phistarmphi_phipt_%d_SE",i_pt_phi);
    h_mMass_SE_1D[KEY_SE] = (TH1F*)((TH1F*)File_SE->Get(KEY_SE.c_str()))->Clone();

    KEY_ME = Form("phi_phistarmphi_phipt_%d_ME",i_pt_phi);
    h_mMass_ME_1D[KEY_ME] = (TH1F*)((TH1F*)File_ME->Get(KEY_ME.c_str()))->Clone();

    KEY_SM = Form("phi_phistarmphi_phipt_%d",i_pt_phi);
    h_mMass_SM_1D[KEY_SM] = (TH1F*)h_mMass_SE_1D[KEY_SE]->Clone(KEY_SM.c_str());
    h_mMass_SM_1D[KEY_SM]->SetTitle(KEY_SM.c_str());
    h_mMass_SM_1D[KEY_SM]->Add(h_mMass_ME_1D[KEY_ME],-1.0);     
    
  }
  TGraMap g_mPhi_pT; 
  TGraMap g_mRcPhi_pT; 
  TGraMap g_mMcPhi_pT; 
  TGraMap g_mPhi_pT_ratio; 
  TGraMap g_mPhi_pT_mcratio; 

  string KEY_pT_ratio = Form("phiratio_pt");
  g_mPhi_pT_ratio[KEY_pT_ratio] = new TGraphAsymmErrors();

  string KEY_pT_mcratio = Form("mcphiratio_pt");
  g_mPhi_pT_mcratio[KEY_pT_mcratio] = new TGraphAsymmErrors();

  string KEY_pT = Form("phi_pt");
  g_mPhi_pT[KEY_pT] = new TGraphAsymmErrors();

  string KEY_RC_pT = Form("rcphi_pt");
  g_mRcPhi_pT[KEY_RC_pT] = new TGraphAsymmErrors();

  string KEY_MC_pT = Form("mcphi_pt");
  g_mMcPhi_pT[KEY_MC_pT] = new TGraphAsymmErrors();

  double scaling_yields[3] = {0.0};
  double scaling_rcyields[3] = {0.0};
  double scaling_mcyields[3] = {0.0};
  
  int scalingidx[12] = {0,0,0,0,0,0,1,1,1,2,2,2};

  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++) // pt bin 
  {
    double total_yield = 0.0;
    double rctotal_yield = 0.0;
    double mctotal_yield = 0.0;
    double total_errorerror = 0.0;
    double rctotal_errorerror = 0.0;
    double mctotal_errorerror = 0.0;
    for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
    {
      string KEY_SM = Form("phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      double error;
      double yield = h_mMass_SM[KEY_SM]->IntegralAndError(0,-1,0,-1,error);
      scaling_yields[scalingidx[i_pt_phi]] += yield;
      cout << "pt = " << i_pt_phi << ", scalingidx = " << scalingidx[i_pt_phi] << endl;
      total_yield += yield;
      total_errorerror += (error*error);

      string KEY_RC = Form("rcphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      double rcerror;
      double rcyield = h_mMass_RC[KEY_RC]->IntegralAndError(0,-1,0,-1,rcerror);
      scaling_rcyields[scalingidx[i_pt_phi]] += rcyield;
      rctotal_yield += rcyield;
      rctotal_errorerror += (rcerror*rcerror);

      string KEY_MC = Form("mcphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      double mcerror;
      double mcyield = h_mMass_MC[KEY_MC]->IntegralAndError(0,-1,0,-1,mcerror);
      scaling_mcyields[scalingidx[i_pt_phi]] += mcyield;
      mctotal_yield += mcyield;
      mctotal_errorerror += (mcerror*mcerror);
    }
    double pt_mean = (data_constants::phi_pt_low[i_pt_phi]+data_constants::phi_pt_high[i_pt_phi])/2.0;
    double total_error = TMath::Sqrt(total_errorerror);
    double rctotal_error = TMath::Sqrt(rctotal_errorerror);
    double mctotal_error = TMath::Sqrt(mctotal_errorerror);
    g_mPhi_pT[KEY_pT]->SetPoint(i_pt_phi,pt_mean,total_yield);
    g_mPhi_pT[KEY_pT]->SetPointError(i_pt_phi,0.0,0.0,total_error,total_error);
    g_mRcPhi_pT[KEY_RC_pT]->SetPoint(i_pt_phi,pt_mean,rctotal_yield);
    g_mRcPhi_pT[KEY_RC_pT]->SetPointError(i_pt_phi,0.0,0.0,rctotal_error,rctotal_error);
    g_mMcPhi_pT[KEY_MC_pT]->SetPoint(i_pt_phi,pt_mean,mctotal_yield);
    g_mMcPhi_pT[KEY_MC_pT]->SetPointError(i_pt_phi,0.0,0.0,mctotal_error,mctotal_error);

  }  
  // Scale RC to the data
  //g_mRcPhi_pT[KEY_RC_pT]->Scale(scaling_yields/scaling_rcyields,"y");
  //g_mMcPhi_pT[KEY_MC_pT]->Scale(scaling_yields/scaling_mcyields,"y");

  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++) // pt bin 
  {
    double pt, yield, yielderror;
    double rcyield, rcyielderror;
    double mcyield, mcyielderror;
 
    g_mPhi_pT[KEY_pT]->GetPoint(i_pt_phi,pt,yield);
    yielderror = g_mPhi_pT[KEY_pT]->GetErrorYhigh(i_pt_phi);
    g_mRcPhi_pT[KEY_RC_pT]->GetPoint(i_pt_phi,pt,rcyield);
    rcyielderror = g_mRcPhi_pT[KEY_RC_pT]->GetErrorYhigh(i_pt_phi);
    g_mMcPhi_pT[KEY_MC_pT]->GetPoint(i_pt_phi,pt,mcyield);
    mcyielderror = g_mMcPhi_pT[KEY_MC_pT]->GetErrorYhigh(i_pt_phi);



    double rcscale = scaling_yields[scalingidx[i_pt_phi]]/scaling_rcyields[scalingidx[i_pt_phi]];
    double mcscale = scaling_yields[scalingidx[i_pt_phi]]/scaling_mcyields[scalingidx[i_pt_phi]];
 
    rcyield *= rcscale;
    mcyield *= mcscale;
    rcyielderror *= rcscale;
    mcyielderror *= mcscale;

    g_mRcPhi_pT[KEY_RC_pT]->SetPoint(i_pt_phi,pt,rcyield);
    g_mRcPhi_pT[KEY_RC_pT]->SetPointError(i_pt_phi,0.0,0.0,rcyielderror,rcyielderror);
    g_mMcPhi_pT[KEY_MC_pT]->SetPoint(i_pt_phi,pt,mcyield);
    g_mMcPhi_pT[KEY_MC_pT]->SetPointError(i_pt_phi,0.0,0.0,mcyielderror,mcyielderror);

    cout << "i_pt_phi = " << i_pt_phi << ", yield = " << yield << ", rcyield = " << rcyield << endl;
    cout << "i_pt_phi = " << i_pt_phi << ", yield = " << yield << ", mcyield = " << mcyield << endl;
    cout << "i_pt_phi = " << i_pt_phi << ", yielderror = " << yielderror << ", rcyielderror = " << rcyielderror << endl;
    cout << "i_pt_phi = " << i_pt_phi << ", yielderror = " << yielderror << ", mcyielderror = " << mcyielderror << endl;
    
    double ratio, mcratio, ratioerror, mcratioerror;
    if(rcyield == 0) 
    {
      ratio = 1;
      ratioerror = 0;
    }
    if(rcyield != 0)
    {
      ratio = yield/(rcyield);
      ratioerror = ratio*TMath::Sqrt(yielderror*yielderror/yield/yield + rcyielderror*rcyielderror/rcyield/rcyield);
    }

    if(mcyield == 0) 
    {
      mcratio = 1;
      mcratioerror = 0;
    }
    if(mcyield != 0)
    {
      mcratio = yield/(mcyield);
      mcratioerror = mcratio*TMath::Sqrt(yielderror*yielderror/yield/yield + mcyielderror*mcyielderror/mcyield/mcyield);
    }

    g_mPhi_pT_ratio[KEY_pT_ratio]->SetPoint(i_pt_phi,pt,ratio);
    g_mPhi_pT_ratio[KEY_pT_ratio]->SetPointError(i_pt_phi,0.0,0.0,ratioerror,ratioerror);
    g_mPhi_pT_mcratio[KEY_pT_mcratio]->SetPoint(i_pt_phi,pt,mcratio);
    g_mPhi_pT_mcratio[KEY_pT_mcratio]->SetPointError(i_pt_phi,0.0,0.0,mcratioerror,mcratioerror);
    cout << "i_pt_phi = " << i_pt_phi << ", ratio = " << ratio << ", ratioerror = " << ratioerror << endl;
    cout << "i_pt_phi = " << i_pt_phi << ", mcratio = " << mcratio << ", mcratioerror = " << mcratioerror << endl;
  } 

  // y vs phi
  TH2FMap h_mMass_phiy;   
  TH2FMap h_mRcMass_phiy;   
  TH2FMap h_mMcMass_phiy;   
  TH2FMap h_mMass_phiy_ratio, h_mMass_phiy_mcratio;   
  TH1FMap h_mMass_phi, h_mMass_y;   
  TH1FMap h_mRcMass_phi, h_mRcMass_y;
  TH1FMap h_mMcMass_phi, h_mMcMass_y;
  TH1FMap h_mMass_phi_ratio, h_mMass_y_ratio;
  TH1FMap h_mMass_phi_mcratio, h_mMass_y_mcratio;

  // phi* vs phi
  TH2FMap h_mMass_phistarphi;   
  TH2FMap h_mRcMass_phistarphi;   
  TH2FMap h_mMass_phistarphi_ratio;   
  TH2FMap h_mMcMass_phistarphi;   
  TH2FMap h_mMass_phistarphi_mcratio;   
  TH1FMap h_mMass_phistar;
  TH1FMap h_mRcMass_phistar;
  TH1FMap h_mMass_phistar_ratio;
  TH1FMap h_mMcMass_phistar;
  TH1FMap h_mMass_phistar_mcratio;

  // phi*-phi
  TH1FMap h_mMass_phistarmphi;
  TH1FMap h_mRcMass_phistarmphi;
  TH1FMap h_mMass_phistarmphi_ratio;
  TH1FMap h_mMcMass_phistarmphi;
  TH1FMap h_mMass_phistarmphi_mcratio;

  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++) // pt bin 
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
    {
      string KEY_SM = Form("phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_phiy = Form("phi_phiy_phipt_%d",i_pt_phi);
      if(i_phi == 0) h_mMass_phiy[KEY_phiy] = (TH2F*)h_mMass_SM[KEY_SM]->Clone(KEY_phiy.c_str()); // First bins
      else h_mMass_phiy[KEY_phiy]->Add(h_mMass_SM[KEY_SM],1.0); // All other bins

      string KEY_RC = Form("rcphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_RC_phiy = Form("rcphi_phiy_phipt_%d",i_pt_phi);
      if(i_phi == 0) h_mRcMass_phiy[KEY_RC_phiy] = (TH2F*)h_mMass_RC[KEY_RC]->Clone(KEY_RC_phiy.c_str()); // First bins
      else h_mRcMass_phiy[KEY_RC_phiy]->Add(h_mMass_RC[KEY_RC],1.0); // All other bins

      string KEY_MC = Form("mcphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_MC_phiy = Form("mcphi_phiy_phipt_%d",i_pt_phi);
      if(i_phi == 0) h_mMcMass_phiy[KEY_MC_phiy] = (TH2F*)h_mMass_MC[KEY_MC]->Clone(KEY_MC_phiy.c_str()); // First bins
      else h_mMcMass_phiy[KEY_MC_phiy]->Add(h_mMass_MC[KEY_MC],1.0); // All other bins

      string KEY_phi = Form("phi_phi_phipt_%d",i_pt_phi);
      string KEY_y   = Form("phi_y_phipt_%d",i_pt_phi);
      string KEY_SM_phi = Form("phi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_SM_y   = Form("phi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mMass_phi[KEY_SM_phi] = (TH1F*) (h_mMass_SM[KEY_SM]->ProjectionX(KEY_SM_phi.c_str(),0,-1,"e"));
      h_mMass_y[KEY_SM_y]     = (TH1F*) (h_mMass_SM[KEY_SM]->ProjectionY(KEY_SM_y.c_str(),0,-1,"e"));
      if(i_phi == 9) // Last bins
      {
        h_mMass_phi[KEY_phi] = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionX(KEY_phi.c_str(),0,-1,"e"));
        h_mMass_y[KEY_y]     = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionY(KEY_y.c_str(),0,-1,"e"));
      } 

      string KEY_RC_phi = Form("rcphi_phi_phipt_%d",i_pt_phi);
      string KEY_RC_y   = Form("rcphi_y_phipt_%d",i_pt_phi);
      if(i_phi == 9) // Last bins
      {
        h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
        h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));
      } 
      KEY_RC_phi = Form("rcphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_RC_y   = Form("rcphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mMass_RC[KEY_RC]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
      h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mMass_RC[KEY_RC]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));

      string KEY_MC_phi = Form("mcphi_phi_phipt_%d",i_pt_phi);
      string KEY_MC_y   = Form("mcphi_y_phipt_%d",i_pt_phi);
      if(i_phi == 9) // Last bins
      {
        h_mMcMass_phi[KEY_MC_phi] = (TH1F*) (h_mMcMass_phiy[KEY_MC_phiy]->ProjectionX(KEY_MC_phi.c_str(),0,-1,"e"));
        h_mMcMass_y[KEY_MC_y]     = (TH1F*) (h_mMcMass_phiy[KEY_MC_phiy]->ProjectionY(KEY_MC_y.c_str(),0,-1,"e"));
      } 
      KEY_MC_phi = Form("mcphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_MC_y   = Form("mcphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mMcMass_phi[KEY_MC_phi] = (TH1F*) (h_mMass_MC[KEY_MC]->ProjectionX(KEY_MC_phi.c_str(),0,-1,"e"));
      h_mMcMass_y[KEY_MC_y]     = (TH1F*) (h_mMass_MC[KEY_MC]->ProjectionY(KEY_MC_y.c_str(),0,-1,"e"));
    }
  }

  cout << "Completed first phiy loop" << endl;

  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++) // pt bin 
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++)
    {
      /////////// DELTA //////////////////////
      string KEY_phiy_ratio = Form("phiratio_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_phiy_mcratio = Form("mcphiratio_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_phiy = Form("phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_RC_phiy = Form("rcphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_MC_phiy = Form("mcphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
       
      h_mMass_phiy_ratio[KEY_phiy_ratio] = (TH2F*)h_mMass_SM[KEY_phiy]->Clone(KEY_phiy_ratio.c_str()); // All other bins
      h_mMass_phiy_mcratio[KEY_phiy_mcratio] = (TH2F*)h_mMass_SM[KEY_phiy]->Clone(KEY_phiy_mcratio.c_str()); // All other bins
      float phiy_yield   = h_mMass_SM[KEY_phiy]->Integral(0,-1,0,-1);
      //cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
      float phiy_rcyield = h_mMass_RC[KEY_RC_phiy]->Integral(0,-1,0,-1);
      float phiy_mcyield = h_mMass_MC[KEY_MC_phiy]->Integral(0,-1,0,-1);

      h_mMass_RC[KEY_RC_phiy]->Scale(phiy_yield/phiy_rcyield);
      h_mMass_MC[KEY_MC_phiy]->Scale(phiy_yield/phiy_mcyield);
      h_mMass_phiy_ratio[KEY_phiy_ratio]->Divide(h_mMass_RC[KEY_RC_phiy]);
      h_mMass_phiy_mcratio[KEY_phiy_mcratio]->Divide(h_mMass_MC[KEY_MC_phiy]);
   
      string KEY_phi = Form("phi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_y   = Form("phi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_RC_phi = Form("rcphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_RC_y   = Form("rcphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_MC_phi = Form("mcphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_MC_y   = Form("mcphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_phi_ratio = Form("phiratio_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_phi_mcratio = Form("mcphiratio_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_y_ratio   = Form("phiratio_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_y_mcratio   = Form("mcphiratio_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mMass_phi_ratio[KEY_phi_ratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_ratio.c_str()); // All other bins
      h_mMass_phi_mcratio[KEY_phi_mcratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_mcratio.c_str()); // All other bins
      float phi_yield   = h_mMass_phi[KEY_phi]->Integral(0,-1);
      float phi_rcyield = h_mRcMass_phi[KEY_RC_phi]->Integral(0,-1);
      float phi_mcyield = h_mMcMass_phi[KEY_MC_phi]->Integral(0,-1);

      //cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
      h_mRcMass_phi[KEY_RC_phi]->Scale(phi_yield/phi_rcyield);
      h_mMcMass_phi[KEY_MC_phi]->Scale(phi_yield/phi_mcyield);
      h_mMass_phi_ratio[KEY_phi_ratio]->Divide(h_mRcMass_phi[KEY_RC_phi]);
      h_mMass_phi_mcratio[KEY_phi_mcratio]->Divide(h_mMcMass_phi[KEY_MC_phi]);

      h_mMass_y_ratio[KEY_y_ratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_ratio.c_str()); // All other bins
      h_mMass_y_mcratio[KEY_y_mcratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_mcratio.c_str()); // All other bins
      float y_yield   = h_mMass_y[KEY_y]->Integral(0,-1);
      float y_rcyield = h_mRcMass_y[KEY_RC_y]->Integral(0,-1);
      float y_mcyield = h_mMcMass_y[KEY_MC_y]->Integral(0,-1);

      //cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
      h_mRcMass_y[KEY_RC_y]->Scale(y_yield/y_rcyield);
      h_mMcMass_y[KEY_MC_y]->Scale(y_yield/y_mcyield);
      h_mMass_y_ratio[KEY_y_ratio]->Divide(h_mRcMass_y[KEY_RC_y]);
      h_mMass_y_mcratio[KEY_y_mcratio]->Divide(h_mMcMass_y[KEY_MC_y]);
    }
  }
  cout << "Completed second phiy loop" << endl;

  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++) // pt bin 
  {
    string KEY_phiy_ratio = Form("phiratio_phiy_phipt_%d",i_pt_phi);
    string KEY_phiy = Form("phi_phiy_phipt_%d",i_pt_phi);
    string KEY_RC_phiy = Form("rcphi_phiy_phipt_%d",i_pt_phi);
    string KEY_phiy_mcratio = Form("mcphiratio_phiy_phipt_%d",i_pt_phi);
    string KEY_MC_phiy = Form("mcphi_phiy_phipt_%d",i_pt_phi);

    h_mMass_phiy_ratio[KEY_phiy_ratio] = (TH2F*)h_mMass_phiy[KEY_phiy]->Clone(KEY_phiy_ratio.c_str()); // All other bins
    double phiy_yield   = h_mMass_phiy[KEY_phiy]->Integral(1,h_mMass_phiy[KEY_phiy]->GetNbinsX(),1,h_mMass_phiy[KEY_phiy]->GetNbinsY());
    //cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    double phiy_rcyield = h_mRcMass_phiy[KEY_RC_phiy]->Integral(1,h_mRcMass_phiy[KEY_RC_phiy]->GetNbinsX(),1,h_mRcMass_phiy[KEY_RC_phiy]->GetNbinsY());

    h_mRcMass_phiy[KEY_RC_phiy]->Scale(phiy_yield/phiy_rcyield);
    h_mMass_phiy_ratio[KEY_phiy_ratio]->Divide(h_mRcMass_phiy[KEY_RC_phiy]);

    h_mMass_phiy_mcratio[KEY_phiy_mcratio] = (TH2F*)h_mMass_phiy[KEY_phiy]->Clone(KEY_phiy_mcratio.c_str()); // All other bins
    //cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    double phiy_mcyield = h_mMcMass_phiy[KEY_MC_phiy]->Integral(1,h_mMcMass_phiy[KEY_MC_phiy]->GetNbinsX(),1,h_mMcMass_phiy[KEY_MC_phiy]->GetNbinsY());

    h_mMcMass_phiy[KEY_MC_phiy]->Scale(phiy_yield/phiy_mcyield);
    h_mMass_phiy_mcratio[KEY_phiy_mcratio]->Divide(h_mMcMass_phiy[KEY_MC_phiy]);

    string KEY_phi = Form("phi_phi_phipt_%d",i_pt_phi);
    string KEY_y   = Form("phi_y_phipt_%d",i_pt_phi);
    string KEY_RC_phi = Form("rcphi_phi_phipt_%d",i_pt_phi);
    string KEY_RC_y   = Form("rcphi_y_phipt_%d",i_pt_phi);
    string KEY_phi_ratio = Form("phiratio_phi_phipt_%d",i_pt_phi);
    string KEY_y_ratio   = Form("phiratio_y_phipt_%d",i_pt_phi);
    string KEY_MC_phi = Form("mcphi_phi_phipt_%d",i_pt_phi);
    string KEY_MC_y   = Form("mcphi_y_phipt_%d",i_pt_phi);
    string KEY_phi_mcratio = Form("mcphiratio_phi_phipt_%d",i_pt_phi);
    string KEY_y_mcratio   = Form("mcphiratio_y_phipt_%d",i_pt_phi);

    h_mMass_phi_ratio[KEY_phi_ratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_ratio.c_str()); // All other bins
    double phi_yield   = h_mMass_phi[KEY_phi]->Integral(1,h_mMass_phi[KEY_phi]->GetNbinsX());
    double phi_rcyield = h_mRcMass_phi[KEY_RC_phi]->Integral(1,h_mRcMass_phi[KEY_RC_phi]->GetNbinsX());

    //cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
    h_mRcMass_phi[KEY_RC_phi]->Scale(phi_yield/phi_rcyield);
    h_mMass_phi_ratio[KEY_phi_ratio]->Divide(h_mRcMass_phi[KEY_RC_phi]);

    h_mMass_y_ratio[KEY_y_ratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_ratio.c_str()); // All other bins
    double y_yield   = h_mMass_y[KEY_y]->Integral(1,h_mMass_y[KEY_y]->GetNbinsX());
    double y_rcyield = h_mRcMass_y[KEY_RC_y]->Integral(1,h_mRcMass_y[KEY_RC_y]->GetNbinsX());

    //cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
    h_mRcMass_y[KEY_RC_y]->Scale(y_yield/y_rcyield);
    h_mMass_y_ratio[KEY_y_ratio]->Divide(h_mRcMass_y[KEY_RC_y]);

    h_mMass_phi_mcratio[KEY_phi_mcratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_mcratio.c_str()); // All other bins
    double phi_mcyield = h_mMcMass_phi[KEY_MC_phi]->Integral(1,h_mMcMass_phi[KEY_MC_phi]->GetNbinsX());

    //cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
    h_mMcMass_phi[KEY_MC_phi]->Scale(phi_yield/phi_mcyield);
    h_mMass_phi_mcratio[KEY_phi_mcratio]->Divide(h_mMcMass_phi[KEY_MC_phi]);

    h_mMass_y_mcratio[KEY_y_mcratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_mcratio.c_str()); // All other bins
    double y_mcyield = h_mMcMass_y[KEY_MC_y]->Integral(1,h_mMcMass_y[KEY_MC_y]->GetNbinsX());

    //cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
    h_mMcMass_y[KEY_MC_y]->Scale(y_yield/y_mcyield);
    h_mMass_y_mcratio[KEY_y_mcratio]->Divide(h_mMcMass_y[KEY_MC_y]);
  } 
  cout << "Completed third phiy loop" << endl;


  //// phi* vs phi

  float yields_largerbins[4] = {0.0};
  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++) // pt bin 
  {
    string KEY_phistarphi_ratio = Form("phiratio_phistarphi_phipt_%d",i_pt_phi);
    string KEY_phistarphi = Form("phi_phistarphi_phipt_%d",i_pt_phi);
    string KEY_RC_phistarphi = Form("rcphi_phistarphi_phipt_%d",i_pt_phi);
    string KEY_phistarphi_mcratio = Form("mcphiratio_phistarphi_phipt_%d",i_pt_phi);
    string KEY_MC_phistarphi = Form("mcphi_phistarphi_phipt_%d",i_pt_phi);

    h_mMass_phistarphi_ratio[KEY_phistarphi_ratio] = (TH2F*)h_mMass_SM[KEY_phistarphi]->Clone(KEY_phistarphi_ratio.c_str()); // All other bins
    double phistarphi_yield   = h_mMass_SM[KEY_phistarphi]->Integral(0,-1,0,-1);
    //cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phistarphi_yield << endl;
    double phistarphi_rcyield = h_mMass_RC[KEY_RC_phistarphi]->Integral(0,-1,0,-1);
    yields_largerbins[int(i_pt_phi)/3] += phistarphi_yield;
 
    h_mMass_RC[KEY_RC_phistarphi]->Scale(phistarphi_yield/phistarphi_rcyield);
    h_mMass_phistarphi_ratio[KEY_phistarphi_ratio]->Divide(h_mMass_RC[KEY_RC_phistarphi]);

    h_mMass_phistarphi_mcratio[KEY_phistarphi_mcratio] = (TH2F*)h_mMass_SM[KEY_phistarphi]->Clone(KEY_phistarphi_mcratio.c_str()); // All other bins
    //cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phistarphi_yield << endl;
    double phistarphi_mcyield = h_mMass_MC[KEY_MC_phistarphi]->Integral(0,-1,0,-1);

    h_mMass_MC[KEY_MC_phistarphi]->Scale(phistarphi_yield/phistarphi_mcyield);
    h_mMass_phistarphi_mcratio[KEY_phistarphi_mcratio]->Divide(h_mMass_MC[KEY_MC_phistarphi]);

    string KEY_phistar   = Form("phi_phistar_phipt_%d",i_pt_phi);
    string KEY_RC_phistar   = Form("rcphi_phistar_phipt_%d",i_pt_phi);
    string KEY_phistar_ratio   = Form("phiratio_phistar_phipt_%d",i_pt_phi);
    string KEY_MC_phistar   = Form("mcphi_phistar_phipt_%d",i_pt_phi);
    string KEY_phistar_mcratio   = Form("mcphiratio_phistar_phipt_%d",i_pt_phi);
     
    h_mMass_phistar[KEY_phistar] = (TH1F*) h_mMass_SM[KEY_phistarphi]->ProjectionY(KEY_phistar.c_str(),0,-1,"e");
    h_mMass_phistar_ratio[KEY_phistar_ratio] = (TH1F*)h_mMass_phistar[KEY_phistar]->Clone(KEY_phistar_ratio.c_str()); // All other bins
    double phistar_yield   = h_mMass_phistar[KEY_phistar]->Integral(1,h_mMass_phistar[KEY_phistar]->GetNbinsX());
    h_mRcMass_phistar[KEY_RC_phistar] = (TH1F*) h_mMass_RC[KEY_RC_phistarphi]->ProjectionY(KEY_RC_phistar.c_str(),0,-1,"e");
    double phistar_rcyield = h_mRcMass_phistar[KEY_RC_phistar]->Integral(1,h_mRcMass_phistar[KEY_RC_phistar]->GetNbinsX());

    //cout << "Yield in pt bin y " << i_pt_phi << " is = " << phistar_yield << endl;
    h_mRcMass_phistar[KEY_RC_phistar]->Scale(phistar_yield/phistar_rcyield);
    h_mMass_phistar_ratio[KEY_phistar_ratio]->Divide(h_mRcMass_phistar[KEY_RC_phistar]);

    h_mMcMass_phistar[KEY_MC_phistar] = (TH1F*) h_mMass_MC[KEY_MC_phistarphi]->ProjectionY(KEY_MC_phistar.c_str(),0,-1,"e");
    h_mMass_phistar_mcratio[KEY_phistar_mcratio] = (TH1F*)h_mMass_phistar[KEY_phistar]->Clone(KEY_phistar_mcratio.c_str()); // All other bins
    double phistar_mcyield = h_mMcMass_phistar[KEY_MC_phistar]->Integral(1,h_mMcMass_phistar[KEY_MC_phistar]->GetNbinsX());

    //cout << "Yield in pt bin y " << i_pt_phi << " is = " << phistar_yield << endl;
    h_mMcMass_phistar[KEY_MC_phistar]->Scale(phistar_yield/phistar_mcyield);
    h_mMass_phistar_mcratio[KEY_phistar_mcratio]->Divide(h_mMcMass_phistar[KEY_MC_phistar]);
  } 
  cout << "Completed phistar loop" << endl;
 
  for(int ipt = 0; ipt < 4; ipt++)
  {
    cout << "ipt = " << ipt << ", yield of phi-mesons = " << yields_largerbins[ipt] << endl;
  }

  // phi*-phi
  for(Int_t i_pt_phi = 0; i_pt_phi < 12; i_pt_phi++) // pt bin 
  {
    string KEY_phistarmphi_ratio = Form("phiratio_phistarmphi_phipt_%d",i_pt_phi);
    string KEY_phistarmphi = Form("phi_phistarmphi_phipt_%d",i_pt_phi);
    string KEY_RC_phistarmphi = Form("rcphi_phistarmphi_phipt_%d",i_pt_phi);
    string KEY_phistarmphi_mcratio = Form("mcphiratio_phistarmphi_phipt_%d",i_pt_phi);
    string KEY_MC_phistarmphi = Form("mcphi_phistarmphi_phipt_%d",i_pt_phi);

    h_mMass_phistarmphi_ratio[KEY_phistarmphi_ratio] = (TH1F*)h_mMass_SM_1D[KEY_phistarmphi]->Clone(KEY_phistarmphi_ratio.c_str()); // All other bins
    double phistarmphi_yield   = h_mMass_SM_1D[KEY_phistarmphi]->Integral(0,-1);
    //cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phistarmphi_yield << endl;
    double phistarmphi_rcyield = h_mMass_RC_1D[KEY_RC_phistarmphi]->Integral(0,-1);

    h_mMass_RC_1D[KEY_RC_phistarmphi]->Scale(phistarmphi_yield/phistarmphi_rcyield);
    h_mMass_phistarmphi_ratio[KEY_phistarmphi_ratio]->Divide(h_mMass_RC_1D[KEY_RC_phistarmphi]);

    h_mMass_phistarmphi_mcratio[KEY_phistarmphi_mcratio] = (TH1F*)h_mMass_SM_1D[KEY_phistarmphi]->Clone(KEY_phistarmphi_mcratio.c_str()); // All other bins
    //cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phistarmphi_yield << endl;
    double phistarmphi_mcyield = h_mMass_MC_1D[KEY_MC_phistarmphi]->Integral(0,-1);

    h_mMass_MC_1D[KEY_MC_phistarmphi]->Scale(phistarmphi_yield/phistarmphi_mcyield);
    h_mMass_phistarmphi_mcratio[KEY_phistarmphi_mcratio]->Divide(h_mMass_MC_1D[KEY_MC_phistarmphi]);
  } 
  cout << "Completed phistar-phi loop" << endl;

  

  TCanvas *c1 = new TCanvas("c1","c1",10,10,1600,1200);     
  c1->Divide(4,3);
  for(int i = 0; i < 12; i++)
  {
    c1->cd(i+1);
    c1->cd(i+1)->SetLeftMargin(0.15);
    c1->cd(i+1)->SetBottomMargin(0.15);
    c1->cd(i+1)->SetTicks(1,1);
    c1->cd(i+1)->SetGrid(0,0);  
  } 

  TCanvas *craw = new TCanvas("craw","craw",10,10,1500,600);     
  craw->Divide(5,2);
  for(int i = 0; i < 10; i++)
  {
    craw->cd(i+1);
    craw->cd(i+1)->SetLeftMargin(0.15);
    craw->cd(i+1)->SetBottomMargin(0.15);
    craw->cd(i+1)->SetTicks(1,1);
    craw->cd(i+1)->SetGrid(0,0);  
  } 

  std::string phiphilabel[4] = {"#phi","#phi_{K+}-#phi_{K-}","#phi_{#phi}-#phi_{K+}","#phi_{#phi}-#phi_{K-}"};
  std::string phietalabel[4] = {"y","#eta_{K+}-#eta_{K-}","#eta_{#phi}-#eta_{K+}","#eta_{#phi}-#eta_{K-}"};
  std::string philabel[4] = {"","","plusphi","minusphi"};
  std::string phiextralabel[4] = {"","_delta","_delta","_delta"};
  

  for(int iset = 0; iset < 1; iset++)
  {
    std::string outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_SE.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    std::string outputstart = Form("%s[",outputname.c_str());
    std::string outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_SE = Form("phi%s%s_phipt_%d_cos2phistarphi_%d_SE",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_SE[KEY_SE]->SetTitle(Form("SE %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_SE[KEY_SE]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_SE[KEY_SE]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

        h_mMass_SE[KEY_SE]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_ME.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("phi%s%s_phipt_%d_cos2phistarphi_%d_ME",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_ME[KEY_ME]->SetTitle(Form("ME %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_ME[KEY_ME]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_ME[KEY_ME]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

        h_mMass_ME[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_SM.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("phi%s%s_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_SM[KEY_ME]->SetTitle(Form("SE-ME %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_SM[KEY_ME]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_SM[KEY_ME]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

        h_mMass_SM[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_IndividualCos2PhiStarPhiBins_RC.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("rcphi%s%s_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_RC[KEY_ME]->SetTitle(Form("RC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_RC[KEY_ME]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_RC[KEY_ME]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

        h_mMass_RC[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phi_IndividualCos2PhiStarPhiBins_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      TLegend *legphi = new TLegend(0.825,0.55,0.975,0.75);
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_phi = Form("phi%s%s_phi_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        string KEY_RC_phi = Form("rcphi%s%s_phi_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_phi[KEY_phi]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phi[KEY_phi]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_phi[KEY_phi]->GetYaxis()->SetTitle("Counts");
        h_mMass_phi[KEY_phi]->SetMarkerStyle(20);
        h_mMass_phi[KEY_phi]->SetMarkerColor(kOrange+7);
        h_mMass_phi[KEY_phi]->SetLineColor(kOrange+7);

        int min = h_mMass_phi[KEY_phi]->GetMinimum();
        int max = h_mMass_phi[KEY_phi]->GetMaximum();
        if(h_mRcMass_phi[KEY_RC_phi]->GetMinimum() < min) min = h_mRcMass_phi[KEY_RC_phi]->GetMinimum();
        if(h_mRcMass_phi[KEY_RC_phi]->GetMaximum() > max) max = h_mRcMass_phi[KEY_RC_phi]->GetMaximum();
        h_mMass_phi[KEY_phi]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

        h_mMass_phi[KEY_phi]->Draw("pE");

        h_mRcMass_phi[KEY_RC_phi]->SetMarkerStyle(24);
        h_mRcMass_phi[KEY_RC_phi]->SetMarkerColor(kBlack);
        h_mRcMass_phi[KEY_RC_phi]->SetLineColor(kBlack);
        h_mRcMass_phi[KEY_RC_phi]->Draw("pE same");

        if(i_phi == 0)
        {
          legphi->AddEntry(h_mMass_phi[KEY_phi],"Data","p");
          legphi->AddEntry(h_mRcMass_phi[KEY_RC_phi],"RC","p");
        }
        legphi->Draw("same");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phi_IndividualCos2PhiStarPhiBins_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      TLegend *legphi = new TLegend(0.825,0.55,0.975,0.75);
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_phi = Form("phi%s%s_phi_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        string KEY_MC_phi = Form("mcphi%s%s_phi_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_phi[KEY_phi]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phi[KEY_phi]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_phi[KEY_phi]->GetYaxis()->SetTitle("Counts");
        h_mMass_phi[KEY_phi]->SetMarkerStyle(20);
        h_mMass_phi[KEY_phi]->SetMarkerColor(kOrange+7);
        h_mMass_phi[KEY_phi]->SetLineColor(kOrange+7);

        int min = h_mMass_phi[KEY_phi]->GetMinimum();
        int max = h_mMass_phi[KEY_phi]->GetMaximum();
        if(h_mMcMass_phi[KEY_MC_phi]->GetMinimum() < min) min = h_mMcMass_phi[KEY_MC_phi]->GetMinimum();
        if(h_mMcMass_phi[KEY_MC_phi]->GetMaximum() > max) max = h_mMcMass_phi[KEY_MC_phi]->GetMaximum();
        h_mMass_phi[KEY_phi]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

        h_mMass_phi[KEY_phi]->Draw("pE");

        h_mMcMass_phi[KEY_MC_phi]->SetMarkerStyle(24);
        h_mMcMass_phi[KEY_MC_phi]->SetMarkerColor(kBlack);
        h_mMcMass_phi[KEY_MC_phi]->SetLineColor(kBlack);
        h_mMcMass_phi[KEY_MC_phi]->Draw("pE same");

        if(i_phi == 0)
        {
          legphi->AddEntry(h_mMass_phi[KEY_phi],"Data","p");
          legphi->AddEntry(h_mMcMass_phi[KEY_MC_phi],"MC","p");
        }
        legphi->Draw("same");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_y_IndividualCos2PhiStarPhiBins_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      TLegend *legy = new TLegend(0.825,0.55,0.975,0.75);
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_y = Form("phi%s%s_y_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        string KEY_RC_y = Form("rcphi%s%s_y_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_y[KEY_y]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_y[KEY_y]->GetXaxis()->SetTitle(phietalabel[iset].c_str());
        h_mMass_y[KEY_y]->GetYaxis()->SetTitle("Counts");
        h_mMass_y[KEY_y]->SetMarkerStyle(20);
        h_mMass_y[KEY_y]->SetMarkerColor(kOrange+7);
        h_mMass_y[KEY_y]->SetLineColor(kOrange+7);

        int min = h_mMass_y[KEY_y]->GetMinimum();
        int max = h_mMass_y[KEY_y]->GetMaximum();
        if(h_mRcMass_y[KEY_RC_y]->GetMinimum() < min) min = h_mRcMass_y[KEY_RC_y]->GetMinimum();
        if(h_mRcMass_y[KEY_RC_y]->GetMaximum() > max) max = h_mRcMass_y[KEY_RC_y]->GetMaximum();
        h_mMass_y[KEY_y]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

        h_mMass_y[KEY_y]->Draw("pE");

        h_mRcMass_y[KEY_RC_y]->SetMarkerStyle(24);
        h_mRcMass_y[KEY_RC_y]->SetMarkerColor(kBlack);
        h_mRcMass_y[KEY_RC_y]->SetLineColor(kBlack);
        h_mRcMass_y[KEY_RC_y]->Draw("pE same");

        if(i_phi == 0)
        {
          legy->AddEntry(h_mMass_y[KEY_y],"Data","p");
          legy->AddEntry(h_mRcMass_y[KEY_RC_y],"RC","p");
        }
        legy->Draw("same");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_y_IndividualCos2PhiStarPhiBins_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      TLegend *legy = new TLegend(0.825,0.55,0.975,0.75);
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_y = Form("phi%s%s_y_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        string KEY_MC_y = Form("mcphi%s%s_y_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_y[KEY_y]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_y[KEY_y]->GetXaxis()->SetTitle(phietalabel[iset].c_str());
        h_mMass_y[KEY_y]->GetYaxis()->SetTitle("Counts");
        h_mMass_y[KEY_y]->SetMarkerStyle(20);
        h_mMass_y[KEY_y]->SetMarkerColor(kOrange+7);
        h_mMass_y[KEY_y]->SetLineColor(kOrange+7);

        int min = h_mMass_y[KEY_y]->GetMinimum();
        int max = h_mMass_y[KEY_y]->GetMaximum();
        if(h_mMcMass_y[KEY_MC_y]->GetMinimum() < min) min = h_mMcMass_y[KEY_MC_y]->GetMinimum();
        if(h_mMcMass_y[KEY_MC_y]->GetMaximum() > max) max = h_mMcMass_y[KEY_MC_y]->GetMaximum();
        h_mMass_y[KEY_y]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

        h_mMass_y[KEY_y]->Draw("pE");

        h_mMcMass_y[KEY_MC_y]->SetMarkerStyle(24);
        h_mMcMass_y[KEY_MC_y]->SetMarkerColor(kBlack);
        h_mMcMass_y[KEY_MC_y]->SetLineColor(kBlack);
        h_mMcMass_y[KEY_MC_y]->Draw("pE same");

        if(i_phi == 0)
        {
          legy->AddEntry(h_mMass_y[KEY_y],"Data","p");
          legy->AddEntry(h_mMcMass_y[KEY_MC_y],"MC","p");
        }
        legy->Draw("same");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phi_IndividualCos2PhiStarPhiBins_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("phiratio%s%s_phi_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_phi_ratio[KEY_ME]->SetTitle(Form("Data/RC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phi_ratio[KEY_ME]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_phi_ratio[KEY_ME]->GetYaxis()->SetTitle("Data/RC");

        h_mMass_phi_ratio[KEY_ME]->Draw("pE");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_y_IndividualCos2PhiStarPhiBins_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("phiratio%s%s_y_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_y_ratio[KEY_ME]->SetTitle(Form("Data/RC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_y_ratio[KEY_ME]->GetXaxis()->SetTitle(phietalabel[iset].c_str());
        h_mMass_y_ratio[KEY_ME]->GetYaxis()->SetTitle("Data/RC");

        h_mMass_y_ratio[KEY_ME]->Draw("pE");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_IndividualCos2PhiStarPhiBins_MC.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("mcphi%s%s_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_MC[KEY_ME]->SetTitle(Form("MC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_MC[KEY_ME]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_MC[KEY_ME]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

        h_mMass_MC[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phi_IndividualCos2PhiStarPhiBins_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("mcphiratio%s%s_phi_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_phi_mcratio[KEY_ME]->SetTitle(Form("Data/MC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phi_mcratio[KEY_ME]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_phi_mcratio[KEY_ME]->GetYaxis()->SetTitle("Data/MC");

        h_mMass_phi_mcratio[KEY_ME]->Draw("pE");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_y_IndividualCos2PhiStarPhiBins_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("mcphiratio%s%s_y_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_y_mcratio[KEY_ME]->SetTitle(Form("Data/MC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_y_mcratio[KEY_ME]->GetXaxis()->SetTitle(phietalabel[iset].c_str());
        h_mMass_y_mcratio[KEY_ME]->GetYaxis()->SetTitle("Data/MC");

        h_mMass_y_mcratio[KEY_ME]->Draw("pE");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_IndividualCos2PhiStarPhiBins_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("phiratio%s%s_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_phiy_ratio[KEY_ME]->SetTitle(Form("Data/RC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phiy_ratio[KEY_ME]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_phiy_ratio[KEY_ME]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

        h_mMass_phiy_ratio[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_IndividualCos2PhiStarPhiBins_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 12; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("mcphiratio%s%s_phipt_%d_cos2phistarphi_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i,i_phi);
        h_mMass_phiy_mcratio[KEY_ME]->SetTitle(Form("Data/MC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phiy_mcratio[KEY_ME]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
        h_mMass_phiy_mcratio[KEY_ME]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

        h_mMass_phiy_mcratio[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());
  }


  /////////////////////////////// No Cos2PhiStarPhi separation //////////////////////////////
  for(int iset = 0; iset < 1; iset++)
  {
    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phistarphi_ratio = Form("phi%sratio%s_phistarphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phistarphi_ratio[KEY_phistarphi_ratio]->SetTitle(Form("Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phistarphi_ratio[KEY_phistarphi_ratio]->GetXaxis()->SetTitle("#phi");
      h_mMass_phistarphi_ratio[KEY_phistarphi_ratio]->GetYaxis()->SetTitle("#phi*");

      h_mMass_phistarphi_ratio[KEY_phistarphi_ratio]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistarphi_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phistarphi = Form("phi%s%s_phistarphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_SM[KEY_phistarphi]->SetTitle(Form("Data %1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_SM[KEY_phistarphi]->GetXaxis()->SetTitle("#phi");
      h_mMass_SM[KEY_phistarphi]->GetYaxis()->SetTitle("phi*");

      h_mMass_SM[KEY_phistarphi]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistarphi_data.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_RC_phistarphi = Form("rcphi%s%s_phistarphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_RC[KEY_RC_phistarphi]->SetTitle(Form("RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_RC[KEY_RC_phistarphi]->GetXaxis()->SetTitle("phi");
      h_mMass_RC[KEY_RC_phistarphi]->GetYaxis()->SetTitle("phi*");

      h_mMass_RC[KEY_RC_phistarphi]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistarphi_rc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phiy_ratio = Form("phi%sratio%s_phiy_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phiy_ratio[KEY_phiy_ratio]->SetTitle(Form("Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phiy_ratio[KEY_phiy_ratio]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
      h_mMass_phiy_ratio[KEY_phiy_ratio]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

      h_mMass_phiy_ratio[KEY_phiy_ratio]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phiy = Form("phi%s%s_phiy_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phiy[KEY_phiy]->SetTitle(Form("Data %1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phiy[KEY_phiy]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
      h_mMass_phiy[KEY_phiy]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

      h_mMass_phiy[KEY_phiy]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_data.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_RC_phiy = Form("rcphi%s%s_phiy_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mRcMass_phiy[KEY_RC_phiy]->SetTitle(Form("RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mRcMass_phiy[KEY_RC_phiy]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
      h_mRcMass_phiy[KEY_RC_phiy]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

      h_mRcMass_phiy[KEY_RC_phiy]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_rc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phistar_ratio = Form("phi%sratio%s_phistar_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phistar_ratio[KEY_phistar_ratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phistar_ratio[KEY_phistar_ratio]->GetXaxis()->SetTitle("phi*");
      h_mMass_phistar_ratio[KEY_phistar_ratio]->GetYaxis()->SetTitle("Data/RC");
      h_mMass_phistar_ratio[KEY_phistar_ratio]->SetMarkerStyle(20);

      h_mMass_phistar_ratio[KEY_phistar_ratio]->Draw("pE");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistar_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    TLegend *legphistar = new TLegend(0.825,0.55,0.975,0.75);
    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phistar = Form("phi%s%s_phistar_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      string KEY_RC_phistar = Form("rcphi%s%s_phistar_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phistar[KEY_phistar]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phistar[KEY_phistar]->GetXaxis()->SetTitle("#phi*");
      h_mMass_phistar[KEY_phistar]->GetYaxis()->SetTitle("Counts");
      h_mMass_phistar[KEY_phistar]->SetMarkerStyle(20);
      h_mMass_phistar[KEY_phistar]->SetMarkerColor(kOrange+7);
      h_mMass_phistar[KEY_phistar]->SetLineColor(kOrange+7);

      int min = h_mMass_phistar[KEY_phistar]->GetMinimum();
      int max = h_mMass_phistar[KEY_phistar]->GetMaximum();
      if(h_mRcMass_phistar[KEY_RC_phistar]->GetMinimum() < min) min = h_mRcMass_phistar[KEY_RC_phistar]->GetMinimum();
      if(h_mRcMass_phistar[KEY_RC_phistar]->GetMaximum() > max) max = h_mRcMass_phistar[KEY_RC_phistar]->GetMaximum();
      h_mMass_phistar[KEY_phistar]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

      h_mMass_phistar[KEY_phistar]->Draw("pE");

      h_mRcMass_phistar[KEY_RC_phistar]->SetMarkerStyle(24);
      h_mRcMass_phistar[KEY_RC_phistar]->SetMarkerColor(kBlack);
      h_mRcMass_phistar[KEY_RC_phistar]->SetLineColor(kBlack);
      h_mRcMass_phistar[KEY_RC_phistar]->Draw("pE same");

      if(i == 0)
      {
        legphistar->AddEntry(h_mMass_phistar[KEY_phistar],"Data","p");
        legphistar->AddEntry(h_mRcMass_phistar[KEY_RC_phistar],"RC","p");
      }
      legphistar->Draw("same");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistar_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phistarmphi_ratio = Form("phi%sratio%s_phistarmphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phistarmphi_ratio[KEY_phistarmphi_ratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phistarmphi_ratio[KEY_phistarmphi_ratio]->GetXaxis()->SetTitle("#phi*-#phi");
      h_mMass_phistarmphi_ratio[KEY_phistarmphi_ratio]->GetYaxis()->SetTitle("Data/RC");
      h_mMass_phistarmphi_ratio[KEY_phistarmphi_ratio]->SetMarkerStyle(20);

      h_mMass_phistarmphi_ratio[KEY_phistarmphi_ratio]->Draw("pE");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistarmphi_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    TLegend *legphistarmphi = new TLegend(0.825,0.55,0.975,0.75);
    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phistarmphi = Form("phi%s%s_phistarmphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      string KEY_RC_phistarmphi = Form("rcphi%s%s_phistarmphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_SM_1D[KEY_phistarmphi]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_SM_1D[KEY_phistarmphi]->GetXaxis()->SetTitle("#phi*-#phi");
      h_mMass_SM_1D[KEY_phistarmphi]->GetYaxis()->SetTitle("Counts");
      h_mMass_SM_1D[KEY_phistarmphi]->SetMarkerStyle(20);
      h_mMass_SM_1D[KEY_phistarmphi]->SetMarkerColor(kOrange+7);
      h_mMass_SM_1D[KEY_phistarmphi]->SetLineColor(kOrange+7);

      int min = h_mMass_SM_1D[KEY_phistarmphi]->GetMinimum();
      int max = h_mMass_SM_1D[KEY_phistarmphi]->GetMaximum();
      if(h_mMass_RC_1D[KEY_RC_phistarmphi]->GetMinimum() < min) min = h_mMass_RC_1D[KEY_RC_phistarmphi]->GetMinimum();
      if(h_mMass_RC_1D[KEY_RC_phistarmphi]->GetMaximum() > max) max = h_mMass_RC_1D[KEY_RC_phistarmphi]->GetMaximum();
      h_mMass_SM_1D[KEY_phistarmphi]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

      h_mMass_SM_1D[KEY_phistarmphi]->Draw("pE");

      h_mMass_RC_1D[KEY_RC_phistarmphi]->SetMarkerStyle(24);
      h_mMass_RC_1D[KEY_RC_phistarmphi]->SetMarkerColor(kBlack);
      h_mMass_RC_1D[KEY_RC_phistarmphi]->SetLineColor(kBlack);
      h_mMass_RC_1D[KEY_RC_phistarmphi]->Draw("pE same");

      if(i == 0)
      {
        legphistarmphi->AddEntry(h_mMass_SM_1D[KEY_phistarmphi],"Data","p");
        legphistarmphi->AddEntry(h_mMass_RC_1D[KEY_RC_phistarmphi],"RC","p");
      }
      legphistarmphi->Draw("same");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistarmphi_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phi_ratio = Form("phi%sratio%s_phi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phi_ratio[KEY_phi_ratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phi_ratio[KEY_phi_ratio]->GetXaxis()->SetTitle("phi*-phi");
      h_mMass_phi_ratio[KEY_phi_ratio]->GetYaxis()->SetTitle("Data/RC");
      h_mMass_phi_ratio[KEY_phi_ratio]->SetMarkerStyle(20);

      h_mMass_phi_ratio[KEY_phi_ratio]->Draw("pE");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phi_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    TLegend *legphi = new TLegend(0.825,0.55,0.975,0.75);
    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phi = Form("phi%s%s_phi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      string KEY_RC_phi = Form("rcphi%s%s_phi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phi[KEY_phi]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phi[KEY_phi]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
      h_mMass_phi[KEY_phi]->GetYaxis()->SetTitle("Counts");
      h_mMass_phi[KEY_phi]->SetMarkerStyle(20);
      h_mMass_phi[KEY_phi]->SetMarkerColor(kOrange+7);
      h_mMass_phi[KEY_phi]->SetLineColor(kOrange+7);

      int min = h_mMass_phi[KEY_phi]->GetMinimum();
      int max = h_mMass_phi[KEY_phi]->GetMaximum();
      if(h_mRcMass_phi[KEY_RC_phi]->GetMinimum() < min) min = h_mRcMass_phi[KEY_RC_phi]->GetMinimum();
      if(h_mRcMass_phi[KEY_RC_phi]->GetMaximum() > max) max = h_mRcMass_phi[KEY_RC_phi]->GetMaximum();
      h_mMass_phi[KEY_phi]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

      h_mMass_phi[KEY_phi]->Draw("pE");

      h_mRcMass_phi[KEY_RC_phi]->SetMarkerStyle(24);
      h_mRcMass_phi[KEY_RC_phi]->SetMarkerColor(kBlack);
      h_mRcMass_phi[KEY_RC_phi]->SetLineColor(kBlack);
      h_mRcMass_phi[KEY_RC_phi]->Draw("pE same");

      if(i == 0)
      {
        legphi->AddEntry(h_mMass_phi[KEY_phi],"Data","p");
        legphi->AddEntry(h_mRcMass_phi[KEY_RC_phi],"RC","p");
      }
      legphi->Draw("same");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phi_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_y_ratio = Form("phi%sratio%s_y_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_y_ratio[KEY_y_ratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_y_ratio[KEY_y_ratio]->GetXaxis()->SetTitle(phietalabel[iset].c_str());
      h_mMass_y_ratio[KEY_y_ratio]->GetYaxis()->SetTitle("Data/RC");
      h_mMass_y_ratio[KEY_y_ratio]->SetMarkerStyle(20);

      h_mMass_y_ratio[KEY_y_ratio]->Draw("pE");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_y_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    TLegend *legy = new TLegend(0.825,0.55,0.975,0.75);
    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_y = Form("phi%s%s_y_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      string KEY_RC_y = Form("rcphi%s%s_y_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_y[KEY_y]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_y[KEY_y]->GetXaxis()->SetTitle(phietalabel[iset].c_str());
      h_mMass_y[KEY_y]->GetYaxis()->SetTitle("Counts");
      h_mMass_y[KEY_y]->SetMarkerStyle(20);
      h_mMass_y[KEY_y]->SetMarkerColor(kOrange+7);
      h_mMass_y[KEY_y]->SetLineColor(kOrange+7);

      int min = h_mMass_y[KEY_y]->GetMinimum();
      int max = h_mMass_y[KEY_y]->GetMaximum();
      if(h_mRcMass_y[KEY_RC_y]->GetMinimum() < min) min = h_mRcMass_y[KEY_RC_y]->GetMinimum();
      if(h_mRcMass_y[KEY_RC_y]->GetMaximum() > max) max = h_mRcMass_y[KEY_RC_y]->GetMaximum();
      h_mMass_y[KEY_y]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

      h_mMass_y[KEY_y]->Draw("pE");

      h_mRcMass_y[KEY_RC_y]->SetMarkerStyle(24);
      h_mRcMass_y[KEY_RC_y]->SetMarkerColor(kBlack);
      h_mRcMass_y[KEY_RC_y]->SetLineColor(kBlack);
      h_mRcMass_y[KEY_RC_y]->Draw("pE same");

      if(i == 0)
      {
        legy->AddEntry(h_mMass_y[KEY_y],"Data","p");
        legy->AddEntry(h_mRcMass_y[KEY_RC_y],"RC","p");
      }
      legy->Draw("same");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_y_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));
   
    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phiy_mcratio = Form("mcphi%sratio%s_phiy_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phiy_mcratio[KEY_phiy_mcratio]->SetTitle(Form("Data/MC %1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phiy_mcratio[KEY_phiy_mcratio]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
      h_mMass_phiy_mcratio[KEY_phiy_mcratio]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

      h_mMass_phiy_mcratio[KEY_phiy_mcratio]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phi_mcratio = Form("mcphi%sratio%s_phi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phi_mcratio[KEY_phi_mcratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phi_mcratio[KEY_phi_mcratio]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
      h_mMass_phi_mcratio[KEY_phi_mcratio]->GetYaxis()->SetTitle("Data/MC");
      h_mMass_phi_mcratio[KEY_phi_mcratio]->SetMarkerStyle(20);

      h_mMass_phi_mcratio[KEY_phi_mcratio]->Draw("pE");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phi_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    TLegend *legphimc = new TLegend(0.825,0.55,0.975,0.75);
    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phi = Form("phi%s%s_phi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      string KEY_MC_phi = Form("mcphi%s%s_phi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phi[KEY_phi]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phi[KEY_phi]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
      h_mMass_phi[KEY_phi]->GetYaxis()->SetTitle("Counts");
      h_mMass_phi[KEY_phi]->SetMarkerStyle(20);
      h_mMass_phi[KEY_phi]->SetMarkerColor(kOrange+7);
      h_mMass_phi[KEY_phi]->SetLineColor(kOrange+7);

      int min = h_mMass_phi[KEY_phi]->GetMinimum();
      int max = h_mMass_phi[KEY_phi]->GetMaximum();
      if(h_mMcMass_phi[KEY_MC_phi]->GetMinimum() < min) min = h_mMcMass_phi[KEY_MC_phi]->GetMinimum();
      if(h_mMcMass_phi[KEY_MC_phi]->GetMaximum() > max) max = h_mMcMass_phi[KEY_MC_phi]->GetMaximum();
      h_mMass_phi[KEY_phi]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

      h_mMass_phi[KEY_phi]->Draw("pE");

      h_mMcMass_phi[KEY_MC_phi]->SetMarkerStyle(24);
      h_mMcMass_phi[KEY_MC_phi]->SetMarkerColor(kBlack);
      h_mMcMass_phi[KEY_MC_phi]->SetLineColor(kBlack);
      h_mMcMass_phi[KEY_MC_phi]->Draw("pE same");

      if(i == 0)
      {
        legphimc->AddEntry(h_mMass_phi[KEY_phi],"Data","p");
        legphimc->AddEntry(h_mMcMass_phi[KEY_MC_phi],"MC","p");
      }
      legphimc->Draw("same");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phi_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_y_mcratio = Form("mcphi%sratio%s_y_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_y_mcratio[KEY_y_mcratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_y_mcratio[KEY_y_mcratio]->GetXaxis()->SetTitle(phietalabel[iset].c_str());
      h_mMass_y_mcratio[KEY_y_mcratio]->GetYaxis()->SetTitle("Data/MC");
      h_mMass_y_mcratio[KEY_y_mcratio]->SetMarkerStyle(20);

      h_mMass_y_mcratio[KEY_y_mcratio]->Draw("pE");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_y_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phistar_mcratio = Form("mcphi%sratio%s_phistar_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phistar_mcratio[KEY_phistar_mcratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phistar_mcratio[KEY_phistar_mcratio]->GetXaxis()->SetTitle("#phi*");
      h_mMass_phistar_mcratio[KEY_phistar_mcratio]->GetYaxis()->SetTitle("Data/MC");
      h_mMass_phistar_mcratio[KEY_phistar_mcratio]->SetMarkerStyle(20);

      h_mMass_phistar_mcratio[KEY_phistar_mcratio]->Draw("pE");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistar_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phistarmphi_mcratio = Form("mcphi%sratio%s_phistarmphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phistarmphi_mcratio[KEY_phistarmphi_mcratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phistarmphi_mcratio[KEY_phistarmphi_mcratio]->GetXaxis()->SetTitle("#phi*-#phi");
      h_mMass_phistarmphi_mcratio[KEY_phistarmphi_mcratio]->GetYaxis()->SetTitle("Data/MC");
      h_mMass_phistarmphi_mcratio[KEY_phistarmphi_mcratio]->SetMarkerStyle(20);

      h_mMass_phistarmphi_mcratio[KEY_phistarmphi_mcratio]->Draw("pE");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistarmphi_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    TLegend *legymc = new TLegend(0.825,0.55,0.975,0.75);
    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_y = Form("phi%s%s_y_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      string KEY_MC_y = Form("mcphi%s%s_y_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_y[KEY_y]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_y[KEY_y]->GetXaxis()->SetTitle(phietalabel[iset].c_str());
      h_mMass_y[KEY_y]->GetYaxis()->SetTitle("Counts");
      h_mMass_y[KEY_y]->SetMarkerStyle(20);
      h_mMass_y[KEY_y]->SetMarkerColor(kOrange+7);
      h_mMass_y[KEY_y]->SetLineColor(kOrange+7);

      int min = h_mMass_y[KEY_y]->GetMinimum();
      int max = h_mMass_y[KEY_y]->GetMaximum();
      if(h_mMcMass_y[KEY_MC_y]->GetMinimum() < min) min = h_mMcMass_y[KEY_MC_y]->GetMinimum();
      if(h_mMcMass_y[KEY_MC_y]->GetMaximum() > max) max = h_mMcMass_y[KEY_MC_y]->GetMaximum();
      h_mMass_y[KEY_y]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

      h_mMass_y[KEY_y]->Draw("pE");

      h_mMcMass_y[KEY_MC_y]->SetMarkerStyle(24);
      h_mMcMass_y[KEY_MC_y]->SetMarkerColor(kBlack);
      h_mMcMass_y[KEY_MC_y]->SetLineColor(kBlack);
      h_mMcMass_y[KEY_MC_y]->Draw("pE same");

      if(i == 0)
      {
        legymc->AddEntry(h_mMass_y[KEY_y],"Data","p");
        legymc->AddEntry(h_mMcMass_y[KEY_MC_y],"MC","p");
      }
      legymc->Draw("same");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_y_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    TLegend *legphistarmphimc = new TLegend(0.825,0.55,0.975,0.75);
    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phistarmphi = Form("phi%s%s_phistarmphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      string KEY_MC_phistarmphi = Form("mcphi%s%s_phistarmphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_SM_1D[KEY_phistarmphi]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_SM_1D[KEY_phistarmphi]->GetXaxis()->SetTitle("#phi*-#phi");
      h_mMass_SM_1D[KEY_phistarmphi]->GetYaxis()->SetTitle("Counts");
      h_mMass_SM_1D[KEY_phistarmphi]->SetMarkerStyle(20);
      h_mMass_SM_1D[KEY_phistarmphi]->SetMarkerColor(kOrange+7);
      h_mMass_SM_1D[KEY_phistarmphi]->SetLineColor(kOrange+7);

      int min = h_mMass_SM_1D[KEY_phistarmphi]->GetMinimum();
      int max = h_mMass_SM_1D[KEY_phistarmphi]->GetMaximum();
      if(h_mMass_MC_1D[KEY_MC_phistarmphi]->GetMinimum() < min) min = h_mMass_MC_1D[KEY_MC_phistarmphi]->GetMinimum();
      if(h_mMass_MC_1D[KEY_MC_phistarmphi]->GetMaximum() > max) max = h_mMass_MC_1D[KEY_MC_phistarmphi]->GetMaximum();
      h_mMass_SM_1D[KEY_phistarmphi]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

      h_mMass_SM_1D[KEY_phistarmphi]->Draw("pE");

      h_mMass_MC_1D[KEY_MC_phistarmphi]->SetMarkerStyle(24);
      h_mMass_MC_1D[KEY_MC_phistarmphi]->SetMarkerColor(kBlack);
      h_mMass_MC_1D[KEY_MC_phistarmphi]->SetLineColor(kBlack);
      h_mMass_MC_1D[KEY_MC_phistarmphi]->Draw("pE same");

      if(i == 0)
      {
        legphistarmphimc->AddEntry(h_mMass_SM_1D[KEY_phistarmphi],"Data","p");
        legphistarmphimc->AddEntry(h_mMass_MC_1D[KEY_MC_phistarmphi],"MC","p");
      }
      legphistarmphimc->Draw("same");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistarmphi_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    TLegend *legphistarmc = new TLegend(0.825,0.55,0.975,0.75);
    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_phistar = Form("phi%s%s_phistar_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      string KEY_MC_phistar = Form("mcphi%s%s_phistar_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_phistar[KEY_phistar]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_phistar[KEY_phistar]->GetXaxis()->SetTitle("#phi*");
      h_mMass_phistar[KEY_phistar]->GetYaxis()->SetTitle("Counts");
      h_mMass_phistar[KEY_phistar]->SetMarkerStyle(20);
      h_mMass_phistar[KEY_phistar]->SetMarkerColor(kOrange+7);
      h_mMass_phistar[KEY_phistar]->SetLineColor(kOrange+7);

      int min = h_mMass_phistar[KEY_phistar]->GetMinimum();
      int max = h_mMass_phistar[KEY_phistar]->GetMaximum();
      if(h_mMcMass_phistar[KEY_MC_phistar]->GetMinimum() < min) min = h_mMcMass_phistar[KEY_MC_phistar]->GetMinimum();
      if(h_mMcMass_phistar[KEY_MC_phistar]->GetMaximum() > max) max = h_mMcMass_phistar[KEY_MC_phistar]->GetMaximum();
      h_mMass_phistar[KEY_phistar]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

      h_mMass_phistar[KEY_phistar]->Draw("pE");

      h_mMcMass_phistar[KEY_MC_phistar]->SetMarkerStyle(24);
      h_mMcMass_phistar[KEY_MC_phistar]->SetMarkerColor(kBlack);
      h_mMcMass_phistar[KEY_MC_phistar]->SetLineColor(kBlack);
      h_mMcMass_phistar[KEY_MC_phistar]->Draw("pE same");

      if(i == 0)
      {
        legphistarmc->AddEntry(h_mMass_phistar[KEY_phistar],"Data","p");
        legphistarmc->AddEntry(h_mMcMass_phistar[KEY_MC_phistar],"MC","p");
      }
      legphistarmc->Draw("same");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistar_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_MC_phiy = Form("mcphi%s%s_phiy_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMcMass_phiy[KEY_MC_phiy]->SetTitle(Form("MC %1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMcMass_phiy[KEY_MC_phiy]->GetXaxis()->SetTitle(phiphilabel[iset].c_str());
      h_mMcMass_phiy[KEY_MC_phiy]->GetYaxis()->SetTitle(phietalabel[iset].c_str());

      h_mMcMass_phiy[KEY_MC_phiy]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phiy_mc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    for(int i = 0; i < 12; i++)
    {
      c1->cd(i+1);
      string KEY_MC_phistarphi = Form("mcphi%s%s_phistarphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
      h_mMass_MC[KEY_MC_phistarphi]->SetTitle(Form("MC %1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
      h_mMass_MC[KEY_MC_phistarphi]->GetXaxis()->SetTitle("#phi");
      h_mMass_MC[KEY_MC_phistarphi]->GetYaxis()->SetTitle("#phi*");

      h_mMass_MC[KEY_MC_phistarphi]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistarphi_mc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    //for(int i = 0; i < 12; i++)
    //{
    //  c1->cd(i+1);
    //  string KEY_MC_phistarmphi = Form("mcphi%s%s_phistarmphi_phipt_%d",philabel[iset].c_str(),phiextralabel[iset].c_str(),i);
    //  h_mMass_MC_1D[KEY_MC_phistarmphi]->SetTitle(Form("MC %1.1f<p_{T}<%1.1f, 20-60 Centrality",data_constants::phi_pt_low[i],data_constants::phi_pt_high[i]));
    //  h_mMass_MC_1D[KEY_MC_phistarmphi]->GetXaxis()->SetTitle("#phi*-#phi");
    //  h_mMass_MC_1D[KEY_MC_phistarmphi]->GetYaxis()->SetTitle("Counts");

    //  h_mMass_MC_1D[KEY_MC_phistarmphi]->Draw("colz");
    //}
    //c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_%s%s_phistarmphi_mc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),philabel[iset].c_str(),phiextralabel[iset].c_str()));

    TCanvas *cpt = new TCanvas("cpt","cpt",10,10,400,400);     
    for(int i = 0; i < 1; i++)
    {
      cpt->cd();
      cpt->cd()->SetLeftMargin(0.15);
      cpt->cd()->SetBottomMargin(0.15);
      cpt->cd()->SetTicks(1,1);
      cpt->cd()->SetGrid(0,0);  
    } 
    cpt->cd();
    string KEY_pT_ratio = Form("phiratio_pt");
    g_mPhi_pT_ratio[KEY_pT_ratio]->SetTitle(Form("Data/RC 20-60 Centrality"));
    g_mPhi_pT_ratio[KEY_pT_ratio]->GetXaxis()->SetTitle("p_{T}");
    g_mPhi_pT_ratio[KEY_pT_ratio]->GetYaxis()->SetTitle("Data/RC");
    g_mPhi_pT_ratio[KEY_pT_ratio]->SetMarkerStyle(20);

    g_mPhi_pT_ratio[KEY_pT_ratio]->Draw("APE");
    cpt->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_pT_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

    TLegend *legpt = new TLegend(0.6,0.7,0.8,0.9);
    cpt->cd();
    cpt->SetLogy();
    string KEY_pT = Form("phi_pt");
    string KEY_RC_pT = Form("rcphi_pt");
    g_mPhi_pT[KEY_pT]->SetTitle(Form("20-60 Centrality"));
    g_mPhi_pT[KEY_pT]->GetXaxis()->SetTitle("p_{T}");
    g_mPhi_pT[KEY_pT]->GetYaxis()->SetTitle("Counts");
    g_mPhi_pT[KEY_pT]->SetMarkerStyle(20);
    g_mPhi_pT[KEY_pT]->SetMarkerColor(kOrange+7);
    g_mPhi_pT[KEY_pT]->SetLineColor(kOrange+7);

    g_mPhi_pT[KEY_pT]->Draw("APE");

    g_mRcPhi_pT[KEY_RC_pT]->SetMarkerStyle(24);
    g_mRcPhi_pT[KEY_RC_pT]->SetMarkerColor(kBlack);
    g_mRcPhi_pT[KEY_RC_pT]->SetLineColor(kBlack);
    g_mRcPhi_pT[KEY_RC_pT]->Draw("pE same");

    legpt->AddEntry(g_mPhi_pT[KEY_pT],"Data","p");
    legpt->AddEntry(g_mRcPhi_pT[KEY_RC_pT],"RC","p");
    legpt->Draw("same");
    cpt->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_pT_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

    cpt->cd();
    cpt->SetLogy(0);
    string KEY_pT_mcratio = Form("mcphiratio_pt");
    g_mPhi_pT_mcratio[KEY_pT_mcratio]->SetTitle(Form("Data/MC 20-60 Centrality"));
    g_mPhi_pT_mcratio[KEY_pT_mcratio]->GetXaxis()->SetTitle("p_{T}");
    g_mPhi_pT_mcratio[KEY_pT_mcratio]->GetYaxis()->SetTitle("Data/MC");
    g_mPhi_pT_mcratio[KEY_pT_mcratio]->SetMarkerStyle(20);
    g_mPhi_pT_mcratio[KEY_pT_mcratio]->Draw("APE");
    cpt->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_pT_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

    TLegend *legptmc = new TLegend(0.6,0.7,0.8,0.9);
    cpt->cd();
    cpt->SetLogy();
    string KEY_MC_pT = Form("mcphi_pt");
    g_mPhi_pT[KEY_pT]->SetTitle(Form("20-60 Centrality"));
    g_mPhi_pT[KEY_pT]->GetXaxis()->SetTitle("p_{T}");
    g_mPhi_pT[KEY_pT]->GetYaxis()->SetTitle("Counts");
    g_mPhi_pT[KEY_pT]->SetMarkerStyle(20);
    g_mPhi_pT[KEY_pT]->SetMarkerColor(kOrange+7);
    g_mPhi_pT[KEY_pT]->SetLineColor(kOrange+7);

    g_mPhi_pT[KEY_pT]->Draw("APE");

    g_mMcPhi_pT[KEY_MC_pT]->SetMarkerStyle(24);
    g_mMcPhi_pT[KEY_MC_pT]->SetMarkerColor(kBlack);
    g_mMcPhi_pT[KEY_MC_pT]->SetLineColor(kBlack);
    g_mMcPhi_pT[KEY_MC_pT]->Draw("pE same");

    legptmc->AddEntry(g_mPhi_pT[KEY_pT],"Data","p");
    legptmc->AddEntry(g_mMcPhi_pT[KEY_MC_pT],"MC","p");
    legptmc->Draw("same");
    cpt->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_PhiDistribution_pT_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));
  }
}

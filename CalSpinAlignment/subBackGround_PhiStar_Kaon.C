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

void subBackGround_PhiStar_Kaon(int energy = 4, int pid = 0, int year = 0, string date = "20240426", bool random3D = false, int order = 2, string etamode = "eta1_eta1", int deltaonly = 0)
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

  //string InPutFile_RC = "effaccfiles/Phi/19GeV/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240326_Acc_pT0p1_TPC_TOF_kaons_20240327/Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0.root";  
  string folder = "CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240515_fixed2Dfunc_FinerPhiBins_Spectra_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0";
  string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0.root",folder.c_str());  
  TFile *File_RC = TFile::Open(InPutFile_RC.c_str());

  cout << "Loaded All Files" << endl;

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH2FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;
  TH2FMap h_mMass_MC;
  TH2FMap h_mMass_RC;
  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++) // pt bin 
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
    {
      for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
      {
        string KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_RC[KEY_RC] = (TH2F*)((TH2F*)File_RC->Get(KEY_RC.c_str()))->Clone();

        KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_RC[KEY_RC] = (TH2F*)((TH2F*)File_RC->Get(KEY_RC.c_str()))->Clone();

        string KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_MC[KEY_MC] = (TH2F*)((TH2F*)File_RC->Get(KEY_MC.c_str()))->Clone();

        KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_MC[KEY_MC] = (TH2F*)((TH2F*)File_RC->Get(KEY_MC.c_str()))->Clone();

        string KEY_SE = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_SE",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();

        string KEY_ME = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_ME",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();

        string KEY_SM = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
        h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
        h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     

        KEY_SE = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta_SE",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();

        KEY_ME = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta_ME",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();

        KEY_SM = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
        h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
        h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     
      }
      string KEY_RC = Form("rckaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_RC[KEY_RC] = (TH2F*)((TH2F*)File_RC->Get(KEY_RC.c_str()))->Clone();

      KEY_RC = Form("rckaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_RC[KEY_RC] = (TH2F*)((TH2F*)File_RC->Get(KEY_RC.c_str()))->Clone();

      KEY_RC = Form("rckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_RC[KEY_RC] = (TH2F*)((TH2F*)File_RC->Get(KEY_RC.c_str()))->Clone();

      KEY_RC = Form("rckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_RC[KEY_RC] = (TH2F*)((TH2F*)File_RC->Get(KEY_RC.c_str()))->Clone();

      string KEY_MC = Form("mckaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_MC[KEY_MC] = (TH2F*)((TH2F*)File_RC->Get(KEY_MC.c_str()))->Clone();

      KEY_MC = Form("mckaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_MC[KEY_MC] = (TH2F*)((TH2F*)File_RC->Get(KEY_MC.c_str()))->Clone();

      KEY_MC = Form("mckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_MC[KEY_MC] = (TH2F*)((TH2F*)File_RC->Get(KEY_MC.c_str()))->Clone();

      KEY_MC = Form("mckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_MC[KEY_MC] = (TH2F*)((TH2F*)File_RC->Get(KEY_MC.c_str()))->Clone();

      string KEY_SE = Form("kaon_delta_phipt_%d_cos2phistarphi_%d_SE",i_pt_phi,i_phi);
      h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();

      string KEY_ME = Form("kaon_delta_phipt_%d_cos2phistarphi_%d_ME",i_pt_phi,i_phi);
      h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();

      string KEY_SM = Form("kaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
      h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
      h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     

      KEY_SE = Form("kaon_deltarphi_phipt_%d_cos2phistarphi_%d_SE",i_pt_phi,i_phi);
      h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();

      KEY_ME = Form("kaon_deltarphi_phipt_%d_cos2phistarphi_%d_ME",i_pt_phi,i_phi);
      h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();

      KEY_SM = Form("kaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
      h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
      h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     

      KEY_ME = Form("kaonplusphi_delta_phipt_%d_cos2phistarphi_%d_ME",i_pt_phi,i_phi);
      h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();

      KEY_SE = Form("kaonplusphi_delta_phipt_%d_cos2phistarphi_%d_SE",i_pt_phi,i_phi);
      h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();

      KEY_SM = Form("kaonplusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
      h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
      h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     

      KEY_SE = Form("kaonminusphi_delta_phipt_%d_cos2phistarphi_%d_SE",i_pt_phi,i_phi);
      h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();

      KEY_ME = Form("kaonminusphi_delta_phipt_%d_cos2phistarphi_%d_ME",i_pt_phi,i_phi);
      h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();

      KEY_SM = Form("kaonminusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
      h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
      h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     
    }
  }

  TGraMap g_mKaon_pT; 
  TGraMap g_mRcKaon_pT; 
  TGraMap g_mMcKaon_pT; 
  TGraMap g_mKaon_pT_ratio; 
  TGraMap g_mKaon_pT_mcratio; 
  TH2FMap h_mMass_phiy;   
  TH2FMap h_mRcMass_phiy;   
  TH2FMap h_mMcMass_phiy;   
  TH2FMap h_mMass_phieta;   
  TH2FMap h_mRcMass_phieta;   
  TH2FMap h_mMcMass_phieta;   
  TH2FMap h_mMass_phiy_ratio, h_mMass_phiy_mcratio, h_mMass_phieta_ratio, h_mMass_phieta_mcratio;   
  TH1FMap h_mMass_phi, h_mMass_y, h_mMass_eta;   
  TH1FMap h_mRcMass_phi, h_mRcMass_y, h_mRcMass_eta;
  TH1FMap h_mMcMass_phi, h_mMcMass_y, h_mMcMass_eta;
  TH1FMap h_mMass_phi_ratio, h_mMass_y_ratio, h_mMass_eta_ratio;
  TH1FMap h_mMass_phi_mcratio, h_mMass_y_mcratio, h_mMass_eta_mcratio;

  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++) // pt bin 
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
    {
      for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
      {
        string KEY_SM = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        string KEY_phiy = Form("kaon_phiy_phipt_%d",i_pt_phi);
        if(i_phi == 0 && i_pt_kaon == 0) h_mMass_phiy[KEY_phiy] = (TH2F*)h_mMass_SM[KEY_SM]->Clone(KEY_phiy.c_str()); // First bins
        else h_mMass_phiy[KEY_phiy]->Add(h_mMass_SM[KEY_SM],1.0); // All other bins

        string KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        string KEY_RC_phiy = Form("rckaon_phiy_phipt_%d",i_pt_phi);
        if(i_phi == 0 && i_pt_kaon == 0) h_mRcMass_phiy[KEY_RC_phiy] = (TH2F*)h_mMass_RC[KEY_RC]->Clone(KEY_RC_phiy.c_str()); // First bins
        else h_mRcMass_phiy[KEY_RC_phiy]->Add(h_mMass_RC[KEY_RC],1.0); // All other bins

        string KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        string KEY_MC_phiy = Form("mckaon_phiy_phipt_%d",i_pt_phi);
        if(i_phi == 0 && i_pt_kaon == 0) h_mMcMass_phiy[KEY_MC_phiy] = (TH2F*)h_mMass_MC[KEY_MC]->Clone(KEY_MC_phiy.c_str()); // First bins
        else h_mMcMass_phiy[KEY_MC_phiy]->Add(h_mMass_MC[KEY_MC],1.0); // All other bins

        KEY_SM = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        string KEY_phieta = Form("kaon_phieta_phipt_%d",i_pt_phi);
        if(i_phi == 0 && i_pt_kaon == 0) h_mMass_phieta[KEY_phieta] = (TH2F*)h_mMass_SM[KEY_SM]->Clone(KEY_phieta.c_str()); // First bins
        else h_mMass_phieta[KEY_phieta]->Add(h_mMass_SM[KEY_SM],1.0); // All other bins

        KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        string KEY_RC_phieta = Form("rckaon_phieta_phipt_%d",i_pt_phi);
        if(i_phi == 0 && i_pt_kaon == 0) h_mRcMass_phieta[KEY_RC_phieta] = (TH2F*)h_mMass_RC[KEY_RC]->Clone(KEY_RC_phieta.c_str()); // First bins
        else h_mRcMass_phieta[KEY_RC_phieta]->Add(h_mMass_RC[KEY_RC],1.0); // All other bins

        KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        string KEY_MC_phieta = Form("mckaon_phieta_phipt_%d",i_pt_phi);
        if(i_phi == 0 && i_pt_kaon == 0) h_mMcMass_phieta[KEY_MC_phieta] = (TH2F*)h_mMass_MC[KEY_MC]->Clone(KEY_MC_phieta.c_str()); // First bins
        else h_mMcMass_phieta[KEY_MC_phieta]->Add(h_mMass_MC[KEY_MC],1.0); // All other bins

        string KEY_phi = Form("kaon_phi_phipt_%d",i_pt_phi);
        string KEY_y   = Form("kaon_y_phipt_%d",i_pt_phi);
        string KEY_eta   = Form("kaon_eta_phipt_%d",i_pt_phi);
        if(i_phi == 9 && i_pt_kaon == data_constants::kaon_pt_bins-1) // Last bins
        {
          h_mMass_phi[KEY_phi] = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionX(KEY_phi.c_str(),0,-1,"e"));
          h_mMass_y[KEY_y]     = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionY(KEY_y.c_str(),0,-1,"e"));
          h_mMass_eta[KEY_eta] = (TH1F*) (h_mMass_phieta[KEY_phieta]->ProjectionY(KEY_eta.c_str(),0,-1,"e"));
        } 

        string KEY_RC_phi = Form("rckaon_phi_phipt_%d",i_pt_phi);
        string KEY_RC_y   = Form("rckaon_y_phipt_%d",i_pt_phi);
        string KEY_RC_eta   = Form("rckaon_eta_phipt_%d",i_pt_phi);
        if(i_phi == 9 && i_pt_kaon == data_constants::kaon_pt_bins-1) // Last bins
        {
          h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
          h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));
          h_mRcMass_eta[KEY_RC_eta] = (TH1F*) (h_mRcMass_phieta[KEY_RC_phieta]->ProjectionY(KEY_RC_eta.c_str(),0,-1,"e"));
        } 

        string KEY_MC_phi = Form("mckaon_phi_phipt_%d",i_pt_phi);
        string KEY_MC_y   = Form("mckaon_y_phipt_%d",i_pt_phi);
        string KEY_MC_eta   = Form("mckaon_eta_phipt_%d",i_pt_phi);
        if(i_phi == 9 && i_pt_kaon == data_constants::kaon_pt_bins-1) // Last bins
        {
          h_mMcMass_phi[KEY_MC_phi] = (TH1F*) (h_mMcMass_phiy[KEY_MC_phiy]->ProjectionX(KEY_MC_phi.c_str(),0,-1,"e"));
          h_mMcMass_y[KEY_MC_y]     = (TH1F*) (h_mMcMass_phiy[KEY_MC_phiy]->ProjectionY(KEY_MC_y.c_str(),0,-1,"e"));
          h_mMcMass_eta[KEY_MC_eta] = (TH1F*) (h_mMcMass_phieta[KEY_MC_phieta]->ProjectionY(KEY_MC_eta.c_str(),0,-1,"e"));
        } 
      }
      string KEY_SM = Form("kaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_phiy = Form("kaon_delta_phiy_phipt_%d",i_pt_phi);
      //string KEY_phiy_c = Form("kaon_delta_phiy_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      if(i_phi == 0) h_mMass_phiy[KEY_phiy] = (TH2F*)h_mMass_SM[KEY_SM]->Clone(KEY_phiy.c_str()); // First bins
      else h_mMass_phiy[KEY_phiy]->Add(h_mMass_SM[KEY_SM],1.0); // All other bins

      string KEY_RC = Form("rckaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_RC_phiy = Form("rckaon_delta_phiy_phipt_%d",i_pt_phi);
      //string KEY_RC_phiy_c = Form("rckaon_delta_phiy_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      if(i_phi == 0) h_mRcMass_phiy[KEY_RC_phiy] = (TH2F*)h_mMass_RC[KEY_RC]->Clone(KEY_RC_phiy.c_str()); // First bins
      else h_mRcMass_phiy[KEY_RC_phiy]->Add(h_mMass_RC[KEY_RC],1.0); // All other bins

      string KEY_MC = Form("mckaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_MC_phiy = Form("mckaon_delta_phiy_phipt_%d",i_pt_phi);
      if(i_phi == 0) h_mMcMass_phiy[KEY_MC_phiy] = (TH2F*)h_mMass_MC[KEY_MC]->Clone(KEY_MC_phiy.c_str()); // First bins
      else h_mMcMass_phiy[KEY_MC_phiy]->Add(h_mMass_MC[KEY_MC],1.0); // All other bins

      string KEY_phi = Form("kaon_delta_phi_phipt_%d",i_pt_phi);
      string KEY_y   = Form("kaon_delta_y_phipt_%d",i_pt_phi);
      string KEY_SM_phi = Form("kaon_delta_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_SM_y   = Form("kaon_delta_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mMass_phi[KEY_SM_phi] = (TH1F*) (h_mMass_SM[KEY_SM]->ProjectionX(KEY_SM_phi.c_str(),0,-1,"e"));
      h_mMass_y[KEY_SM_y]     = (TH1F*) (h_mMass_SM[KEY_SM]->ProjectionY(KEY_SM_y.c_str(),0,-1,"e"));
      if(i_phi == 9) // Last bins
      {
        h_mMass_phi[KEY_phi] = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionX(KEY_phi.c_str(),0,-1,"e"));
        h_mMass_y[KEY_y]     = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionY(KEY_y.c_str(),0,-1,"e"));
      } 

      string KEY_RC_phi = Form("rckaon_delta_phi_phipt_%d",i_pt_phi);
      string KEY_RC_y   = Form("rckaon_delta_y_phipt_%d",i_pt_phi);
      if(i_phi == 9) // Last bins
      {
        h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
        h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));
      } 
      KEY_RC_phi = Form("rckaon_delta_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_RC_y   = Form("rckaon_delta_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mMass_RC[KEY_RC]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
      h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mMass_RC[KEY_RC]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));

      string KEY_MC_phi = Form("mckaon_delta_phi_phipt_%d",i_pt_phi);
      string KEY_MC_y   = Form("mckaon_delta_y_phipt_%d",i_pt_phi);
      if(i_phi == 9) // Last bins
      {
        h_mMcMass_phi[KEY_MC_phi] = (TH1F*) (h_mMcMass_phiy[KEY_MC_phiy]->ProjectionX(KEY_MC_phi.c_str(),0,-1,"e"));
        h_mMcMass_y[KEY_MC_y]     = (TH1F*) (h_mMcMass_phiy[KEY_MC_phiy]->ProjectionY(KEY_MC_y.c_str(),0,-1,"e"));
      } 
      KEY_MC_phi = Form("mckaon_delta_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_MC_y   = Form("mckaon_delta_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mMcMass_phi[KEY_MC_phi] = (TH1F*) (h_mMass_MC[KEY_MC]->ProjectionX(KEY_MC_phi.c_str(),0,-1,"e"));
      h_mMcMass_y[KEY_MC_y]     = (TH1F*) (h_mMass_MC[KEY_MC]->ProjectionY(KEY_MC_y.c_str(),0,-1,"e"));
      
      ///////////Delta R
      KEY_SM = Form("kaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_phiy = Form("kaon_deltarphi_phiy_phipt_%d",i_pt_phi);
      //KEY_phiy_c = Form("kaon_deltarphi_phiy_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      if(i_phi == 0) h_mMass_phiy[KEY_phiy] = (TH2F*)h_mMass_SM[KEY_SM]->Clone(KEY_phiy.c_str()); // First bins
      else h_mMass_phiy[KEY_phiy]->Add(h_mMass_SM[KEY_SM],1.0); // All other bins

      KEY_RC = Form("rckaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_RC_phiy = Form("rckaon_deltarphi_phiy_phipt_%d",i_pt_phi);
      //KEY_RC_phiy_c = Form("rckaon_deltarphi_phiy_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      if(i_phi == 0) h_mRcMass_phiy[KEY_RC_phiy] = (TH2F*)h_mMass_RC[KEY_RC]->Clone(KEY_RC_phiy.c_str()); // First bins
      else h_mRcMass_phiy[KEY_RC_phiy]->Add(h_mMass_RC[KEY_RC],1.0); // All other bins

      KEY_MC = Form("mckaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_MC_phiy = Form("mckaon_deltarphi_phiy_phipt_%d",i_pt_phi);
      if(i_phi == 0) h_mMcMass_phiy[KEY_MC_phiy] = (TH2F*)h_mMass_MC[KEY_MC]->Clone(KEY_MC_phiy.c_str()); // First bins
      else h_mMcMass_phiy[KEY_MC_phiy]->Add(h_mMass_MC[KEY_MC],1.0); // All other bins

      KEY_phi = Form("kaon_deltarphi_phi_phipt_%d",i_pt_phi);
      KEY_y   = Form("kaon_deltarphi_y_phipt_%d",i_pt_phi);
      KEY_SM_phi = Form("kaon_deltarphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_SM_y   = Form("kaon_deltarphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mMass_phi[KEY_SM_phi] = (TH1F*) (h_mMass_SM[KEY_SM]->ProjectionX(KEY_SM_phi.c_str(),0,-1,"e"));
      h_mMass_y[KEY_SM_y]     = (TH1F*) (h_mMass_SM[KEY_SM]->ProjectionY(KEY_SM_y.c_str(),0,-1,"e"));
      if(i_phi == 9) // Last bins
      {
        h_mMass_phi[KEY_phi] = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionX(KEY_phi.c_str(),0,-1,"e"));
        h_mMass_y[KEY_y]     = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionY(KEY_y.c_str(),0,-1,"e"));
      } 

      KEY_RC_phi = Form("rckaon_deltarphi_phi_phipt_%d",i_pt_phi);
      KEY_RC_y   = Form("rckaon_deltarphi_y_phipt_%d",i_pt_phi);
      if(i_phi == 9) // Last bins
      {
        h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
        h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));
      } 
      KEY_RC_phi = Form("rckaon_deltarphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_RC_y   = Form("rckaon_deltarphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mMass_RC[KEY_RC]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
      h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mMass_RC[KEY_RC]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));

      KEY_MC_phi = Form("mckaon_deltarphi_phi_phipt_%d",i_pt_phi);
      KEY_MC_y   = Form("mckaon_deltarphi_y_phipt_%d",i_pt_phi);
      if(i_phi == 9) // Last bins
      {
        h_mMcMass_phi[KEY_MC_phi] = (TH1F*) (h_mMcMass_phiy[KEY_MC_phiy]->ProjectionX(KEY_MC_phi.c_str(),0,-1,"e"));
        h_mMcMass_y[KEY_MC_y]     = (TH1F*) (h_mMcMass_phiy[KEY_MC_phiy]->ProjectionY(KEY_MC_y.c_str(),0,-1,"e"));
      } 
      KEY_MC_phi = Form("mckaon_deltarphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_MC_y   = Form("mckaon_deltarphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mMcMass_phi[KEY_MC_phi] = (TH1F*) (h_mMass_MC[KEY_MC]->ProjectionX(KEY_MC_phi.c_str(),0,-1,"e"));
      h_mMcMass_y[KEY_MC_y]     = (TH1F*) (h_mMass_MC[KEY_MC]->ProjectionY(KEY_MC_y.c_str(),0,-1,"e"));


      //KEY_SM = Form("kaonplusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      //KEY_phiy = Form("kaonplusphi_delta_phiy_phipt_%d",i_pt_phi);
      //if(i_phi == 0) h_mMass_phiy[KEY_phiy] = (TH2F*)h_mMass_SM[KEY_SM]->Clone(KEY_phiy.c_str()); // First bins
      //else h_mMass_phiy[KEY_phiy]->Add(h_mMass_SM[KEY_SM],1.0); // All other bins

      //KEY_RC = Form("rckaonplusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      //KEY_RC_phiy = Form("rckaonplusphi_delta_phiy_phipt_%d",i_pt_phi);
      //if(i_phi == 0) h_mRcMass_phiy[KEY_RC_phiy] = (TH2F*)h_mMass_RC[KEY_RC]->Clone(KEY_RC_phiy.c_str()); // First bins
      //else h_mRcMass_phiy[KEY_RC_phiy]->Add(h_mMass_RC[KEY_RC],1.0); // All other bins

      //KEY_phi = Form("kaonplusphi_delta_phi_phipt_%d",i_pt_phi);
      //KEY_y   = Form("kaonplusphi_delta_y_phipt_%d",i_pt_phi);
      //if(i_phi == 9) // Last bins
      //{
      //  h_mMass_phi[KEY_phi] = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionX(KEY_phi.c_str(),0,-1,"e"));
      //  h_mMass_y[KEY_y]     = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionY(KEY_y.c_str(),0,-1,"e"));
      //} 

      //KEY_RC_phi = Form("rckaonplusphi_delta_phi_phipt_%d",i_pt_phi);
      //KEY_RC_y   = Form("rckaonplusphi_delta_y_phipt_%d",i_pt_phi);
      //if(i_phi == 9) // Last bins
      //{
      //  h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
      //  h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));
      //} 

      //KEY_SM = Form("kaonminusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      //KEY_phiy = Form("kaonminusphi_delta_phiy_phipt_%d",i_pt_phi);
      //if(i_phi == 0) h_mMass_phiy[KEY_phiy] = (TH2F*)h_mMass_SM[KEY_SM]->Clone(KEY_phiy.c_str()); // First bins
      //else h_mMass_phiy[KEY_phiy]->Add(h_mMass_SM[KEY_SM],1.0); // All other bins

      //KEY_RC = Form("rckaonminusphi_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      //KEY_RC_phiy = Form("rckaonminusphi_delta_phiy_phipt_%d",i_pt_phi);
      //if(i_phi == 0) h_mRcMass_phiy[KEY_RC_phiy] = (TH2F*)h_mMass_RC[KEY_RC]->Clone(KEY_RC_phiy.c_str()); // First bins
      //else h_mRcMass_phiy[KEY_RC_phiy]->Add(h_mMass_RC[KEY_RC],1.0); // All other bins

      //KEY_phi = Form("kaonminusphi_delta_phi_phipt_%d",i_pt_phi);
      //KEY_y   = Form("kaonminusphi_delta_y_phipt_%d",i_pt_phi);
      //if(i_phi == 9) // Last bins
      //{
      //  h_mMass_phi[KEY_phi] = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionX(KEY_phi.c_str(),0,-1,"e"));
      //  h_mMass_y[KEY_y]     = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionY(KEY_y.c_str(),0,-1,"e"));
      //} 

      //KEY_RC_phi = Form("rckaonminusphi_delta_phi_phipt_%d",i_pt_phi);
      //KEY_RC_y   = Form("rckaonminusphi_delta_y_phipt_%d",i_pt_phi);
      //if(i_phi == 9) // Last bins
      //{
      //  h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
      //  h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));
      //} 
    }
  }

  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++) // pt bin 
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++)
    {
      /////////// DELTA //////////////////////
      string KEY_phiy_ratio = Form("kaonratio_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_phiy_mcratio = Form("mckaonratio_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_phiy = Form("kaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_RC_phiy = Form("rckaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_MC_phiy = Form("mckaon_delta_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
       
      h_mMass_phiy_ratio[KEY_phiy_ratio] = (TH2F*)h_mMass_SM[KEY_phiy]->Clone(KEY_phiy_ratio.c_str()); // All other bins
      h_mMass_phiy_mcratio[KEY_phiy_mcratio] = (TH2F*)h_mMass_SM[KEY_phiy]->Clone(KEY_phiy_mcratio.c_str()); // All other bins
      float phiy_yield   = h_mMass_SM[KEY_phiy]->Integral(0,-1,0,-1);
      cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
      float phiy_rcyield = h_mMass_RC[KEY_RC_phiy]->Integral(0,-1,0,-1);
      float phiy_mcyield = h_mMass_MC[KEY_MC_phiy]->Integral(0,-1,0,-1);

      h_mMass_RC[KEY_RC_phiy]->Scale(phiy_yield/phiy_rcyield);
      h_mMass_MC[KEY_MC_phiy]->Scale(phiy_yield/phiy_mcyield);
      h_mMass_phiy_ratio[KEY_phiy_ratio]->Divide(h_mMass_RC[KEY_RC_phiy]);
      h_mMass_phiy_mcratio[KEY_phiy_mcratio]->Divide(h_mMass_MC[KEY_MC_phiy]);
   
      string KEY_phi = Form("kaon_delta_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_y   = Form("kaon_delta_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_RC_phi = Form("rckaon_delta_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_RC_y   = Form("rckaon_delta_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_MC_phi = Form("mckaon_delta_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_MC_y   = Form("mckaon_delta_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_phi_ratio = Form("kaonratio_delta_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_phi_mcratio = Form("mckaonratio_delta_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_y_ratio   = Form("kaonratio_delta_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      string KEY_y_mcratio   = Form("mckaonratio_delta_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mMass_phi_ratio[KEY_phi_ratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_ratio.c_str()); // All other bins
      h_mMass_phi_mcratio[KEY_phi_mcratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_mcratio.c_str()); // All other bins
      float phi_yield   = h_mMass_phi[KEY_phi]->Integral(0,-1);
      float phi_rcyield = h_mRcMass_phi[KEY_RC_phi]->Integral(0,-1);
      float phi_mcyield = h_mMcMass_phi[KEY_MC_phi]->Integral(0,-1);

      cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
      h_mRcMass_phi[KEY_RC_phi]->Scale(phi_yield/phi_rcyield);
      h_mMcMass_phi[KEY_MC_phi]->Scale(phi_yield/phi_mcyield);
      h_mMass_phi_ratio[KEY_phi_ratio]->Divide(h_mRcMass_phi[KEY_RC_phi]);
      h_mMass_phi_mcratio[KEY_phi_mcratio]->Divide(h_mMcMass_phi[KEY_MC_phi]);

      h_mMass_y_ratio[KEY_y_ratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_ratio.c_str()); // All other bins
      h_mMass_y_mcratio[KEY_y_mcratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_mcratio.c_str()); // All other bins
      float y_yield   = h_mMass_y[KEY_y]->Integral(0,-1);
      float y_rcyield = h_mRcMass_y[KEY_RC_y]->Integral(0,-1);
      float y_mcyield = h_mMcMass_y[KEY_MC_y]->Integral(0,-1);

      cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
      h_mRcMass_y[KEY_RC_y]->Scale(y_yield/y_rcyield);
      h_mMcMass_y[KEY_MC_y]->Scale(y_yield/y_mcyield);
      h_mMass_y_ratio[KEY_y_ratio]->Divide(h_mRcMass_y[KEY_RC_y]);
      h_mMass_y_mcratio[KEY_y_mcratio]->Divide(h_mMcMass_y[KEY_MC_y]);

      /////Delta R
      KEY_phiy_ratio = Form("kaonratio_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_phiy_mcratio = Form("mckaonratio_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_phiy = Form("kaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_RC_phiy = Form("rckaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_MC_phiy = Form("mckaon_deltarphi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
       
      h_mMass_phiy_ratio[KEY_phiy_ratio] = (TH2F*)h_mMass_SM[KEY_phiy]->Clone(KEY_phiy_ratio.c_str()); // All other bins
      h_mMass_phiy_mcratio[KEY_phiy_mcratio] = (TH2F*)h_mMass_SM[KEY_phiy]->Clone(KEY_phiy_mcratio.c_str()); // All other bins
      phiy_yield   = h_mMass_SM[KEY_phiy]->Integral(0,-1,0,-1);
      cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
      phiy_rcyield = h_mMass_RC[KEY_RC_phiy]->Integral(0,-1,0,-1);
      phiy_mcyield = h_mMass_MC[KEY_MC_phiy]->Integral(0,-1,0,-1);

      h_mMass_RC[KEY_RC_phiy]->Scale(phiy_yield/phiy_rcyield);
      h_mMass_MC[KEY_MC_phiy]->Scale(phiy_yield/phiy_mcyield);
      h_mMass_phiy_ratio[KEY_phiy_ratio]->Divide(h_mMass_RC[KEY_RC_phiy]);
      h_mMass_phiy_mcratio[KEY_phiy_mcratio]->Divide(h_mMass_MC[KEY_MC_phiy]);
   
      KEY_phi = Form("kaon_deltarphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_y   = Form("kaon_deltarphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_RC_phi = Form("rckaon_deltarphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_RC_y   = Form("rckaon_deltarphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_MC_phi = Form("mckaon_deltarphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_MC_y   = Form("mckaon_deltarphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_phi_ratio = Form("kaonratio_deltarphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_phi_mcratio = Form("mckaonratio_deltarphi_phi_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_y_ratio   = Form("kaonratio_deltarphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);
      KEY_y_mcratio   = Form("mckaonratio_deltarphi_y_phipt_%d_cos2phistarphi_%d",i_pt_phi,i_phi);

      h_mMass_phi_ratio[KEY_phi_ratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_ratio.c_str()); // All other bins
      h_mMass_phi_mcratio[KEY_phi_mcratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_mcratio.c_str()); // All other bins
      phi_yield   = h_mMass_phi[KEY_phi]->Integral(0,-1);
      phi_rcyield = h_mRcMass_phi[KEY_RC_phi]->Integral(0,-1);
      phi_mcyield = h_mMcMass_phi[KEY_MC_phi]->Integral(0,-1);

      cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
      h_mRcMass_phi[KEY_RC_phi]->Scale(phi_yield/phi_rcyield);
      h_mMcMass_phi[KEY_MC_phi]->Scale(phi_yield/phi_mcyield);
      h_mMass_phi_ratio[KEY_phi_ratio]->Divide(h_mRcMass_phi[KEY_RC_phi]);
      h_mMass_phi_mcratio[KEY_phi_mcratio]->Divide(h_mMcMass_phi[KEY_MC_phi]);

      h_mMass_y_ratio[KEY_y_ratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_ratio.c_str()); // All other bins
      h_mMass_y_mcratio[KEY_y_mcratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_mcratio.c_str()); // All other bins
      y_yield   = h_mMass_y[KEY_y]->Integral(0,-1);
      y_rcyield = h_mRcMass_y[KEY_RC_y]->Integral(0,-1);
      y_mcyield = h_mMcMass_y[KEY_MC_y]->Integral(0,-1);

      cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
      h_mRcMass_y[KEY_RC_y]->Scale(y_yield/y_rcyield);
      h_mMcMass_y[KEY_MC_y]->Scale(y_yield/y_mcyield);
      h_mMass_y_ratio[KEY_y_ratio]->Divide(h_mRcMass_y[KEY_RC_y]);
      h_mMass_y_mcratio[KEY_y_mcratio]->Divide(h_mMcMass_y[KEY_MC_y]);
    }
  }

  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++) // pt bin 
  {
    string KEY_phiy_ratio = Form("kaonratio_phiy_phipt_%d",i_pt_phi);
    string KEY_phieta_ratio = Form("kaonratio_phieta_phipt_%d",i_pt_phi);
    string KEY_phiy = Form("kaon_phiy_phipt_%d",i_pt_phi);
    string KEY_phieta = Form("kaon_phieta_phipt_%d",i_pt_phi);
    string KEY_RC_phiy = Form("rckaon_phiy_phipt_%d",i_pt_phi);
    string KEY_RC_phieta = Form("rckaon_phieta_phipt_%d",i_pt_phi);
    string KEY_phiy_mcratio = Form("mckaonratio_phiy_phipt_%d",i_pt_phi);
    string KEY_phieta_mcratio = Form("mckaonratio_phieta_phipt_%d",i_pt_phi);
    string KEY_MC_phiy = Form("mckaon_phiy_phipt_%d",i_pt_phi);
    string KEY_MC_phieta = Form("mckaon_phieta_phipt_%d",i_pt_phi);
     
    h_mMass_phiy_ratio[KEY_phiy_ratio] = (TH2F*)h_mMass_phiy[KEY_phiy]->Clone(KEY_phiy_ratio.c_str()); // All other bins
    double phiy_yield   = h_mMass_phiy[KEY_phiy]->Integral(0,-1,0,-1);
    cout << "NOTE THIS YIELD IS PHI YIELD Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield/2 << endl;
    double phiy_rcyield = h_mRcMass_phiy[KEY_RC_phiy]->Integral(0,-1,0,-1);

    h_mRcMass_phiy[KEY_RC_phiy]->Scale(phiy_yield/phiy_rcyield);
    h_mMass_phiy_ratio[KEY_phiy_ratio]->Divide(h_mRcMass_phiy[KEY_RC_phiy]);

    h_mMass_phieta_ratio[KEY_phieta_ratio] = (TH2F*)h_mMass_phieta[KEY_phieta]->Clone(KEY_phieta_ratio.c_str()); // All other bins
    double phieta_yield   = h_mMass_phieta[KEY_phieta]->Integral(1,h_mMass_phieta[KEY_phieta]->GetNbinsX(),1,h_mMass_phieta[KEY_phieta]->GetNbinsY());
    cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    double phieta_rcyield = h_mRcMass_phieta[KEY_RC_phieta]->Integral(1,h_mRcMass_phieta[KEY_RC_phieta]->GetNbinsX(),1,h_mRcMass_phieta[KEY_RC_phieta]->GetNbinsY());

    h_mMass_phiy_mcratio[KEY_phiy_mcratio] = (TH2F*)h_mMass_phiy[KEY_phiy]->Clone(KEY_phiy_mcratio.c_str()); // All other bins
    cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    double phiy_mcyield = h_mMcMass_phiy[KEY_MC_phiy]->Integral(1,h_mMcMass_phiy[KEY_MC_phiy]->GetNbinsX(),1,h_mMcMass_phiy[KEY_MC_phiy]->GetNbinsY());

    h_mMcMass_phiy[KEY_MC_phiy]->Scale(phiy_yield/phiy_mcyield);
    h_mMass_phiy_mcratio[KEY_phiy_mcratio]->Divide(h_mMcMass_phiy[KEY_MC_phiy]);

    h_mMass_phieta_mcratio[KEY_phieta_mcratio] = (TH2F*)h_mMass_phieta[KEY_phieta]->Clone(KEY_phieta_mcratio.c_str()); // All other bins
    cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    double phieta_mcyield = h_mMcMass_phieta[KEY_MC_phieta]->Integral(1,h_mMcMass_phieta[KEY_MC_phieta]->GetNbinsX(),1,h_mMcMass_phieta[KEY_MC_phieta]->GetNbinsY());

    h_mMcMass_phieta[KEY_MC_phieta]->Scale(phieta_yield/phieta_mcyield);
    h_mMass_phieta_mcratio[KEY_phieta_mcratio]->Divide(h_mMcMass_phieta[KEY_MC_phieta]);
   
    string KEY_phi = Form("kaon_phi_phipt_%d",i_pt_phi);
    string KEY_y   = Form("kaon_y_phipt_%d",i_pt_phi);
    string KEY_eta = Form("kaon_eta_phipt_%d",i_pt_phi);
    string KEY_RC_phi = Form("rckaon_phi_phipt_%d",i_pt_phi);
    string KEY_RC_y   = Form("rckaon_y_phipt_%d",i_pt_phi);
    string KEY_RC_eta = Form("rckaon_eta_phipt_%d",i_pt_phi);
    string KEY_phi_ratio = Form("kaonratio_phi_phipt_%d",i_pt_phi);
    string KEY_y_ratio   = Form("kaonratio_y_phipt_%d",i_pt_phi);
    string KEY_eta_ratio = Form("kaonratio_eta_phipt_%d",i_pt_phi);
    string KEY_MC_phi = Form("mckaon_phi_phipt_%d",i_pt_phi);
    string KEY_MC_y   = Form("mckaon_y_phipt_%d",i_pt_phi);
    string KEY_MC_eta = Form("mckaon_eta_phipt_%d",i_pt_phi);
    string KEY_phi_mcratio = Form("mckaonratio_phi_phipt_%d",i_pt_phi);
    string KEY_y_mcratio   = Form("mckaonratio_y_phipt_%d",i_pt_phi);
    string KEY_eta_mcratio = Form("mckaonratio_eta_phipt_%d",i_pt_phi);

    h_mMass_phi_ratio[KEY_phi_ratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_ratio.c_str()); // All other bins
    double phi_yield   = h_mMass_phi[KEY_phi]->Integral(1,h_mMass_phi[KEY_phi]->GetNbinsX());
    double phi_rcyield = h_mRcMass_phi[KEY_RC_phi]->Integral(1,h_mRcMass_phi[KEY_RC_phi]->GetNbinsX());

    cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
    h_mRcMass_phi[KEY_RC_phi]->Scale(phi_yield/phi_rcyield);
    h_mMass_phi_ratio[KEY_phi_ratio]->Divide(h_mRcMass_phi[KEY_RC_phi]);

    h_mMass_y_ratio[KEY_y_ratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_ratio.c_str()); // All other bins
    double y_yield   = h_mMass_y[KEY_y]->Integral(1,h_mMass_y[KEY_y]->GetNbinsX());
    double y_rcyield = h_mRcMass_y[KEY_RC_y]->Integral(1,h_mRcMass_y[KEY_RC_y]->GetNbinsX());

    cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
    h_mRcMass_y[KEY_RC_y]->Scale(y_yield/y_rcyield);
    h_mMass_y_ratio[KEY_y_ratio]->Divide(h_mRcMass_y[KEY_RC_y]);

    h_mMass_eta_ratio[KEY_eta_ratio] = (TH1F*)h_mMass_eta[KEY_eta]->Clone(KEY_eta_ratio.c_str()); // All other bins
    double eta_yield   = h_mMass_eta[KEY_eta]->Integral(1,h_mMass_eta[KEY_eta]->GetNbinsX());
    double eta_rcyield = h_mRcMass_eta[KEY_RC_eta]->Integral(1,h_mRcMass_eta[KEY_RC_eta]->GetNbinsX());

    cout << "Yield in pt bin y " << i_pt_phi << " is = " << eta_yield << endl;
    h_mRcMass_eta[KEY_RC_eta]->Scale(eta_yield/eta_rcyield);
    h_mMass_eta_ratio[KEY_eta_ratio]->Divide(h_mRcMass_eta[KEY_RC_eta]);

    h_mMass_phi_mcratio[KEY_phi_mcratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_mcratio.c_str()); // All other bins
    double phi_mcyield = h_mMcMass_phi[KEY_MC_phi]->Integral(1,h_mMcMass_phi[KEY_MC_phi]->GetNbinsX());

    cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
    h_mMcMass_phi[KEY_MC_phi]->Scale(phi_yield/phi_mcyield);
    h_mMass_phi_mcratio[KEY_phi_mcratio]->Divide(h_mMcMass_phi[KEY_MC_phi]);

    h_mMass_y_mcratio[KEY_y_mcratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_mcratio.c_str()); // All other bins
    double y_mcyield = h_mMcMass_y[KEY_MC_y]->Integral(1,h_mMcMass_y[KEY_MC_y]->GetNbinsX());

    cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
    h_mMcMass_y[KEY_MC_y]->Scale(y_yield/y_mcyield);
    h_mMass_y_mcratio[KEY_y_mcratio]->Divide(h_mMcMass_y[KEY_MC_y]);

    h_mMass_eta_mcratio[KEY_eta_mcratio] = (TH1F*)h_mMass_eta[KEY_eta]->Clone(KEY_eta_mcratio.c_str()); // All other bins
    double eta_mcyield = h_mMcMass_eta[KEY_MC_eta]->Integral(1,h_mMcMass_eta[KEY_MC_eta]->GetNbinsX());

    cout << "Yield in pt bin y " << i_pt_phi << " is = " << eta_yield << endl;
    h_mMcMass_eta[KEY_MC_eta]->Scale(eta_yield/eta_mcyield);
    h_mMass_eta_mcratio[KEY_eta_mcratio]->Divide(h_mMcMass_eta[KEY_MC_eta]);

    /////////// DELTA //////////////////////
    KEY_phiy_ratio = Form("kaonratio_delta_phiy_phipt_%d",i_pt_phi);
    KEY_phiy_mcratio = Form("mckaonratio_delta_phiy_phipt_%d",i_pt_phi);
    KEY_phiy = Form("kaon_delta_phiy_phipt_%d",i_pt_phi);
    KEY_RC_phiy = Form("rckaon_delta_phiy_phipt_%d",i_pt_phi);
    KEY_MC_phiy = Form("mckaon_delta_phiy_phipt_%d",i_pt_phi);
     
    h_mMass_phiy_ratio[KEY_phiy_ratio] = (TH2F*)h_mMass_phiy[KEY_phiy]->Clone(KEY_phiy_ratio.c_str()); // All other bins
    phiy_yield   = h_mMass_phiy[KEY_phiy]->Integral(0,-1,0,-1);
    cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    phiy_rcyield = h_mRcMass_phiy[KEY_RC_phiy]->Integral(0,-1,0,-1);

    h_mRcMass_phiy[KEY_RC_phiy]->Scale(phiy_yield/phiy_rcyield);
    h_mMass_phiy_ratio[KEY_phiy_ratio]->Divide(h_mRcMass_phiy[KEY_RC_phiy]);

    h_mMass_phiy_mcratio[KEY_phiy_mcratio] = (TH2F*)h_mMass_phiy[KEY_phiy]->Clone(KEY_phiy_mcratio.c_str()); // All other bins
    cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    phiy_mcyield = h_mMcMass_phiy[KEY_MC_phiy]->Integral(0,-1,0,-1);

    h_mMcMass_phiy[KEY_MC_phiy]->Scale(phiy_yield/phiy_mcyield);
    h_mMass_phiy_mcratio[KEY_phiy_mcratio]->Divide(h_mMcMass_phiy[KEY_MC_phiy]);
   
    KEY_phi = Form("kaon_delta_phi_phipt_%d",i_pt_phi);
    KEY_y   = Form("kaon_delta_y_phipt_%d",i_pt_phi);
    KEY_RC_phi = Form("rckaon_delta_phi_phipt_%d",i_pt_phi);
    KEY_RC_y   = Form("rckaon_delta_y_phipt_%d",i_pt_phi);
    KEY_MC_phi = Form("mckaon_delta_phi_phipt_%d",i_pt_phi);
    KEY_MC_y   = Form("mckaon_delta_y_phipt_%d",i_pt_phi);
    KEY_phi_ratio = Form("kaonratio_delta_phi_phipt_%d",i_pt_phi);
    KEY_y_ratio   = Form("kaonratio_delta_y_phipt_%d",i_pt_phi);
    KEY_phi_mcratio = Form("mckaonratio_delta_phi_phipt_%d",i_pt_phi);
    KEY_y_mcratio   = Form("mckaonratio_delta_y_phipt_%d",i_pt_phi);

    h_mMass_phi_ratio[KEY_phi_ratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_ratio.c_str()); // All other bins
    phi_yield   = h_mMass_phi[KEY_phi]->Integral(0,-1);
    phi_rcyield = h_mRcMass_phi[KEY_RC_phi]->Integral(0,-1);

    cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
    h_mRcMass_phi[KEY_RC_phi]->Scale(phi_yield/phi_rcyield);
    h_mMass_phi_ratio[KEY_phi_ratio]->Divide(h_mRcMass_phi[KEY_RC_phi]);

    h_mMass_phi_mcratio[KEY_phi_mcratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_mcratio.c_str()); // All other bins
    phi_yield   = h_mMass_phi[KEY_phi]->Integral(0,-1);
    phi_mcyield = h_mMcMass_phi[KEY_MC_phi]->Integral(0,-1);

    cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
    h_mMcMass_phi[KEY_MC_phi]->Scale(phi_yield/phi_mcyield);
    h_mMass_phi_mcratio[KEY_phi_mcratio]->Divide(h_mMcMass_phi[KEY_MC_phi]);

    h_mMass_y_ratio[KEY_y_ratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_ratio.c_str()); // All other bins
    y_yield   = h_mMass_y[KEY_y]->Integral(0,-1);
    y_rcyield = h_mRcMass_y[KEY_RC_y]->Integral(0,-1);

    cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
    h_mRcMass_y[KEY_RC_y]->Scale(y_yield/y_rcyield);
    h_mMass_y_ratio[KEY_y_ratio]->Divide(h_mRcMass_y[KEY_RC_y]);

    h_mMass_y_mcratio[KEY_y_mcratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_mcratio.c_str()); // All other bins
    y_yield   = h_mMass_y[KEY_y]->Integral(0,-1);
    y_mcyield = h_mMcMass_y[KEY_MC_y]->Integral(0,-1);

    cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
    h_mMcMass_y[KEY_MC_y]->Scale(y_yield/y_mcyield);
    h_mMass_y_mcratio[KEY_y_mcratio]->Divide(h_mMcMass_y[KEY_MC_y]);

    /////////// DELTA R //////////////////////
    KEY_phiy_ratio = Form("kaonratio_deltarphi_phiy_phipt_%d",i_pt_phi);
    KEY_phiy_mcratio = Form("mckaonratio_deltarphi_phiy_phipt_%d",i_pt_phi);
    KEY_phiy = Form("kaon_deltarphi_phiy_phipt_%d",i_pt_phi);
    KEY_RC_phiy = Form("rckaon_deltarphi_phiy_phipt_%d",i_pt_phi);
    KEY_MC_phiy = Form("mckaon_deltarphi_phiy_phipt_%d",i_pt_phi);
     
    h_mMass_phiy_ratio[KEY_phiy_ratio] = (TH2F*)h_mMass_phiy[KEY_phiy]->Clone(KEY_phiy_ratio.c_str()); // All other bins
    phiy_yield   = h_mMass_phiy[KEY_phiy]->Integral(0,-1,0,-1);
    cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    phiy_rcyield = h_mRcMass_phiy[KEY_RC_phiy]->Integral(0,-1,0,-1);

    h_mRcMass_phiy[KEY_RC_phiy]->Scale(phiy_yield/phiy_rcyield);
    h_mMass_phiy_ratio[KEY_phiy_ratio]->Divide(h_mRcMass_phiy[KEY_RC_phiy]);

    h_mMass_phiy_mcratio[KEY_phiy_mcratio] = (TH2F*)h_mMass_phiy[KEY_phiy]->Clone(KEY_phiy_mcratio.c_str()); // All other bins
    cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    phiy_mcyield = h_mMcMass_phiy[KEY_MC_phiy]->Integral(0,-1,0,-1);

    h_mMcMass_phiy[KEY_MC_phiy]->Scale(phiy_yield/phiy_mcyield);
    h_mMass_phiy_mcratio[KEY_phiy_mcratio]->Divide(h_mMcMass_phiy[KEY_MC_phiy]);
   
    KEY_phi = Form("kaon_deltarphi_phi_phipt_%d",i_pt_phi);
    KEY_y   = Form("kaon_deltarphi_y_phipt_%d",i_pt_phi);
    KEY_RC_phi = Form("rckaon_deltarphi_phi_phipt_%d",i_pt_phi);
    KEY_RC_y   = Form("rckaon_deltarphi_y_phipt_%d",i_pt_phi);
    KEY_MC_phi = Form("mckaon_deltarphi_phi_phipt_%d",i_pt_phi);
    KEY_MC_y   = Form("mckaon_deltarphi_y_phipt_%d",i_pt_phi);
    KEY_phi_ratio = Form("kaonratio_deltarphi_phi_phipt_%d",i_pt_phi);
    KEY_y_ratio   = Form("kaonratio_deltarphi_y_phipt_%d",i_pt_phi);
    KEY_phi_mcratio = Form("mckaonratio_deltarphi_phi_phipt_%d",i_pt_phi);
    KEY_y_mcratio   = Form("mckaonratio_deltarphi_y_phipt_%d",i_pt_phi);

    h_mMass_phi_ratio[KEY_phi_ratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_ratio.c_str()); // All other bins
    phi_yield   = h_mMass_phi[KEY_phi]->Integral(0,-1);
    phi_rcyield = h_mRcMass_phi[KEY_RC_phi]->Integral(0,-1);

    cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
    h_mRcMass_phi[KEY_RC_phi]->Scale(phi_yield/phi_rcyield);
    h_mMass_phi_ratio[KEY_phi_ratio]->Divide(h_mRcMass_phi[KEY_RC_phi]);

    h_mMass_phi_mcratio[KEY_phi_mcratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_mcratio.c_str()); // All other bins
    phi_yield   = h_mMass_phi[KEY_phi]->Integral(0,-1);
    phi_mcyield = h_mMcMass_phi[KEY_MC_phi]->Integral(0,-1);

    cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
    h_mMcMass_phi[KEY_MC_phi]->Scale(phi_yield/phi_mcyield);
    h_mMass_phi_mcratio[KEY_phi_mcratio]->Divide(h_mMcMass_phi[KEY_MC_phi]);

    h_mMass_y_ratio[KEY_y_ratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_ratio.c_str()); // All other bins
    y_yield   = h_mMass_y[KEY_y]->Integral(0,-1);
    y_rcyield = h_mRcMass_y[KEY_RC_y]->Integral(0,-1);

    cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
    h_mRcMass_y[KEY_RC_y]->Scale(y_yield/y_rcyield);
    h_mMass_y_ratio[KEY_y_ratio]->Divide(h_mRcMass_y[KEY_RC_y]);

    h_mMass_y_mcratio[KEY_y_mcratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_mcratio.c_str()); // All other bins
    y_yield   = h_mMass_y[KEY_y]->Integral(0,-1);
    y_mcyield = h_mMcMass_y[KEY_MC_y]->Integral(0,-1);

    cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
    h_mMcMass_y[KEY_MC_y]->Scale(y_yield/y_mcyield);
    h_mMass_y_mcratio[KEY_y_mcratio]->Divide(h_mMcMass_y[KEY_MC_y]);


    /////////// DELTA phi-K+//////////////////////
    //KEY_phiy_ratio = Form("kaonplusphiratio_delta_phiy_phipt_%d",i_pt_phi);
    //KEY_phiy = Form("kaonplusphi_delta_phiy_phipt_%d",i_pt_phi);
    //KEY_RC_phiy = Form("rckaonplusphi_delta_phiy_phipt_%d",i_pt_phi);
    // 
    //h_mMass_phiy_ratio[KEY_phiy_ratio] = (TH2F*)h_mMass_phiy[KEY_phiy]->Clone(KEY_phiy_ratio.c_str()); // All other bins
    //phiy_yield   = h_mMass_phiy[KEY_phiy]->Integral(0,-1,0,-1);
    //cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    //phiy_rcyield = h_mRcMass_phiy[KEY_RC_phiy]->Integral(0,-1,0,-1);

    //h_mRcMass_phiy[KEY_RC_phiy]->Scale(phiy_yield/phiy_rcyield);
    //h_mMass_phiy_ratio[KEY_phiy_ratio]->Divide(h_mRcMass_phiy[KEY_RC_phiy]);
   
    //KEY_phi = Form("kaonplusphi_delta_phi_phipt_%d",i_pt_phi);
    //KEY_y   = Form("kaonplusphi_delta_y_phipt_%d",i_pt_phi);
    //KEY_RC_phi = Form("rckaonplusphi_delta_phi_phipt_%d",i_pt_phi);
    //KEY_RC_y   = Form("rckaonplusphi_delta_y_phipt_%d",i_pt_phi);
    //KEY_phi_ratio = Form("kaonplusphiratio_delta_phi_phipt_%d",i_pt_phi);
    //KEY_y_ratio   = Form("kaonplusphiratio_delta_y_phipt_%d",i_pt_phi);

    //h_mMass_phi_ratio[KEY_phi_ratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_ratio.c_str()); // All other bins
    //phi_yield   = h_mMass_phi[KEY_phi]->Integral(0,-1);
    //phi_rcyield = h_mRcMass_phi[KEY_RC_phi]->Integral(0,-1);

    //cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
    //h_mRcMass_phi[KEY_RC_phi]->Scale(phi_yield/phi_rcyield);
    //h_mMass_phi_ratio[KEY_phi_ratio]->Divide(h_mRcMass_phi[KEY_RC_phi]);

    //h_mMass_y_ratio[KEY_y_ratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_ratio.c_str()); // All other bins
    //y_yield   = h_mMass_y[KEY_y]->Integral(0,-1);
    //y_rcyield = h_mRcMass_y[KEY_RC_y]->Integral(0,-1);

    //cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
    //h_mRcMass_y[KEY_RC_y]->Scale(y_yield/y_rcyield);
    //h_mMass_y_ratio[KEY_y_ratio]->Divide(h_mRcMass_y[KEY_RC_y]);

    /////////// DELTA phi-K-//////////////////////
    //KEY_phiy_ratio = Form("kaonminusphiratio_delta_phiy_phipt_%d",i_pt_phi);
    //KEY_phiy = Form("kaonminusphi_delta_phiy_phipt_%d",i_pt_phi);
    //KEY_RC_phiy = Form("rckaonminusphi_delta_phiy_phipt_%d",i_pt_phi);
    // 
    //h_mMass_phiy_ratio[KEY_phiy_ratio] = (TH2F*)h_mMass_phiy[KEY_phiy]->Clone(KEY_phiy_ratio.c_str()); // All other bins
    //phiy_yield   = h_mMass_phiy[KEY_phiy]->Integral(0,-1,0,-1);
    //cout << "Yield in pt bin phi and y " << i_pt_phi << " is = " << phiy_yield << endl;
    //phiy_rcyield = h_mRcMass_phiy[KEY_RC_phiy]->Integral(0,-1,0,-1);

    //h_mRcMass_phiy[KEY_RC_phiy]->Scale(phiy_yield/phiy_rcyield);
    //h_mMass_phiy_ratio[KEY_phiy_ratio]->Divide(h_mRcMass_phiy[KEY_RC_phiy]);
   
    //KEY_phi = Form("kaonminusphi_delta_phi_phipt_%d",i_pt_phi);
    //KEY_y   = Form("kaonminusphi_delta_y_phipt_%d",i_pt_phi);
    //KEY_RC_phi = Form("rckaonminusphi_delta_phi_phipt_%d",i_pt_phi);
    //KEY_RC_y   = Form("rckaonminusphi_delta_y_phipt_%d",i_pt_phi);
    //KEY_phi_ratio = Form("kaonminusphiratio_delta_phi_phipt_%d",i_pt_phi);
    //KEY_y_ratio   = Form("kaonminusphiratio_delta_y_phipt_%d",i_pt_phi);

    //h_mMass_phi_ratio[KEY_phi_ratio] = (TH1F*)h_mMass_phi[KEY_phi]->Clone(KEY_phi_ratio.c_str()); // All other bins
    //phi_yield   = h_mMass_phi[KEY_phi]->Integral(0,-1);
    //phi_rcyield = h_mRcMass_phi[KEY_RC_phi]->Integral(0,-1);

    //cout << "Yield in pt bin phi " << i_pt_phi << " is = " << phi_yield << endl;
    //h_mRcMass_phi[KEY_RC_phi]->Scale(phi_yield/phi_rcyield);
    //h_mMass_phi_ratio[KEY_phi_ratio]->Divide(h_mRcMass_phi[KEY_RC_phi]);

    //h_mMass_y_ratio[KEY_y_ratio] = (TH1F*)h_mMass_y[KEY_y]->Clone(KEY_y_ratio.c_str()); // All other bins
    //y_yield   = h_mMass_y[KEY_y]->Integral(0,-1);
    //y_rcyield = h_mRcMass_y[KEY_RC_y]->Integral(0,-1);

    //cout << "Yield in pt bin y " << i_pt_phi << " is = " << y_yield << endl;
    //h_mRcMass_y[KEY_RC_y]->Scale(y_yield/y_rcyield);
    //h_mMass_y_ratio[KEY_y_ratio]->Divide(h_mRcMass_y[KEY_RC_y]);
  } 

  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++) // pt bin 
  {
    string KEY_pT_ratio = Form("kaonratio_pt_phipt_%d",i_pt_phi);
    g_mKaon_pT_ratio[KEY_pT_ratio] = new TGraphAsymmErrors();

    string KEY_pT = Form("kaon_pt_phipt_%d",i_pt_phi);
    g_mKaon_pT[KEY_pT] = new TGraphAsymmErrors();

    string KEY_RC_pT = Form("rckaon_pt_phipt_%d",i_pt_phi);
    g_mRcKaon_pT[KEY_RC_pT] = new TGraphAsymmErrors();

    double scaling_yields = 0.0;
    double scaling_rcyields = 0.0;

    for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
    {
      double total_yield = 0.0;
      double total_errorerror = 0.0;
      double rctotal_yield = 0.0;
      double rctotal_errorerror = 0.0;
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        string KEY_SM = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        double error;
        double yield = h_mMass_SM[KEY_SM]->IntegralAndError(0,-1,0,-1,error);
        scaling_yields += yield;
        total_yield += yield;
        total_errorerror += (error*error);

        string KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        double rcerror;
        double rcyield = h_mMass_RC[KEY_RC]->IntegralAndError(0,-1,0,-1,rcerror);
        scaling_rcyields += rcyield;
        rctotal_yield += rcyield;
        rctotal_errorerror += (rcerror*rcerror);
      }
      double pt_mean = (data_constants::kaon_pt_low[i_pt_kaon]+data_constants::kaon_pt_high[i_pt_kaon])/2.0;
      double total_error = TMath::Sqrt(total_errorerror);
      double rctotal_error = TMath::Sqrt(rctotal_errorerror);
      g_mKaon_pT[KEY_pT]->SetPoint(i_pt_kaon,pt_mean,total_yield);
      g_mKaon_pT[KEY_pT]->SetPointError(i_pt_kaon,0.0,0.0,total_error,total_error);
      g_mRcKaon_pT[KEY_RC_pT]->SetPoint(i_pt_kaon,pt_mean,rctotal_yield);
      g_mRcKaon_pT[KEY_RC_pT]->SetPointError(i_pt_kaon,0.0,0.0,rctotal_error,rctotal_error);

    }
    // Scale RC to the data
    g_mRcKaon_pT[KEY_RC_pT]->Scale(scaling_yields/scaling_rcyields,"y");

    for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
    {
      double pt, yield, yielderror;
      double rcyield, rcyielderror;
 
      g_mKaon_pT[KEY_pT]->GetPoint(i_pt_kaon,pt,yield);
      yielderror = g_mKaon_pT[KEY_pT]->GetErrorYhigh(i_pt_kaon);
      g_mRcKaon_pT[KEY_RC_pT]->GetPoint(i_pt_kaon,pt,rcyield);
      rcyielderror = g_mRcKaon_pT[KEY_RC_pT]->GetErrorYhigh(i_pt_kaon);

      cout << "i_pt_kaon = " << i_pt_kaon << ", yielderror = " << yielderror << ", rcyielderror = " << rcyielderror << endl;
      
      double ratio, ratioerror;
      if(rcyield == 0) 
      {
        ratio = 1;
        ratioerror = 0;
      }
      else
      {
        ratio = yield/(rcyield);
        ratioerror = ratio*TMath::Sqrt(yielderror*yielderror/yield/yield + rcyielderror*rcyielderror/rcyield/rcyield);
      }
      g_mKaon_pT_ratio[KEY_pT_ratio]->SetPoint(i_pt_kaon,pt,ratio);
      g_mKaon_pT_ratio[KEY_pT_ratio]->SetPointError(i_pt_kaon,0.0,0.0,ratioerror,ratioerror);
      cout << "i_pt_kaon = " << i_pt_kaon << ", ratio = " << ratio << ", ratioerror = " << ratioerror << endl;
    } 
  }

  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++) // pt bin 
  {
    string KEY_pT_mcratio = Form("kaonratio_pt_phipt_%d",i_pt_phi);
    g_mKaon_pT_mcratio[KEY_pT_mcratio] = new TGraphAsymmErrors();

    string KEY_pT = Form("kaon_pt_phipt_%d",i_pt_phi);
    g_mKaon_pT[KEY_pT] = new TGraphAsymmErrors();

    string KEY_MC_pT = Form("mckaon_pt_phipt_%d",i_pt_phi);
    g_mMcKaon_pT[KEY_MC_pT] = new TGraphAsymmErrors();

    double scaling_yields = 0.0;
    double scaling_mcyields = 0.0;

    for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
    {
      double total_yield = 0.0;
      double total_errorerror = 0.0;
      double mctotal_yield = 0.0;
      double mctotal_errorerror = 0.0;
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        string KEY_SM = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        double error;
        double yield = h_mMass_SM[KEY_SM]->IntegralAndError(0,-1,0,-1,error);
        scaling_yields += yield;
        total_yield += yield;
        total_errorerror += (error*error);

        string KEY_MC = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        double mcerror;
        double mcyield = h_mMass_MC[KEY_MC]->IntegralAndError(0,-1,0,-1,mcerror);
        scaling_mcyields += mcyield;
        mctotal_yield += mcyield;
        mctotal_errorerror += (mcerror*mcerror);
      }
      double pt_mean = (data_constants::kaon_pt_low[i_pt_kaon]+data_constants::kaon_pt_high[i_pt_kaon])/2.0;
      double total_error = TMath::Sqrt(total_errorerror);
      double mctotal_error = TMath::Sqrt(mctotal_errorerror);
      g_mKaon_pT[KEY_pT]->SetPoint(i_pt_kaon,pt_mean,total_yield);
      g_mKaon_pT[KEY_pT]->SetPointError(i_pt_kaon,0.0,0.0,total_error,total_error);
      g_mMcKaon_pT[KEY_MC_pT]->SetPoint(i_pt_kaon,pt_mean,mctotal_yield);
      g_mMcKaon_pT[KEY_MC_pT]->SetPointError(i_pt_kaon,0.0,0.0,mctotal_error,mctotal_error);

    }
    // Scale MC to the data
    g_mMcKaon_pT[KEY_MC_pT]->Scale(scaling_yields/scaling_mcyields,"y");

    for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
    {
      double pt, yield, yielderror;
      double mcyield, mcyielderror;
 
      g_mKaon_pT[KEY_pT]->GetPoint(i_pt_kaon,pt,yield);
      yielderror = g_mKaon_pT[KEY_pT]->GetErrorYhigh(i_pt_kaon);
      g_mMcKaon_pT[KEY_MC_pT]->GetPoint(i_pt_kaon,pt,mcyield);
      mcyielderror = g_mMcKaon_pT[KEY_MC_pT]->GetErrorYhigh(i_pt_kaon);

      cout << "i_pt_kaon = " << i_pt_kaon << ", yielderror = " << yielderror << ", mcyielderror = " << mcyielderror << endl;
      
      double ratio, ratioerror;
      if(mcyield == 0) 
      {
        ratio = 1;
        ratioerror = 0;
      }
      else
      {
        ratio = yield/(mcyield);
        ratioerror = ratio*TMath::Sqrt(yielderror*yielderror/yield/yield + mcyielderror*mcyielderror/mcyield/mcyield);
      }
      g_mKaon_pT_mcratio[KEY_pT_mcratio]->SetPoint(i_pt_kaon,pt,ratio);
      g_mKaon_pT_mcratio[KEY_pT_mcratio]->SetPointError(i_pt_kaon,0.0,0.0,ratioerror,ratioerror);
      cout << "i_pt_kaon = " << i_pt_kaon << ", ratio = " << ratio << ", ratioerror = " << ratioerror << endl;
    } 
  }
  



  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);     
  c1->Divide(2,2);
  for(int i = 0; i < 4; i++)
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

  std::string kaonphilabel[5] = {"#phi","#phi_{K+}-#phi_{K-}","#phi_{#DeltaR}","#phi_{#phi}-#phi_{K+}","#phi_{#phi}-#phi_{K-}"};
  std::string kaonetalabel[5] = {"y","#eta_{K+}-#eta_{K-}","#DeltaR=#sqrt{(#Delta#eta)^{2}+(#Delta#phi)^{2}}","#eta_{#phi}-#eta_{K+}","#eta_{#phi}-#eta_{K-}"};
  std::string kaonlabel[5] = {"","","","plusphi","minusphi"};
  std::string kaonextralabel[5] = {"","_delta","_deltarphi","_delta","_delta"};
  
  if(!deltaonly)
  {
    for(int iset = 0; iset < 1; iset++)
    {
      std::string outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_SE.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
      std::string outputstart = Form("%s[",outputname.c_str());
      std::string outputstop  = Form("%s]",outputname.c_str());
      craw->Print(outputstart.c_str());

      for(int i = 0; i < 4; i++)
      {
        for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
        {
          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
          {
            craw->cd(i_phi+1);
            string KEY_SE = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_SE",i,i_phi,i_pt_kaon);
            h_mMass_SE[KEY_SE]->SetTitle(Form("SE %1.1f<p_{T,#phi}<%1.1f, %1.1f<p_{T,K}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],data_constants::kaon_pt_low[i_pt_kaon],data_constants::kaon_pt_high[i_pt_kaon],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
            h_mMass_SE[KEY_SE]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
            h_mMass_SE[KEY_SE]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

            h_mMass_SE[KEY_SE]->Draw("colz");
          }
          craw->Update();
          craw->Print(outputname.c_str());
        }
      }
      craw->Print(outputstop.c_str());

      outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_ME.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
      outputstart = Form("%s[",outputname.c_str());
      outputstop  = Form("%s]",outputname.c_str());
      craw->Print(outputstart.c_str());

      for(int i = 0; i < 4; i++)
      {
        for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
        {
          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
          {
            craw->cd(i_phi+1);
            string KEY_ME = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_ME",i,i_phi,i_pt_kaon);
            h_mMass_ME[KEY_ME]->SetTitle(Form("ME %1.1f<p_{T,#phi}<%1.1f, %1.1f<p_{T,K}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],data_constants::kaon_pt_low[i_pt_kaon],data_constants::kaon_pt_high[i_pt_kaon],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
            h_mMass_ME[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
            h_mMass_ME[KEY_ME]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

            h_mMass_ME[KEY_ME]->Draw("colz");
          }
          craw->Update();
          craw->Print(outputname.c_str());
        }
      }
      craw->Print(outputstop.c_str());

      outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_SM.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
      outputstart = Form("%s[",outputname.c_str());
      outputstop  = Form("%s]",outputname.c_str());
      craw->Print(outputstart.c_str());

      for(int i = 0; i < 4; i++)
      {
        for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
        {
          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
          {
            craw->cd(i_phi+1);
            string KEY_ME = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i,i_phi,i_pt_kaon);
            h_mMass_SM[KEY_ME]->SetTitle(Form("SE-ME %1.1f<p_{T,#phi}<%1.1f, %1.1f<p_{T,K}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],data_constants::kaon_pt_low[i_pt_kaon],data_constants::kaon_pt_high[i_pt_kaon],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
            h_mMass_SM[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
            h_mMass_SM[KEY_ME]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

            h_mMass_SM[KEY_ME]->Draw("colz");
          }
          craw->Update();
          craw->Print(outputname.c_str());
        }
      }
      craw->Print(outputstop.c_str());

      outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_MC.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
      outputstart = Form("%s[",outputname.c_str());
      outputstop  = Form("%s]",outputname.c_str());
      craw->Print(outputstart.c_str());

      for(int i = 0; i < 4; i++)
      {
        for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
        {
          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
          {
            craw->cd(i_phi+1);
            string KEY_ME = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i,i_phi,i_pt_kaon);
            h_mMass_MC[KEY_ME]->SetTitle(Form("MC %1.1f<p_{T,#phi}<%1.1f, %1.1f<p_{T,K}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],data_constants::kaon_pt_low[i_pt_kaon],data_constants::kaon_pt_high[i_pt_kaon],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
            h_mMass_MC[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
            h_mMass_MC[KEY_ME]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

            h_mMass_MC[KEY_ME]->Draw("colz");
          }
          craw->Update();
          craw->Print(outputname.c_str());
        }
      }
      craw->Print(outputstop.c_str());

      outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_RC.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
      outputstart = Form("%s[",outputname.c_str());
      outputstop  = Form("%s]",outputname.c_str());
      craw->Print(outputstart.c_str());

      for(int i = 0; i < 4; i++)
      {
        for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
        {
          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
          {
            craw->cd(i_phi+1);
            string KEY_ME = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i,i_phi,i_pt_kaon);
            h_mMass_RC[KEY_ME]->SetTitle(Form("RC %1.1f<p_{T,#phi}<%1.1f, %1.1f<p_{T,K}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],data_constants::kaon_pt_low[i_pt_kaon],data_constants::kaon_pt_high[i_pt_kaon],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
            h_mMass_RC[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
            h_mMass_RC[KEY_ME]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

            h_mMass_RC[KEY_ME]->Draw("colz");
          }
          craw->Update();
          craw->Print(outputname.c_str());
        }
      }
      craw->Print(outputstop.c_str());

      outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phieta_SE.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
      outputstart = Form("%s[",outputname.c_str());
      outputstop  = Form("%s]",outputname.c_str());
      craw->Print(outputstart.c_str());

      for(int i = 0; i < 4; i++)
      {
        for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
        {
          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
          {
            craw->cd(i_phi+1);
            string KEY_SE = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta_SE",i,i_phi,i_pt_kaon);
            h_mMass_SE[KEY_SE]->SetTitle(Form("SE %1.1f<p_{T,#phi}<%1.1f, %1.1f<p_{T,K}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],data_constants::kaon_pt_low[i_pt_kaon],data_constants::kaon_pt_high[i_pt_kaon],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
            h_mMass_SE[KEY_SE]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
            h_mMass_SE[KEY_SE]->GetYaxis()->SetTitle("#eta");

            h_mMass_SE[KEY_SE]->Draw("colz");
          }
          craw->Update();
          craw->Print(outputname.c_str());
        }
      }
      craw->Print(outputstop.c_str());

      outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phieta_MC.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
      outputstart = Form("%s[",outputname.c_str());
      outputstop  = Form("%s]",outputname.c_str());
      craw->Print(outputstart.c_str());

      for(int i = 0; i < 4; i++)
      {
        for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
        {
          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
          {
            craw->cd(i_phi+1);
            string KEY_SE = Form("mckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i,i_phi,i_pt_kaon);
            h_mMass_MC[KEY_SE]->SetTitle(Form("MC %1.1f<p_{T,#phi}<%1.1f, %1.1f<p_{T,K}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],data_constants::kaon_pt_low[i_pt_kaon],data_constants::kaon_pt_high[i_pt_kaon],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
            h_mMass_MC[KEY_SE]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
            h_mMass_MC[KEY_SE]->GetYaxis()->SetTitle("#eta");

            h_mMass_MC[KEY_SE]->Draw("colz");
          }
          craw->Update();
          craw->Print(outputname.c_str());
        }
      }
      craw->Print(outputstop.c_str());

      outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phieta_RC.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
      outputstart = Form("%s[",outputname.c_str());
      outputstop  = Form("%s]",outputname.c_str());
      craw->Print(outputstart.c_str());

      for(int i = 0; i < 4; i++)
      {
        for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
        {
          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
          {
            craw->cd(i_phi+1);
            string KEY_SE = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i,i_phi,i_pt_kaon);
            h_mMass_RC[KEY_SE]->SetTitle(Form("RC %1.1f<p_{T,#phi}<%1.1f, %1.1f<p_{T,K}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],data_constants::kaon_pt_low[i_pt_kaon],data_constants::kaon_pt_high[i_pt_kaon],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
            h_mMass_RC[KEY_SE]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
            h_mMass_RC[KEY_SE]->GetYaxis()->SetTitle("#eta");

            h_mMass_RC[KEY_SE]->Draw("colz");
          }
          craw->Update();
          craw->Print(outputname.c_str());
        }
      }
      craw->Print(outputstop.c_str());

      outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phieta_ME.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
      outputstart = Form("%s[",outputname.c_str());
      outputstop  = Form("%s]",outputname.c_str());
      craw->Print(outputstart.c_str());

      for(int i = 0; i < 4; i++)
      {
        for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
        {
          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
          {
            craw->cd(i_phi+1);
            string KEY_ME = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta_ME",i,i_phi,i_pt_kaon);
            h_mMass_ME[KEY_ME]->SetTitle(Form("ME %1.1f<p_{T,#phi}<%1.1f, %1.1f<p_{T,K}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],data_constants::kaon_pt_low[i_pt_kaon],data_constants::kaon_pt_high[i_pt_kaon],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
            h_mMass_ME[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
            h_mMass_ME[KEY_ME]->GetYaxis()->SetTitle("#eta");

            h_mMass_ME[KEY_ME]->Draw("colz");
          }
          craw->Update();
          craw->Print(outputname.c_str());
        }
      }
      craw->Print(outputstop.c_str());

      outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phieta_SM.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
      outputstart = Form("%s[",outputname.c_str());
      outputstop  = Form("%s]",outputname.c_str());
      craw->Print(outputstart.c_str());

      for(int i = 0; i < 4; i++)
      {
        for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
        {
          for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
          {
            craw->cd(i_phi+1);
            string KEY_ME = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i,i_phi,i_pt_kaon);
            h_mMass_SM[KEY_ME]->SetTitle(Form("SE-ME %1.1f<p_{T,#phi}<%1.1f, %1.1f<p_{T,K}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],data_constants::kaon_pt_low[i_pt_kaon],data_constants::kaon_pt_high[i_pt_kaon],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
            h_mMass_SM[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
            h_mMass_SM[KEY_ME]->GetYaxis()->SetTitle("#eta");

            h_mMass_SM[KEY_ME]->Draw("colz");
          }
          craw->Update();
          craw->Print(outputname.c_str());
        }
      }
      craw->Print(outputstop.c_str());
    }
  }

  for(int iset = 1; iset < 3; iset++)
  {
    std::string outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_SE.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    std::string outputstart = Form("%s[",outputname.c_str());
    std::string outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_SE = Form("kaon%s%s_phipt_%d_cos2phistarphi_%d_SE",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_SE[KEY_SE]->SetTitle(Form("SE %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_SE[KEY_SE]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_SE[KEY_SE]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

        h_mMass_SE[KEY_SE]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_ME.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("kaon%s%s_phipt_%d_cos2phistarphi_%d_ME",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_ME[KEY_ME]->SetTitle(Form("ME %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_ME[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_ME[KEY_ME]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

        h_mMass_ME[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_SM.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("kaon%s%s_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_SM[KEY_ME]->SetTitle(Form("SE-ME %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_SM[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_SM[KEY_ME]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

        h_mMass_SM[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_IndividualCos2PhiStarPhiBins_RC.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("rckaon%s%s_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_RC[KEY_ME]->SetTitle(Form("RC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_RC[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_RC[KEY_ME]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

        h_mMass_RC[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phi_IndividualCos2PhiStarPhiBins_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      TLegend *legphi = new TLegend(0.825,0.55,0.975,0.75);
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_phi = Form("kaon%s%s_phi_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        string KEY_RC_phi = Form("rckaon%s%s_phi_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_phi[KEY_phi]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phi[KEY_phi]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
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

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phi_IndividualCos2PhiStarPhiBins_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      TLegend *legphi = new TLegend(0.825,0.55,0.975,0.75);
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_phi = Form("kaon%s%s_phi_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        string KEY_MC_phi = Form("mckaon%s%s_phi_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_phi[KEY_phi]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phi[KEY_phi]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
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

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_y_IndividualCos2PhiStarPhiBins_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      TLegend *legy = new TLegend(0.825,0.55,0.975,0.75);
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_y = Form("kaon%s%s_y_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        string KEY_RC_y = Form("rckaon%s%s_y_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_y[KEY_y]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_y[KEY_y]->GetXaxis()->SetTitle(kaonetalabel[iset].c_str());
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

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_y_IndividualCos2PhiStarPhiBins_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      TLegend *legy = new TLegend(0.825,0.55,0.975,0.75);
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_y = Form("kaon%s%s_y_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        string KEY_MC_y = Form("mckaon%s%s_y_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_y[KEY_y]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_y[KEY_y]->GetXaxis()->SetTitle(kaonetalabel[iset].c_str());
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

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phi_IndividualCos2PhiStarPhiBins_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("kaonratio%s%s_phi_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_phi_ratio[KEY_ME]->SetTitle(Form("Data/RC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phi_ratio[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_phi_ratio[KEY_ME]->GetYaxis()->SetTitle("Data/RC");

        h_mMass_phi_ratio[KEY_ME]->Draw("pE");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_y_IndividualCos2PhiStarPhiBins_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("kaonratio%s%s_y_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_y_ratio[KEY_ME]->SetTitle(Form("Data/RC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_y_ratio[KEY_ME]->GetXaxis()->SetTitle(kaonetalabel[iset].c_str());
        h_mMass_y_ratio[KEY_ME]->GetYaxis()->SetTitle("Data/RC");

        h_mMass_y_ratio[KEY_ME]->Draw("pE");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_IndividualCos2PhiStarPhiBins_MC.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("mckaon%s%s_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_MC[KEY_ME]->SetTitle(Form("MC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_MC[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_MC[KEY_ME]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

        h_mMass_MC[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phi_IndividualCos2PhiStarPhiBins_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("mckaonratio%s%s_phi_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_phi_mcratio[KEY_ME]->SetTitle(Form("Data/MC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phi_mcratio[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_phi_mcratio[KEY_ME]->GetYaxis()->SetTitle("Data/MC");

        h_mMass_phi_mcratio[KEY_ME]->Draw("pE");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_y_IndividualCos2PhiStarPhiBins_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("mckaonratio%s%s_y_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_y_mcratio[KEY_ME]->SetTitle(Form("Data/MC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_y_mcratio[KEY_ME]->GetXaxis()->SetTitle(kaonetalabel[iset].c_str());
        h_mMass_y_mcratio[KEY_ME]->GetYaxis()->SetTitle("Data/MC");

        h_mMass_y_mcratio[KEY_ME]->Draw("pE");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_IndividualCos2PhiStarPhiBins_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("kaonratio%s%s_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_phiy_ratio[KEY_ME]->SetTitle(Form("Data/RC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phiy_ratio[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_phiy_ratio[KEY_ME]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

        h_mMass_phiy_ratio[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());

    outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_IndividualCos2PhiStarPhiBins_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str());
    outputstart = Form("%s[",outputname.c_str());
    outputstop  = Form("%s]",outputname.c_str());
    craw->Print(outputstart.c_str());

    for(int i = 0; i < 4; i++)
    {
      for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
      {
        craw->cd(i_phi+1);
        string KEY_ME = Form("mckaonratio%s%s_phipt_%d_cos2phistarphi_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i,i_phi);
        h_mMass_phiy_mcratio[KEY_ME]->SetTitle(Form("Data/MC %1.1f<p_{T,#phi}<%1.1f, %1.1f<cos(2#phi*-2phi)<%1.1f, 20-60 Cent",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2],(float(i_phi)-5.)/5.,(float(i_phi)-4.)/5.));
        h_mMass_phiy_mcratio[KEY_ME]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_phiy_mcratio[KEY_ME]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

        h_mMass_phiy_mcratio[KEY_ME]->Draw("colz");
      }
      craw->Update();
      craw->Print(outputname.c_str());
    }
    craw->Print(outputstop.c_str());
  }


  /////////////////////////////// No Cos2PhiStarPhi separation //////////////////////////////
  for(int iset = 0; iset < 3; iset++)
  {
    for(int i = 0; i < 4; i++)
    {
      c1->cd(i+1);
      string KEY_phiy_ratio = Form("kaon%sratio%s_phiy_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
      h_mMass_phiy_ratio[KEY_phiy_ratio]->SetTitle(Form("Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mMass_phiy_ratio[KEY_phiy_ratio]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
      h_mMass_phiy_ratio[KEY_phiy_ratio]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

      h_mMass_phiy_ratio[KEY_phiy_ratio]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

    if(iset == 0)
    {
      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_phieta_ratio = Form("kaon%sratio%s_phieta_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mMass_phieta_ratio[KEY_phieta_ratio]->SetTitle(Form("Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mMass_phieta_ratio[KEY_phieta_ratio]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_phieta_ratio[KEY_phieta_ratio]->GetYaxis()->SetTitle("#eta");

        h_mMass_phieta_ratio[KEY_phieta_ratio]->Draw("colz");
      }
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phieta_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_phieta = Form("kaon%s%s_phieta_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mMass_phieta[KEY_phieta]->SetTitle(Form("Data %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mMass_phieta[KEY_phieta]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_phieta[KEY_phieta]->GetYaxis()->SetTitle("#eta");

        h_mMass_phieta[KEY_phieta]->Draw("colz");
      }
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phieta_data.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_RC_phieta = Form("rckaon%s%s_phieta_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mRcMass_phieta[KEY_RC_phieta]->SetTitle(Form("RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mRcMass_phieta[KEY_RC_phieta]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mRcMass_phieta[KEY_RC_phieta]->GetYaxis()->SetTitle("#eta");

        h_mRcMass_phieta[KEY_RC_phieta]->Draw("colz");
      }
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phieta_rc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_eta_ratio = Form("kaon%sratio%s_eta_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mMass_eta_ratio[KEY_eta_ratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mMass_eta_ratio[KEY_eta_ratio]->GetXaxis()->SetTitle("#eta");
        h_mMass_eta_ratio[KEY_eta_ratio]->GetYaxis()->SetTitle("Data/RC");
        h_mMass_eta_ratio[KEY_eta_ratio]->SetMarkerStyle(20);

        h_mMass_eta_ratio[KEY_eta_ratio]->Draw("pE");
      }
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_eta_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

      TLegend *legphi = new TLegend(0.825,0.55,0.975,0.75);
      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_eta = Form("kaon%s%s_eta_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        string KEY_RC_eta = Form("rckaon%s%s_eta_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mMass_eta[KEY_eta]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mMass_eta[KEY_eta]->GetXaxis()->SetTitle("#eta");
        h_mMass_eta[KEY_eta]->GetYaxis()->SetTitle("Counts");
        h_mMass_eta[KEY_eta]->SetMarkerStyle(20);
        h_mMass_eta[KEY_eta]->SetMarkerColor(kOrange+7);
        h_mMass_eta[KEY_eta]->SetLineColor(kOrange+7);

        int min = h_mMass_eta[KEY_eta]->GetMinimum();
        int max = h_mMass_eta[KEY_eta]->GetMaximum();
        if(h_mRcMass_eta[KEY_RC_eta]->GetMinimum() < min) min = h_mRcMass_eta[KEY_RC_eta]->GetMinimum();
        if(h_mRcMass_eta[KEY_RC_eta]->GetMaximum() > max) max = h_mRcMass_eta[KEY_RC_eta]->GetMaximum();
        h_mMass_eta[KEY_eta]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

        h_mMass_eta[KEY_eta]->Draw("pE");

        h_mRcMass_eta[KEY_RC_eta]->SetMarkerStyle(24);
        h_mRcMass_eta[KEY_RC_eta]->SetMarkerColor(kBlack);
        h_mRcMass_eta[KEY_RC_eta]->SetLineColor(kBlack);
        h_mRcMass_eta[KEY_RC_eta]->Draw("pE same");

        if(i == 0)
        {
          legphi->AddEntry(h_mMass_eta[KEY_eta],"Data","p");
          legphi->AddEntry(h_mRcMass_eta[KEY_RC_eta],"RC","p");
        }
        legphi->Draw("same");
      }
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_eta_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));
    }

   
    for(int i = 0; i < 4; i++)
    {
      c1->cd(i+1);
      string KEY_phiy = Form("kaon%s%s_phiy_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
      h_mMass_phiy[KEY_phiy]->SetTitle(Form("Data %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mMass_phiy[KEY_phiy]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
      h_mMass_phiy[KEY_phiy]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

      h_mMass_phiy[KEY_phiy]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_data.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

    for(int i = 0; i < 4; i++)
    {
      c1->cd(i+1);
      string KEY_RC_phiy = Form("rckaon%s%s_phiy_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
      h_mRcMass_phiy[KEY_RC_phiy]->SetTitle(Form("RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mRcMass_phiy[KEY_RC_phiy]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
      h_mRcMass_phiy[KEY_RC_phiy]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

      h_mRcMass_phiy[KEY_RC_phiy]->Draw("colz");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_rc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

    for(int i = 0; i < 4; i++)
    {
      c1->cd(i+1);
      string KEY_phi_ratio = Form("kaon%sratio%s_phi_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
      h_mMass_phi_ratio[KEY_phi_ratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mMass_phi_ratio[KEY_phi_ratio]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
      h_mMass_phi_ratio[KEY_phi_ratio]->GetYaxis()->SetTitle("Data/RC");
      h_mMass_phi_ratio[KEY_phi_ratio]->SetMarkerStyle(20);

      h_mMass_phi_ratio[KEY_phi_ratio]->Draw("pE");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phi_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

    TLegend *legphi = new TLegend(0.825,0.55,0.975,0.75);
    for(int i = 0; i < 4; i++)
    {
      c1->cd(i+1);
      string KEY_phi = Form("kaon%s%s_phi_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
      string KEY_RC_phi = Form("rckaon%s%s_phi_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
      h_mMass_phi[KEY_phi]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mMass_phi[KEY_phi]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
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
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phi_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

    for(int i = 0; i < 4; i++)
    {
      c1->cd(i+1);
      string KEY_y_ratio = Form("kaon%sratio%s_y_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
      h_mMass_y_ratio[KEY_y_ratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mMass_y_ratio[KEY_y_ratio]->GetXaxis()->SetTitle(kaonetalabel[iset].c_str());
      h_mMass_y_ratio[KEY_y_ratio]->GetYaxis()->SetTitle("Data/RC");
      h_mMass_y_ratio[KEY_y_ratio]->SetMarkerStyle(20);

      h_mMass_y_ratio[KEY_y_ratio]->Draw("pE");
    }
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_y_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

    TLegend *legy = new TLegend(0.825,0.55,0.975,0.75);
    for(int i = 0; i < 4; i++)
    {
      c1->cd(i+1);
      string KEY_y = Form("kaon%s%s_y_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
      string KEY_RC_y = Form("rckaon%s%s_y_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
      h_mMass_y[KEY_y]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mMass_y[KEY_y]->GetXaxis()->SetTitle(kaonetalabel[iset].c_str());
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
    c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_y_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));
   
    //if(iset == 1)
    //{
      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_phiy_mcratio = Form("mckaon%sratio%s_phiy_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mMass_phiy_mcratio[KEY_phiy_mcratio]->SetTitle(Form("Data/MC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mMass_phiy_mcratio[KEY_phiy_mcratio]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_phiy_mcratio[KEY_phiy_mcratio]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

        h_mMass_phiy_mcratio[KEY_phiy_mcratio]->Draw("colz");
      }
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_phi_mcratio = Form("mckaon%sratio%s_phi_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mMass_phi_mcratio[KEY_phi_mcratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mMass_phi_mcratio[KEY_phi_mcratio]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMass_phi_mcratio[KEY_phi_mcratio]->GetYaxis()->SetTitle("Data/MC");
        h_mMass_phi_mcratio[KEY_phi_mcratio]->SetMarkerStyle(20);

        h_mMass_phi_mcratio[KEY_phi_mcratio]->Draw("pE");
      }
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phi_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

      TLegend *legphimc = new TLegend(0.825,0.55,0.975,0.75);
      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_phi = Form("kaon%s%s_phi_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        string KEY_MC_phi = Form("mckaon%s%s_phi_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mMass_phi[KEY_phi]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mMass_phi[KEY_phi]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
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
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phi_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_y_mcratio = Form("mckaon%sratio%s_y_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mMass_y_mcratio[KEY_y_mcratio]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mMass_y_mcratio[KEY_y_mcratio]->GetXaxis()->SetTitle(kaonetalabel[iset].c_str());
        h_mMass_y_mcratio[KEY_y_mcratio]->GetYaxis()->SetTitle("Data/MC");
        h_mMass_y_mcratio[KEY_y_mcratio]->SetMarkerStyle(20);

        h_mMass_y_mcratio[KEY_y_mcratio]->Draw("pE");
      }
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_y_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

      TLegend *legymc = new TLegend(0.825,0.55,0.975,0.75);
      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_y = Form("kaon%s%s_y_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        string KEY_MC_y = Form("mckaon%s%s_y_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mMass_y[KEY_y]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mMass_y[KEY_y]->GetXaxis()->SetTitle(kaonetalabel[iset].c_str());
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
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_y_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));

      for(int i = 0; i < 4; i++)
      {
        c1->cd(i+1);
        string KEY_MC_phiy = Form("mckaon%s%s_phiy_phipt_%d",kaonlabel[iset].c_str(),kaonextralabel[iset].c_str(),i);
        h_mMcMass_phiy[KEY_MC_phiy]->SetTitle(Form("MC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
        h_mMcMass_phiy[KEY_MC_phiy]->GetXaxis()->SetTitle(kaonphilabel[iset].c_str());
        h_mMcMass_phiy[KEY_MC_phiy]->GetYaxis()->SetTitle(kaonetalabel[iset].c_str());

        h_mMcMass_phiy[KEY_MC_phiy]->Draw("colz");
      }
      c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_%s%s_phiy_mc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),kaonlabel[iset].c_str(),kaonextralabel[iset].c_str()));
    }
  //}

  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_pT_ratio = Form("kaonratio_pt_phipt_%d",i);
    g_mKaon_pT_ratio[KEY_pT_ratio]->SetTitle(Form("Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    g_mKaon_pT_ratio[KEY_pT_ratio]->GetXaxis()->SetTitle("p_{T}");
    g_mKaon_pT_ratio[KEY_pT_ratio]->GetYaxis()->SetTitle("Data/RC");
    g_mKaon_pT_ratio[KEY_pT_ratio]->SetMarkerStyle(20);
    //g_mKaon_pT_ratio[KEY_pT_ratio]->SetMarkerColor(kOrange+7);
    //g_mKaon_pT_ratio[KEY_pT_ratio]->SetLineColor(kOrange+7);

    //int min = g_mKaon_pT[KEY_pT]->GetMinimum();
    //int max = g_mKaon_pT[KEY_pT]->GetMaximum();
    //if(g_mRcKaon_pT[KEY_RC_pT]->GetMinimum() < min) min = g_mRcKaon_pT[KEY_RC_pT]->GetMinimum();
    //if(g_mRcKaon_pT[KEY_RC_pT]->GetMaximum() > max) max = g_mRcKaon_pT[KEY_RC_pT]->GetMaximum();
    //g_mKaon_pT[KEY_pT]->SetMinimum(min*0.9);
    //g_mKaon_pT[KEY_pT]->SetMaximum(max*1.1);

    g_mKaon_pT_ratio[KEY_pT_ratio]->Draw("APE");

    //g_mKaon_pT[KEY_RC_pT]->SetMarkerStyle(24);
    //g_mKaon_pT[KEY_RC_pT]->SetMarkerColor(kBlack);
    //g_mKaon_pT[KEY_RC_pT]->SetLineColor(kBlack);
    //g_mKaon_pT[KEY_RC_pT]->Draw("pE same");

    //if(i == 0)
    //{
    //  legpt->AddEntry(h_mMass_y[KEY_y],"Data","p");
    //  legpt->AddEntry(h_mRcMass_y[KEY_RC_y],"RC","p");
    //}
    //legpt->Draw("same");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pT_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

  TLegend *legpt = new TLegend(0.6,0.7,0.8,0.9);
  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_pT = Form("kaon_pt_phipt_%d",i);
    string KEY_RC_pT = Form("rckaon_pt_phipt_%d",i);
    g_mKaon_pT[KEY_pT]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    g_mKaon_pT[KEY_pT]->GetXaxis()->SetTitle("p_{T}");
    g_mKaon_pT[KEY_pT]->GetYaxis()->SetTitle("Counts");
    g_mKaon_pT[KEY_pT]->SetMarkerStyle(20);
    g_mKaon_pT[KEY_pT]->SetMarkerColor(kOrange+7);
    g_mKaon_pT[KEY_pT]->SetLineColor(kOrange+7);

    //int min = g_mKaon_pT[KEY_pT]->GetMinimum();
    //int max = g_mKaon_pT[KEY_pT]->GetMaximum();
    //if(g_mRcKaon_pT[KEY_RC_pT]->GetMinimum() < min) min = g_mRcKaon_pT[KEY_RC_pT]->GetMinimum();
    //if(g_mRcKaon_pT[KEY_RC_pT]->GetMaximum() > max) max = g_mRcKaon_pT[KEY_RC_pT]->GetMaximum();
    //g_mKaon_pT[KEY_pT]->SetMinimum(min*0.9);
    //g_mKaon_pT[KEY_pT]->SetMaximum(max*1.1);

    g_mKaon_pT[KEY_pT]->Draw("APE");

    g_mRcKaon_pT[KEY_RC_pT]->SetMarkerStyle(24);
    g_mRcKaon_pT[KEY_RC_pT]->SetMarkerColor(kBlack);
    g_mRcKaon_pT[KEY_RC_pT]->SetLineColor(kBlack);
    g_mRcKaon_pT[KEY_RC_pT]->Draw("pE same");

    if(i == 0)
    {
      legpt->AddEntry(g_mKaon_pT[KEY_pT],"Data","p");
      legpt->AddEntry(g_mRcKaon_pT[KEY_RC_pT],"RC","p");
    }
    legpt->Draw("same");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pT_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_pT_mcratio = Form("kaonratio_pt_phipt_%d",i);
    g_mKaon_pT_mcratio[KEY_pT_mcratio]->SetTitle(Form("Data/MC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    g_mKaon_pT_mcratio[KEY_pT_mcratio]->GetXaxis()->SetTitle("p_{T}");
    g_mKaon_pT_mcratio[KEY_pT_mcratio]->GetYaxis()->SetTitle("Data/MC");
    g_mKaon_pT_mcratio[KEY_pT_mcratio]->SetMarkerStyle(20);
    //g_mKaon_pT_mcratio[KEY_pT_mcratio]->SetMarkerColor(kOrange+7);
    //g_mKaon_pT_mcratio[KEY_pT_mcratio]->SetLineColor(kOrange+7);

    //int min = g_mKaon_pT[KEY_pT]->GetMinimum();
    //int max = g_mKaon_pT[KEY_pT]->GetMaximum();
    //if(g_mMcKaon_pT[KEY_MC_pT]->GetMinimum() < min) min = g_mMcKaon_pT[KEY_MC_pT]->GetMinimum();
    //if(g_mMcKaon_pT[KEY_MC_pT]->GetMaximum() > max) max = g_mMcKaon_pT[KEY_MC_pT]->GetMaximum();
    //g_mKaon_pT[KEY_pT]->SetMinimum(min*0.9);
    //g_mKaon_pT[KEY_pT]->SetMaximum(max*1.1);

    g_mKaon_pT_mcratio[KEY_pT_mcratio]->Draw("APE");

    //g_mKaon_pT[KEY_MC_pT]->SetMarkerStyle(24);
    //g_mKaon_pT[KEY_MC_pT]->SetMarkerColor(kBlack);
    //g_mKaon_pT[KEY_MC_pT]->SetLineColor(kBlack);
    //g_mKaon_pT[KEY_MC_pT]->Draw("pE same");

    //if(i == 0)
    //{
    //  legpt->AddEntry(h_mMass_y[KEY_y],"Data","p");
    //  legpt->AddEntry(h_mMcMass_y[KEY_MC_y],"MC","p");
    //}
    //legpt->Draw("same");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pT_datamcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

  TLegend *legptmc = new TLegend(0.6,0.7,0.8,0.9);
  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_pT = Form("kaon_pt_phipt_%d",i);
    string KEY_MC_pT = Form("mckaon_pt_phipt_%d",i);
    g_mKaon_pT[KEY_pT]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    g_mKaon_pT[KEY_pT]->GetXaxis()->SetTitle("p_{T}");
    g_mKaon_pT[KEY_pT]->GetYaxis()->SetTitle("Counts");
    g_mKaon_pT[KEY_pT]->SetMarkerStyle(20);
    g_mKaon_pT[KEY_pT]->SetMarkerColor(kOrange+7);
    g_mKaon_pT[KEY_pT]->SetLineColor(kOrange+7);

    //int min = g_mKaon_pT[KEY_pT]->GetMinimum();
    //int max = g_mKaon_pT[KEY_pT]->GetMaximum();
    //if(g_mMcKaon_pT[KEY_MC_pT]->GetMinimum() < min) min = g_mMcKaon_pT[KEY_MC_pT]->GetMinimum();
    //if(g_mMcKaon_pT[KEY_MC_pT]->GetMaximum() > max) max = g_mMcKaon_pT[KEY_MC_pT]->GetMaximum();
    //g_mKaon_pT[KEY_pT]->SetMinimum(min*0.9);
    //g_mKaon_pT[KEY_pT]->SetMaximum(max*1.1);

    g_mKaon_pT[KEY_pT]->Draw("APE");

    g_mMcKaon_pT[KEY_MC_pT]->SetMarkerStyle(24);
    g_mMcKaon_pT[KEY_MC_pT]->SetMarkerColor(kBlack);
    g_mMcKaon_pT[KEY_MC_pT]->SetLineColor(kBlack);
    g_mMcKaon_pT[KEY_MC_pT]->Draw("pE same");

    if(i == 0)
    {
      legptmc->AddEntry(g_mKaon_pT[KEY_pT],"Data","p");
      legptmc->AddEntry(g_mMcKaon_pT[KEY_MC_pT],"MC","p");
    }
    legptmc->Draw("same");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pT_datamc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

  // write background subtracted histograms to output file
  //string outputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_KaonDistributions_%s.root",vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  //TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  //File_OutPut->cd();
  //for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++) // pt bin 
  //{
  //  string KEY_phiy_ratio = Form("kaonratio_phiy_phipt_%d",i_pt_phi);
  //  h_mMass_phiy_ratio[KEY_phiy_ratio]->Write();

  //  string KEY_phiy = Form("kaon_phiy_phipt_%d",i_pt_phi);
  //  h_mMass_phiy[KEY_phiy]->Write();

  //  string KEY_RC_phiy = Form("rckaon_phiy_phipt_%d",i_pt_phi);
  //  h_mRcMass_phiy[KEY_RC_phiy]->Write();

  //  string KEY_phi_ratio = Form("kaonratio_phi_phipt_%d",i_pt_phi);
  //  h_mMass_phi_ratio[KEY_phi_ratio]->Write();

  //  string KEY_phi = Form("kaon_phi_phipt_%d",i_pt_phi);
  //  h_mMass_phi[KEY_phi]->Write();

  //  string KEY_y_ratio   = Form("kaonratio_y_phipt_%d",i_pt_phi);
  //  h_mMass_y_ratio[KEY_y_ratio]->Write();    

  //  string KEY_y   = Form("kaon_y_phipt_%d",i_pt_phi);
  //  h_mMass_y[KEY_y]->Write();    

  //  string KEY_RC_phi = Form("rckaon_phi_phipt_%d",i_pt_phi);
  //  h_mRcMass_phi[KEY_RC_phi]->Write(); 

  //  string KEY_RC_y   = Form("rckaon_y_phipt_%d",i_pt_phi);
  //  h_mRcMass_y[KEY_RC_y]->Write();     

  //  string KEY_pT = Form("kaon_pt_phipt_%d",i_pt_phi);
  //  g_mKaon_pT[KEY_pT]->SetName(KEY_pT.c_str());
  //  g_mKaon_pT[KEY_pT]->Write();

  //  string KEY_RC_pT = Form("rckaon_pt_phipt_%d",i_pt_phi);
  //  g_mRcKaon_pT[KEY_RC_pT]->SetName(KEY_RC_pT.c_str());  
  //  g_mRcKaon_pT[KEY_RC_pT]->Write();  
  //}
  //File_OutPut->Close();
}

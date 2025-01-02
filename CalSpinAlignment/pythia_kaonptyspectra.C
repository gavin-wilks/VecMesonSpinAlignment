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

void pythia_kaonptyspectra(int energy = 4, int pid = 0, int year = 0, string date = "20240828", bool random3D = false, int order = 2, string etamode = "eta1_eta1", int deltaonly = 0, string sim = "Embed")
{
  //original was 20240715
  std::string EP[2] = {"","2nd"};
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  //string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_deltar_fixedagain.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_cent2060_kaonspectra_basicinfo_fpoly.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  //string InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(order == 1) InPutFile_SE = Form("../data/Yields_Phi_SE_%s_%s_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(random3D) InPutFile_SE = Form("../data/3DRandom/Yields_Phi_SE_%s_%s_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_deltar_fixedagain.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_cent2060_kaonspectra_basicinfo_fpoly.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  //string InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(order == 1) InPutFile_ME = Form("../data/Yields_Phi_ME_%s_%s_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str(),etamode.c_str());
  if(random3D) InPutFile_ME = Form("../data/3DRandom/Yields_Phi_ME_%s_%s_3DRandom.root",vmsa::mBeamEnergy[energy].c_str(),date.c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());

  //string folder = "PhiEmbeddingWithWeights_PIDEff";

  //string folder = "PhiEmbeddingWithWeights_PIDEff_20240725";
  //string folder = "PhiEmbeddingWithWeights_PIDEff_20240726";
  //string folder = "Pythia_NoWeights_20240731";
  //string folder = "KaonOnlyReweight";
  //string folder = "Pythia_20240808_EtaCutFirst";
  //string folder = "PhiEmbedding_NoHelicity_20240822EP_phiweight1";
  //string folder = "PhiEmbedding_20240822_NOMCphicuts";
  //string folder = "PhiEmbedding_20240822_MCphicuts";
  //string folder = "PhiEmbedding_20240822EP_phiswitch0";
  //string folder = "PhiEmbedding_WithWeights_20240809";
  //string folder = "PhiEmbedding_20240822_nsigcut";
  string folder = "PhiEmbedding_weights_nomccuts_rapidity";
  //if(yspectra == 68)
  //{
  //  spectra = "PhiEmbedding_20240821_phiweight0";
  //  filename = "PhiEmbedding_20240821_phiswitch0.root"; 
  //}
  //if(yspectra == 69)
  //{
  //  spectra = "PhiEmbedding_20240821_phiweight1";
  //  filename = "PhiEmbedding_20240821_phiswitch1.root"; 
  //}
  //if(yspectra == 70)
  //{
  //  spectra = "PhiEmbedding_NoHelicity_20240822EP_phiweight1";
  //  filename = "Embedding_19GeV_20240820.root"; 
  //}
  //if(yspectra == 71)
  //{
  //  spectra = "PhiEmbedding_Helicity2D_20240822EP_phiweight1";
  //  filename = "PhiEmbedding_20240822EP_phiswitch1.root"; 
  //}
  //string folder = "Pythia_PhiMesonWeights_20240731_Retry";
  //string folder = "Pythia_PhiMesonWeights_KaonWeights_20240731";
  //string folder = "Pythia_PhiMesonWeights_20240801_Helicity2D";
  //string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbed_Weights_PIDEff_1.root",folder.c_str());  
  //string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbed_WithWeights_PIDEff_20240725_FixedMC_0.root",folder.c_str());  
  //string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",folder.c_str());  
  //string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbedding_20240822EP_phiswitch1.root",folder.c_str());  
  //string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbedding_20240822_NOMCphicuts.root",folder.c_str());  
  //string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbedding_20240822_MCphicuts.root",folder.c_str());  
  //string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbedding_20240822_FinerToFEta_phiswitch0.root",folder.c_str());  
  //string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbedding_WithWeights_20240809.root",folder.c_str());  
  //string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbedding_20240822_nsigcut.root",folder.c_str());  
  string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbedding_weights_nomccuts_rapidity.root",folder.c_str());  
  //string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbedding_20240822_FinerToFEta_phiswitch1.root",folder.c_str());  
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
        string KEY_SE = Form("kaonplus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_SE",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();

        string KEY_ME = Form("kaonplus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_ME",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();

        string KEY_SM = Form("kaonplus_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
        h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
        h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     

        KEY_SE = Form("kaonminus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_SE",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();

        KEY_ME = Form("kaonminus_phipt_%d_cos2phistarphi_%d_kaonpt_%d_ME",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();

        KEY_SM = Form("kaonminus_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
        h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
        h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     

        //KEY_SE = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta_SE",i_pt_phi,i_phi,i_pt_kaon);
        //h_mMass_SE[KEY_SE] = (TH2F*)((TH2F*)File_SE->Get(KEY_SE.c_str()))->Clone();

        //KEY_ME = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta_ME",i_pt_phi,i_phi,i_pt_kaon);
        //h_mMass_ME[KEY_ME] = (TH2F*)((TH2F*)File_ME->Get(KEY_ME.c_str()))->Clone();

        //KEY_SM = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d_eta",i_pt_phi,i_phi,i_pt_kaon);
        //h_mMass_SM[KEY_SM] = (TH2F*)h_mMass_SE[KEY_SE]->Clone(KEY_SM.c_str());
        //h_mMass_SM[KEY_SM]->SetTitle(KEY_SM.c_str());
        //h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);     
      }
    }
  }
  string KEY_RC = Form("rc9_kminus_cent9");
  h_mMass_RC[KEY_RC] = (TH2F*)((TH2F*)File_RC->Get(KEY_RC.c_str()))->Clone();

  //string KEY_MC = Form("mc_kminus_cent9");
  //h_mMass_MC[KEY_MC] = (TH2F*)((TH2F*)File_RC->Get(KEY_MC.c_str()))->Clone();

  KEY_RC = Form("rc9_kplus_cent9");
  h_mMass_RC[KEY_RC] = (TH2F*)((TH2F*)File_RC->Get(KEY_RC.c_str()))->Clone();

  //KEY_MC = Form("mc_kplus_cent9");
  //h_mMass_MC[KEY_MC] = (TH2F*)((TH2F*)File_RC->Get(KEY_MC.c_str()))->Clone();

  TH2FMap h_mMass_phiy;   
  TH2FMap h_mRcMass_phiy;   
  TH1FMap h_mMass_y;   
  TH1FMap h_mRcMass_y;

  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++) // pt bin 
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
    {
      for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
      {
        string KEY_SM = Form("kaonplus_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        //string KEY_phiy = Form("kaonplus_phiy_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        string KEY_phiy = Form("kaonplus_phiy_kaonpt_%d",i_pt_kaon);
        if(i_phi == 0 && i_pt_phi == 0) h_mMass_phiy[KEY_phiy] = (TH2F*)h_mMass_SM[KEY_SM]->Clone(KEY_phiy.c_str()); // First bins
        else h_mMass_phiy[KEY_phiy]->Add(h_mMass_SM[KEY_SM],1.0); // All other bins

        //string KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        //string KEY_RC_phiy = Form("rckaon_phiy_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        //if(i_phi == 0) h_mRcMass_phiy[KEY_RC_phiy] = (TH2F*)h_mMass_RC[KEY_RC]->Clone(KEY_RC_phiy.c_str()); // First bins
        //else h_mRcMass_phiy[KEY_RC_phiy]->Add(h_mMass_RC[KEY_RC],1.0); // All other bins

        //string KEY_phi = Form("kaon_phi_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        //string KEY_y   = Form("kaonplus_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        string KEY_y   = Form("kaonplus_y_kaonpt_%d",i_pt_kaon);
        if(i_phi == 9 && i_pt_phi == 3) // Last bins
        {
          //h_mMass_phi[KEY_phi] = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionX(KEY_phi.c_str(),0,-1,"e"));
          h_mMass_y[KEY_y]     = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionY(KEY_y.c_str(),0,-1,"e"));
        } 

        //string KEY_RC_phi = Form("rckaon_phi_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        //string KEY_RC_y   = Form("rckaon_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        //if(i_phi == 9) // Last bins
        //{
        //  //h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
        //  h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));
        //} 
      }
    }
  }
  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++) // pt bin 
  {
    for(Int_t i_phi = 0; i_phi < 10; i_phi++) // phi-psi bin
    {
      for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
      {
        string KEY_SM = Form("kaonminus_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        //string KEY_phiy = Form("kaonminus_phiy_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        string KEY_phiy = Form("kaonminus_phiy_kaonpt_%d",i_pt_kaon);
        if(i_phi == 0 && i_pt_phi == 0) h_mMass_phiy[KEY_phiy] = (TH2F*)h_mMass_SM[KEY_SM]->Clone(KEY_phiy.c_str()); // First bins
        else h_mMass_phiy[KEY_phiy]->Add(h_mMass_SM[KEY_SM],1.0); // All other bins

        //string KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        //string KEY_RC_phiy = Form("rckaon_phiy_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        //if(i_phi == 0) h_mRcMass_phiy[KEY_RC_phiy] = (TH2F*)h_mMass_RC[KEY_RC]->Clone(KEY_RC_phiy.c_str()); // First bins
        //else h_mRcMass_phiy[KEY_RC_phiy]->Add(h_mMass_RC[KEY_RC],1.0); // All other bins

        //string KEY_phi = Form("kaon_phi_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        //string KEY_y   = Form("kaonminus_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        string KEY_y   = Form("kaonminus_y_kaonpt_%d",i_pt_kaon);
        if(i_phi == 9 && i_pt_phi == 3) // Last bins
        {
          //h_mMass_phi[KEY_phi] = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionX(KEY_phi.c_str(),0,-1,"e"));
          h_mMass_y[KEY_y]     = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionY(KEY_y.c_str(),0,-1,"e"));
        } 

        //string KEY_RC_phi = Form("rckaon_phi_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        //string KEY_RC_y   = Form("rckaon_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        //if(i_phi == 9) // Last bins
        //{
        //  //h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
        //  h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));
        //} 
      }
    }
  }

  const int rebinconstant = 3;
  const int pt_total_rebinned = data_constants::kaon_pt_bins/rebinconstant;
  cout << "pt_total_rebinned = " << pt_total_rebinned << endl;
  const int pt_total_rebinned_p1 = data_constants::kaon_pt_bins/rebinconstant+1;
  cout << "pt_total_rebinned_p1 = " << pt_total_rebinned_p1 << endl;
  double kaon_pt_low_rebin[pt_total_rebinned_p1] = {0.0};

  for(int ipt = 0; ipt < data_constants::kaon_pt_bins; ipt++)
  {
    cout << "ipt = " << ipt << ", range = " << data_constants::kaon_pt_low[ipt] << "," << data_constants::kaon_pt_high[ipt] << endl;
  }

  for(int ipt = 0; ipt < pt_total_rebinned; ipt++)
  {
    kaon_pt_low_rebin[ipt] = data_constants::kaon_pt_low[ipt*rebinconstant];
  }
  kaon_pt_low_rebin[pt_total_rebinned] = data_constants::kaon_pt_high[data_constants::kaon_pt_bins-1];

  for(int ipt = 0; ipt < pt_total_rebinned_p1; ipt++)
  {
    cout << "ipt = " << ipt << ", pt = " << kaon_pt_low_rebin[ipt] << endl;;
  }

  TH2FMap h_mMass_pty;
  TH2FMap h_mRcMass_pty;
  TH2FMap h_mRatio_pty;
  TH1FMap h_mMass_pt;
  TH1FMap h_mRcMass_pt;
  TH1FMap h_mRatio_pt;
  TH1FMap h_mRatio_y;

  cout << "BEFORE CREATING FILE" << endl;  
  string outputfile = "../output/kaon_ptyratio_19GeV_NoPhiWeights.root";
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  cout << "AFTER CREATING FILE" << endl;  

  for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
  {
    int rebin_pt_kaon = i_pt_kaon/rebinconstant;
    //string KEY_y_rebin = Form("kaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,rebin_pt_kaon);
    string KEY_y_rebin = Form("kaonplusrebin_y_kaonpt_%d",rebin_pt_kaon);
    //string KEY_RC_y_rebin = Form("rckaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,rebin_pt_kaon);

    //string KEY_y = Form("kaon_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
    string KEY_y = Form("kaonplus_y_kaonpt_%d",i_pt_kaon);
    //string KEY_RC_y = Form("rckaon_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
    if(i_pt_kaon%rebinconstant == 0)
    {
      h_mMass_y[KEY_y_rebin]      = (TH1F*) h_mMass_y[KEY_y]->Clone(KEY_y_rebin.c_str());
      //h_mRcMass_y[KEY_RC_y_rebin] = (TH1F*) h_mRcMass_y[KEY_RC_y]->Clone(KEY_RC_y_rebin.c_str());
    }
    else
    {
      h_mMass_y[KEY_y_rebin]->Add(h_mMass_y[KEY_y],1.0);     
      //h_mRcMass_y[KEY_RC_y_rebin]->Add(h_mRcMass_y[KEY_RC_y],1.0);
    }

    KEY_y_rebin = Form("kaonminusrebin_y_kaonpt_%d",rebin_pt_kaon);
    //string KEY_RC_y_rebin = Form("rckaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,rebin_pt_kaon);

    //string KEY_y = Form("kaon_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
    KEY_y = Form("kaonminus_y_kaonpt_%d",i_pt_kaon);
    //string KEY_RC_y = Form("rckaon_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
    if(i_pt_kaon%rebinconstant == 0)
    {
      h_mMass_y[KEY_y_rebin]      = (TH1F*) h_mMass_y[KEY_y]->Clone(KEY_y_rebin.c_str());
      //h_mRcMass_y[KEY_RC_y_rebin] = (TH1F*) h_mRcMass_y[KEY_RC_y]->Clone(KEY_RC_y_rebin.c_str());
    }
    else
    {
      h_mMass_y[KEY_y_rebin]->Add(h_mMass_y[KEY_y],1.0);     
      //h_mRcMass_y[KEY_RC_y_rebin]->Add(h_mRcMass_y[KEY_RC_y],1.0);
    }
  }
  
  //string KEY_pty    = Form("pty_phipt_%d",i_pt_phi);   
  //string KEY_RC_pty = Form("rc_pty_phipt_%d",i_pt_phi);   
  //string KEY_ratio_pty = Form("ratio_pty_phipt_%d",i_pt_phi);   
  string KEY_pty    = Form("kplus_pty");   
  string KEY_RC_pty = Form("kplus_rc_pty");   
  string KEY_ratio_pty = Form("kplus_ratio_pty");   
  h_mMass_pty[KEY_pty] = new TH2F(KEY_pty.c_str(),KEY_pty.c_str(),pt_total_rebinned,kaon_pt_low_rebin,data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
  h_mRcMass_pty[KEY_RC_pty] = new TH2F(KEY_RC_pty.c_str(),KEY_RC_pty.c_str(),pt_total_rebinned,kaon_pt_low_rebin,data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);

  for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins/rebinconstant; i_pt_kaon++)
  {
    //string KEY_y_rebin = Form("kaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
    //string KEY_RC_y_rebin = Form("rckaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
    string KEY_y_rebin = Form("kaonplusrebin_y_kaonpt_%d",i_pt_kaon);
    //string KEY_RC_y_rebin = Form("rckaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
    for(int ibin = 1; ibin <= h_mMass_y[KEY_y_rebin]->GetNbinsX(); ibin++)
    {
      double val = h_mMass_y[KEY_y_rebin]->GetBinContent(ibin);
      double err = h_mMass_y[KEY_y_rebin]->GetBinError(ibin);
      //double valRC = h_mRcMass_y[KEY_RC_y_rebin]->GetBinContent(ibin);
      //double errRC = h_mRcMass_y[KEY_RC_y_rebin]->GetBinError(ibin);

      //cout << "bin = " << ibin << ", bincenter = " << h_mMass_y[KEY_y_rebin]->GetBinCenter(ibin) << endl;
      
      h_mMass_pty[KEY_pty]->SetBinContent(i_pt_kaon+1,ibin,val);
      h_mMass_pty[KEY_pty]->SetBinError(i_pt_kaon+1,ibin,err);
      //h_mRcMass_pty[KEY_RC_pty]->SetBinContent(i_pt_kaon+1,ibin,valRC);
      //h_mRcMass_pty[KEY_RC_pty]->SetBinError(i_pt_kaon+1,ibin,errRC);
    }
  }
  h_mMass_pty[KEY_pty]->RebinY(2);
  cout << "Data: " << endl;
  cout << "NbinsX = " << h_mMass_pty[KEY_pty]->GetNbinsX() << endl;
  cout << "NbinsY = " << h_mMass_pty[KEY_pty]->GetNbinsY() << endl;

  h_mMass_RC["rc9_kplus_cent9"]->RebinY(2);
  h_mMass_RC["rc9_kplus_cent9"]->RebinX(rebinconstant);
  cout << "RC: " << endl;
  cout << "NbinsX = " << h_mMass_RC["rc9_kplus_cent9"]->GetNbinsX() << endl;
  cout << "NbinsY = " << h_mMass_RC["rc9_kplus_cent9"]->GetNbinsY() << endl;

  double yields = h_mMass_pty[KEY_pty]->Integral(0,-1,0,-1);
  //double rcyields = h_mrcmass_pty[key_rc_pty]->integral(0,-1,0,-1);
  double rcyields = h_mMass_RC["rc9_kplus_cent9"]->Integral(0,-1,0,-1);
 
  h_mMass_RC["rc9_kplus_cent9"]->Scale(yields/rcyields);

  int n_bins = h_mMass_pty[KEY_pty]->GetNbinsX();
  for (int i = 1; i <= n_bins; ++i)
  {
    double bin_low_edge = h_mMass_pty[KEY_pty]->GetXaxis()->GetBinLowEdge(i);
    double bin_up_edge  = h_mMass_pty[KEY_pty]->GetXaxis()->GetBinLowEdge(i+1);
    double bin_low_edge_rc = h_mMass_RC["rc9_kplus_cent9"]->GetXaxis()->GetBinLowEdge(i);
    double bin_up_edge_rc  = h_mMass_RC["rc9_kplus_cent9"]->GetXaxis()->GetBinLowEdge(i+1);
    if(bin_low_edge != bin_low_edge_rc || bin_up_edge != bin_up_edge_rc)
    {
      cout << "Bin Data" << i << ": [" << bin_low_edge << ", " << bin_up_edge << ")" << endl;
      cout << "Bin RC  " << i << ": [" << bin_low_edge_rc << ", " << bin_up_edge_rc << ")" << endl;
    }
  }
  int n_binsy = h_mMass_pty[KEY_pty]->GetNbinsY();
  for (int i = 1; i <= n_binsy; ++i)
  {
    double bin_low_edge = h_mMass_pty[KEY_pty]->GetYaxis()->GetBinLowEdge(i);
    double bin_up_edge  = h_mMass_pty[KEY_pty]->GetYaxis()->GetBinLowEdge(i+1);
    double bin_low_edge_rc = h_mMass_RC["rc9_kplus_cent9"]->GetYaxis()->GetBinLowEdge(i);
    double bin_up_edge_rc  = h_mMass_RC["rc9_kplus_cent9"]->GetYaxis()->GetBinLowEdge(i+1);
    if(bin_low_edge != bin_low_edge_rc || bin_up_edge != bin_up_edge_rc)
    {
      cout << "Bin Data" << i << ": [" << bin_low_edge << ", " << bin_up_edge << ")" << endl;
      cout << "Bin RC  " << i << ": [" << bin_low_edge_rc << ", " << bin_up_edge_rc << ")" << endl;
    }
  }

  h_mRatio_pty[KEY_ratio_pty] = (TH2F*) h_mMass_pty[KEY_pty]->Clone(KEY_ratio_pty.c_str());
  h_mRatio_pty[KEY_ratio_pty]->Divide(h_mMass_RC["rc9_kplus_cent9"]);
  
  string KEY_pt    = Form("kplus_pt");   
  string KEY_RC_pt = Form("kplus_rc_pt");   
  string KEY_ratio_pt = Form("kplus_ratio_pt");   
  h_mMass_pt[KEY_pt] = (TH1F*) h_mMass_pty[KEY_pty]->ProjectionX(KEY_pt.c_str(),0,-1,"e");
  h_mRcMass_pt[KEY_RC_pt] = (TH1F*) h_mMass_RC["rc9_kplus_cent9"]->ProjectionX(KEY_RC_pt.c_str(),0,-1,"e");
  h_mRatio_pt[KEY_ratio_pt] = (TH1F*) h_mMass_pt[KEY_pt]->Clone(KEY_ratio_pt.c_str());
  h_mRatio_pt[KEY_ratio_pt]->Divide(h_mRcMass_pt[KEY_RC_pt]);

  string KEY_y    = Form("kplus_y");   
  string KEY_RC_y = Form("kplus_rc_y");   
  string KEY_ratio_y = Form("kplus_ratio_y");   
  h_mMass_y[KEY_y] = (TH1F*) h_mMass_pty[KEY_pty]->ProjectionY(KEY_y.c_str(),0,-1,"e");
  h_mRcMass_y[KEY_RC_y] = (TH1F*) h_mMass_RC["rc9_kplus_cent9"]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e");
  h_mRatio_y[KEY_ratio_y] = (TH1F*) h_mMass_y[KEY_y]->Clone(KEY_ratio_y.c_str());
  h_mRatio_y[KEY_ratio_y]->Divide(h_mRcMass_y[KEY_RC_y]);

  h_mRatio_pty[KEY_ratio_pty]->Write();


  KEY_pty    = Form("kminus_pty");   
  KEY_RC_pty = Form("kminus_rc_pty");   
  KEY_ratio_pty = Form("kminus_ratio_pty");   
  h_mMass_pty[KEY_pty] = new TH2F(KEY_pty.c_str(),KEY_pty.c_str(),pt_total_rebinned,kaon_pt_low_rebin,data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
  h_mRcMass_pty[KEY_RC_pty] = new TH2F(KEY_RC_pty.c_str(),KEY_RC_pty.c_str(),pt_total_rebinned,kaon_pt_low_rebin,data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);

  for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins/rebinconstant; i_pt_kaon++)
  {
    string KEY_y_rebin = Form("kaonminusrebin_y_kaonpt_%d",i_pt_kaon);
    //string KEY_RC_y_rebin = Form("rckaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
    for(int ibin = 1; ibin <= h_mMass_y[KEY_y_rebin]->GetNbinsX(); ibin++)
    {
      double val = h_mMass_y[KEY_y_rebin]->GetBinContent(ibin);
      double err = h_mMass_y[KEY_y_rebin]->GetBinError(ibin);
      //double valRC = h_mRcMass_y[KEY_RC_y_rebin]->GetBinContent(ibin);
      //double errRC = h_mRcMass_y[KEY_RC_y_rebin]->GetBinError(ibin);

      //cout << "bin = " << ibin << ", bincenter = " << h_mMass_y[KEY_y_rebin]->GetBinCenter(ibin) << endl;
      
      h_mMass_pty[KEY_pty]->SetBinContent(i_pt_kaon+1,ibin,val);
      h_mMass_pty[KEY_pty]->SetBinError(i_pt_kaon+1,ibin,err);
      //h_mRcMass_pty[KEY_RC_pty]->SetBinContent(i_pt_kaon+1,ibin,valRC);
      //h_mRcMass_pty[KEY_RC_pty]->SetBinError(i_pt_kaon+1,ibin,errRC);
    }
  }

  h_mMass_pty[KEY_pty]->RebinY(2);
  h_mMass_RC["rc9_kminus_cent9"]->RebinY(2);
  h_mMass_RC["rc9_kminus_cent9"]->RebinX(rebinconstant);

  yields = h_mMass_pty[KEY_pty]->Integral(0,-1,0,-1);
  //double rcyields = h_mrcmass_pty[key_rc_pty]->integral(0,-1,0,-1);
  rcyields = h_mMass_RC["rc9_kminus_cent9"]->Integral(0,-1,0,-1);
 
  h_mMass_RC["rc9_kminus_cent9"]->Scale(yields/rcyields);

  h_mRatio_pty[KEY_ratio_pty] = (TH2F*) h_mMass_pty[KEY_pty]->Clone(KEY_ratio_pty.c_str());
  h_mRatio_pty[KEY_ratio_pty]->Divide(h_mMass_RC["rc9_kplus_cent9"]);
  
  KEY_pt    = Form("kminus_pt");   
  KEY_RC_pt = Form("kminus_rc_pt");   
  KEY_ratio_pt = Form("kminus_ratio_pt");   
  h_mMass_pt[KEY_pt] = (TH1F*) h_mMass_pty[KEY_pty]->ProjectionX(KEY_pt.c_str(),0,-1,"e");
  h_mRcMass_pt[KEY_RC_pt] = (TH1F*) h_mMass_RC["rc9_kminus_cent9"]->ProjectionX(KEY_RC_pt.c_str(),0,-1,"e");
  h_mRatio_pt[KEY_ratio_pt] = (TH1F*) h_mMass_pt[KEY_pt]->Clone(KEY_ratio_pt.c_str());
  h_mRatio_pt[KEY_ratio_pt]->Divide(h_mRcMass_pt[KEY_RC_pt]);

  KEY_y    = Form("kminus_y");   
  KEY_RC_y = Form("kminus_rc_y");   
  KEY_ratio_y = Form("kminus_ratio_y");   
  h_mMass_y[KEY_y] = (TH1F*) h_mMass_pty[KEY_pty]->ProjectionY(KEY_y.c_str(),0,-1,"e");
  h_mRcMass_y[KEY_RC_y] = (TH1F*) h_mMass_RC["rc9_kminus_cent9"]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e");
  h_mRatio_y[KEY_ratio_y] = (TH1F*) h_mMass_y[KEY_y]->Clone(KEY_ratio_y.c_str());
  h_mRatio_y[KEY_ratio_y]->Divide(h_mRcMass_y[KEY_RC_y]);

  h_mRatio_pty[KEY_ratio_pty]->Write();

  for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins/rebinconstant; i_pt_kaon++)
  {
    string KEY_pty    = Form("kplus_pty");   
    string KEY_y_rebin_plus_d = Form("dkaonplusrebin_y_kaonpt_%d",i_pt_kaon);
    h_mMass_y[KEY_y_rebin_plus_d] = (TH1F*) h_mMass_pty[KEY_pty]->ProjectionY(KEY_y_rebin_plus_d.c_str(),i_pt_kaon+1,i_pt_kaon+2,"e");
    string KEY_y_rebin_minus_d = Form("dkaonminusrebin_y_kaonpt_%d",i_pt_kaon);
    h_mMass_y[KEY_y_rebin_minus_d] = (TH1F*) h_mMass_pty[KEY_pty]->ProjectionY(KEY_y_rebin_minus_d.c_str(),i_pt_kaon+1,i_pt_kaon+2,"e");

    //string KEY_y_rebin_plus = Form("kaonplusrebin_y_kaonpt_%d",i_pt_kaon);
    //h_mMass_y[KEY_y_rebin_plus]->Rebin(2);
    //string KEY_y_rebin_minus = Form("kaonminusrebin_y_kaonpt_%d",i_pt_kaon);
    //h_mMass_y[KEY_y_rebin_minus]->Rebin(2);

    string KEY_y_rebin_plus_rc = Form("rckaonplusrebin_y_kaonpt_%d",i_pt_kaon);
    h_mRcMass_y[KEY_y_rebin_plus_rc] = (TH1F*) h_mMass_RC["rc9_kplus_cent9"]->ProjectionY(KEY_y_rebin_plus_rc.c_str(),i_pt_kaon+1,i_pt_kaon+2,"e");
    string KEY_y_rebin_minus_rc = Form("rckaonminusrebin_y_kaonpt_%d",i_pt_kaon);
    h_mRcMass_y[KEY_y_rebin_minus_rc] = (TH1F*) h_mMass_RC["rc9_kminus_cent9"]->ProjectionY(KEY_y_rebin_minus_rc.c_str(),i_pt_kaon+1,i_pt_kaon+2,"e");

    double bin_low_edge_rc = h_mMass_RC["rc9_kplus_cent9"]->GetXaxis()->GetBinLowEdge(i_pt_kaon+1);
    double bin_up_edge_rc  = h_mMass_RC["rc9_kplus_cent9"]->GetXaxis()->GetBinLowEdge(i_pt_kaon+2);
    cout << "Bin RC  " << i_pt_kaon << ": [" << bin_low_edge_rc << ", " << bin_up_edge_rc << ")" << endl;

    string KEY_y_rebin_plus_ratio = Form("kaonplusrebin_y_kaonpt_%d_ratio",i_pt_kaon);
    string KEY_y_rebin_minus_ratio = Form("kaonminusrebin_y_kaonpt_%d_ratio",i_pt_kaon);
    h_mRatio_y[KEY_y_rebin_plus_ratio] = (TH1F*) h_mMass_y[KEY_y_rebin_plus_d]->Clone(KEY_y_rebin_plus_ratio.c_str());
    h_mRatio_y[KEY_y_rebin_plus_ratio]->Divide(h_mRcMass_y[KEY_y_rebin_plus_rc]);
    h_mRatio_y[KEY_y_rebin_minus_ratio] = (TH1F*) h_mMass_y[KEY_y_rebin_minus_d]->Clone(KEY_y_rebin_minus_ratio.c_str());
    h_mRatio_y[KEY_y_rebin_minus_ratio]->Divide(h_mRcMass_y[KEY_y_rebin_minus_rc]);
  } 

  TCanvas *cy = new TCanvas("cy","cy",10,10,800,400);
  cy->Divide(2,1);
  for(int i = 0; i < 2; i++)
  {
    cy->cd(i+1);
    cy->cd(i+1)->SetLeftMargin(0.15);
    cy->cd(i+1)->SetBottomMargin(0.15);
    cy->cd(i+1)->SetTicks(1,1);
    cy->cd(i+1)->SetGrid(0,0);  
  } 
  
  string outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_y_individualptbins_ratio_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim.c_str());
  string outputstart = Form("%s[",outputname.c_str()); 
  string outputstop = Form("%s]",outputname.c_str()); 

  cy->Print(outputstart.c_str());
  
  string histname_1[2] = {"plus","minus"};
  string histlabel[2] = {"K^{+}","K^{-}"};
  for(int ipt = 0; ipt < 31; ipt++)
  {
    for(int i = 0; i < 2; i++)
    {
      cy->cd(i+1);
      string KEY_ratio_y = Form("kaon%srebin_y_kaonpt_%d_ratio",histname_1[i].c_str(),ipt);   
      h_mRatio_y[KEY_ratio_y]->SetTitle(Form("%s Data/RC %1.3f<p_{T}<%1.3f, 20-60 Centrality",histlabel[i].c_str(),data_constants::kaon_pt_low[3*ipt],data_constants::kaon_pt_high[3*ipt+2]));
      h_mRatio_y[KEY_ratio_y]->GetXaxis()->SetTitle("y");
      h_mRatio_y[KEY_ratio_y]->GetYaxis()->SetTitle("Data/RC");

      h_mRatio_y[KEY_ratio_y]->Draw("pE");  
    }
    cy->Update();
    cy->Print(outputname.c_str());
  }
  cy->Update();
  cy->Print(outputstop.c_str());

  outputname = Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_y_individualptbins_datarc_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  cy->Print(outputstart.c_str());
  
  for(int ipt = 0; ipt < 31; ipt++)
  {
    for(int i = 0; i < 2; i++)
    {
      cy->cd(i+1);
      string KEY_y = Form("dkaon%srebin_y_kaonpt_%d",histname_1[i].c_str(),ipt);   
      string KEY_rc_y = Form("rckaon%srebin_y_kaonpt_%d",histname_1[i].c_str(),ipt);   
      h_mMass_y[KEY_y]->SetTitle(Form("%s %1.3f<p_{T}<%1.3f, 20-60 Centrality",histlabel[i].c_str(),data_constants::kaon_pt_low[3*ipt],data_constants::kaon_pt_high[3*ipt+2]));
      h_mMass_y[KEY_y]->GetXaxis()->SetTitle("y");
      h_mMass_y[KEY_y]->GetYaxis()->SetTitle("Counts");
     
      h_mMass_y[KEY_y]->SetMarkerStyle(20);
      h_mMass_y[KEY_y]->SetMarkerColor(kOrange+7);
      h_mMass_y[KEY_y]->SetLineColor(kOrange+7);

      h_mRcMass_y[KEY_rc_y]->SetMarkerStyle(24);
      h_mRcMass_y[KEY_rc_y]->SetMarkerColor(kBlue);
      h_mRcMass_y[KEY_rc_y]->SetLineColor(kBlue);
      
      int max = h_mMass_y[KEY_y]->GetMaximum();
      int min = h_mMass_y[KEY_y]->GetMinimum();
      if(h_mRcMass_y[KEY_rc_y]->GetMaximum() > max) max = h_mRcMass_y[KEY_rc_y]->GetMaximum();
      if(h_mRcMass_y[KEY_rc_y]->GetMinimum() < min) max = h_mRcMass_y[KEY_rc_y]->GetMinimum();
      h_mMass_y[KEY_y]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);

      h_mMass_y[KEY_y]->Draw("pE");  
      h_mRcMass_y[KEY_rc_y]->Draw("pE same");

      TLegend *leg = new TLegend(0.4,0.2,0.6,0.4);
      leg->AddEntry(h_mMass_y[KEY_y],"Data","p");
      leg->AddEntry(h_mRcMass_y[KEY_rc_y],"Pythia RC","p");
      leg->Draw("same");
    }
    cy->Update();
    cy->Print(outputname.c_str());
  }
  cy->Update();
  cy->Print(outputstop.c_str());
 
     
  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,400);     
  c1->Divide(2,1);
  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    c1->cd(i+1)->SetLeftMargin(0.15);
    c1->cd(i+1)->SetBottomMargin(0.15);
    c1->cd(i+1)->SetTicks(1,1);
    c1->cd(i+1)->SetGrid(0,0);  
  } 

  string histname[2] = {"kplus","kminus"};
  /////////////////////////////// No Cos2PhiStarPhi separation //////////////////////////////
  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_pty = Form("%s_pty",histname[i].c_str());   
    h_mMass_pty[KEY_pty]->SetTitle(Form("%s Data %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mMass_pty[KEY_pty]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mMass_pty[KEY_pty]->GetYaxis()->SetTitle("y");

    h_mMass_pty[KEY_pty]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pty_data_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim.c_str()));
 
  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_RC_pty = Form("rc9_%s_cent9",histname[i].c_str());   
    h_mMass_RC[KEY_RC_pty]->SetTitle(Form("%s RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mMass_RC[KEY_RC_pty]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mMass_RC[KEY_RC_pty]->GetYaxis()->SetTitle("y");

    h_mMass_RC[KEY_RC_pty]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pty_rc_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim.c_str()));

  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_ratio_pty = Form("%s_ratio_pty",histname[i].c_str());   
    h_mRatio_pty[KEY_ratio_pty]->SetTitle(Form("%s Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mRatio_pty[KEY_ratio_pty]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mRatio_pty[KEY_ratio_pty]->GetXaxis()->SetRangeUser(0.0,1.0);
    h_mRatio_pty[KEY_ratio_pty]->GetYaxis()->SetTitle("y");

    h_mRatio_pty[KEY_ratio_pty]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pty_datarcratio_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim.c_str()));

  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_ratio_pt = Form("%s_ratio_pt",histname[i].c_str());   
    h_mRatio_pt[KEY_ratio_pt]->SetTitle(Form("%s Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mRatio_pt[KEY_ratio_pt]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mRatio_pt[KEY_ratio_pt]->GetYaxis()->SetTitle("Data/RC");

    h_mRatio_pt[KEY_ratio_pt]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_ptonly_datarcratio_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim.c_str()));

  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_ratio_y = Form("%s_ratio_y",histname[i].c_str());   
    h_mRatio_y[KEY_ratio_y]->SetTitle(Form("%s Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mRatio_y[KEY_ratio_y]->GetXaxis()->SetTitle("y");
    h_mRatio_y[KEY_ratio_y]->GetYaxis()->SetTitle("Data/RC");

    h_mRatio_y[KEY_ratio_y]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_yonly_datarcratio_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim.c_str()));

  TLegend *legy = new TLegend(0.825,0.55,0.975,0.75);
  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_y = Form("%s_y",histname[i].c_str());   
    string KEY_RC_y = Form("%s_rc_y",histname[i].c_str());   
    h_mMass_y[KEY_y]->SetTitle(Form("%s %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mMass_y[KEY_y]->GetXaxis()->SetTitle("y");
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
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_yonly_datarc_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim.c_str()));
   
  TLegend *legpt = new TLegend(0.6,0.7,0.8,0.9);
  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    c1->cd(i+1)->SetLogy();
    string KEY_pt = Form("%s_pt",histname[i].c_str());   
    string KEY_RC_pt = Form("%s_rc_pt",histname[i].c_str());   
    h_mMass_pt[KEY_pt]->SetTitle(Form("%s %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mMass_pt[KEY_pt]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mMass_pt[KEY_pt]->GetYaxis()->SetTitle("Counts");
    h_mMass_pt[KEY_pt]->SetMarkerStyle(20);
    h_mMass_pt[KEY_pt]->SetMarkerColor(kOrange+7);
    h_mMass_pt[KEY_pt]->SetLineColor(kOrange+7);

    int min = h_mMass_pt[KEY_pt]->GetMinimum();
    int max = h_mMass_pt[KEY_pt]->GetMaximum();
    if(h_mRcMass_pt[KEY_RC_pt]->GetMinimum() < min) min = h_mRcMass_pt[KEY_RC_pt]->GetMinimum();
    if(h_mRcMass_pt[KEY_RC_pt]->GetMaximum() > max) max = h_mRcMass_pt[KEY_RC_pt]->GetMaximum();
    h_mMass_pt[KEY_pt]->GetYaxis()->SetRangeUser(1/*min*0.9*/,max*2.0);

    h_mMass_pt[KEY_pt]->Draw("pE");

    h_mRcMass_pt[KEY_RC_pt]->SetMarkerStyle(24);
    h_mRcMass_pt[KEY_RC_pt]->SetMarkerColor(kBlack);
    h_mRcMass_pt[KEY_RC_pt]->SetLineColor(kBlack);
    h_mRcMass_pt[KEY_RC_pt]->Draw("pE same");

    if(i == 0)
    {
      legpt->AddEntry(h_mMass_pt[KEY_pt],"Data","p");
      legpt->AddEntry(h_mRcMass_pt[KEY_RC_pt],"RC","p");
    }
    legpt->Draw("same");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_ptonly_datarc_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim.c_str()));
  File_OutPut->Close();
}

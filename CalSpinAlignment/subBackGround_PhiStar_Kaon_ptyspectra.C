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

void subBackGround_PhiStar_Kaon_ptyspectra(int energy = 4, int pid = 0, int year = 0, string date = "20240426", bool random3D = false, int order = 2, string etamode = "eta1_eta1", int deltaonly = 0)
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
  //string folder = "CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240523_fixed2Dfunc_PhiSpectraWeighted_rerho1n10.0_ysigma1000_noalignment";
  //string folder = "CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_PhiWeight_rerho1n10.0_ysigma1000_noalignment";
  //string folder = "CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_PhiAndKaonWeight_rerho1n10.0_ysigma1000_noalignment";
  //string folder = "CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_PhiAndKaonWeightAgain_rerho1n10.0_ysigma1000_noalignment";
  //string folder = "CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_NoWeight_rerho1n10.0_ysigma1000_noalignment";
 /// string folder = "Pythia_PhiEmbedding_PIDEff";
  string folder = "PhiEmbeddingWithWeights_PIDEff_20240726";
  string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",folder.c_str());  
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
    }
  }

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
        string KEY_SM = Form("kaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        string KEY_phiy = Form("kaon_phiy_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        if(i_phi == 0) h_mMass_phiy[KEY_phiy] = (TH2F*)h_mMass_SM[KEY_SM]->Clone(KEY_phiy.c_str()); // First bins
        else h_mMass_phiy[KEY_phiy]->Add(h_mMass_SM[KEY_SM],1.0); // All other bins

        string KEY_RC = Form("rckaon_phipt_%d_cos2phistarphi_%d_kaonpt_%d",i_pt_phi,i_phi,i_pt_kaon);
        string KEY_RC_phiy = Form("rckaon_phiy_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        if(i_phi == 0) h_mRcMass_phiy[KEY_RC_phiy] = (TH2F*)h_mMass_RC[KEY_RC]->Clone(KEY_RC_phiy.c_str()); // First bins
        else h_mRcMass_phiy[KEY_RC_phiy]->Add(h_mMass_RC[KEY_RC],1.0); // All other bins

        //string KEY_phi = Form("kaon_phi_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        string KEY_y   = Form("kaon_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        if(i_phi == 9) // Last bins
        {
          //h_mMass_phi[KEY_phi] = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionX(KEY_phi.c_str(),0,-1,"e"));
          h_mMass_y[KEY_y]     = (TH1F*) (h_mMass_phiy[KEY_phiy]->ProjectionY(KEY_y.c_str(),0,-1,"e"));
        } 

        //string KEY_RC_phi = Form("rckaon_phi_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        string KEY_RC_y   = Form("rckaon_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
        if(i_phi == 9) // Last bins
        {
          //h_mRcMass_phi[KEY_RC_phi] = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionX(KEY_RC_phi.c_str(),0,-1,"e"));
          h_mRcMass_y[KEY_RC_y]     = (TH1F*) (h_mRcMass_phiy[KEY_RC_phiy]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e"));
        } 
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
  string outputfile = "../output/kaon_ptyratio_19GeV.root";
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  cout << "AFTER CREATING FILE" << endl;  

  for(Int_t i_pt_phi = 0; i_pt_phi < 4; i_pt_phi++) // pt bin 
  {
    for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins; i_pt_kaon++)
    {
      int rebin_pt_kaon = i_pt_kaon/rebinconstant;
      string KEY_y_rebin = Form("kaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,rebin_pt_kaon);
      string KEY_RC_y_rebin = Form("rckaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,rebin_pt_kaon);

      string KEY_y = Form("kaon_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
      string KEY_RC_y = Form("rckaon_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
      if(i_pt_kaon%rebinconstant == 0)
      {
        h_mMass_y[KEY_y_rebin]      = (TH1F*) h_mMass_y[KEY_y]->Clone(KEY_y_rebin.c_str());
        h_mRcMass_y[KEY_RC_y_rebin] = (TH1F*) h_mRcMass_y[KEY_RC_y]->Clone(KEY_RC_y_rebin.c_str());
      }
      else
      {
        h_mMass_y[KEY_y_rebin]->Add(h_mMass_y[KEY_y],1.0);     
        h_mRcMass_y[KEY_RC_y_rebin]->Add(h_mRcMass_y[KEY_RC_y],1.0);
      }
    }
  
    string KEY_pty    = Form("pty_phipt_%d",i_pt_phi);   
    string KEY_RC_pty = Form("rc_pty_phipt_%d",i_pt_phi);   
    string KEY_ratio_pty = Form("ratio_pty_phipt_%d",i_pt_phi);   
    h_mMass_pty[KEY_pty] = new TH2F(KEY_pty.c_str(),KEY_pty.c_str(),pt_total_rebinned,kaon_pt_low_rebin,data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);
    h_mRcMass_pty[KEY_RC_pty] = new TH2F(KEY_RC_pty.c_str(),KEY_RC_pty.c_str(),pt_total_rebinned,kaon_pt_low_rebin,data_constants::kaon_rapidity_bins*2,-data_constants::kaon_rapidity_max,data_constants::kaon_rapidity_max);

    for(Int_t i_pt_kaon = 0; i_pt_kaon < data_constants::kaon_pt_bins/rebinconstant; i_pt_kaon++)
    {
      string KEY_y_rebin = Form("kaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
      string KEY_RC_y_rebin = Form("rckaonrebin_y_phipt_%d_kaonpt_%d",i_pt_phi,i_pt_kaon);
      for(int ibin = 1; ibin <= h_mMass_y[KEY_y_rebin]->GetNbinsX(); ibin++)
      {
        double val = h_mMass_y[KEY_y_rebin]->GetBinContent(ibin);
        double err = h_mMass_y[KEY_y_rebin]->GetBinError(ibin);
        double valRC = h_mRcMass_y[KEY_RC_y_rebin]->GetBinContent(ibin);
        double errRC = h_mRcMass_y[KEY_RC_y_rebin]->GetBinError(ibin);

        //cout << "bin = " << ibin << ", bincenter = " << h_mMass_y[KEY_y_rebin]->GetBinCenter(ibin) << endl;
        
        h_mMass_pty[KEY_pty]->SetBinContent(i_pt_kaon+1,ibin,val);
        h_mMass_pty[KEY_pty]->SetBinError(i_pt_kaon+1,ibin,err);
        h_mRcMass_pty[KEY_RC_pty]->SetBinContent(i_pt_kaon+1,ibin,valRC);
        h_mRcMass_pty[KEY_RC_pty]->SetBinError(i_pt_kaon+1,ibin,errRC);
      }
    }

    h_mMass_pty[KEY_pty]->RebinY(2);
    h_mRcMass_pty[KEY_RC_pty]->RebinY(2);
 
    double yields = h_mMass_pty[KEY_pty]->Integral(0,-1,0,-1);
    double rcyields = h_mRcMass_pty[KEY_RC_pty]->Integral(0,-1,0,-1);
 
    h_mRcMass_pty[KEY_RC_pty]->Scale(yields/rcyields);

    h_mRatio_pty[KEY_ratio_pty] = (TH2F*) h_mMass_pty[KEY_pty]->Clone(KEY_ratio_pty.c_str());
    h_mRatio_pty[KEY_ratio_pty]->Divide(h_mRcMass_pty[KEY_RC_pty]);
    
    string KEY_pt    = Form("pt_phipt_%d",i_pt_phi);   
    string KEY_RC_pt = Form("rc_pt_phipt_%d",i_pt_phi);   
    string KEY_ratio_pt = Form("ratio_pt_phipt_%d",i_pt_phi);   
    h_mMass_pt[KEY_pt] = (TH1F*) h_mMass_pty[KEY_pty]->ProjectionX(KEY_pt.c_str(),0,-1,"e");
    h_mRcMass_pt[KEY_RC_pt] = (TH1F*) h_mRcMass_pty[KEY_RC_pty]->ProjectionX(KEY_RC_pt.c_str(),0,-1,"e");
    h_mRatio_pt[KEY_ratio_pt] = (TH1F*) h_mMass_pt[KEY_pt]->Clone(KEY_ratio_pt.c_str());
    h_mRatio_pt[KEY_ratio_pt]->Divide(h_mRcMass_pt[KEY_RC_pt]);

    string KEY_y    = Form("y_phipt_%d",i_pt_phi);   
    string KEY_RC_y = Form("rc_y_phipt_%d",i_pt_phi);   
    string KEY_ratio_y = Form("ratio_y_phipt_%d",i_pt_phi);   
    h_mMass_y[KEY_y] = (TH1F*) h_mMass_pty[KEY_pty]->ProjectionY(KEY_y.c_str(),0,-1,"e");
    h_mRcMass_y[KEY_RC_y] = (TH1F*) h_mRcMass_pty[KEY_RC_pty]->ProjectionY(KEY_RC_y.c_str(),0,-1,"e");
    h_mRatio_y[KEY_ratio_y] = (TH1F*) h_mMass_y[KEY_y]->Clone(KEY_ratio_y.c_str());
    h_mRatio_y[KEY_ratio_y]->Divide(h_mRcMass_y[KEY_RC_y]);

    h_mRatio_pty[KEY_ratio_pty]->Write();
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

  /////////////////////////////// No Cos2PhiStarPhi separation //////////////////////////////
  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_pty = Form("pty_phipt_%d",i);   
    h_mMass_pty[KEY_pty]->SetTitle(Form("Data %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    h_mMass_pty[KEY_pty]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mMass_pty[KEY_pty]->GetYaxis()->SetTitle("y");

    h_mMass_pty[KEY_pty]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pty_data.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));
 
  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_RC_pty = Form("rc_pty_phipt_%d",i);   
    h_mRcMass_pty[KEY_RC_pty]->SetTitle(Form("RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    h_mRcMass_pty[KEY_RC_pty]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mRcMass_pty[KEY_RC_pty]->GetYaxis()->SetTitle("y");

    h_mRcMass_pty[KEY_RC_pty]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pty_rc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_ratio_pty = Form("ratio_pty_phipt_%d",i);   
    h_mRatio_pty[KEY_ratio_pty]->SetTitle(Form("Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    h_mRatio_pty[KEY_ratio_pty]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mRatio_pty[KEY_ratio_pty]->GetYaxis()->SetTitle("y");

    h_mRatio_pty[KEY_ratio_pty]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pty_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_ratio_pt = Form("ratio_pt_phipt_%d",i);   
    h_mRatio_pt[KEY_ratio_pt]->SetTitle(Form("Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    h_mRatio_pt[KEY_ratio_pt]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mRatio_pt[KEY_ratio_pt]->GetYaxis()->SetTitle("Data/RC");

    h_mRatio_pt[KEY_ratio_pt]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_ptonly_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_ratio_y = Form("ratio_y_phipt_%d",i);   
    h_mRatio_y[KEY_ratio_y]->SetTitle(Form("Data/RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    h_mRatio_y[KEY_ratio_y]->GetXaxis()->SetTitle("y");
    h_mRatio_y[KEY_ratio_y]->GetYaxis()->SetTitle("Data/RC");

    h_mRatio_y[KEY_ratio_y]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_yonly_datarcratio.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));

  TLegend *legy = new TLegend(0.825,0.55,0.975,0.75);
  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_y = Form("y_phipt_%d",i);   
    string KEY_RC_y = Form("rc_y_phipt_%d",i);   
    h_mMass_y[KEY_y]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
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
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_yonly_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));
   
  TLegend *legpt = new TLegend(0.6,0.7,0.8,0.9);
  for(int i = 0; i < 4; i++)
  {
    c1->cd(i+1);
    string KEY_pt = Form("pt_phipt_%d",i);   
    string KEY_RC_pt = Form("rc_pt_phipt_%d",i);   
    h_mMass_pt[KEY_pt]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    h_mMass_pt[KEY_pt]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mMass_pt[KEY_pt]->GetYaxis()->SetTitle("Counts");
    h_mMass_pt[KEY_pt]->SetMarkerStyle(20);
    h_mMass_pt[KEY_pt]->SetMarkerColor(kOrange+7);
    h_mMass_pt[KEY_pt]->SetLineColor(kOrange+7);

    int min = h_mMass_pt[KEY_pt]->GetMinimum();
    int max = h_mMass_pt[KEY_pt]->GetMaximum();
    if(h_mRcMass_pt[KEY_RC_pt]->GetMinimum() < min) min = h_mRcMass_pt[KEY_RC_pt]->GetMinimum();
    if(h_mRcMass_pt[KEY_RC_pt]->GetMaximum() > max) max = h_mRcMass_pt[KEY_RC_pt]->GetMaximum();
    h_mMass_pt[KEY_pt]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

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
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_ptonly_datarc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str()));
  File_OutPut->Close();
}

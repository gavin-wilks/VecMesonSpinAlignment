#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "phi_data_constants_19GeV.h"
//#include "../Utility/draw.h"
#include "../../../draw.h"
#include "../Utility/functions.h"
#include "resolution_pt.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void phiembedding_pid(Int_t energy = 4, Int_t pid = 0, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 50, string sim_py = "Rc5", string sim_em = "Rc8")//defaultF = 0 is BESII, defaultF = 1 is BESI
{
 
  std::string spectra = "";
  if(yspectra == 0) spectra = "_NoRapiditySpectra";
  if(yspectra == 1) spectra = "_NoRapiditySpectra_EP";
  if(yspectra == 2) spectra = "_NoRapiditySpectra_EP_noV2";
  if(yspectra == 3) spectra = "_WithRapiditySpectra";
  if(yspectra == 4) spectra = "_WithRapiditySpectra_HalfSigma";
  if(yspectra == 5) spectra = "_NoRapiditySpectra_FixedFirstEP";
  if(yspectra == 6) spectra = "_NoRapiditySpectra_FixedFirstEP_InputRho";
  if(yspectra == 7) spectra = "_NoRapiditySpectra_InputRho0p4_Fixed";
  if(yspectra == 8) spectra = "_NoRapiditySpectra_v2times3";
  if(yspectra == 9) spectra = "_NoRapiditySpectra_FixedFirstEP_AcceptanceOnly";
  if(yspectra == 10) spectra = "_NoRapiditySpectra_FixedFirstEP_ptmin0p5";
  if(yspectra == 11) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins";
  if(yspectra == 12) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_InputRho0p5";
  if(yspectra == 13) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_InputRho0p5_3xV2";
  if(yspectra == 14) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_Eta0p5";
  if(yspectra == 15) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_rhoforeachstudy";
  if(yspectra == 16) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP";
  if(yspectra == 17) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc";
  if(yspectra == 18) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc_pT0p2";
  if(yspectra == 19) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc_pT0p2_TPC";
  if(yspectra == 20) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc_pT0p2_TPC_TOF";
  if(yspectra == 21) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_Acc1Rc_pT0p1_TPC_TOF_kaons_20240327";
  if(yspectra == 22) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_phis_20240423";
  if(yspectra == 23) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_phis_20240426_withRho_deltar";
  if(yspectra == 24) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240508_fixed2Dfunc_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0";
  if(yspectra == 25) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240508_fixed2Dfunc_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0_SpectraWeight";
  if(yspectra == 26) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240523_fixed2Dfunc_PhiSpectra_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0_helicityrho0.3333";
  if(yspectra == 27) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_NoWeight_rerho1n10.0_ysigma1000_noalignment";
  if(yspectra == 28) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_PhiWeight_rerho1n10.0_ysigma1000_noalignment";
  if(yspectra == 29) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_PhiAndKaonWeight_rerho1n10.0_ysigma1000_noalignment";
  if(yspectra == 30) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_PhiAndKaonWeightAgain_rerho1n10.0_ysigma1000_noalignment";
  
  if(yspectra == 35) spectra = "PhiEmbeddingWithWeights";
  if(yspectra == 36) spectra = "PhiEmbeddingWithWeights_PIDEff_20240722";
  if(yspectra == 37) spectra = "PhiEmbeddingNoWeights";
  if(yspectra == 38) spectra = "PhiEmbeddingWithWeights_PIDEff_20240725";
  if(yspectra == 39) spectra = "PhiEmbeddingWithWeights_PIDEff_20240726";
  if(yspectra == 50) spectra = "PhiEmbeddingNoWeights_PIDEff_20240730";

  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  //string inputfile = Form("effaccfiles/%s/%s/%s/Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  //TFile *File_InPut = TFile::Open(inputfile.c_str());

  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_Weights_PIDEff_2.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_WithWeights_PIDEff_20240722_1.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_NoWeights_PIDEff_20240723_1.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_NoWeights_PIDEff_20240724_fourMomentum_1.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_WithWeights_PIDEff_20240725_FixedMC_0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_WithWeights_PIDEff_20240725_FixedMC_FixedAngle_0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_NoWeights_PIDEff_20240730_FixedMC_FixedAngle_RC_0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  TFile *File_InPutRC = TFile::Open(inputfileRC.c_str());

  TH3F *h_mCounts_RC[2];  
  TH1F *h_mCounts_RC_m2[2][data_constants::kaon_pt_bins];  
  TH1F *h_mCounts_RC_nsigma[2][data_constants::kaon_pt_bins];  
  TH1F *h_mCounts_RC_pt[2][4]; // 0: no cut, 1: m2 cut only, 2: nsigma cut only, 3: both cuts  
  
  int cutlowm2[4]  = {0 ,16,0 ,16};
  int cuthighm2[4] = {-1,35,-1,35};
  int cutlowns[4]  = {0 ,0 ,26,26};
  int cuthighns[4] = {-1,-1,75,75};

  cout << inputfileRC << endl;
 
  string particle[2] = {"kplus","kminus"};

  for(int ipar = 0; ipar < 2; ipar++)
  {
    for(int i_cent = 2; i_cent <= 5; i_cent++)
    { 
      string KEY_Embedding = Form("rc8_%s_m2_cent%d",particle[ipar].c_str(),i_cent);
      string KEY_Embedding9 = Form("rc8_%s_m2_cent%d",particle[ipar].c_str(),9);
      cout << KEY_Embedding << endl;
      if(i_cent == 2) h_mCounts_RC[ipar] = (TH3F*)((TH3F*) File_InPutRC->Get(KEY_Embedding.c_str()))->Clone(KEY_Embedding9.c_str());
      else h_mCounts_RC[ipar]->Add((TH3F*) File_InPutRC->Get(KEY_Embedding.c_str()),1.0);

      if(i_cent == 5)
      {
        for(int ipt = 0; ipt < data_constants::kaon_pt_bins; ipt++)
        {
          string KEY = Form("m2_pt%d",ipt);
          h_mCounts_RC_m2[ipar][ipt] = (TH1F*) h_mCounts_RC[ipar]->ProjectionY(KEY.c_str(),1+ipt,2+ipt,0,-1,"e");
          KEY = Form("nsigma_pt%d",ipt);
          h_mCounts_RC_nsigma[ipar][ipt] = (TH1F*) h_mCounts_RC[ipar]->ProjectionZ(KEY.c_str(),1+ipt,2+ipt,0,-1,"e");
        }
        for(int icut = 0; icut < 4; icut++)
        {
          string KEY = Form("pt_cut%d",icut);
          h_mCounts_RC_pt[ipar][icut] = (TH1F*) h_mCounts_RC[ipar]->ProjectionX(KEY.c_str(),cutlowm2[icut],cuthighm2[icut],cutlowns[icut],cuthighns[icut],"e");
          h_mCounts_RC_pt[ipar][icut]->Rebin(3);
          h_mCounts_RC_pt[ipar][icut]->Print();
        }
      }
    }
  }

  TCanvas *c = new TCanvas("c","c",10,10,800,400);
  c->Divide(2,1);

  int marker[4] = {kBlack,kOrange+7,kBlue,kGray+2};
  string cut[4] = {"no cut","m^{2} cut","n#sigma_{K} cut","both cuts"};

  TLegend *leg = new TLegend(0.5,0.2,0.7,0.4);
  for(int ipar = 0; ipar < 2; ipar++)
  {
    c->cd(ipar+1)->SetLeftMargin(0.15);
    c->cd(ipar+1)->SetBottomMargin(0.15);
    c->cd(ipar+1)->SetTicks(1,1);
    c->cd(ipar+1)->SetGrid(0,0);
    c->cd(ipar+1)->SetLogy();
    double max = -1000000000;
    double min = 1000000000;
    for(int icut = 0; icut < 4; icut++)
    {
      double tmax = h_mCounts_RC_pt[ipar][icut]->GetMaximum();
      double tmin = h_mCounts_RC_pt[ipar][icut]->GetMinimum();
      if(tmax > max) max = tmax;
      if(tmin < min) min = tmin;
    }

    for(int icut = 0; icut < 4; icut++)
    {
      h_mCounts_RC_pt[ipar][icut]->SetTitle(Form("%s Centrality %d",particle[ipar].c_str(),9));
      h_mCounts_RC_pt[ipar][icut]->GetXaxis()->SetTitle("p_{T}");
      h_mCounts_RC_pt[ipar][icut]->GetYaxis()->SetTitle("Counts");
      h_mCounts_RC_pt[ipar][icut]->SetMarkerStyle(20);
      h_mCounts_RC_pt[ipar][icut]->SetMarkerColor(marker[icut]);
      h_mCounts_RC_pt[ipar][icut]->SetLineColor(marker[icut]);
 
      cout << "Max = " << max << ", Min = " << min << endl;

      h_mCounts_RC_pt[ipar][icut]->GetYaxis()->SetRangeUser(0.01,max*1.1);

      if(icut == 0) h_mCounts_RC_pt[ipar][icut]->Draw("pE");
      else h_mCounts_RC_pt[ipar][icut]->Draw("pE same");
       
      if(ipar == 0) leg->AddEntry(h_mCounts_RC_pt[ipar][icut],cut[icut].c_str(),"p");
    }
    leg->Draw("same"); 
  }
  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/PythiaEmbedding_ptspectraCuts_%s_Order%d_%s_Cent%d_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(), spectra.c_str(),correction.c_str(),order,etamode.c_str(),9,sim_em.c_str()));  

  string outputname = Form("figures/%s/%s/pTstudy/%s/PythiaEmbedding_m2_%s_Order%d_%s_Cent%d_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(), spectra.c_str(),correction.c_str(),order,etamode.c_str(),9,sim_em.c_str());
  string outputstart = Form("%s[",outputname.c_str());
  string outputstop = Form("%s]",outputname.c_str());
 
  c->Print(outputstart.c_str());

  for(int ipt = 0; ipt < data_constants::kaon_pt_bins; ipt++)
  {
    for(int ipar = 0; ipar < 2; ipar++)
    {
      c->cd(ipar+1);
      h_mCounts_RC_m2[ipar][ipt]->SetTitle(Form("%s Centrality %d, %1.3f<p_{T}<%1.3f",particle[ipar].c_str(),9,data_constants::kaon_pt_low[ipt],data_constants::kaon_pt_high[ipt]));
      h_mCounts_RC_m2[ipar][ipt]->GetXaxis()->SetTitle("m^{2}");
      h_mCounts_RC_m2[ipar][ipt]->GetYaxis()->SetTitle("Counts");
      h_mCounts_RC_m2[ipar][ipt]->SetMarkerStyle(20);
      h_mCounts_RC_m2[ipar][ipt]->Draw("pE");
    }
    c->Update();
    c->Print(outputname.c_str());
  }
  c->Print(outputstop.c_str()); 

  outputname = Form("figures/%s/%s/pTstudy/%s/PythiaEmbedding_nsigmakaon_%s_Order%d_%s_Cent%d_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(), spectra.c_str(),correction.c_str(),order,etamode.c_str(),9,sim_em.c_str());
  outputstart = Form("%s[",outputname.c_str());
  outputstop = Form("%s]",outputname.c_str());
 
  c->Print(outputstart.c_str());

  for(int ipt = 0; ipt < data_constants::kaon_pt_bins; ipt++)
  {
    for(int ipar = 0; ipar < 2; ipar++)
    {
      c->cd(ipar+1);
      h_mCounts_RC_nsigma[ipar][ipt]->SetTitle(Form("%s Centrality %d, %1.3f<p_{T}<%1.3f",particle[ipar].c_str(),9,data_constants::kaon_pt_low[ipt],data_constants::kaon_pt_high[ipt]));
      h_mCounts_RC_nsigma[ipar][ipt]->GetXaxis()->SetTitle("n#sigma_{K}");
      h_mCounts_RC_nsigma[ipar][ipt]->GetYaxis()->SetTitle("Counts");
      h_mCounts_RC_nsigma[ipar][ipt]->SetMarkerStyle(20);
      h_mCounts_RC_nsigma[ipar][ipt]->Draw("pE");
    }
    c->Update();
    c->Print(outputname.c_str());
  }
  c->Print(outputstop.c_str()); 

}


#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
//#include "../Utility/draw.h"
#include "../../../draw.h"
#include "../Utility/functions.h"
#include "resolution_pt.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void compareDataRcPhiStar_PhiEmbed(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 46)//defaultF = 0 is BESII, defaultF = 1 is BESI
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
  if(yspectra == 26) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240520_fixed2Dfunc_FinerBins_helicityandglobalrho00fromdata";
  if(yspectra == 27) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240523_fixed2Dfunc_PhiSpectra_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0_helicityrho0.3333";
  if(yspectra == 28) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240524_fixed2Dfunc_OnlyCosWeighted_rerho1n10.0_ysigma1000_noalignment";
  if(yspectra == 29) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240507_fixed2Dfunc_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0";

  std::string filename = "";
  if(yspectra == 35) 
  {
    spectra = "PhiEmbeddingWithWeights";
    filename = "PhiEmbed_Full_WithWeights_1.root";
  }
  if(yspectra == 36) 
  {
    spectra = "PhiEmbeddingWithWeights_PIDEff";    
    filename = "PhiEmbed_Weights_PIDEff_2.root";
  }
  if(yspectra == 37)
  {
    spectra = "Pythia_PhiEmbedding_PIDEff";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 39)
  {
    spectra = "PhiEmbeddingWithWeights_PIDEff_20240726";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 40)
  {
    spectra = "Pythia_PIDEff_Helicity2D_20240729";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 41)
  {
    spectra = "Pythia_PIDEff_Helicity2Dplusdelta_20240729";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 42)
  {
    spectra = "Pythia_PIDEff_Helicity2Dminusdelta_20240729";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 43)
  {
    spectra = "Pythia_PIDEff_Helicity2Dminusdeltaplusrhodelta_20240729";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 44)
  {
    spectra = "Pythia_PIDEff_20240729";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 45)
  {
    spectra = "Pythia_PIDEff_Helicity2Dminusdeltaplusrhodelta_20240730";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 46)
  {
    spectra = "Pythia_PIDEff_Sigmay1p5_20240731";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
 
  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  string inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_%sPhiPtSys_%s_PolySys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_%sPhiPtSys_%s_PolySys_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1D *h_mCounts[4];  

  //string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra.c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputfileRC = Form("effaccfiles/%s/%s/%s/%s",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),filename.c_str());
  TFile *File_InPutRC = TFile::Open(inputfileRC.c_str());
  TH1D *h_mCounts_RC[4];  

  TH1D *h_mCounts_Ratio[4];  
  
  cout << inputfile << endl;
  cout << inputfileRC << endl;
 
  for(int i_pt = 2; i_pt < 6; ++i_pt) // pt loop
  {
    string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    h_mCounts[i_pt-2] = (TH1D*) File_InPut->Get(KEY_counts.c_str());
    string KEY_counts_RC = Form("h_mRc5EffPhiS_Cent_9_Pt_%d",i_pt);
    h_mCounts_RC[i_pt-2] = (TH1D*) File_InPutRC->Get(KEY_counts_RC.c_str());
    if(yspectra < 6) h_mCounts_RC[i_pt-2]->Rebin(2);

    float data = h_mCounts[i_pt-2]->Integral(1,10);
    float rc   = h_mCounts_RC[i_pt-2]->Integral(1,10);

    h_mCounts_RC[i_pt-2]->Scale(data/rc);

    h_mCounts_Ratio[i_pt-2] = (TH1D*) h_mCounts[i_pt-2]->Clone();
    h_mCounts_Ratio[i_pt-2]->Divide(h_mCounts_RC[i_pt-2]);

  }

  TGraphAsymmErrors *gslope = new TGraphAsymmErrors();
  TGraphAsymmErrors *gintercept = new TGraphAsymmErrors();
  TF1 *fits[4];
  float pt[4] = {1.5,2.1,2.7,3.6};

  TCanvas *c = new TCanvas("c","c",10,10,800,800);
  c->Divide(2,2);
  for(int i = 0; i < 4; i++)
  {
    c->cd(i+1);
    c->cd(i+1)->SetLeftMargin(0.15);
    c->cd(i+1)->SetBottomMargin(0.15);
    c->cd(i+1)->SetTicks(1,1);
    c->cd(i+1)->SetGrid(0,0);
    
    h_mCounts_Ratio[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    h_mCounts_Ratio[i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
    h_mCounts_Ratio[i]->GetYaxis()->SetTitle("(Data)/(RC from Simulation)");
    h_mCounts_Ratio[i]->Draw("pE");
   
    fits[i] = new TF1(Form("fit%d",i),"[0]+[1]*x",-1.0,1.0);
    h_mCounts_Ratio[i]->Fit(fits[i],"NMI");
    fits[i]->SetLineColor(kRed);
    fits[i]->Draw("l same");
    
    gslope->SetPoint(i,pt[i],fits[i]->GetParameter(1)); 
    gslope->SetPointError(i,0.0,0.0,fits[i]->GetParError(1),fits[i]->GetParError(1));
    gintercept->SetPoint(i,pt[i],fits[i]->GetParameter(0)); 
    gintercept->SetPointError(i,0.0,0.0,fits[i]->GetParError(0),fits[i]->GetParError(0));
  }

  c->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_Ratio_%s_Order%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectra.c_str()));  

  TCanvas *cfit = new TCanvas("cfit","cfit",10,10,800,400);
  cfit->Divide(2,1);
  
  cfit->cd(1);
  cfit->cd(1)->SetLeftMargin(0.15);  
  cfit->cd(1)->SetBottomMargin(0.15);
  cfit->cd(1)->SetTicks(1,1);
  cfit->cd(1)->SetGrid(0,0);

  gslope->SetMarkerStyle(20);
  gslope->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gslope->GetYaxis()->SetTitle("d(Data/RC)/d(cos(2#phi*-2#phi))");
  gslope->SetMarkerSize(1.2);
  gslope->Draw("APE"); 

  cfit->cd(2);
  cfit->cd(2)->SetLeftMargin(0.15);  
  cfit->cd(2)->SetBottomMargin(0.15);
  cfit->cd(2)->SetTicks(1,1);
  cfit->cd(2)->SetGrid(0,0);

  gintercept->SetMarkerStyle(20);
  gintercept->SetMarkerSize(1.2);
  gintercept->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gintercept->GetYaxis()->SetTitle("(Data/RC) intercept");
  gintercept->Draw("APE"); 

  cfit->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_FitParams_%s_Order%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectra.c_str()));  

  TLegend *leg = new TLegend(0.4,0.7,0.6,0.9);
  for(int i = 0; i < 4; i++)
  {
    c->cd(i+1);
    c->cd(i+1)->SetLeftMargin(0.15);
    c->cd(i+1)->SetBottomMargin(0.15);
    c->cd(i+1)->SetTicks(1,1);
    c->cd(i+1)->SetGrid(0,0);
    
    h_mCounts[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    h_mCounts[i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
    h_mCounts[i]->GetYaxis()->SetTitle("Counts");
    h_mCounts[i]->SetMarkerStyle(20);
    h_mCounts[i]->SetMarkerColor(kOrange+7);
    h_mCounts[i]->SetLineColor(kOrange+7);

    int min = h_mCounts[i]->GetMinimum();
    int max = h_mCounts[i]->GetMaximum();
    if(h_mCounts_RC[i]->GetMinimum() < min) min = h_mCounts_RC[i]->GetMinimum();
    if(h_mCounts_RC[i]->GetMaximum() > max) max = h_mCounts_RC[i]->GetMaximum();
    h_mCounts[i]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

    h_mCounts[i]->Draw("pE");

    h_mCounts_RC[i]->SetMarkerStyle(24);
    h_mCounts_RC[i]->SetMarkerColor(kBlack);
    h_mCounts_RC[i]->SetLineColor(kBlack);
    h_mCounts_RC[i]->Draw("pE same");

    if(i == 0)
    {
      leg->AddEntry(h_mCounts[i],"Data","p");
      leg->AddEntry(h_mCounts_RC[i],"RC","p");
    }
    leg->Draw("same");
  }

  c->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_Comparison_%s_Order%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectra.c_str()));  


}

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nPt = g_rho->GetN();
  TBox *bSys[nPt];
  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot sys errors
  {
    double pt, rho;
    g_rho->GetPoint(i_pt,pt,rho);
    double err = g_rho->GetErrorYhigh(i_pt);

    bSys[i_pt] = new TBox(pt-0.08,rho-err,pt+0.08,rho+err);
    bSys[i_pt]->SetFillColor(0);
    bSys[i_pt]->SetFillStyle(0);
    bSys[i_pt]->SetLineStyle(1);
    bSys[i_pt]->SetLineWidth(1);
    bSys[i_pt]->SetLineColor(plot_color);
    bSys[i_pt]->Draw("l Same");
  }
}


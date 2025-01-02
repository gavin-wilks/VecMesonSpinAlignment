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

void compareSimV2(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 1)//defaultF = 0 is BESII, defaultF = 1 is BESI
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
  
  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  //string inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_%sPhiPtSys_%s_PolySys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  //if(order == 1) inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_%sPhiPtSys_%s_PolySys_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  //TFile *File_InPut = TFile::Open(inputfile.c_str());
  //TH1D *h_mCounts[4];  

  //string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra.c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra.c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPutRC = TFile::Open(inputfileRC.c_str());
  TH2D *h_mCounts2D_MC[4];  
  TH2D *h_mCounts2D_RC[4];  
  TH1D *h_mCounts_MC[4];  
  TH1D *h_mCounts_RC[4];  

  TH1D *h_mCounts_Ratio[4];  

  TGraphAsymmErrors *g_mMcv2 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_mRcv2 = new TGraphAsymmErrors();
  
  //cout << inputfile << endl;
  cout << inputfileRC << endl;
 
  for(int i_pt = 2; i_pt < 6; ++i_pt) // pt loop
  {
    string KEY_counts_RC = Form("h_mRcEffCosEP_Cent_9_Pt_%d",i_pt);
    h_mCounts2D_RC[i_pt-2] = (TH2D*) File_InPutRC->Get(KEY_counts_RC.c_str());
    h_mCounts_RC[i_pt-2] = h_mCounts2D_RC[i_pt-2]->ProjectionY();
    string KEY_counts_MC = Form("h_mMcEffCosEP_Cent_9_Pt_%d",i_pt);
    h_mCounts2D_MC[i_pt-2] = (TH2D*) File_InPutRC->Get(KEY_counts_MC.c_str());
    h_mCounts_MC[i_pt-2] = h_mCounts2D_MC[i_pt-2]->ProjectionY();


    //float data = h_mCounts[i_pt-2]->Integral(1,10);
    //float rc   = h_mCounts_RC[i_pt-2]->Integral(1,10);

    //h_mCounts_RC[i_pt-2]->Scale(data/rc);

    //h_mCounts_Ratio[i_pt-2] = (TH1D*) h_mCounts[i_pt-2]->Clone();
    //h_mCounts_Ratio[i_pt-2]->Divide(h_mCounts_RC[i_pt-2]);

  }

  TF1 *Rcv2[4];
  TF1 *Mcv2[4];

  double pt[4] = {1.5,2.1,2.7,3.6};

  TCanvas *c = new TCanvas("c","c",600,10,800,800);
  c->Divide(2,2);
  for(int i = 0; i < 4; i++)
  {
    c->cd(i+1);
    c->cd(i+1)->SetLeftMargin(0.15);
    c->cd(i+1)->SetBottomMargin(0.15);
    c->cd(i+1)->SetTicks(1,1);
    c->cd(i+1)->SetGrid(0,0);
    
    Rcv2[i] = new TF1(Form("Rcv2_%d",i),flow,-TMath::Pi()/2.0,TMath::Pi()/2.0,2);
    h_mCounts_RC[i]->Fit(Rcv2[i],"MRI");    
    g_mRcv2->SetPoint(i,pt[i],Rcv2[i]->GetParameter(0));
    g_mRcv2->SetPointError(i,0.0,0.0,Rcv2[i]->GetParError(0),Rcv2[i]->GetParError(0));

    h_mCounts_RC[i]->SetTitle(Form("RC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    h_mCounts_RC[i]->GetXaxis()->SetTitle("#phi-#Psi_{2}");
    h_mCounts_RC[i]->GetYaxis()->SetTitle("RC Counts");
    h_mCounts_RC[i]->Draw("pE");
  }

  c->SaveAs(Form("figures/%s/%s/pTstudy/PhiPsi_Rc_%s_Order%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectra.c_str()));  


  for(int i = 0; i < 4; i++)
  {
    c->cd(i+1);
    c->cd(i+1)->SetLeftMargin(0.15);
    c->cd(i+1)->SetBottomMargin(0.15);
    c->cd(i+1)->SetTicks(1,1);
    c->cd(i+1)->SetGrid(0,0);

    Mcv2[i] = new TF1(Form("Mcv2_%d",i),flow,-TMath::Pi()/2.0,TMath::Pi()/2.0,2);
    h_mCounts_MC[i]->Fit(Mcv2[i],"MRI");    
    g_mMcv2->SetPoint(i,pt[i],Mcv2[i]->GetParameter(0));
    g_mMcv2->SetPointError(i,0.0,0.0,Mcv2[i]->GetParError(0),Mcv2[i]->GetParError(0));
    
    h_mCounts_MC[i]->SetTitle(Form("MC %1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
    h_mCounts_MC[i]->GetXaxis()->SetTitle("#phi-#Psi_{2}");
    h_mCounts_MC[i]->GetYaxis()->SetTitle("MC Counts");
    h_mCounts_MC[i]->Draw("pE");
  }
  c->SaveAs(Form("figures/%s/%s/pTstudy/PhiPsi_Mc_%s_Order%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectra.c_str()));  

  
  TCanvas *c2 = new TCanvas("c2","c2",600,10,400,400);
  c2->cd();
  c2->cd()->SetLeftMargin(0.15);
  c2->cd()->SetBottomMargin(0.15);
  c2->cd()->SetTicks(1,1);
  c2->cd()->SetGrid(0,0);

  TLegend *leg = new TLegend(0.7,0.2,0.9,0.4);
  
  g_mRcv2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_mRcv2->GetYaxis()->SetTitle("Observed v_{2}");
  g_mRcv2->SetMarkerStyle(20);
  g_mRcv2->SetMarkerColor(kOrange+7);

  g_mMcv2->SetMarkerStyle(20);
  g_mMcv2->SetMarkerColor(kBlack);

  g_mRcv2->Draw("APE");
  g_mMcv2->Draw("PE same");

  leg->AddEntry(g_mMcv2,"MC","p");
  leg->AddEntry(g_mRcv2,"RC (with cuts)","p");
  leg->Draw("same");  

  c2->SaveAs(Form("figures/%s/%s/pTstudy/PhiPsi_v2comparison_%s_Order%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectra.c_str()));  

 // TLegend *leg = new TLegend(0.4,0.7,0.6,0.9);
 // for(int i = 0; i < 4; i++)
 // {
 //   c->cd(i+1);
 //   c->cd(i+1)->SetLeftMargin(0.15);
 //   c->cd(i+1)->SetBottomMargin(0.15);
 //   c->cd(i+1)->SetTicks(1,1);
 //   c->cd(i+1)->SetGrid(0,0);
 //   
 //   h_mCounts[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
 //   h_mCounts[i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
 //   h_mCounts[i]->GetYaxis()->SetTitle("Counts");
 //   h_mCounts[i]->SetMarkerStyle(20);
 //   h_mCounts[i]->SetMarkerColor(kOrange+7);
 //   h_mCounts[i]->SetLineColor(kOrange+7);

 //   int min = h_mCounts[i]->GetMinimum();
 //   int max = h_mCounts[i]->GetMaximum();
 //   if(h_mCounts_RC[i]->GetMinimum() < min) min = h_mCounts_RC[i]->GetMinimum();
 //   if(h_mCounts_RC[i]->GetMaximum() > max) max = h_mCounts_RC[i]->GetMaximum();
 //   h_mCounts[i]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

 //   h_mCounts[i]->Draw("pE");

 //   h_mCounts_RC[i]->SetMarkerStyle(24);
 //   h_mCounts_RC[i]->SetMarkerColor(kBlack);
 //   h_mCounts_RC[i]->SetLineColor(kBlack);
 //   h_mCounts_RC[i]->Draw("pE same");

 //   if(i == 0)
 //   {
 //     leg->AddEntry(h_mCounts[i],"Data","p");
 //     leg->AddEntry(h_mCounts_RC[i],"MC","p");
 //   }
 //   leg->Draw("same");
 // }

 // c->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_Comparison_%s_Order%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectra.c_str()));  


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


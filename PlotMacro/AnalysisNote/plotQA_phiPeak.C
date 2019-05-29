#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/draw.h"

void plotQA_phiPeak(int energy = 6)
{
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  const int ptBin = 3;
  const int thetaBin = 3;

  string InPutFile_SE = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Phi/Yields/Yields_SE_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_SE = TFile::Open(InPutFile_SE.c_str());
  string KEY_InPutSE = Form("pt_%d_Centrality_9_CosThetaStar_%d_2nd_Dca_0_Sig_0_%s_SE",ptBin,thetaBin,vmsa::mPID[0].c_str());
  TH1F *h_mInPut_SE = (TH1F*)File_SE->Get(KEY_InPutSE.c_str())->Clone(); 
  string KEY_SE = Form("pt_%d_Centrality_9_CosThetaStar_%d_2nd_Dca_0_Sig_0_%s_SE_QA",ptBin,thetaBin,vmsa::mPID[0].c_str());
  TH1F *h_mMass_SE = (TH1F*)h_mInPut_SE->Clone(KEY_SE.c_str());
  string KEY_SM = Form("pt_%d_Centrality_9_CosThetaStar_%d_2nd_Dca_0_Sig_0_%s_SM_QA",ptBin,thetaBin,vmsa::mPID[0].c_str());
  TH1F *h_mMass_SM = (TH1F*)h_mMass_SE->Clone(KEY_SM.c_str());

  string InPutFile_ME = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Phi/Yields/Yields_ME_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_ME = TFile::Open(InPutFile_ME.c_str());
  string KEY_InPutME = Form("pt_%d_Centrality_9_CosThetaStar_%d_2nd_Dca_0_Sig_0_%s_ME",ptBin,thetaBin,vmsa::mPID[0].c_str());
  TH1F *h_mInPut_ME = (TH1F*)File_ME->Get(KEY_InPutME.c_str())->Clone(); 
  string KEY_ME = Form("pt_%d_Centrality_9_CosThetaStar_%d_2nd_Dca_0_Sig_0_%s_ME_QA",ptBin,thetaBin,vmsa::mPID[0].c_str());
  TH1F *h_mMass_ME = (TH1F*)h_mInPut_ME->Clone(KEY_ME.c_str());


  int Norm_bin_start = h_mMass_SE->FindBin(vmsa::Norm_Start[0][0]);
  int Norm_bin_stop  = h_mMass_SE->FindBin(vmsa::Norm_Stop[0][0]);

  float Inte_SE = h_mMass_SE->Integral(Norm_bin_start,Norm_bin_stop);
  float Inte_ME = h_mMass_ME->Integral(Norm_bin_start,Norm_bin_stop);

  h_mMass_ME->Scale(Inte_SE/Inte_ME);
  h_mMass_SM->Add(h_mMass_ME,-1.0);

  string leg_energy = Form("Au+Au %s & 20-60%%",vmsa::mBeamEnergy[energy].c_str());
  TCanvas *c_peak = new TCanvas("c_peak","c_peak",10,10,800,800);
  c_peak->cd();
  c_peak->cd()->SetLeftMargin(0.15);
  c_peak->cd()->SetBottomMargin(0.15);
  c_peak->cd()->SetTicks(1,1);
  c_peak->cd()->SetGrid(0,0);
  h_mMass_SE->SetTitle(leg_energy.c_str());
  h_mMass_SE->SetStats(0);
  h_mMass_SE->GetXaxis()->SetTitle("M(K^{+},K^{-})(GeV/c^{2})");
  h_mMass_SE->GetXaxis()->CenterTitle();
  h_mMass_SE->GetXaxis()->SetTitleSize(0.06);
  h_mMass_SE->GetXaxis()->SetTitleOffset(0.9);
  h_mMass_SE->GetXaxis()->SetLabelSize(0.04);
  h_mMass_SE->SetNdivisions(505,"X");

  h_mMass_SE->GetYaxis()->SetTitle("Yields");
  h_mMass_SE->GetYaxis()->CenterTitle();
  h_mMass_SE->GetYaxis()->SetTitleOffset(0.97);
  h_mMass_SE->GetYaxis()->SetTitleSize(0.06);
  h_mMass_SE->GetYaxis()->SetLabelSize(0.04);
  h_mMass_SE->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE->GetMaximum(),1.1*h_mMass_SE->GetMaximum());
  h_mMass_SE->SetNdivisions(505,"Y");
  h_mMass_SE->SetMarkerStyle(24);
  h_mMass_SE->SetMarkerSize(1.4);
  h_mMass_SE->DrawCopy("PE");

  h_mMass_ME->SetLineColor(2);
  h_mMass_ME->SetFillColor(2);
  h_mMass_ME->SetFillStyle(3002);
  h_mMass_ME->DrawCopy("h same");

  h_mMass_SM->SetLineColor(4);
  h_mMass_SM->SetFillColor(4);
  h_mMass_SM->SetFillStyle(3004);
  h_mMass_SM->DrawCopy("h same");

  string leg_pt = Form("%1.1f < p_{T} < %1.1f GeV/c",vmsa::pt_low[energy][ptBin],vmsa::pt_up[energy][ptBin]);
  plotTopLegend((char*)leg_pt.c_str(),0.2,0.8,0.04,1,0.0,42,1);
  string leg_cos = Form("%d/7 < cos(#theta*) < %d/7",thetaBin-1,thetaBin);
  plotTopLegend((char*)leg_cos.c_str(),0.22,0.75,0.04,1,0.0,42,1);

  TLegend *leg = new TLegend(0.5,0.4,0.8,0.55);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(h_mMass_SE,"Same Event (Sig+Bg)","p");
  leg->AddEntry(h_mMass_ME,"Mixed Event (Bg)","f");
  leg->AddEntry(h_mMass_SM,"Sig","f");
  leg->Draw("same");

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/RhoSig/c_peak_AuAu%s_Pt%d_Theta%d.eps",vmsa::mBeamEnergy[energy].c_str(),ptBin,thetaBin);
  c_peak->SaveAs(FigName.c_str());
}

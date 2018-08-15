#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
// #include "../../Utility/functions.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/type.h"
#include "../StRoot/StToFMatchMaker/StToFMatchCons.h"
#include "TGraphAsymmErrors.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

TH1D* CalEffError(TH1D *h_TPC, TH1D *h_ToF, std::string HistName);

void plotQA_Tracks(int energy = 6)
{
  int pidQA = tof::pidQA;
  int chargeQA = tof::chargeQA;
  // int centQA = tof::centQA;
  int centQA = 9;
  // int etaQA = tof::etaQA;
  // int phiQA = tof::phiQA;
  int etaQA = 10;
  int phiQA = 12;

  string inputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/file_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH3D *h_mTracks_TPC[2][3][10]; // pt, eta, phi distribution as a function of charge | pid | centrality
  TH3D *h_mTracks_ToF[2][3][10];
  TH1DMap h_mCounts_TPC; // counts for TPC tracks
  TH1DMap h_mCounts_ToF; // counts for ToF tracks

  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName_Trakcs_TPC = Form("h_mTracks_TPC_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mTracks_TPC[i_charge][i_pid][i_cent] = (TH3D*)File_InPut->Get(HistName_Trakcs_TPC.c_str());

	string HistName_Trakcs_ToF = Form("h_mTracks_ToF_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mTracks_ToF[i_charge][i_pid][i_cent] = (TH3D*)File_InPut->Get(HistName_Trakcs_ToF.c_str());

	for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
	{
	  for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
	  {
	    string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mCounts_TPC[HistName_TPC] = (TH1D*)File_InPut->Get(HistName_TPC.c_str());
	    string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mCounts_ToF[HistName_ToF] = (TH1D*)File_InPut->Get(HistName_ToF.c_str());
	  }
	}
      }
    }
  }

  TCanvas *c_Tracks = new TCanvas("c_Tracks","c_Tracks",10,10,1200,600);
  c_Tracks->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_Tracks->cd(i_pad+1);
    c_Tracks->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Tracks->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Tracks->cd(i_pad+1)->SetTicks(1,1);
    c_Tracks->cd(i_pad+1)->SetGrid(0,0);
  }

  c_Tracks->cd(1);
  c_Tracks->cd(1)->SetLogy();
  string HistQA_Trakcs_TPC = Form("h_mTracks_TPC_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[pidQA].c_str(),tof::mCharge[chargeQA].c_str(),centQA,etaQA,phiQA);
  TH1D *h_QA_TPC = (TH1D*)h_mTracks_TPC[chargeQA][pidQA][centQA]->ProjectionX(HistQA_Trakcs_TPC.c_str(),etaQA+1,etaQA+1,phiQA+1,phiQA+1)->Clone();
  h_QA_TPC->SetLineColor(1);
  h_QA_TPC->DrawCopy("hE");

  string HistQA_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[pidQA].c_str(),tof::mCharge[chargeQA].c_str(),centQA,etaQA,phiQA);
  h_mCounts_TPC[HistQA_TPC]->SetMarkerStyle(24);
  h_mCounts_TPC[HistQA_TPC]->SetMarkerSize(1.2);
  h_mCounts_TPC[HistQA_TPC]->SetMarkerColor(2);
  h_mCounts_TPC[HistQA_TPC]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mCounts_TPC[HistQA_TPC]->GetXaxis()->CenterTitle();
  h_mCounts_TPC[HistQA_TPC]->GetYaxis()->SetTitle("Counts");
  h_mCounts_TPC[HistQA_TPC]->GetYaxis()->CenterTitle();
  h_mCounts_TPC[HistQA_TPC]->DrawCopy("pE same");

  c_Tracks->cd(2);
  c_Tracks->cd(2)->SetLogy();
  string HistQA_Trakcs_ToF = Form("h_mTracks_ToF_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[pidQA].c_str(),tof::mCharge[chargeQA].c_str(),centQA,etaQA,phiQA);
  TH1D *h_QA_ToF = (TH1D*)h_mTracks_ToF[chargeQA][pidQA][centQA]->ProjectionX(HistQA_Trakcs_ToF.c_str(),etaQA+1,etaQA+1,phiQA+1,phiQA+1)->Clone();
  h_QA_ToF->SetLineColor(1);
  h_QA_ToF->DrawCopy("hE");

  string HistQA_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[pidQA].c_str(),tof::mCharge[chargeQA].c_str(),centQA,etaQA,phiQA);
  h_mCounts_ToF[HistQA_ToF]->SetMarkerStyle(24);
  h_mCounts_ToF[HistQA_ToF]->SetMarkerSize(1.2);
  h_mCounts_ToF[HistQA_ToF]->SetMarkerColor(2);
  h_mCounts_ToF[HistQA_ToF]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mCounts_ToF[HistQA_ToF]->GetXaxis()->CenterTitle();
  h_mCounts_ToF[HistQA_ToF]->GetYaxis()->SetTitle("Counts");
  h_mCounts_ToF[HistQA_ToF]->GetYaxis()->CenterTitle();
  h_mCounts_ToF[HistQA_ToF]->DrawCopy("pE same");
  c_Tracks->SaveAs("../../figures/c_Tracks.eps");

  TCanvas *c_Efficiency = new TCanvas("c_Efficiency","c_Efficiency",10,10,1200,400);
  c_Efficiency->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_Efficiency->cd(i_pad+1);
    c_Efficiency->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Efficiency->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Efficiency->cd(i_pad+1)->SetTicks(1,1);
    c_Efficiency->cd(i_pad+1)->SetGrid(0,0);
  }

  c_Efficiency->cd(1);
  string HistEff_Tracks = Form("h_mEffTracks_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[pidQA].c_str(),tof::mCharge[chargeQA].c_str(),centQA,etaQA,phiQA);
  TH1D *h_mEff_Tracks = CalEffError(h_QA_TPC,h_QA_ToF,HistEff_Tracks.c_str());
  h_mEff_Tracks->SetTitle("CalEffError");
  h_mEff_Tracks->SetMarkerStyle(20);
  h_mEff_Tracks->SetMarkerSize(0.8);
  h_mEff_Tracks->SetMarkerColor(1);
  h_mEff_Tracks->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mEff_Tracks->GetXaxis()->CenterTitle();
  h_mEff_Tracks->GetYaxis()->SetTitle("Efficiency");
  h_mEff_Tracks->GetYaxis()->SetRangeUser(0.0,1.2);
  h_mEff_Tracks->GetYaxis()->CenterTitle();
  h_mEff_Tracks->Draw("PE");

  c_Efficiency->cd(2);
  string HistEff = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[pidQA].c_str(),tof::mCharge[chargeQA].c_str(),centQA,etaQA,phiQA);
  TH1D *h_mEfficiency = (TH1D*)h_mCounts_ToF[HistQA_ToF]->Clone(HistEff.c_str());
  h_mEfficiency->Reset();
  TH1D *h_ToF_B = (TH1D*)h_mCounts_ToF[HistQA_ToF]->Clone("h_ToF_B");
  TH1D *h_TPC_B = (TH1D*)h_mCounts_TPC[HistQA_TPC]->Clone("h_TPC_B");
  h_mEfficiency->Divide(h_ToF_B,h_TPC_B,1,1,"B");
  h_mEfficiency->SetTitle("TH1::Divide(\"B\")");
  h_mEfficiency->SetMarkerStyle(24);
  h_mEfficiency->SetMarkerSize(1.2);
  h_mEfficiency->SetMarkerColor(2);
  h_mEfficiency->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mEfficiency->GetXaxis()->CenterTitle();
  h_mEfficiency->GetYaxis()->SetTitle("Efficiency");
  h_mEfficiency->GetYaxis()->SetRangeUser(0.0,1.2);
  h_mEfficiency->Draw("pE");

  c_Efficiency->cd(3);
  TH1D *h_ToF = (TH1D*)h_mCounts_ToF[HistQA_ToF]->Clone("h_ToF");
  TH1D *h_TPC = (TH1D*)h_mCounts_TPC[HistQA_TPC]->Clone("h_TPC");
  TH1D *h_frame = (TH1D*)h_mEfficiency->Clone("h_frame");
  h_frame->Reset();
  h_frame->SetTitle("TGraphAsymmErrors::BayesDivide");
  h_frame->Draw("pE");
  TGraphAsymmErrors *g_mEfficiency = new TGraphAsymmErrors();
  g_mEfficiency->BayesDivide(h_ToF,h_TPC);
  g_mEfficiency->SetMarkerStyle(24);
  g_mEfficiency->SetMarkerSize(1.2);
  g_mEfficiency->SetMarkerColor(2);
  g_mEfficiency->Draw("pE same");

  c_Efficiency->SaveAs("../../figures/c_ToFMatchEff_method.eps");
}

TH1D* CalEffError(TH1D *h_TPC, TH1D *h_ToF, std::string HistName)
{
  TH1D* h_ratio = (TH1D*)h_ToF->Clone();
  h_ratio->Divide(h_TPC);
  for(int i_bin = 1; i_bin < h_ratio->GetNbinsX()+1; ++i_bin)
  {
    double n = h_TPC->GetBinContent(i_bin);
    double k = h_ToF->GetBinContent(i_bin);
    double variance = (k+1.0)*(k+2.0)/((n+2.0)*(n+3.0))-(k+1.0)*(k+1.0)/((n+2.0)*(n+2.0));
    double sigma = TMath::Sqrt(variance);
    if(n > 0.0 && k > 0.0) h_ratio->SetBinError(i_bin,sigma);
  }
  h_ratio->SetName(HistName.c_str());

  return h_ratio;
}

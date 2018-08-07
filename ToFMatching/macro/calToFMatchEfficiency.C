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

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

TH1D* CalEffError(TH1D *h_TPC, TH1D *h_ToF, std::string HistName);

void calToFMatchEfficiency(int energy = 6)
{
  string inputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/file_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH3D *h_mTracks_TPC[2][3][9]; // pt, eta, phi distribution as a function of charge | pid | centrality
  TH3D *h_mTracks_ToF[2][3][9];
  TH1DMap h_mCounts_TPC; // counts for TPC tracks
  TH1DMap h_mCounts_ToF; // counts for ToF tracks

  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 9; ++i_cent)
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

#if _PlotQA_
  TCanvas *c_Tracks = new TCanvas("c_Tracks","c_Tracks",10,10,1200,400);
  c_Tracks->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_Tracks->cd(i_pad+1);
    c_Tracks->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Tracks->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Tracks->cd(i_pad+1)->SetTicks(1,1);
    c_Tracks->cd(i_pad+1)->SetGrid(0,0);
  }

  c_Tracks->cd(1);
  c_Tracks->cd(1)->SetLogy();
  string QA_Pt_TPC = Form("h_mPt_TPC_%s%s_Cent_%d",tof::mPID_ToF[0].c_str(),tof::mCharge[0].c_str(),tof::centQA);
  TH1D *h_QApt_TPC = (TH1D*)h_mTracks_TPC[tof::chargeQA][tof::pidQA][tof::centQA]->Project3D("x")->Clone(QA_Pt_TPC.c_str());
  h_QApt_TPC->SetMarkerStyle(20);
  h_QApt_TPC->SetMarkerSize(1.2);
  h_QApt_TPC->SetMarkerColor(kGray+2);
  h_QApt_TPC->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_QApt_TPC->GetXaxis()->CenterTitle();
  h_QApt_TPC->GetYaxis()->SetTitle("Counts");
  h_QApt_TPC->GetYaxis()->CenterTitle();
  h_QApt_TPC->DrawCopy("pE");

  string QA_Pt_ToF = Form("h_mPt_ToF_%s%s_Cent_%d",tof::mPID_ToF[0].c_str(),tof::mCharge[0].c_str(),tof::centQA);
  TH1D *h_QApt_ToF = (TH1D*)h_mTracks_ToF[tof::chargeQA][tof::pidQA][tof::centQA]->Project3D("x")->Clone(QA_Pt_ToF.c_str());
  h_QApt_ToF->SetMarkerStyle(24);
  h_QApt_ToF->SetMarkerSize(1.2);
  h_QApt_ToF->SetMarkerColor(2);
  h_QApt_ToF->DrawCopy("pE same");

  c_Tracks->cd(2);
  string QA_Eta_TPC = Form("h_mEta_TPC_%s%s_Cent_%d",tof::mPID_ToF[0].c_str(),tof::mCharge[0].c_str(),tof::centQA);
  TH1D *h_QAeta_TPC = (TH1D*)h_mTracks_TPC[tof::chargeQA][tof::pidQA][tof::centQA]->Project3D("y")->Clone(QA_Eta_TPC.c_str());
  h_QAeta_TPC->SetMarkerStyle(20);
  h_QAeta_TPC->SetMarkerSize(1.2);
  h_QAeta_TPC->SetMarkerColor(kGray+2);
  h_QAeta_TPC->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_QAeta_TPC->GetXaxis()->CenterTitle();
  h_QAeta_TPC->GetYaxis()->SetTitle("Counts");
  h_QAeta_TPC->GetYaxis()->CenterTitle();
  h_QAeta_TPC->DrawCopy("pE");

  string QA_Eta_ToF = Form("h_mEta_ToF_%s%s_Cent_%d",tof::mPID_ToF[0].c_str(),tof::mCharge[0].c_str(),tof::centQA);
  TH1D *h_QAeta_ToF = (TH1D*)h_mTracks_ToF[tof::chargeQA][tof::pidQA][tof::centQA]->Project3D("x")->Clone(QA_Eta_ToF.c_str());
  h_QAeta_ToF->SetMarkerStyle(24);
  h_QAeta_ToF->SetMarkerSize(1.2);
  h_QAeta_ToF->SetMarkerColor(2);
  h_QAeta_ToF->DrawCopy("pE same");

  c_Tracks->cd(3);
  string QA_Phi_TPC = Form("h_mPhi_TPC_%s%s_Cent_%d",tof::mPID_ToF[0].c_str(),tof::mCharge[0].c_str(),tof::centQA);
  TH1D *h_QAphi_TPC = (TH1D*)h_mTracks_TPC[tof::chargeQA][tof::pidQA][tof::centQA]->Project3D("y")->Clone(QA_Phi_TPC.c_str());
  h_QAphi_TPC->SetMarkerStyle(20);
  h_QAphi_TPC->SetMarkerSize(1.2);
  h_QAphi_TPC->SetMarkerColor(kGray+2);
  h_QAphi_TPC->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_QAphi_TPC->GetXaxis()->CenterTitle();
  h_QAphi_TPC->GetYaxis()->SetTitle("Counts");
  h_QAphi_TPC->GetYaxis()->CenterTitle();
  h_QAphi_TPC->DrawCopy("pE");

  string QA_Phi_ToF = Form("h_mEta_ToF_%s%s_Cent_%d",tof::mPID_ToF[0].c_str(),tof::mCharge[0].c_str(),tof::centQA);
  TH1D *h_QAphi_ToF = (TH1D*)h_mTracks_ToF[tof::chargeQA][tof::pidQA][tof::centQA]->Project3D("x")->Clone(QA_Phi_ToF.c_str());
  h_QAphi_ToF->SetMarkerStyle(24);
  h_QAphi_ToF->SetMarkerSize(1.2);
  h_QAphi_ToF->SetMarkerColor(2);
  h_QAphi_ToF->DrawCopy("pE same");
#endif
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

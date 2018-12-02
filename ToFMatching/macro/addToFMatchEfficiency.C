#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/type.h"
#include "../StRoot/StToFMatchMaker/StToFMatchCons.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

void addToFMatchEfficiency(int energy = 6)
{
  string inputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/file_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH3D *h_mTracks_TPC[2][3][10]; // pt, eta, phi distribution as a function of charge | pid | centrality
  TH3D *h_mTracks_ToF[2][3][10];
  TH1DMap h_mCounts_TPC; // counts for TPC tracks
  TH1DMap h_mCounts_ToF; // counts for ToF tracks

  TH1D *h_FrameEta_ToF = (TH1D*)File_InPut->Get("h_FrameEta_ToF");
  TH1D *h_FramePhi_ToF = (TH1D*)File_InPut->Get("h_FramePhi_ToF");

  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName_Trakcs_TPC = Form("h_mTracks_TPC_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mTracks_TPC[i_charge][i_pid][i_cent] = (TH3D*)File_InPut->Get(HistName_Trakcs_TPC.c_str());

	string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_TPC[HistName_TPC] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionX(HistName_TPC.c_str())->Clone();

	string HistName_Trakcs_ToF = Form("h_mTracks_ToF_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mTracks_ToF[i_charge][i_pid][i_cent] = (TH3D*)File_InPut->Get(HistName_Trakcs_ToF.c_str());

	string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_ToF[HistName_ToF] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionX(HistName_ToF.c_str())->Clone();
      }
    }
  }

  //---------------efficiency calculation------------------------
  TH1DMap h_mEfficiency; // ToF matching efficiency
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);

	string HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency[HistName] = (TH1D*)h_mCounts_ToF[HistName_ToF]->Clone(HistName.c_str());
	h_mEfficiency[HistName]->SetTitle(HistName.c_str());
	h_mEfficiency[HistName]->Reset();
	h_mEfficiency[HistName]->Divide(h_mCounts_ToF[HistName_ToF],h_mCounts_TPC[HistName_TPC],1,1,"B");
      }
    }
  }

#if _PlotQA_
  int pidQA = tof::pidQA;
  int chargeQA = tof::chargeQA;
  int centQA = 9;

  TCanvas *c_Efficiency = new TCanvas("c_Efficiency","c_Efficiency",10,10,800,800);
  c_Efficiency->cd();
  c_Efficiency->cd()->SetLeftMargin(0.15);
  c_Efficiency->cd()->SetBottomMargin(0.15);
  c_Efficiency->cd()->SetTicks(1,1);
  c_Efficiency->cd()->SetGrid(0,0);

  string HistEff = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[pidQA].c_str(),tof::mCharge[chargeQA].c_str(),centQA);
  // h_mEfficiency[HistEff]->SetTitle("TH1::Divide(\"B\")");
  h_mEfficiency[HistEff]->SetMarkerStyle(24);
  h_mEfficiency[HistEff]->SetMarkerSize(1.2);
  h_mEfficiency[HistEff]->SetMarkerColor(2);
  h_mEfficiency[HistEff]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mEfficiency[HistEff]->GetXaxis()->CenterTitle();
  h_mEfficiency[HistEff]->GetYaxis()->SetTitle("Efficiency");
  h_mEfficiency[HistEff]->GetYaxis()->SetRangeUser(0.0,1.2);
  h_mEfficiency[HistEff]->Draw("pE");
#endif 

  //---------------output------------------------
  string outputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/centEff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency[HistName]->Write();
      }
    }
  }
  h_FrameEta_ToF->Write();
  h_FramePhi_ToF->Write();
  File_OutPut->Close();
}

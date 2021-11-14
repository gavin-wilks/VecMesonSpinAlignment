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
#include "../../ToFMatching/StRoot/StToFMatchMaker/StToFMatchCons.h"

void plotQA_TofEta(int energy = 1)
{
  // string inputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/file_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/file_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  // calculate efficiency vs centrality
  TH3D *h_mTracks_TPC[2][3][10]; // pt, eta, phi distribution as a function of charge | pid | centrality
  TH3D *h_mTracks_ToF[2][3][10];
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
      }
    }
  }
  // calculate efficiencya vs. eta
  TH1DMap h_mCounts_TPC_eta; // efficiency vs. cent & eta
  TH1DMap h_mCounts_ToF_eta;
  TH1DMap h_mEfficiency_eta;
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_Eta",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_TPC_eta[HistName_TPC] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionY(HistName_TPC.c_str())->Clone();

	string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_Eta",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_ToF_eta[HistName_ToF] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionY(HistName_ToF.c_str())->Clone();

	string HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_eta[HistName] = (TH1D*)h_mCounts_ToF_eta[HistName_ToF]->Clone(HistName.c_str());
	h_mEfficiency_eta[HistName]->SetTitle(HistName.c_str());
	h_mEfficiency_eta[HistName]->Reset();
	h_mEfficiency_eta[HistName]->Divide(h_mCounts_ToF_eta[HistName_ToF],h_mCounts_TPC_eta[HistName_TPC],1,1,"B");
      }
    }
  }


  string HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta",tof::mPID_ToF[0].c_str(),tof::mCharge[0].c_str(),2);
  h_mEfficiency_eta[HistName]->Draw();
}

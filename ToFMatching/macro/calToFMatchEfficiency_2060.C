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
#include "../../Utility/functions.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
//#include "../../Utility/type.h"
#include "../StRoot/StToFMatchMaker/StToFMatchCons.h"


typedef std::map<std::string,TH1D*> TH1DMap;

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

void calToFMatchEfficiency_2060(int energy = 4)
{
  // string inputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/file_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  std::string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/file_%s_ToFMatching.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH1D *h_FrameEta_ToF = (TH1D*)File_InPut->Get("h_FrameEta_ToF");
  TH1D *h_FramePhi_ToF = (TH1D*)File_InPut->Get("h_FramePhi_ToF");

  TH3D *h_mTracks_TPC[2][3][10]; // pt, eta, phi distribution as a function of charge | pid | centrality
  TH3D *h_mTracks_ToF[2][3][10];
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
      }
      for(int i_cent = 2; i_cent <= 5; ++i_cent) // make centrality 9 to be 20-60%
      {
	if(i_cent == 2)
	{
	  string HistName_Trakcs_TPC = Form("h_mTracks_TPC_%s%s_Cent_9",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str());
	  h_mTracks_TPC[i_charge][i_pid][9] = (TH3D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->Clone(HistName_Trakcs_TPC.c_str());
	  h_mTracks_TPC[i_charge][i_pid][9]->Sumw2();

	  string HistName_Trakcs_ToF = Form("h_mTracks_ToF_%s%s_Cent_9",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str());
	  h_mTracks_ToF[i_charge][i_pid][9] = (TH3D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->Clone(HistName_Trakcs_ToF.c_str());
	  h_mTracks_ToF[i_charge][i_pid][9]->Sumw2();
	}
	else
	{
	  h_mTracks_TPC[i_charge][i_pid][9]->Add(h_mTracks_TPC[i_charge][i_pid][i_cent]);
	  h_mTracks_ToF[i_charge][i_pid][9]->Add(h_mTracks_ToF[i_charge][i_pid][i_cent]);
	}
      }
    }
  }


  // calculate efficiency vs centrality
  TH1DMap h_mCounts_TPC_cent; // efficiency vs. cent
  TH1DMap h_mCounts_ToF_cent;
  TH1DMap h_mEfficiency_cent;
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_TPC_cent[HistName_TPC] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionX(HistName_TPC.c_str())->Clone();

	string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mCounts_ToF_cent[HistName_ToF] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionX(HistName_ToF.c_str())->Clone();

	string HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_cent[HistName] = (TH1D*)h_mCounts_ToF_cent[HistName_ToF]->Clone(HistName.c_str());
	h_mEfficiency_cent[HistName]->SetTitle(HistName.c_str());
	h_mEfficiency_cent[HistName]->Reset();
	h_mEfficiency_cent[HistName]->Divide(h_mCounts_ToF_cent[HistName_ToF],h_mCounts_TPC_cent[HistName_TPC],1,1,"B");
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
	for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
	{
	  string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta);
	  h_mCounts_TPC_eta[HistName_TPC] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionX(HistName_TPC.c_str(),i_eta+1,i_eta+1)->Clone();

	  string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta);
	  h_mCounts_ToF_eta[HistName_ToF] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionX(HistName_ToF.c_str(),i_eta+1,i_eta+1)->Clone();

	  string HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta);
	  h_mEfficiency_eta[HistName] = (TH1D*)h_mCounts_ToF_eta[HistName_ToF]->Clone(HistName.c_str());
	  h_mEfficiency_eta[HistName]->SetTitle(HistName.c_str());
	  h_mEfficiency_eta[HistName]->Reset();
	  h_mEfficiency_eta[HistName]->Divide(h_mCounts_ToF_eta[HistName_ToF],h_mCounts_TPC_eta[HistName_TPC],1,1,"B");
	}
      }
    }
  }

  TH1DMap h_mCounts_TPC; // efficiency vs. cent & eta & phi
  TH1DMap h_mCounts_ToF;
  TH1DMap h_mEfficiency;
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
	{
	  for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
	  {
	    string HistName_TPC = Form("h_mCounts_TPC_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mCounts_TPC[HistName_TPC] = (TH1D*)h_mTracks_TPC[i_charge][i_pid][i_cent]->ProjectionX(HistName_TPC.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1)->Clone();

	    string HistName_ToF = Form("h_mCounts_ToF_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mCounts_ToF[HistName_ToF] = (TH1D*)h_mTracks_ToF[i_charge][i_pid][i_cent]->ProjectionX(HistName_ToF.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1)->Clone();

	    string HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mEfficiency[HistName] = (TH1D*)h_mCounts_ToF[HistName_ToF]->Clone(HistName.c_str());
	    h_mEfficiency[HistName]->SetTitle(HistName.c_str());
	    h_mEfficiency[HistName]->Reset();
	    h_mEfficiency[HistName]->Divide(h_mCounts_ToF[HistName_ToF],h_mCounts_TPC[HistName_TPC],1,1,"B");
	    cout << "Calculating => " << HistName.c_str() << endl;
	  }
	}
      }
    }
  }


  //---------------output------------------------
  // string outputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string outputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/Data/Efficiency/ToF/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_cent[HistName]->Write();
	for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
	{
	  HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta);
	  h_mEfficiency_eta[HistName]->Write();
	  for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
	  {
	    HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mEfficiency[HistName]->Write();
	  }
	}
      }
    }
  }
  h_FrameEta_ToF->Write();
  h_FramePhi_ToF->Write();
  File_OutPut->Close();
}

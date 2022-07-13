#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include <TRandom3.h>
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/draw.h"
#include "../../Utility/type.h"
#include "../StRoot/StToFMatchMaker/StToFMatchCons.h"

using namespace std;

void comTofMatchEfficiency(int energy = 4, int pid = 0, int charge = 0)
{
  gStyle->SetOptDate(0);
  gRandom->SetSeed();

  int NCentralityMax = 9;
  int NEtaMax = 10;
  int NPhiMax = 12;

  /*string input_0080 = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut_0080 = TFile::Open(input_0080.c_str());

  TH1DMap h_mEfficiency_0080;
  File_InPut_0080->cd();
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_0080[HistName] = (TH1D*)File_InPut_0080->Get(HistName.c_str());
	for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
	{
	  HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta);
	  h_mEfficiency_0080[HistName] = (TH1D*)File_InPut_0080->Get(HistName.c_str());
	  for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
	  {
	    HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mEfficiency_0080[HistName] = (TH1D*)File_InPut_0080->Get(HistName.c_str());
	  }
	}
      }
    }
  }*/


  string input_2060 = Form("../../output/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut_2060 = TFile::Open(input_2060.c_str());

  TH1DMap h_mEfficiency_2060;
  File_InPut_2060->cd();
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = 0; i_cent < 10; ++i_cent)
      {
	string HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_2060[HistName] = (TH1D*)File_InPut_2060->Get(HistName.c_str());
	for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
	{
	  HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta);
	  h_mEfficiency_2060[HistName] = (TH1D*)File_InPut_2060->Get(HistName.c_str());
	  for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
	  {
	    HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mEfficiency_2060[HistName] = (TH1D*)File_InPut_2060->Get(HistName.c_str());
	  }
	}
      }
    }
  }

  // int const cent_bin = floor(NCentralityMax * gRandom->Rndm());
  int const cent_bin = 9;
  int const eta_bin = floor(NEtaMax * gRandom->Rndm());
  cout << "cent_bin = " << cent_bin << ", eta_bin = " << eta_bin << endl;

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1600,1200);
  c_play->Divide(4,3);
  for(int i_pad = 0; i_pad < 12; ++i_pad)
  {
    c_play->cd(i_pad+1)->SetLeftMargin(0.15);
    c_play->cd(i_pad+1)->SetBottomMargin(0.15);
    c_play->cd(i_pad+1)->SetGrid(0,0);
    c_play->cd(i_pad+1)->SetTicks(1,1);
    string HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[pid].c_str(),tof::mCharge[charge].c_str(),cent_bin,eta_bin,i_pad);

    h_mEfficiency_2060[HistName]->SetMarkerStyle(24);
    h_mEfficiency_2060[HistName]->SetMarkerSize(1.2);
    h_mEfficiency_2060[HistName]->SetMarkerColor(2);
    h_mEfficiency_2060[HistName]->SetLineColor(2);
    h_mEfficiency_2060[HistName]->Draw("pE");

   // h_mEfficiency_0080[HistName]->SetMarkerStyle(20);
    //h_mEfficiency_0080[HistName]->SetMarkerSize(1.0);
    //h_mEfficiency_0080[HistName]->SetMarkerColor(1);
    //h_mEfficiency_0080[HistName]->SetLineColor(1);
    //h_mEfficiency_0080[HistName]->Draw("pE same");
  }
}

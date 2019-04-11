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

using namespace std;

void comTpcTrackEfficiency(int energy = 6, int pid = 0, int year = 0)
{
  gStyle->SetOptDate(0);
  gRandom->SetSeed();

  int NCentralityMax = 9;
  int NEtaMax = 10;
  int NPhiMax = 12;

  string innput_0080 = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str());
  TFile *File_InPut_0080 = TFile::Open(innput_0080.c_str());

  TH1DMap h_mEfficiency_0080;
  File_InPut_0080->cd();
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistNameEff;

    HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEfficiency_0080[HistNameEff] = (TH1D*)File_InPut_0080->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;

    HistNameEff = Form("h_mEffEta_Cent_%d",i_cent);
    h_mEfficiency_0080[HistNameEff] = (TH1D*)File_InPut_0080->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;

    HistNameEff = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEfficiency_0080[HistNameEff] = (TH1D*)File_InPut_0080->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;
    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mEfficiency_0080[HistNameEff] = (TH1D*)File_InPut_0080->Get(HistNameEff.c_str());
	cout << "read in => " << HistNameEff.c_str() << endl;
      }
    }
  }

  string innput_2060 = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mParType[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mYear[year].c_str());
  TFile *File_InPut_2060 = TFile::Open(innput_2060.c_str());

  TH1DMap h_mEfficiency_2060;
  File_InPut_2060->cd();
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistNameEff;

    HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEfficiency_2060[HistNameEff] = (TH1D*)File_InPut_2060->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;

    HistNameEff = Form("h_mEffEta_Cent_%d",i_cent);
    h_mEfficiency_2060[HistNameEff] = (TH1D*)File_InPut_2060->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;

    HistNameEff = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEfficiency_2060[HistNameEff] = (TH1D*)File_InPut_2060->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;
    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
	HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mEfficiency_2060[HistNameEff] = (TH1D*)File_InPut_2060->Get(HistNameEff.c_str());
	cout << "read in => " << HistNameEff.c_str() << endl;
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
    string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",cent_bin,eta_bin,i_pad);

    h_mEfficiency_2060[HistNameEff]->SetMarkerStyle(24);
    h_mEfficiency_2060[HistNameEff]->SetMarkerSize(1.2);
    h_mEfficiency_2060[HistNameEff]->SetMarkerColor(2);
    h_mEfficiency_2060[HistNameEff]->SetLineColor(2);
    h_mEfficiency_2060[HistNameEff]->Draw("pE");

    h_mEfficiency_0080[HistNameEff]->SetMarkerStyle(20);
    h_mEfficiency_0080[HistNameEff]->SetMarkerSize(1.0);
    h_mEfficiency_0080[HistNameEff]->SetMarkerColor(1);
    h_mEfficiency_0080[HistNameEff]->SetLineColor(1);
    h_mEfficiency_0080[HistNameEff]->Draw("pE same");
  }
}

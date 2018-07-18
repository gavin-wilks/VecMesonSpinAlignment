#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/draw.h"

using namespace std;

void plotKaonEff(const int pid = 0, const string surffix = "pT80_phi2")
{
  gStyle->SetOptDate(0);
  const int Centrality[3] = {9,0,1};
  const int Color[3] = {1,2,4};
  string InPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/SpinAlignment/Embedding/%s/Efficiency/Eff_200GeV_StMcEvent_run11_pr_%s.root",vmsa::mParType[pid].c_str(),surffix.c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TH1D *h_mEffPt[3];
  TH1D *h_mEffEta[3];
  TH1D *h_mEffPhi[3];

  for(int i_cent = 0; i_cent < 3; ++i_cent)
  {
    string HistPt = Form("h_mEffPt_Cent_%d",Centrality[i_cent]);
    h_mEffPt[i_cent] = (TH1D*)File_InPut->Get(HistPt.c_str());
    h_mEffPt[i_cent]->SetLineColor(Color[i_cent]);

    string HistEta = Form("h_mEffEta_Cent_%d",Centrality[i_cent]);
    h_mEffEta[i_cent] = (TH1D*)File_InPut->Get(HistEta.c_str());
    h_mEffEta[i_cent]->SetLineColor(Color[i_cent]);

    string HistPhi = Form("h_mEffPhi_Cent_%d",Centrality[i_cent]);
    h_mEffPhi[i_cent] = (TH1D*)File_InPut->Get(HistPhi.c_str());
    h_mEffPhi[i_cent]->SetLineColor(Color[i_cent]);
  }

  TCanvas *c_KaonEff = new TCanvas("c_KaonEff","c_KaonEff",10,10,1500,500);
  c_KaonEff->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_KaonEff->cd(i_pad+1)->SetLeftMargin(0.15);
    c_KaonEff->cd(i_pad+1)->SetBottomMargin(0.15);
    c_KaonEff->cd(i_pad+1)->SetGrid(0,0);
    c_KaonEff->cd(i_pad+1)->SetTicks(1,1);
  }

  c_KaonEff->cd(1);
  h_mEffPt[0]->GetYaxis()->SetRangeUser(0.0,1.2);
  h_mEffPt[0]->Draw("pE");
  h_mEffPt[1]->Draw("pE same");
  h_mEffPt[2]->Draw("pE same");

  c_KaonEff->cd(2);
  h_mEffEta[0]->GetYaxis()->SetRangeUser(0.0,1.2);
  h_mEffEta[0]->Draw("pE");
  h_mEffEta[1]->Draw("pE same");
  h_mEffEta[2]->Draw("pE same");

  c_KaonEff->cd(3);
  h_mEffPhi[0]->GetYaxis()->SetRangeUser(0.0,1.2);
  h_mEffPhi[0]->Draw("pE");
  h_mEffPhi[1]->Draw("pE same");
  h_mEffPhi[2]->Draw("pE same");
  string fig_KaonEff = Form("../figures/KaonEfficiencyQA/c_%sEff_%s.pdf",vmsa::mParType[pid].c_str(),surffix.c_str());
  c_KaonEff->SaveAs(fig_KaonEff.c_str());

  TH1D *h_EffDiff[3][10][12];
  TCanvas *c_EffDiff[10];
  for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
  {
    string canvas = Form("c_EffDiff_%d",i_eta);
    // c_EffDiff[i_eta] = new TCanvas(canvas.c_str(),canvas.c_str(),10,10,1200,600);
    // c_EffDiff[i_eta]->Divide(2,1);
    c_EffDiff[i_eta] = new TCanvas(canvas.c_str(),canvas.c_str(),10,10,1200,900);
    c_EffDiff[i_eta]->Divide(4,3);
    for(int i_phi = 0; i_phi < 12; ++i_phi)
    {
      c_EffDiff[i_eta]->cd(i_phi+1)->SetLeftMargin(0.15);
      c_EffDiff[i_eta]->cd(i_phi+1)->SetBottomMargin(0.15);
      c_EffDiff[i_eta]->cd(i_phi+1)->SetGrid(0,0);
      c_EffDiff[i_eta]->cd(i_phi+1)->SetTicks(1,1);
      for(int i_cent = 0; i_cent < 3; ++i_cent)
      {
	string HistDiff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",Centrality[i_cent],i_eta,i_phi);
	h_EffDiff[i_cent][i_eta][i_phi] = (TH1D*)File_InPut->Get(HistDiff.c_str());
	h_EffDiff[i_cent][i_eta][i_phi]->SetLineColor(Color[i_cent]);
      }
      h_EffDiff[0][i_eta][i_phi]->GetYaxis()->SetRangeUser(0.0,1.2);
      h_EffDiff[0][i_eta][i_phi]->Draw("pE");
      h_EffDiff[1][i_eta][i_phi]->Draw("pE same");
      h_EffDiff[2][i_eta][i_phi]->Draw("pE same");
    }
    string fig_KaonDiff = Form("../figures/KaonEfficiencyQA/c_%sDiff_Eta_%d_%s.pdf",vmsa::mParType[pid].c_str(),i_eta,surffix.c_str());
    c_EffDiff[i_eta]->SaveAs(fig_KaonDiff.c_str());
  }
}

#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/draw.h"

using namespace std;

void plotKaonEff(int energy = 6, int cent = 9)
{
  gStyle->SetOptDate(0);
  // string inputKplus = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/Kplus/Efficiency/Eff_%s_StMcEvent_run11_pr_2060_EP.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputKplus = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/Kplus/Efficiency/Eff_%s_StMcEvent_run11_pr_EP.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kplus = TFile::Open(inputKplus.c_str());
  TH1D *h_mEffKplus[10];
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent)
  {
    string HistName = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEffKplus[i_cent] = (TH1D*)File_Kplus->Get(HistName.c_str());
  }

  // string inputKminus = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/Kminus/Efficiency/Eff_%s_StMcEvent_run11_pr_2060_EP.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputKminus = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/Kminus/Efficiency/Eff_%s_StMcEvent_run11_pr_EP.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_Kminus = TFile::Open(inputKminus.c_str());
  TH1D *h_mEffKminus[10];
  for(int i_cent = vmsa::Cent_start; i_cent < vmsa::Cent_stop; ++i_cent)
  {
    string HistName = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEffKminus[i_cent] = (TH1D*)File_Kminus->Get(HistName.c_str());
  }

  TCanvas *c_eff = new TCanvas("c_eff","c_eff",10,10,800,800);
  c_eff->SetLeftMargin(0.15);
  c_eff->SetBottomMargin(0.15);
  c_eff->SetGrid(0,0);
  c_eff->SetTicks(1,1);
  TH1D *h_frame = new TH1D("h_frame","h_frame",100,-TMath::Pi(),TMath::Pi());
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10);
    h_frame->SetBinError(i_bin+1,1);
  }
  string legEnergy = Form("AuAu %s 20%%-60%%",vmsa::mBeamEnergy[energy].c_str());
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetTitle("#phi-#Psi_{2}");
  h_frame->GetXaxis()->CenterTitle();
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetNdivisions(505);
  h_frame->GetXaxis()->SetRangeUser(-0.5*TMath::Pi(),0.5*TMath::Pi());

  h_frame->GetYaxis()->SetTitle("efficiency");
  h_frame->GetYaxis()->SetTitleSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->SetNdivisions(505);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetRangeUser(0.76,0.79);
  h_frame->Draw("pE");

  h_mEffKplus[cent]->SetMarkerStyle(20);
  h_mEffKplus[cent]->SetMarkerColor(kGray+2);
  h_mEffKplus[cent]->SetMarkerSize(1.4);
  h_mEffKplus[cent]->SetLineColor(kGray+2);
  h_mEffKplus[cent]->Draw("pE same");

  h_mEffKminus[cent]->SetMarkerStyle(24);
  h_mEffKminus[cent]->SetMarkerColor(2);
  h_mEffKminus[cent]->SetMarkerSize(1.4);
  h_mEffKminus[cent]->SetLineColor(2);
  h_mEffKminus[cent]->Draw("pE same");

  TLegend *leg = new TLegend(0.2,0.7,0.5,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  // leg->SetNColumns(2);
  leg->AddEntry(h_mEffKplus[cent],"K^{+}","p");
  leg->AddEntry(h_mEffKminus[cent],"K^{-}","p");
  leg->Draw("same");

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/phiMeson/c_effKaonEPAuAu%s.eps",vmsa::mBeamEnergy[energy].c_str());
  c_eff->SaveAs(FigName.c_str());
}

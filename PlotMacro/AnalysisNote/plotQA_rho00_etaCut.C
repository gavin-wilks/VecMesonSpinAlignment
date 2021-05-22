#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include "../../Utility/draw.h"

void plotQA_rho00_etaCut()
{
  string input_eta1 = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu27GeV/Phi/rho00/eta_1/RawRhoPtSys.root";
  TFile *File_eta1 = TFile::Open(input_eta1.c_str());
  string GraphName = "rhoRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Count";
  TGraphAsymmErrors *g_rho_eta1 = (TGraphAsymmErrors*)File_eta1->Get(GraphName.c_str());
  TH1F *h_frame = (TH1F*)File_eta1->Get("h_frame");

  string input_eta08 = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu27GeV/Phi/rho00/eta_08/RawRhoPtSys.root";
  TFile *File_eta08 = TFile::Open(input_eta08.c_str());
  TGraphAsymmErrors *g_rho_eta08 = (TGraphAsymmErrors*)File_eta08->Get(GraphName.c_str());

  string input_eta06 = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu27GeV/Phi/rho00/eta_06/RawRhoPtSys.root";
  TFile *File_eta06 = TFile::Open(input_eta06.c_str());
  TGraphAsymmErrors *g_rho_eta06 = (TGraphAsymmErrors*)File_eta06->Get(GraphName.c_str());

  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  string leg_energy = "AuAu27GeV & 20-60%";
  h_frame->SetTitle(leg_energy.c_str());
  h_frame->GetYaxis()->SetRangeUser(0.25,0.4);
  h_frame->Draw("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  g_rho_eta1->SetMarkerStyle(20);
  g_rho_eta1->SetMarkerColor(kGray+2);
  g_rho_eta1->SetMarkerSize(1.2);
  g_rho_eta1->Draw("pE same");

  g_rho_eta08->SetMarkerStyle(24);
  g_rho_eta08->SetMarkerColor(2);
  g_rho_eta08->SetMarkerSize(1.2);
  g_rho_eta08->Draw("pE same");

  g_rho_eta06->SetMarkerStyle(24);
  g_rho_eta06->SetMarkerColor(4);
  g_rho_eta06->SetMarkerSize(1.2);
  g_rho_eta06->Draw("pE same");

  TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(g_rho_eta1,"|#eta| < 1.0","p");
  leg->AddEntry(g_rho_eta08,"|#eta| < 0.8","p");
  leg->AddEntry(g_rho_eta06,"|#eta| < 0.6","p");
  leg->Draw("same");

  string FigName = "/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/RhoSig/c_rhoQA_etaCut_AuAu27GeV.eps";
  c_rho->SaveAs(FigName.c_str());
}

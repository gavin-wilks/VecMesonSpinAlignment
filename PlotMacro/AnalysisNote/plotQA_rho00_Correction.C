#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include "../../Utility/draw.h"

using namespace std;

void plotQA_rho00_Correction()
{
  string inputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AnalysisNote_1st/rho00_step_200GeV_Run14.root";
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1D *h_rawRho = (TH1D*)File_InPut->Get("raw_rho00");
  TH1D *h_effRho = (TH1D*)File_InPut->Get("EffCor_rho00");
  TH1D *h_resRho = (TH1D*)File_InPut->Get("ResCor_rho00");
  TH1D *h_etaRho = (TH1D*)File_InPut->Get("EtaCor_rho00");

  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);

  h_rawRho->SetTitle("#rho_{00} Correstion Steps");
  h_rawRho->SetStats(0);
  h_rawRho->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_rawRho->GetYaxis()->SetTitle("#rho_{00}");
  h_rawRho->GetYaxis()->SetRangeUser(0.31,0.35);
  h_rawRho->GetYaxis()->SetNdivisions(505);
  h_rawRho->GetYaxis()->SetLabelSize(0.04);
  h_rawRho->GetYaxis()->SetTitleOffset(1.5);
  h_rawRho->SetMarkerStyle(24);
  h_rawRho->SetMarkerColor(kGray+3);
  h_rawRho->SetMarkerSize(1.2);
  h_rawRho->SetLineColor(kGray+3);
  h_rawRho->Draw("pEX0");
  PlotLine(1.2,5.4,1.0/3.0,1.0/3.0,1,2,2);

  h_effRho->SetMarkerStyle(25);
  h_effRho->SetMarkerColor(kAzure-3);
  h_effRho->SetMarkerSize(1.2);
  h_effRho->SetLineColor(kAzure-3);
  h_effRho->Draw("pEX0 same");

  // h_resRho->SetMarkerStyle(26);
  // h_resRho->SetMarkerColor(kRed+2);
  // h_resRho->SetMarkerSize(1.6);
  // h_resRho->SetLineColor(kRed+2);
  // h_resRho->Draw("pEX0 same");

  h_etaRho->SetMarkerStyle(30);
  h_etaRho->SetMarkerColor(kRed+2);
  h_etaRho->SetMarkerSize(1.6);
  h_etaRho->SetLineColor(kRed+2);
  h_etaRho->Draw("pEX0 same");

  TLegend *leg = new TLegend(0.4,0.7,0.8,0.85);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_rawRho,"raw #rho_{00}","p");
  leg->AddEntry(h_effRho,"efficiency corrected #rho_{00}","p");
  leg->AddEntry(h_etaRho,"final #rho_{00}","p");
  leg->Draw("same");

  c_rho->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/phiMeson/c_rhoCorrectionSteps.eps");
}

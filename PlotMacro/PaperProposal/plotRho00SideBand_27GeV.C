#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"

using namespace std;

void plotRho00SideBand_27GeV()
{
  double invMass_27GeV = 1.02;
  double rhoData_27GeV = 0.360235;
  double errData_27GeV = 0.0025245586;
  TGraphErrors *g_rhoData_27GeV = new TGraphErrors();
  g_rhoData_27GeV->SetPoint(0,invMass_27GeV,rhoData_27GeV);
  g_rhoData_27GeV->SetPointError(0,0.01,errData_27GeV);

  double invMass[3] = {1.00,1.02,1.04};
  double errX[3] = {0.01,0.01,0.01};

  // double rhoData[3] = {0.347066,0.337032,0.317486};
  // double errData[3] = {0.00290068,0.00274057,0.00270384};
  double rhoData[3] = {0.344993,0.337498,0.319533};
  double errData[3] = {0.00173669,0.00168262,0.00166905};
  TGraphErrors *g_rhoData = new TGraphErrors();
  for(int i_point = 0; i_point < 3; ++i_point)
  {
    g_rhoData->SetPoint(i_point,invMass[i_point],rhoData[i_point]);
    g_rhoData->SetPointError(i_point,0.0,errData[i_point]);
  }

  double rhoSim[3] = {0.345257,0.336874,0.328237};
  double errSim[3] = {0.00348644,0.00281394,0.00259895};
  TGraphErrors *g_rhoSim = new TGraphErrors();
  for(int i_point = 0; i_point < 3; ++i_point)
  {
    g_rhoSim->SetPoint(i_point,invMass[i_point]-0.002,rhoSim[i_point]);
    g_rhoSim->SetPointError(i_point,0.0,errSim[i_point]);
  }

  double rhoSimPt[3] = {0.335348,0.333532,0.331146};
  double errSimPt[3] = {0.0014886,0.00118673,0.00107318};
  TGraphErrors *g_rhoSimPt = new TGraphErrors();
  for(int i_point = 0; i_point < 3; ++i_point)
  {
    g_rhoSimPt->SetPoint(i_point,invMass[i_point]+0.002,rhoSimPt[i_point]);
    g_rhoSimPt->SetPointError(i_point,0.0,errSimPt[i_point]);
  }

  double rhoSimV2[3] = {0.337031,0.329063,0.320837};
  double errSimV2[3] = {0.00111772,0.000905791,0.000832921};
  TGraphErrors *g_rhoSimV2 = new TGraphErrors();
  for(int i_point = 0; i_point < 3; ++i_point)
  {
    g_rhoSimV2->SetPoint(i_point,invMass[i_point]+0.004,rhoSimV2[i_point]);
    g_rhoSimV2->SetPointError(i_point,0.0,errSimV2[i_point]);
  }

  TH1F *h_frame = new TH1F("h_frame","h_frame",200,0.98,1.08);
  for(int i_bin = 0; i_bin < 200; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  h_frame->SetTitle("Au+Au 27GeV (Run11+Run18)");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.98,1.06);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("M(K^{+},K^{-})(GeV/c^{2})");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.29,0.37);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.98,1.06,1.0/3.0,1.0/3.0,1,3,2);

  g_rhoData_27GeV->SetMarkerStyle(29);
  g_rhoData_27GeV->SetMarkerColor(2);
  g_rhoData_27GeV->SetLineColor(2);
  g_rhoData_27GeV->SetMarkerSize(1.8);
  g_rhoData_27GeV->Draw("pE same");

  g_rhoData->SetMarkerStyle(30);
  g_rhoData->SetMarkerColor(2);
  g_rhoData->SetLineColor(2);
  g_rhoData->SetMarkerSize(1.4);
  g_rhoData->Draw("pE same");

  g_rhoSim->SetMarkerStyle(20);
  g_rhoSim->SetMarkerColor(kGray+2);
  g_rhoSim->SetLineColor(kGray+2);
  g_rhoSim->SetMarkerSize(1.4);
  g_rhoSim->Draw("pE same");

  g_rhoSimPt->SetMarkerStyle(24);
  g_rhoSimPt->SetMarkerColor(kGray+2);
  g_rhoSimPt->SetLineColor(kGray+2);
  g_rhoSimPt->SetMarkerSize(1.4);
  // g_rhoSimPt->Draw("pE same");

  g_rhoSimV2->SetMarkerStyle(25);
  g_rhoSimV2->SetMarkerColor(4);
  g_rhoSimV2->SetLineColor(4);
  g_rhoSimV2->SetMarkerSize(1.4);
  // g_rhoSimV2->Draw("pE same");

  PlotLine(0.99,0.99,0.29,0.37,1,3,2);
  PlotLine(1.01,1.01,0.29,0.37,1,3,2);
  PlotLine(1.03,1.03,0.29,0.37,1,3,2);
  PlotLine(1.05,1.05,0.29,0.37,1,3,2);

  TLegend *leg = new TLegend(0.2,0.2,0.7,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(g_rhoData_27GeV,"STAR Signal","P");
  leg->AddEntry(g_rhoData,"STAR Side Band (Sig+Bkg)","P");
  leg->AddEntry(g_rhoSim,"Simulation with Kaon pairs from Data","P");
  // leg->AddEntry(g_rhoSimPt,"Simulation with published #phi spectra","P");
  // leg->AddEntry(g_rhoSimV2,"Simulation with v_{2} = 0","P");
  leg->Draw("same");

  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_SideBand_27GeV.eps");
}

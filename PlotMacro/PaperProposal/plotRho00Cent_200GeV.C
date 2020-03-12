#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../Utility/draw.h"
// #include "../../Utility/StSpinAlignmentCons.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);

void plotRho00Cent_200GeV()
{
  string input_run11 = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/RunDependence/rhoCent_200GeV_Run11_LXaxis.root";
  TFile *File_Run11 = TFile::Open(input_run11.c_str());
  TGraphAsymmErrors *g_rhoCent_200GeV_Run11_stat = (TGraphAsymmErrors*)File_Run11->Get("g_rhoCent_200GeV_Run11_2nd_stat");
  TGraphAsymmErrors *g_rhoCent_200GeV_Run11_sys = (TGraphAsymmErrors*)File_Run11->Get("g_rhoCent_200GeV_Run11_2nd_sys");

  string input_run14 = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/RunDependence/rhoCent_200GeV_Run14_LXaxis.root";
  TFile *File_Run14 = TFile::Open(input_run14.c_str());
  TGraphAsymmErrors *g_rhoCent_200GeV_Run14_stat = (TGraphAsymmErrors*)File_Run14->Get("g_rhoCent_200GeV_Run14_2nd_stat");
  TGraphAsymmErrors *g_rhoCent_200GeV_Run14_sys = (TGraphAsymmErrors*)File_Run14->Get("g_rhoCent_200GeV_Run14_2nd_sys");

  TH1F *h_frame = new TH1F("h_frame","h_frame",100,0,100);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
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
  // c_rho00->cd()->SetLogx();
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,80.0);
  h_frame->GetXaxis()->SetNdivisions(510,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("Centrality (%)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.22,0.44);
  // if(energy == 3) h_frame->GetYaxis()->SetRangeUser(0.301,0.44);
  // if(energy == 6) h_frame->GetYaxis()->SetRangeUser(0.22,0.4);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,80.0,1.0/3.0,1.0/3.0,1,3,2);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoCent_200GeV_Run11_stat,30,kGray+2,1.4);
  plotSysErrors(g_rhoCent_200GeV_Run11_sys,kGray+2);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoCent_200GeV_Run14_stat,29,kRed,1.8);
  plotSysErrors(g_rhoCent_200GeV_Run14_sys,kRed);

  Draw_TGAE_Point_new_Symbol(10,0.42,0.0,0.0,0.0,0.0,30,kGray+2,1.4);
  plotTopLegend((char*)"Run11 (2^{nd}-order EP)",13,0.4175,0.04,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(10,0.40,0.0,0.0,0.0,0.0,29,kRed,1.8);
  plotTopLegend((char*)"Run14 (2^{nd}-order EP)",13,0.3975,0.04,1,0.0,42,0);

  string leg_energy = "AuAu 200 GeV (|#eta| < 1)";
  plotTopLegend((char*)leg_energy.c_str(),0.38,0.25,0.04,1,0.0,42,1);
  plotTopLegend((char*)"#phi-meson (1.2 < p_{T}< 5.4 GeV/c)",0.38,0.20,0.04,1,0.0,42,1);

  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_rhoRunSys_200GeV.eps");
  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_rhoRunSys_200GeV.png");
}

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color)
{
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);

    PlotLine(energy-1,energy+1,rho+err,rho+err,plot_color,2,1);
    PlotLine(energy-1,energy-1,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(energy+1,energy+1,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(energy-1,energy+1,rho-err,rho-err,plot_color,2,1);
    PlotLine(energy-1,energy-1,rho-err+0.001,rho-err,plot_color,2,1);
    PlotLine(energy+1,energy+1,rho-err+0.001,rho-err,plot_color,2,1);
  }
}

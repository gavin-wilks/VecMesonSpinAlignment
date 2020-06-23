#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../../Utility/draw.h"
#include "../../../Utility/StSpinAlignmentCons.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);

void plotFig5_Rho00EnergyPhiOnly()
{
  gStyle->SetOptDate(0);
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed;

  const float size_marker = 1.4;
  const float size_font = 0.03;

  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *g_rhoPhi_1st_stat = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_1stEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_1st_sys  = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_1stEP_energy_sys");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys  = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_sys");

  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,0,1000);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
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
  c_rho00->cd()->SetLogx();
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(9.0,240.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.301,0.40);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  // h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(9.0,240.0,1.0/3.0,1.0/3.0,1,3,2);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_1st_stat,style_phi_1st,color_phi_1st,size_marker-0.2);
  plotSysErrors(g_rhoPhi_1st_sys,color_phi_1st);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_2nd_stat,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotSysErrors(g_rhoPhi_2nd_sys,color_phi_2nd);

  // plot legend
  Draw_TGAE_Point_new_Symbol(80,0.39,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
  plotTopLegend((char*)"1^{st}-order EP",90,0.3885,0.04,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(80,0.38,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotTopLegend((char*)"2^{nd}-order EP",90,0.3785,0.04,1,0.0,42,0);

  plotTopLegend((char*)"Au+Au (20-60\% & |y| < 1)",0.45,0.30,0.04,1,0.0,42,1);
  plotTopLegend((char*)"#phi-meson (1.2 < p_{T}< 5.4 GeV/c)",0.4,0.25,0.04,1,0.0,42,1);

#if 0
  //-----------plot mean rho_00------------------
  TFile *File_Input_mean = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/rho00_mean_BES.root");
  TGraphAsymmErrors *g_mean_rho_1st_stat = (TGraphAsymmErrors*)File_Input_mean->Get("rho00_1stEP_mean_stat");
  TGraphAsymmErrors *g_mean_rho_1st_sys  = (TGraphAsymmErrors*)File_Input_mean->Get("rho00_1stEP_mean_sys");
  TGraphAsymmErrors *g_mean_rho_2nd_stat = (TGraphAsymmErrors*)File_Input_mean->Get("rho00_2ndEP_mean_stat");
  TGraphAsymmErrors *g_mean_rho_2nd_sys  = (TGraphAsymmErrors*)File_Input_mean->Get("rho00_2ndEP_mean_sys");
  if(isPlotMean)
  {
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mean_rho_1st_stat,24,kAzure+2,1.4);
    for(int i_energy = 0; i_energy < 1; ++i_energy) // plot sys errors
    {
      double energy, rho;
      g_mean_rho_1st_sys->GetPoint(i_energy,energy,rho);
      double err = g_mean_rho_1st_sys->GetErrorYhigh(i_energy);

      PlotLine(energy*0.95,energy*1.05,rho+err,rho+err,kAzure+2,2,1);
      PlotLine(energy*0.95,energy*0.95,rho+err-0.001,rho+err,kAzure+2,2,1);
      PlotLine(energy*1.05,energy*1.05,rho+err-0.001,rho+err,kAzure+2,2,1);
      PlotLine(energy*0.95,energy*1.05,rho-err,rho-err,kAzure+2,2,1);
      PlotLine(energy*0.95,energy*0.95,rho-err+0.001,rho-err,kAzure+2,2,1);
      PlotLine(energy*1.05,energy*1.05,rho-err+0.001,rho-err,kAzure+2,2,1);
    }

    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mean_rho_2nd_stat,30,kRed,1.8);
    for(int i_energy = 0; i_energy < 1; ++i_energy) // plot sys errors
    {
      double energy, rho;
      g_mean_rho_2nd_sys->GetPoint(i_energy,energy,rho);
      double err = g_mean_rho_2nd_sys->GetErrorYhigh(i_energy);

      PlotLine(energy*0.95,energy*1.05,rho+err,rho+err,kRed,2,1);
      PlotLine(energy*0.95,energy*0.95,rho+err-0.001,rho+err,kRed,2,1);
      PlotLine(energy*1.05,energy*1.05,rho+err-0.001,rho+err,kRed,2,1);
      PlotLine(energy*0.95,energy*1.05,rho-err,rho-err,kRed,2,1);
      PlotLine(energy*0.95,energy*0.95,rho-err+0.001,rho-err,kRed,2,1);
      PlotLine(energy*1.05,energy*1.05,rho-err+0.001,rho-err,kRed,2,1);
    }

  }

  //-----------plot mean rho_00------------------

  //-----------plot theory rho_00------------------
  TFile *File_Input_Theory = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Theory.root");
  TGraph *g_Theory_120 = (TGraph*)File_Input_Theory->Get("g_Theory_120");
  TGraph *g_Theory_128 = (TGraph*)File_Input_Theory->Get("g_Theory_128");
  TGraph *g_Theory_136 = (TGraph*)File_Input_Theory->Get("g_Theory_136");

  g_Theory_136->SetLineColor(2);
  g_Theory_136->SetLineStyle(1);
  g_Theory_136->SetLineWidth(4);
  g_Theory_136->Draw("l Same");

  g_Theory_128->SetLineColor(6);
  g_Theory_128->SetLineStyle(3);
  g_Theory_128->SetLineWidth(4);
  g_Theory_128->Draw("l Same");

  g_Theory_120->SetLineColor(4);
  g_Theory_120->SetLineStyle(9);
  g_Theory_120->SetLineWidth(4);
  g_Theory_120->Draw("l Same");
  //-----------plot theory rho_00------------------
#endif

  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig5_rho00EnergyPhiOnly.eps");
  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig5_rho00EnergyPhiOnly.png");
}

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color)
{
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);

    PlotLine(energy*0.95,energy*1.05,rho+err,rho+err,plot_color,2,1);
    PlotLine(energy*0.95,energy*0.95,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(energy*1.05,energy*1.05,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(energy*0.95,energy*1.05,rho-err,rho-err,plot_color,2,1);
    PlotLine(energy*0.95,energy*0.95,rho-err+0.001,rho-err,plot_color,2,1);
    PlotLine(energy*1.05,energy*1.05,rho-err+0.001,rho-err,plot_color,2,1);
  }
}

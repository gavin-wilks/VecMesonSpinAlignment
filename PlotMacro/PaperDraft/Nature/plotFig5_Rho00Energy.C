#include <string>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TStyle.h>
#include <TF1.h>

#include "../../../Utility/draw.h"
#include "../../../Utility/StSpinAlignmentCons.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);

void plotFig5_Rho00Energy()
{
  gStyle->SetOptDate(0);
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed;
  const int style_phi_ALICE = 30;
  const int color_phi_ALICE = kRed;

  const int style_Kstr = 20;
  const int color_Kstr = kAzure+2;
  const int style_Kstr_ALICE = 24;
  const int color_Kstr_ALICE = kAzure+2;

  const float size_marker = 1.4;
  const float size_font = 0.03;
  
  //----------------------------------------------------------
  // phi-meson STAR
  //beam-energy dependence of phi-meson rho00 from STAR, pT: 1.2 - 5.4 GeV/c, 20-60%
  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *g_rhoPhi_1st_stat = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_1stEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_1st_sys  = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_1stEP_energy_sys");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys  = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_sys");
  //----------------------------------------------------------

  //----------------------------------------------------------
  TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_sNN.root");
  // K* STAR
  //beam-energy dependence of kstar rho00 from STAR, pT: 1.0 - 1.5 GeV/c, 20-60%
  TGraphAsymmErrors *g_rhoKstar_stat       = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_stat");
  TGraphAsymmErrors *g_rhoKstar_sys        = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_sys");
  //data points from ALICE
  //K* ALICE pT: 0.8 - 1.2 GeV/c
  TGraphAsymmErrors *g_rhoKstar_ALICE_stat = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_ALICE_stat");
  TGraphAsymmErrors *g_rhoKstar_ALICE_sys  = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_ALICE_sys");
  //Phi ALICE pT: 1.0 - 5.0 GeV/c
  TGraphAsymmErrors *g_rhoPhi_ALICE_stat   = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoPhi_ALICE_stat");
  TGraphAsymmErrors *g_rhoPhi_ALICE_sys    = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoPhi_ALICE_sys");
  //----------------------------------------------------------

  TH1F *h_frame = new TH1F("h_frame","h_frame",5000,0,5000);
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
  h_frame->GetXaxis()->SetRangeUser(8.0,4096.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.14,0.42);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  // h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(8.0,4096.0,1.0/3.0,1.0/3.0,1,3,2);

  // phi-meson STAR
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_1st_stat,style_phi_1st,color_phi_1st,size_marker-0.2);
  plotSysErrors(g_rhoPhi_1st_sys,color_phi_1st);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_2nd_stat,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotSysErrors(g_rhoPhi_2nd_sys,color_phi_2nd);

  // K* STAR
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKstar_stat,style_Kstr,color_Kstr,size_marker);
  plotSysErrors(g_rhoKstar_sys,color_Kstr);

  // phi-meson ALICE
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_ALICE_stat,style_phi_ALICE,color_phi_ALICE,size_marker);
  plotSysErrors(g_rhoPhi_ALICE_sys,color_phi_ALICE);

  // K* ALICE
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKstar_ALICE_stat,style_Kstr_ALICE,color_Kstr_ALICE,size_marker);
  plotSysErrors(g_rhoKstar_ALICE_sys,color_Kstr_ALICE);

  plotTopLegend((char*)"Au+Au 20-60%",0.25,0.85,size_font+0.01,1,0.0,42,1);
  plotTopLegend((char*)"Pb+Pb 10-50%",0.65,0.85,size_font+0.01,1,0.0,42,1);
  // plotTopLegend((char*)"ALICE (Ref)",0.65,0.85,size_font+0.01,1,0.0,42,1);

  Draw_TGAE_Point_new_Symbol(35,0.22,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotTopLegend((char*)"#phi (2^{nd}-Order EP)",40,0.2165,size_font,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(35,0.20,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
  plotTopLegend((char*)"#phi (1^{st}-Order EP)",40,0.1965,size_font,1,0.0,42,0);
  plotTopLegend((char*)"1.2 < p_{T} < 5.4 GeV/c",35,0.1765,size_font,1,0.0,42,0);
  plotTopLegend((char*)"|y| < 1.0",35,0.1565,size_font,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(450,0.22,0.0,0.0,0.0,0.0,style_Kstr,color_Kstr,size_marker);
  plotTopLegend((char*)"K^{*0} (2^{nd}-Order EP)",550,0.2165,size_font,1,0.0,42,0);
  plotTopLegend((char*)"1.0 < p_{T} < 1.5 GeV/c",450,0.1965,size_font,1,0.0,42,0);
  plotTopLegend((char*)"|y| < 0.5",450,0.1765,size_font,1,0.0,42,0);

  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig5_rho00Energy.eps");
  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig5_rho00Energy.png");
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

#include <string>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TStyle.h>
#include <TF1.h>

#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void genRho00MuB()
{
  gStyle->SetOptDate(0);
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int style_phi_ALICE = 30;
  const int color_phi_ALICE = kGray+1;

  const int style_Kstr = 20;
  const int color_Kstr = kAzure-9;
  const int style_Kstr_ALICE = 24;
  const int color_Kstr_ALICE = kGray+1;

  const float size_marker = 1.4;
  const float size_font = 0.035;

  const float beamEnergy[9] = {7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 54.4, 62.4, 200.0};
  const float muB[9]        = {420, 315, 260, 205, 155, 115, 83, 70, 20};
  std::map<float, float> map_muB;
  for(int iEnergy = 0; iEnergy < 9; ++iEnergy)
  {
    float energy = beamEnergy[iEnergy];
    map_muB[energy] = muB[iEnergy]; 
  }
  
  //----------------------------------------------------------
  // phi-meson STAR
  //beam-energy dependence of phi-meson rho00 from STAR, pT: 1.2 - 5.4 GeV/c, 20-60%
  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/BESII/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat_Laxis = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys_Laxis = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_sys");

  TGraphAsymmErrors *g_rhoMuBPhi_stat = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rhoMuBPhi_sys = new TGraphAsymmErrors();
  for(int iEnergy = 0; iEnergy < g_rhoPhi_2nd_stat_Laxis->GetN(); ++iEnergy)
  {
    double energy, rho;
    g_rhoPhi_2nd_stat_Laxis->GetPoint(iEnergy,energy,rho);
    double err_stat = g_rhoPhi_2nd_stat_Laxis->GetErrorYhigh(iEnergy);
    double err_sys = g_rhoPhi_2nd_sys_Laxis->GetErrorYhigh(iEnergy);
    // cout<< "iEnergy: " << iEnergy << " Energy " << energy << " GeV & muB " << map_muB[energy] << ": " << rho << " #pm " << err_stat << " #pm " << err_sys << endl;
    g_rhoMuBPhi_stat->SetPoint(iEnergy,map_muB[energy],rho);
    g_rhoMuBPhi_stat->SetPointError(iEnergy,0.0,0.0,err_stat,err_stat);
    g_rhoMuBPhi_sys->SetPoint(iEnergy,map_muB[energy],rho);
    g_rhoMuBPhi_sys->SetPointError(iEnergy,0.0,0.0,err_sys,err_sys);
  }
  //----------------------------------------------------------

  //----------------------------------------------------------
  // K* STAR
  //beam-energy dependence of kstar rho00 from STAR, pT: 1.0 - 1.5 GeV/c, 20-60%
  TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_sNN_June9_2021.root");
  TGraphAsymmErrors *g_rhoKstar_stat       = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_stat");
  TGraphAsymmErrors *g_rhoKstar_sys        = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_sys");

  TGraphAsymmErrors *g_rhoMuBKstar_stat = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rhoMuBKstar_sys = new TGraphAsymmErrors();
  for(int iEnergy = 0; iEnergy < g_rhoKstar_stat->GetN(); ++iEnergy)
  {
    double energy, rho;
    g_rhoKstar_stat->GetPoint(iEnergy,energy,rho);
    double err_stat = g_rhoKstar_stat->GetErrorYhigh(iEnergy);
    double err_sys = g_rhoKstar_sys->GetErrorYhigh(iEnergy);
    // cout<< "iEnergy: " << iEnergy << " Energy " << energy << " GeV & muB " << map_muB[energy] << ": " << rho << " #pm " << err_stat << " #pm " << err_sys << endl;
    g_rhoMuBKstar_stat->SetPoint(iEnergy,map_muB[energy],rho);
    g_rhoMuBKstar_stat->SetPointError(iEnergy,0.0,0.0,err_stat,err_stat);
    g_rhoMuBKstar_sys->SetPoint(iEnergy,map_muB[energy],rho);
    g_rhoMuBKstar_sys->SetPointError(iEnergy,0.0,0.0,err_sys,err_sys);
  }
  //----------------------------------------------------------

  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  // c_rho00->cd()->SetLogx();

  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,0,1000);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,500.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#mu_{B} (MeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.27,0.43);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,500.0,1.0/3.0,1.0/3.0,1,2,2);

  // K* STAR
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoMuBKstar_stat,style_Kstr,color_Kstr,size_marker-0.4);
  // plotSysErrors(g_rhoKstar_sys,color_Kstr+2);
  plotSysErrorsBox(g_rhoMuBKstar_sys,color_Kstr+2);

  // phi-meson STAR
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoMuBPhi_stat,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  // plotSysErrors(g_rhoPhi_2nd_sys_Laxis,color_phi_2nd);
  plotSysErrorsBox(g_rhoMuBPhi_sys,color_phi_2nd);

  Draw_TGAE_Point_new_Symbol(40,0.42,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotTopLegend((char*)"#phi    (|y| < 1.0 & 1.2 < p_{T} < 5.4 GeV/c)",48,0.4175,size_font,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(40,0.405,0.0,0.0,0.0,0.0,style_Kstr,color_Kstr,size_marker-0.4);
  plotTopLegend((char*)"K^{*0} (|y| < 1.0 & 1.0 < p_{T} < 5.0 GeV/c)",48,0.4025,size_font,1,0.0,42,0);

  // c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/TestMuB/rho00_muB.eps");
  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/TestMuB/VecMeson/rho00MuB.root","RECREATE");
  File_OutPut->cd();
  g_rhoMuBPhi_stat->SetName("g_rhoMuBPhi_stat");
  g_rhoMuBPhi_stat->Write();
  g_rhoMuBPhi_sys->SetName("g_rhoMuBPhi_sys");
  g_rhoMuBPhi_sys->Write();
  g_rhoMuBKstar_stat->SetName("g_rhoMuBKstar_stat");
  g_rhoMuBKstar_stat->Write();
  g_rhoMuBKstar_sys->SetName("g_rhoMuBKstar_sys");
  g_rhoMuBKstar_sys->Write();
  File_OutPut->Close();
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

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nEnergy = g_rho->GetN();
  TBox *bSys[nEnergy];
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);

    // bSys[i_energy] = new TBox(energy/1.08,rho-err,energy*1.08,rho+err);
    bSys[i_energy] = new TBox(energy-4.5,rho-err,energy+4.5,rho+err);
    bSys[i_energy]->SetFillColor(0);
    bSys[i_energy]->SetFillStyle(0);
    bSys[i_energy]->SetLineStyle(1);
    bSys[i_energy]->SetLineWidth(1);
    bSys[i_energy]->SetLineColor(plot_color);
    bSys[i_energy]->Draw("l Same");
  }
}

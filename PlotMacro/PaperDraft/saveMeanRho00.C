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

void saveMeanRho00()
{
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
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(9.0,240.0,1.0/3.0,1.0/3.0,1,3,2);

  float mean_energy = 30.0;
  float rho_mean_1st = 0.356459;
  float err_stat_1st = 0.00531003;
  float err_sys_1st = 0.00373923;

  TGraphAsymmErrors *g_mean_rho_1st_stat = new TGraphAsymmErrors();
  g_mean_rho_1st_stat->SetPoint(0,mean_energy+2,rho_mean_1st);
  g_mean_rho_1st_stat->SetPointError(0,0.0,0.0,err_stat_1st,err_stat_1st);

  TGraphAsymmErrors *g_mean_rho_1st_sys = new TGraphAsymmErrors();
  g_mean_rho_1st_sys->SetPoint(0,mean_energy+2,rho_mean_1st);
  g_mean_rho_1st_sys->SetPointError(0,0.0,0.0,err_sys_1st,err_sys_1st);

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

  float rho_mean_2nd = 0.349424;
  float err_stat_2nd = 0.00169799;
  float err_sys_2nd = 0.000872537;
  TGraphAsymmErrors *g_mean_rho_2nd_stat = new TGraphAsymmErrors();
  g_mean_rho_2nd_stat->SetPoint(0,mean_energy,rho_mean_2nd);
  g_mean_rho_2nd_stat->SetPointError(0,0.0,0.0,err_stat_2nd,err_stat_2nd);

  TGraphAsymmErrors *g_mean_rho_2nd_sys = new TGraphAsymmErrors();
  g_mean_rho_2nd_sys->SetPoint(0,mean_energy,rho_mean_2nd);
  g_mean_rho_2nd_sys->SetPointError(0,0.0,0.0,err_sys_2nd,err_sys_2nd);

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

  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/rho00_mean_BES.root","RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  g_mean_rho_1st_stat->SetName("rho00_1stEP_mean_stat");
  g_mean_rho_1st_stat->Write();
  g_mean_rho_1st_sys->SetName("rho00_1stEP_mean_sys");
  g_mean_rho_1st_sys->Write();

  g_mean_rho_2nd_stat->SetName("rho00_2ndEP_mean_stat");
  g_mean_rho_2nd_stat->Write();
  g_mean_rho_2nd_sys->SetName("rho00_2ndEP_mean_sys");
  g_mean_rho_2nd_sys->Write();
  File_OutPut->Close();
}

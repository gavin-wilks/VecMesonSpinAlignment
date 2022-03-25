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

void genPHMuB()
{
  const int style_L = 29;
  const int color_L = kRed-4;
  const int style_LBar = 20;
  const int color_LBar = kAzure-9;

  const float size_marker = 1.4;
  const float size_font = 0.035;

  gStyle->SetOptDate(0);
  const int nEnergy = 9;
  const float beamEnergy[nEnergy] = {3.0, 7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 200.0};
  const float muB[nEnergy]        = {750, 420, 315, 260, 205, 155, 115, 70, 20};

  float PH[nEnergy] = {4.90823,1.8,1.19,1.17,0.84,0.93,0.45,1.2,0.277};
  float errPH_statHigh[nEnergy] = {0.81389,0.6,0.35,0.4,0.27,0.25,0.4,1.0,0.04};
  float errPH_statLow[nEnergy] = {-0.81389,-0.6,-0.35,-0.4,-0.27,-0.25,-0.4,-1.0,-0.04};
  float errPH_sysHigh[nEnergy] = {0.15485,0,0,0,0,0,0,0,0.039};
  float errPH_sysLow[nEnergy] = {-0.15485,-0.18,-0.18,-0.27,-0.18,-0.18,-0.18,0,-0.049};

  float PHbar[nEnergy] = {0,7.7,1.59,2.01,1.34,1.1,0.83,1.5,0.24};
  float errPHbar_statHigh[nEnergy] = {0,3.2,1.1,1.1,0.5,0.4,0.5,1.4,0.045};
  float errPHbar_statLow[nEnergy]  = {0,-3.2,-1.1,-1.1,-0.5,-0.4,-0.5,-1.4,-0.045};
  float errPHbar_sysHigh[nEnergy] = {0,0,0,0.35,0,0,0,0,0.061};
  float errPHbar_sysLow[nEnergy]  = {0,-0.9,-0.13,-0.13,-0.13,-0.13,-0.13,0,-0.045};
  
  TGraphAsymmErrors *g_PHMuB_stat = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_PHMuB_sys = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_PHbarMuB_stat = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_PHbarMuB_sys = new TGraphAsymmErrors();
  for(int iEnergy = 0; iEnergy < nEnergy; ++iEnergy)
  {
    g_PHMuB_stat->SetPoint(iEnergy,muB[iEnergy],PH[iEnergy]);
    g_PHMuB_stat->SetPointEYhigh(iEnergy,errPH_statHigh[iEnergy]);
    g_PHMuB_stat->SetPointEYlow(iEnergy,-1.0*errPH_statLow[iEnergy]);
    g_PHMuB_sys->SetPoint(iEnergy,muB[iEnergy],PH[iEnergy]);
    g_PHMuB_sys->SetPointEYhigh(iEnergy,errPH_sysHigh[iEnergy]);
    g_PHMuB_sys->SetPointEYlow(iEnergy,-1.0*errPH_sysLow[iEnergy]);

    g_PHbarMuB_stat->SetPoint(iEnergy,muB[iEnergy],PHbar[iEnergy]);
    g_PHbarMuB_stat->SetPointEYhigh(iEnergy,errPHbar_statHigh[iEnergy]);
    g_PHbarMuB_stat->SetPointEYlow(iEnergy,-1.0*errPHbar_statLow[iEnergy]);
    g_PHbarMuB_sys->SetPoint(iEnergy,muB[iEnergy],PHbar[iEnergy]);
    g_PHbarMuB_sys->SetPointEYhigh(iEnergy,errPHbar_sysHigh[iEnergy]);
    g_PHbarMuB_sys->SetPointEYlow(iEnergy,-1.0*errPHbar_sysLow[iEnergy]);
  }
  g_PHbarMuB_stat->RemovePoint(0);
  g_PHbarMuB_sys->RemovePoint(0);
  //----------------------------------------------------------

  TCanvas *c_PH = new TCanvas("c_PH","c_PH",10,10,800,800);
  c_PH->cd();
  c_PH->cd()->SetLeftMargin(0.15);
  c_PH->cd()->SetBottomMargin(0.15);
  c_PH->cd()->SetTicks(1,1);
  c_PH->cd()->SetGrid(0,0);
  // c_PH->cd()->SetLogx();

  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,0,1000);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,800.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#mu_{B} (MeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(-0.5,12);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#bar{P}_{H} (%)");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,800.0,0.0,0.0,1,2,2);

  // g_PHbarMuB_sys->Draw("pE same");

  // Lambda STAR
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_PHMuB_stat,style_L,color_L,size_marker+0.2);
  plotSysErrorsBox(g_PHMuB_sys,color_L+2);

  // antiLambda STAR
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_PHbarMuB_stat,style_LBar,color_LBar,size_marker-0.2);
  plotSysErrorsBox(g_PHbarMuB_sys,color_LBar+2);

  Draw_TGAE_Point_new_Symbol(50,11,0.0,0.0,0.0,0.0,style_L,color_L,size_marker+0.2);
  plotTopLegend((char*)"#Lambda ",65,10.8,size_font,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(50,10,0.0,0.0,0.0,0.0,style_LBar,color_LBar,size_marker-0.2);
  plotTopLegend((char*)"#bar{#Lambda}",65,9.8,size_font,1,0.0,42,0);

  // c_PH->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/TestMuB/PH_muB.eps");
  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/TestMuB/Lambda/PHMuB.root","RECREATE");
  File_OutPut->cd();
  g_PHMuB_stat->SetName("g_PHMuB_stat");
  g_PHMuB_stat->Write();
  g_PHMuB_sys->SetName("g_PHMuB_sys");
  g_PHMuB_sys->Write();
  g_PHbarMuB_stat->SetName("g_PHbarMuB_stat");
  g_PHbarMuB_stat->Write();
  g_PHbarMuB_sys->SetName("g_PHbarMuB_sys");
  g_PHbarMuB_sys->Write();
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
    double errHigh = g_rho->GetErrorYhigh(i_energy);
    double errLow  = g_rho->GetErrorYlow(i_energy);

    // bSys[i_energy] = new TBox(energy/1.08,rho-err,energy*1.08,rho+err);
    bSys[i_energy] = new TBox(energy-7.5,rho-errLow,energy+7.5,rho+errHigh);
    bSys[i_energy]->SetFillColor(0);
    bSys[i_energy]->SetFillStyle(0);
    bSys[i_energy]->SetLineStyle(1);
    bSys[i_energy]->SetLineWidth(1);
    bSys[i_energy]->SetLineColor(plot_color);
    bSys[i_energy]->Draw("l Same");
  }
}

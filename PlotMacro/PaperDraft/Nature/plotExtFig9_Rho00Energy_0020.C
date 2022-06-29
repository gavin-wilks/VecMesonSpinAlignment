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
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

double rho00_theory(double *x_var, double *par)
{
  double s12 = x_var[0];
  double scurrent = par[0]; // strangeness current
  double ms = par[1]; // mass of s-quark MeV

  double mphi = 1020.0; // MeV
  double c1 = 300.0/(pow(200.0,1.0/3.0) * pow((-0.4+0.39*log(200.0*200.0)),1.0/3.0)); // 300.0/Teff[200.0]
  double gphi = 2.0*sqrt(2.0);

  double Teff = pow(s12,1.0/3.0) * pow((-0.4+0.39*log(s12*s12)),1.0/3.0);

  double denom_phifield = 27.0*pow(ms,4.0)*pow(mphi,4.0)*pow(c1,2.0)*pow(Teff,2.0);
  double numer_phifield = scurrent*pow(197.0,8.0)*1.8*1.0e+5; // <p^2>_phi = 0.18 GeV^2

  double rho00 = 1.0/3.0 + numer_phifield/denom_phifield;

  return rho00;
}


void plotExtFig9_Rho00Energy_0020()
{
  gStyle->SetOptDate(0);
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int style_phi_ALICE = 30;
  const int color_phi_ALICE = kGray+1;
  const int colorDiff_phi = 0;

  const int style_Kstr = 20;
  const int color_Kstr = kAzure-9;
  const int style_Kstr_ALICE = 24;
  const int color_Kstr_ALICE = kGray+1;
  const int colorDiff_Kstr = 2;

  const float size_marker = 1.4;
  const float size_font = 0.035;
  
  //----------------------------------------------------------
  // phi-meson STAR
  //beam-energy dependence of phi-meson rho00 from STAR, pT: 1.2 - 5.4 GeV/c, 0-20%
  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/NewF_JHChen/rhoCent_0020_Laxis.root");
  TGraphAsymmErrors *g_rhoPhi_1st_stat = (TGraphAsymmErrors*)File_InputPhi->Get("g_rho1st_0020_stat");
  TGraphAsymmErrors *g_rhoPhi_1st_sys = (TGraphAsymmErrors*)File_InputPhi->Get("g_rho1st_0020_sys");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat = (TGraphAsymmErrors*)File_InputPhi->Get("g_rho2nd_0020_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys = (TGraphAsymmErrors*)File_InputPhi->Get("g_rho2nd_0020_sys");
  //----------------------------------------------------------

  //----------------------------------------------------------
  TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_Cent_0020.root");
  // K* STAR
  //beam-energy dependence of kstar rho00 from STAR, pT: 1.0 - 1.5 GeV/c, 0-20%
  TGraphAsymmErrors *g_rhoKstar_stat       = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoCentKstar_0020_stat");
  TGraphAsymmErrors *g_rhoKstar_sys        = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoCentKstar_0020_sys");
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
  h_frame->GetXaxis()->SetRangeUser(7.0,4996.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.24,0.42);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  // h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(7.0,4996.0,1.0/3.0,1.0/3.0,1,2,2);

  // K* STAR
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKstar_stat,style_Kstr,color_Kstr,colorDiff_Kstr,size_marker-0.4);
  plotSysErrorsBox(g_rhoKstar_sys,color_Kstr+2);

  // phi-meson STAR
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_1st_stat,style_phi_1st,color_phi_1st,colorDiff_phi,size_marker-0.2);
  plotSysErrorsBox(g_rhoPhi_1st_sys,color_phi_1st);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_2nd_stat,style_phi_2nd,color_phi_2nd,colorDiff_phi,size_marker+0.2);
  plotSysErrorsBox(g_rhoPhi_2nd_sys,color_phi_2nd);

  Draw_TGAE_Point_new_Symbol(10,0.41,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
  plotTopLegend((char*)"#phi",12,0.4075,size_font,1,0.0,42,0);
  plotTopLegend((char*)"1^{st}-order EP (1.2 < p_{T} < 5.4 GeV/c)",19,0.4075,size_font,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(10,0.395,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotTopLegend((char*)"#phi",12,0.3925,size_font,1,0.0,42,0);
  plotTopLegend((char*)"2^{nd}-order EP (1.2 < p_{T} < 5.4 GeV/c)",19,0.3925,size_font,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(10,0.380,0.0,0.0,0.0,0.0,style_Kstr,color_Kstr,size_marker-0.4);
  plotTopLegend((char*)"K^{*0}",12,0.3775,size_font,1,0.0,42,0);
  plotTopLegend((char*)"2^{nd}-order EP (1.0 < p_{T} < 5.0 GeV/c)",19,0.3775,size_font,1,0.0,42,0);

  plotTopLegend((char*)"Au+Au (0-20\% & |y| < 1.0)",0.50,0.20,size_font,1,0.0,42,1);

  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/NatureSubmission/NewF_JHChen/extFig9_rho00Energy_0020.eps");
  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/NatureSubmission/NewF_JHChen/extFig9_rho00Energy_0020.png");
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

    bSys[i_energy] = new TBox(energy/1.08,rho-err,energy*1.08,rho+err);
    // bSys[i_energy] = new TBox(energy-1.5,rho-err,energy+1.5,rho+err);
    bSys[i_energy]->SetFillColor(0);
    bSys[i_energy]->SetFillStyle(0);
    bSys[i_energy]->SetLineStyle(1);
    bSys[i_energy]->SetLineWidth(1);
    bSys[i_energy]->SetLineColor(plot_color);
    bSys[i_energy]->Draw("l Same");
  }
}

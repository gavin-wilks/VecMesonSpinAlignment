#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../Utility/draw.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"

using namespace std;

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

void plot2ndRho00EnergyTheoryEP()
{
  bool isPlotMean = false;
  double font_size = 0.035;

  gStyle->SetOptDate(0);
  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/BESII/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *g_rho_2nd_stat_Laxis = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rho_2nd_sys_Laxis  = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_energy_sys");

  TFile *File_Input_EP = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/BESII/rho00_stat_sys_Xaxis.root");
  TGraphAsymmErrors *g_rho_2nd_stat_Xaxis = (TGraphAsymmErrors*)File_Input_EP->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rho_2nd_sys_Xaxis = (TGraphAsymmErrors*)File_Input_EP->Get("rho00_2ndEP_energy_sys");

  // plot frame
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

  h_frame->GetYaxis()->SetRangeUser(0.291,0.40);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(9.0,240.0,1.0/3.0,1.0/3.0,1,3,2);

  // plot 2nd EP rho00 in normal direction
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_2nd_stat_Laxis,29,kRed,1.8);
  for(int i_energy = 0; i_energy < g_rho_2nd_sys_Laxis->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho_2nd_sys_Laxis->GetPoint(i_energy,energy,rho);
    double err = g_rho_2nd_sys_Laxis->GetErrorYhigh(i_energy);

    PlotLine(energy*0.95,energy*1.05,rho+err,rho+err,kRed+2,2,1);
    PlotLine(energy*0.95,energy*0.95,rho+err-0.001,rho+err,kRed+2,2,1);
    PlotLine(energy*1.05,energy*1.05,rho+err-0.001,rho+err,kRed+2,2,1);
    PlotLine(energy*0.95,energy*1.05,rho-err,rho-err,kRed+2,2,1);
    PlotLine(energy*0.95,energy*0.95,rho-err+0.001,rho-err,kRed+2,2,1);
    PlotLine(energy*1.05,energy*1.05,rho-err+0.001,rho-err,kRed+2,2,1);
  }
  // Draw_TGAE_Point_new_Symbol(80,0.38,0.0,0.0,0.0,0.0,29,kRed,1.8);
  // // plotTopLegend((char*)"2^{nd} EP Normal",90,0.3790,font_size,1,0.0,42,0);
  // plotTopLegend((char*)"Out-of-Plane #rho_{00}",90,0.3790,font_size,1,0.0,42,0);

  // plot 2nd EP rho00 in tangent direction
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rho_2nd_stat_Xaxis,20,kAzure+2,1.4);
  for(int i_energy = 0; i_energy < g_rho_2nd_sys_Xaxis->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho_2nd_sys_Xaxis->GetPoint(i_energy,energy,rho);
    double err = g_rho_2nd_sys_Xaxis->GetErrorYhigh(i_energy);

    PlotLine(energy*0.95,energy*1.05,rho+err,rho+err,kAzure+2,2,1);
    PlotLine(energy*0.95,energy*0.95,rho+err-0.001,rho+err,kAzure+2,2,1);
    PlotLine(energy*1.05,energy*1.05,rho+err-0.001,rho+err,kAzure+2,2,1);
    PlotLine(energy*0.95,energy*1.05,rho-err,rho-err,kAzure+2,2,1);
    PlotLine(energy*0.95,energy*0.95,rho-err+0.001,rho-err,kAzure+2,2,1);
    PlotLine(energy*1.05,energy*1.05,rho-err+0.001,rho-err,kAzure+2,2,1);
  }
  // Draw_TGAE_Point_new_Symbol(80,0.39,0.0,0.0,0.0,0.0,20,kAzure+2,1.4);
  // // plotTopLegend((char*)"2^{nd} EP Tangent",90,0.3890,font_size,1,0.0,42,0);
  // plotTopLegend((char*)"In-Plane #rho_{00}",90,0.3890,font_size,1,0.0,42,0);

  // plot theory curve
  double scurrent[3] = {400.0,600.0,1000.0};
  double ms = 450.0;
  int Style[3] = {kDashed,kDotted,kSolid};
  int Color[3] = {kBlue,kMagenta,kRed};
  TF1 *f_rho00_theory[3];
  for(int i_line = 0; i_line <3; ++i_line)
  {
    string FuncName = Form("f_rho00_theory_%d",i_line);
    f_rho00_theory[i_line] = new TF1(FuncName.c_str(),rho00_theory,20,200,2);
    f_rho00_theory[i_line]->FixParameter(0,scurrent[i_line]);
    f_rho00_theory[i_line]->FixParameter(1,ms);
    f_rho00_theory[i_line]->SetLineColor(Color[i_line]);
    f_rho00_theory[i_line]->SetLineStyle(Style[i_line]);
    f_rho00_theory[i_line]->SetLineWidth(4);
    f_rho00_theory[i_line]->Draw("l same");
  }

  // plot legend

  // normal
  Draw_TGAE_Point_new_Symbol(45,0.385,0.0,0.0,0.0,0.0,29,kRed,1.8);
  plotTopLegend((char*)"Out-of-Plane (2^{nd} EP)",50,0.3835,font_size+0.005,1,0.0,42,0);

  // tangent
  Draw_TGAE_Point_new_Symbol(45,0.375,0.0,0.0,0.0,0.0,20,kAzure+2,1.4);
  plotTopLegend((char*)"In-Plane (2^{nd} EP)",50,0.3735,font_size+0.005,1,0.0,42,0);

  // theory
  for(int i_line = 0; i_line <3; ++i_line)
  {
    string leg_scurrent = Form("g^{4}_{#phi}(#partialj^{(i)}_{s}/#partialt)^{2} = %1.0f fm^{-4} [*]",scurrent[i_line]);
    plotTopLegend((char*)leg_scurrent.c_str(),50,0.315-i_line*0.009,font_size,1,0.0,42,0);
    PlotLine(35,48,0.3157-i_line*0.009,0.3157-i_line*0.009,Color[i_line],4,Style[i_line]);
  }

  // experimental information
  // plotTopLegend((char*)"Au+Au (20-60\% & |#eta| < 1)",15,0.3890,font_size,1,0.0,42,0);
  // plotTopLegend((char*)"#phi-meson (1.2 < p_{T}< 5.4 GeV/c)",14,0.3790,font_size,1,0.0,42,0);
  plotTopLegend((char*)"Au+Au 20-60\%",10,0.392,font_size,1,0.0,42,0);
  // plotTopLegend((char*)"#phi (1.2 < p_{T}< 5.4 GeV/c & |#eta| < 1)",30,0.392,font_size,1,0.0,42,0);
  plotTopLegend((char*)"1.2 < p_{T}< 5.4 GeV/c",30,0.392,font_size,1,0.0,42,0);
  plotTopLegend((char*)"|#eta| < 1",140,0.392,font_size,1,0.0,42,0);

  // plotTopLegend((char*)"#rho_{00} = 1/3",100,0.328,0.04,1,0.0,42,0);

  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/BESII/c_rhoSys_energy_Theory_2ndEP_BESII.eps");
}

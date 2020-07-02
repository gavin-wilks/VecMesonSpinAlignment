#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../../Utility/draw.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);

double rho00_theory(double *x_var, double *par)
{
  double s12 = x_var[0];
  double scurrent = par[0]; // strangeness current
  double ms = par[1]; // mass of s-quark MeV

  double mphi = 1020.0; // MeV
  double c1 = 300.0/(pow(200.0,1.0/3.0) * pow((-0.4+0.39*log(200.0*200.0)),1.0/3.0)); // 300.0/Teff[200.0]
  double gphi = 2.0*sqrt(2.0); // included in scurrent

  double Teff = pow(s12,1.0/3.0) * pow((-0.4+0.39*log(s12*s12)),1.0/3.0);

  double denom_phifield = 27.0*pow(ms,4.0)*pow(mphi,4.0)*pow(c1,2.0)*pow(Teff,2.0);
  double numer_phifield = scurrent*pow(197.0,8.0)*1.8*1.0e+5; // <p^2>_phi = 0.18 GeV^2

  double rho00 = 1.0/3.0 + numer_phifield/denom_phifield;

  return rho00;
}

void plotFig6_Rho00EnergyEP()
{
  gStyle->SetOptDate(0);
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed;
  const int style_phi_EP  = 33;
  const int color_phi_EP  = 9;

  const float size_marker = 1.4;
  const float size_font = 0.04;

  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat_Laxis = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys_Laxis  = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_energy_sys");
  TGraphAsymmErrors *g_rhoPhi_2nd_fit_Laxis  = new TGraphAsymmErrors();
  for(int i_energy = 0; i_energy < g_rhoPhi_2nd_stat_Laxis->GetN(); ++i_energy) // combine stat & sys for fit
  {
    double energy, rho;
    g_rhoPhi_2nd_stat_Laxis->GetPoint(i_energy,energy,rho);
    double err_stat = g_rhoPhi_2nd_stat_Laxis->GetErrorYhigh(i_energy);
    double err_sys = g_rhoPhi_2nd_sys_Laxis->GetErrorYhigh(i_energy);
    double err_fit = TMath::Sqrt(err_stat*err_stat+err_sys*err_sys);

    g_rhoPhi_2nd_fit_Laxis->SetPoint(i_energy,energy,rho);
    g_rhoPhi_2nd_fit_Laxis->SetPointError(i_energy,0.0,0.0,err_fit,err_fit);
  }

  TFile *File_Input_EP = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rho00_stat_sys_Xaxis.root");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat_Xaxis = (TGraphAsymmErrors*)File_Input_EP->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys_Xaxis = (TGraphAsymmErrors*)File_Input_EP->Get("rho00_2ndEP_energy_sys");
  TGraphAsymmErrors *g_rhoPhi_2nd_fit_Xaxis  = new TGraphAsymmErrors();
  for(int i_energy = 0; i_energy < g_rhoPhi_2nd_stat_Xaxis->GetN(); ++i_energy) // combine stat & sys for fit
  {
    double energy, rho;
    g_rhoPhi_2nd_stat_Xaxis->GetPoint(i_energy,energy,rho);
    double err_stat = g_rhoPhi_2nd_stat_Xaxis->GetErrorYhigh(i_energy);
    double err_sys = g_rhoPhi_2nd_sys_Xaxis->GetErrorYhigh(i_energy);
    double err_fit = TMath::Sqrt(err_stat*err_stat+err_sys*err_sys);

    g_rhoPhi_2nd_fit_Xaxis->SetPoint(i_energy,energy,rho);
    g_rhoPhi_2nd_fit_Xaxis->SetPointError(i_energy,0.0,0.0,err_fit,err_fit);
  }

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

  h_frame->GetYaxis()->SetRangeUser(0.291,0.42);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(9.0,240.0,1.0/3.0,1.0/3.0,1,3,2);

  // plot 2nd EP rho00 in normal direction
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_2nd_stat_Laxis,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotSysErrors(g_rhoPhi_2nd_sys_Laxis,color_phi_2nd);

  // plot 2nd EP rho00 in tangent direction
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_2nd_stat_Xaxis,style_phi_EP,color_phi_EP,size_marker+0.4);
  plotSysErrors(g_rhoPhi_2nd_sys_Xaxis,color_phi_EP);

  // plot theory curve
  double scurrent[3] = {400.0,600.0,1000.0};
  double ms = 450.0;
  int Style[3] = {kDashed,kDotted,kSolid};
  int Color[3] = {kBlue,kMagenta,kRed};
  TF1 *f_rho00_Laxis = new TF1("f_rho00_Laxis",rho00_theory,1,201,2);
  f_rho00_Laxis->SetParameter(0,1000.0);
  f_rho00_Laxis->FixParameter(1,ms);
  f_rho00_Laxis->SetLineColor(kRed);
  f_rho00_Laxis->SetLineStyle(kDashed);
  f_rho00_Laxis->SetLineWidth(4);
  f_rho00_Laxis->SetRange(19.0,200.0);
  g_rhoPhi_2nd_fit_Laxis->Fit(f_rho00_Laxis,"MNR");
  f_rho00_Laxis->Draw("l same");
  double chi2_Laxis = f_rho00_Laxis->GetChisquare();
  int ndf_Laxis = f_rho00_Laxis->GetNDF();
  double chi2_ndf_Laxis = chi2_Laxis/(double)ndf_Laxis;
  double p_Laxis = TMath::Prob(chi2_Laxis,ndf_Laxis);
  cout << "chi2 for Laxis: " << chi2_Laxis << endl;
  cout << "ndf for Laxis: " << ndf_Laxis << endl;
  cout << "chi2/ndf for Laxis: " << chi2_ndf_Laxis  << ", p_Laxis: " << p_Laxis << endl;
  cout << "C^{y}_{s} = " << f_rho00_Laxis->GetParameter(0) << " +/- " << f_rho00_Laxis->GetParError(0) << endl;

  TF1 *f_rho00_Xaxis = new TF1("f_rho00_Xaxis",rho00_theory,1,201,2);
  f_rho00_Xaxis->SetParameter(0,200.0);
  f_rho00_Xaxis->FixParameter(1,ms);
  f_rho00_Xaxis->SetLineColor(kBlue);
  f_rho00_Xaxis->SetLineStyle(kDashed);
  f_rho00_Xaxis->SetLineWidth(4);
  f_rho00_Xaxis->SetRange(19.0,200.0);
  g_rhoPhi_2nd_fit_Xaxis->Fit(f_rho00_Xaxis,"MNR");
  f_rho00_Xaxis->Draw("l same");
  double chi2_Xaxis = f_rho00_Xaxis->GetChisquare();
  int ndf_Xaxis = f_rho00_Xaxis->GetNDF();
  double chi2_ndf_Xaxis = chi2_Xaxis/(double)ndf_Xaxis;
  double p_Xaxis = TMath::Prob(chi2_Xaxis,ndf_Xaxis);
  cout << "chi2 for Xaxis: " << chi2_Xaxis << endl;
  cout << "ndf for Xaxis: " << ndf_Xaxis << endl;
  cout << "chi2/ndf for Xaxis: " << chi2_ndf_Xaxis << ", p_Xaxis: " << p_Xaxis << endl;
  cout << "C^{y}_{s} = " << f_rho00_Xaxis->GetParameter(0) << " +/- " << f_rho00_Xaxis->GetParError(0) << endl;

  // plot legend

  plotTopLegend((char*)"Au+Au 20-60%",10,0.41,size_font,1,0.0,42,0);
  plotTopLegend((char*)"|y| < 1 & 1.2 < p_{T}< 5.4 GeV/c",30,0.41,size_font,1,0.0,42,0);

  // normal
  Draw_TGAE_Point_new_Symbol(20,0.40,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotTopLegend((char*)"#phi (Out-of-Plane & 2^{nd}-Order EP)",23,0.3985,size_font,1,0.0,42,0);

  // tangent
  Draw_TGAE_Point_new_Symbol(20,0.39,0.0,0.0,0.0,0.0,style_phi_EP,color_phi_EP,size_marker+0.4);
  plotTopLegend((char*)"#phi (In-Plane & 2^{nd}-Order EP)",23,0.3885,size_font,1,0.0,42,0);

  // theory
  string leg_current_Laxis = Form("C^{(y)}_{s} = %1.0f #pm %1.0f fm^{-8}",f_rho00_Laxis->GetParameter(0),f_rho00_Laxis->GetParError(0));
  string leg_chi2_Laxis = Form("#chi^{2}/ndf: %1.1f",chi2_ndf_Laxis);
  string leg_p_Laxis = Form("p-value: %1.3f", p_Laxis);
  string leg_stat_Laxis = Form("#chi^{2}/ndf: %1.1f & p-value: %1.3f",chi2_ndf_Laxis,p_Laxis);

  string leg_current_Xaxis = Form("C^{(y)}_{s} = %1.0f #pm %1.0f fm^{-8}",f_rho00_Xaxis->GetParameter(0),f_rho00_Xaxis->GetParError(0));
  string leg_chi2_Xaxis = Form("#chi^{2}/ndf: %1.1f",chi2_ndf_Xaxis);
  string leg_p_Xaxis = Form("p-value: %1.3f", p_Xaxis);
  string leg_stat_Xaxis = Form("#chi^{2}/ndf: %1.1f & p-value: %1.3f",chi2_ndf_Xaxis,p_Xaxis);

  TLegend *leg = new TLegend(0.45,0.17,0.85,0.32);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(f_rho00_Laxis,leg_current_Laxis.c_str(),"l");
  // leg->AddEntry((TObject*)0,leg_stat_Laxis.c_str(),"");
  // leg->AddEntry((TObject*)0,leg_chi2_Laxis.c_str(),"");
  // leg->AddEntry((TObject*)0,leg_p_Laxis.c_str(),"");
  leg->AddEntry(f_rho00_Xaxis,leg_current_Xaxis.c_str(),"l");
  // leg->AddEntry((TObject*)0,leg_stat_Xaxis.c_str(),"");
  // leg->AddEntry((TObject*)0,leg_chi2_Xaxis.c_str(),"");
  // leg->AddEntry((TObject*)0,leg_p_Xaxis.c_str(),"");
  leg->Draw("same");

  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig6_rho00EnergyEP.eps");
  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig6_rho00EnergyEP.png");
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

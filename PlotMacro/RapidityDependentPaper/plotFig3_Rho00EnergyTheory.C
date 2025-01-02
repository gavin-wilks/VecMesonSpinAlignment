#include <string>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TStyle.h>
#include <TF1.h>

#include "draw.h"
//#include "../../../Utility/StSpinAlignmentCons.h"

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


void plotFig3_Rho00EnergyTheory()
{
  gStyle->SetOptDate(0);
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kGray+1;//kRed-4;
  const int style_phi_ALICE = 30;
  const int color_phi_ALICE = kGray+1;
  const int colorDiff_phi = 0;
  const int cphiII = kRed-4;
  const int stylephiII = 29;

  const int style_Kstr = 20;
  const int color_Kstr = kGray+1;//kAzure-9;
  const int style_Kstr_ALICE = 24;
  const int color_Kstr_ALICE = kGray+1;
  const int colorDiff_Kstr = 2;

  const float size_marker = 1.4;
  const float size_font = 0.025;//0.035-0.005;
  
  //----------------------------------------------------------
  // phi-meson STAR
  //beam-energy dependence of phi-meson rho00 from STAR, pT: 1.2 - 5.4 GeV/c, 20-60%
  TFile *File_InputPhi = TFile::Open("../VectorMesonSpinAlignment/VecMesonSpinAlignment/data/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat_Laxis = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys_Laxis = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_sys");
  TGraphAsymmErrors *g_rho_2nd_fit_Laxis  = new TGraphAsymmErrors();
  for(int i_energy = 0; i_energy < g_rhoPhi_2nd_stat_Laxis->GetN(); ++i_energy) // combine stat & sys for fit
  {
    double energy, rho;
    g_rhoPhi_2nd_stat_Laxis->GetPoint(i_energy,energy,rho);
    double err_stat = g_rhoPhi_2nd_stat_Laxis->GetErrorYhigh(i_energy);
    double err_sys = g_rhoPhi_2nd_sys_Laxis->GetErrorYhigh(i_energy);
    double err_fit = TMath::Sqrt(err_stat*err_stat+err_sys*err_sys);

    g_rho_2nd_fit_Laxis->SetPoint(i_energy,energy,rho);
    g_rho_2nd_fit_Laxis->SetPointError(i_energy,0.0,0.0,err_fit,err_fit);
  }

  // TFile *File_InputPhi_EP = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/BESII/rho00_stat_sys_Xaxis.root");
  // TGraphAsymmErrors *g_rhoPhi_2nd_stat_Xaxis = (TGraphAsymmErrors*)File_InputPhi_EP->Get("rho00_2ndEP_energy_stat");
  // TGraphAsymmErrors *g_rhoPhi_2nd_sys_Xaxis = (TGraphAsymmErrors*)File_InputPhi_EP->Get("rho00_2ndEP_energy_sys");
  //----------------------------------------------------------

  //----------------------------------------------------------
  TFile *File_InputKstar = TFile::Open("./data_Kstar_rho00_sNN_June9_2021.root");
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
  h_frame->GetXaxis()->SetRangeUser(2.0,4996.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.28,0.4);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  // h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(2.0,4996.0,1.0/3.0,1.0/3.0,1,2,2);

  // K* STAR
  //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKstar_stat,style_Kstr,color_Kstr,colorDiff_Kstr,size_marker-0.4);
  // plotSysErrors(g_rhoKstar_sys,color_Kstr+2);
  //plotSysErrorsBox(g_rhoKstar_sys,color_Kstr+2);

  // phi-meson STAR BESII
  //TGraphAsymmErrors *g_rhoKStar_stat_BESII = new TGraphAsymmErrors();
  //TGraphAsymmErrors *g_rhoKStar_sys_BESII = new TGraphAsymmErrors();
  //double energy, pt, val, stat, sys;
  //g_rhoKstar_stat->GetPoint(2,energy,val);
  //stat = g_rhoKstar_stat->GetErrorYhigh(2);
  //sys  = g_rhoKstar_sys->GetErrorYhigh(2);
  ////g_rhoKStar_stat_BESII->SetPoint(0,energy,1./3.);
  ////g_rhoKStar_stat_BESII->SetPointError(0,0.0,0.0,stat*TMath::Sqrt(36./460.),stat*TMath::Sqrt(36./460.));
  ////g_rhoKStar_sys_BESII->SetPoint(0,energy,1./3.);
  ////g_rhoKStar_sys_BESII->SetPointError(0,0.0,0.0,sys*TMath::Sqrt(36./460.),sys*TMath::Sqrt(36./460.));
  ////Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKStar_stat_BESII,style_Kstr,kGreen,colorDiff_Kstr,size_marker-0.4);
  //// plotSysErrors(g_rhoPhi_2nd_sys_Laxis,color_phi_2nd);
  //plotSysErrorsBox(g_rhoKStar_sys_BESII,color_Kstr+2);

  //int boxcolor = kAzure;
  //TBox *bK19;
  //{
  //  double x1BesII = (19.6 - 2.0)/1.08;
  //  double y1BesII = 1.0/3.0 - stat*TMath::Sqrt(36./478.);
  //  double x2BesII = (19.6 - 2.0)*1.08;
  //  double y2BesII = 1.0/3.0 + stat*TMath::Sqrt(36./478.);
  //
  //  cout << "yerr = " << stat*TMath::Sqrt(36./478.) << endl;
  // 
  //  bK19 = new TBox(x1BesII,y1BesII,x2BesII,y2BesII);
  //  bK19->SetFillColor(boxcolor);
  //  bK19->SetFillColorAlpha(boxcolor,0.75);
  //  bK19->SetFillStyle(3001);
  //  bK19->SetLineStyle(2);
  //  bK19->SetLineColor(boxcolor);
  //  bK19->SetLineWidth(1);
  //  bK19->Draw("l Same");
  //}

  // K* ALICE
  //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKstar_ALICE_stat,style_Kstr_ALICE,color_Kstr_ALICE,0,size_marker-0.4);
  // plotSysErrors(g_rhoKstar_ALICE_sys,color_Kstr_ALICE);

  // phi-meson STAR BESI
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_2nd_stat_Laxis,style_phi_2nd,color_phi_2nd,colorDiff_phi,size_marker+0.2);
  // plotSysErrors(g_rhoPhi_2nd_sys_Laxis,color_phi_2nd);
  plotSysErrorsBox(g_rhoPhi_2nd_sys_Laxis,color_phi_2nd);

  // phi-meson STAR BESII
  TGraphAsymmErrors *g_rhoPhi_stat_BESII = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_rhoPhi_sys_BESII = new TGraphAsymmErrors();
 // g_rhoPhi_stat_BESII->SetPoint(0,19.6,0.3622);
 // g_rhoPhi_stat_BESII->SetPointError(0,0.0,0.0,0.0026,0.0026);
 // g_rhoPhi_sys_BESII->SetPoint(0,19.6,0.3622);
 // g_rhoPhi_sys_BESII->SetPointError(0,0.0,0.0,0.0049,0.0049);
  g_rhoPhi_stat_BESII->SetPoint(0,19.6,0.3510);
  g_rhoPhi_stat_BESII->SetPointError(0,0.0,0.0,0.0023,0.0023);
  g_rhoPhi_sys_BESII->SetPoint(0,19.6,0.3510);
  g_rhoPhi_sys_BESII->SetPointError(0,0.0,0.0,0.0013,0.0013);
  g_rhoPhi_stat_BESII->SetPoint(1,14.6,0.3467);
  g_rhoPhi_stat_BESII->SetPointError(1,0.0,0.0,0.0037,0.0037);
  g_rhoPhi_sys_BESII->SetPoint(1,14.6,0.3467);
  g_rhoPhi_sys_BESII->SetPointError(1,0.0,0.0,0.0018,0.0018);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_stat_BESII,stylephiII,cphiII,colorDiff_phi,size_marker+0.2);
  // plotSysErrors(g_rhoPhi_2nd_sys_Laxis,color_phi_2nd);
  plotSysErrorsBox(g_rhoPhi_sys_BESII,cphiII);

  // phi-meson ALICE
  //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_ALICE_stat,style_phi_ALICE,color_phi_ALICE,0,size_marker);
  // plotSysErrors(g_rhoPhi_ALICE_sys,color_phi_ALICE);

  // plot theory curve
  double ms = 450.0;
  int Style[3] = {kDashed,kDotted,kSolid};
  int Color[3] = {kBlue,kMagenta,kRed};
  TF1 *f_rho00_Laxis = new TF1("f_rho00_Laxis",rho00_theory,1,201,2);
  f_rho00_Laxis->SetParameter(0,1000.0);
  f_rho00_Laxis->FixParameter(1,ms);
  f_rho00_Laxis->SetRange(19.0,200.0);
  g_rho_2nd_fit_Laxis->Fit(f_rho00_Laxis,"MNR");
  double chi2_Laxis = f_rho00_Laxis->GetChisquare();
  int ndf_Laxis = f_rho00_Laxis->GetNDF();
  double chi2_ndf_Laxis = chi2_Laxis/(double)ndf_Laxis;
  double p_Laxis = TMath::Prob(chi2_Laxis,ndf_Laxis);
  cout << "chi2 for Laxis: " << chi2_Laxis << endl;
  cout << "ndf for Laxis: " << ndf_Laxis << endl;
  cout << "chi2/ndf for Laxis: " << chi2_ndf_Laxis  << ", p_Laxis: " << p_Laxis << endl;
  cout << "C^{y}_{s} = " << f_rho00_Laxis->GetParameter(0) << " +/- " << f_rho00_Laxis->GetParError(0) << endl;

  f_rho00_Laxis->SetLineColor(kBlack);
  f_rho00_Laxis->SetLineStyle(1);
  f_rho00_Laxis->SetLineWidth(4);
  f_rho00_Laxis->Draw("l same");

  TF1 *f_rho00_plot = new TF1("f_rho00_plot",rho00_theory,1,4000,2);
  f_rho00_plot->FixParameter(0,f_rho00_Laxis->GetParameter(0));
  f_rho00_plot->FixParameter(1,ms);
  f_rho00_plot->SetLineStyle(2);
  f_rho00_plot->SetLineWidth(4);
  f_rho00_plot->SetLineColor(kBlack);
  f_rho00_plot->SetRange(200.0,3000.0);
  f_rho00_plot->Draw("l same");

  // plot Legend
  //plotTopLegend((char*)"Au+Au",0.34,0.85,size_font,1,0.0,42,1);
  //plotTopLegend((char*)"20\% - 60\% Centrality",0.25,0.82,size_font,1,0.0,42,1);
  // plotTopLegend((char*)"Pb+Pb",0.72,0.85,size_font,1,0.0,42,1);
  // plotTopLegend((char*)"10\% - 50\% Centrality",0.63,0.82,size_font,1,0.0,42,1);
  plotTopLegend((char*)"Au+Au (20\% - 60\% Centrality)",10,0.29,size_font,1,0.0,42,0);
  //plotTopLegend((char*)"open: Pb+Pb (10\% - 50\% Centrality)",0.20,0.20,size_font,1,0.0,42,1);

  Draw_TGAE_Point_new_Symbol(40,0.379,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotTopLegend((char*)"#phi BES-I (|y| < 1.0 & 1.2 < p_{T} < 5.4 GeV/c)^{1}",48,0.378,size_font,1,0.0,42,0);

  //Draw_TGAE_Point_new_Symbol(40,0.395,0.0,0.0,0.0,0.0,style_Kstr,color_Kstr,size_marker-0.4);
  //plotTopLegend((char*)"K^{*0} (|y| < 1.0 & 1.0 < p_{T} < 5.0 GeV/c)",48,0.3925,size_font,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(40,0.371,0.0,0.0,0.0,0.0,stylephiII,cphiII,size_marker+0.2);
  plotTopLegend((char*)"#phi BES-II (|y| < 1.0 & 1.2 < p_{T} < 5.4 GeV/c)",48,0.37,size_font,1,0.0,42,0);

  //plotTopLegend((char*)"K^{*0} BES-II Stat. Error Projection",120,0.2925,size_font,1,0.0,42,0);

  plotTopLegend((char*)"STAR Preliminary",10,0.39,size_font,kRed,0.0,42,0);

  //TBox *bK19l;
  //{
  //  double x1BesII = 100/1.08;
  //  double y1BesII = 0.293;
  //  double x2BesII = 100*1.08;
  //  double y2BesII = 0.298;
  //  bK19l = new TBox(x1BesII,y1BesII,x2BesII,y2BesII);
  //  bK19l->SetFillColor(boxcolor);
  //  bK19l->SetFillColorAlpha(boxcolor,0.75);
  //  bK19l->SetFillStyle(3001);
  //  bK19l->SetLineStyle(2);
  //  bK19l->SetLineColor(boxcolor);
  //  bK19l->SetLineWidth(1);
  //  bK19l->Draw("l Same");
  //}

  // theory
  //string leg_current_Laxis = Form("G^{(y)}_{s} = %1.0f #pm %1.0f fm^{-8} (Fit to BES-I only)",f_rho00_Laxis->GetParameter(0),f_rho00_Laxis->GetParError(0));
  string leg_current_Laxis = Form("G^{(y)}_{s} = %1.2f #pm %1.2f m^{4}_{#pi} (Fit to BES-I only)",4.64,0.73);
  // string leg_chi2_Laxis = Form("#chi^{2}/ndf: %1.1f",chi2_ndf_Laxis);
  // string leg_p_Laxis = Form("p-value: %1.3f", p_Laxis);
  // string leg_stat_Laxis = Form("#chi^{2}/ndf: %1.1f & p-value: %1.3f",chi2_ndf_Laxis,p_Laxis);
  // TLegend *leg = new TLegend(0.335,0.18,0.675,0.23);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(10);
  // leg->SetFillStyle(0);
  // leg->AddEntry(f_rho00_Laxis,leg_current_Laxis.c_str(),"l");
  // leg->AddEntry((TObject*)0,leg_stat_Laxis.c_str(),"");
  // leg->Draw("same");
  PlotLine(35,45,0.363,0.363,kBlack,3,1);
  plotTopLegend((char*)leg_current_Laxis.c_str(),48,0.362,size_font,1,0.0,42,0);

  c_rho00->SaveAs("./rho00EnergyTheory.eps");
  c_rho00->SaveAs("./rho00EnergyTheory.png");
  c_rho00->SaveAs("./rho00EnergyTheory.pdf");
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

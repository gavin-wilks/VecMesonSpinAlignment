#include <string>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TStyle.h>
#include <TF1.h>

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

void genRho00MuB_model()
{
  const float beamEnergy[5] = {19.6, 27.0, 39.0, 62.4, 200.0};
  const float muB[5]        = { 205,   155,  115,   70,   20};
  
  //----------------------------------------------------------
  // phi-meson STAR
  //beam-energy dependence of phi-meson rho00 from STAR, pT: 1.2 - 5.4 GeV/c, 20-60%
  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/BESII/rho00_stat_sys_Laxis.root");
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
  //----------------------------------------------------------

  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  c_rho00->cd()->SetLogx();

  TH1F *h_frame = new TH1F("h_frame","h_frame",500,0,500);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(7.0,300.0);
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

  g_rho_2nd_fit_Laxis->Draw("pE same");

  TF1 *f_rho00_Laxis = new TF1("f_rho00_Laxis",rho00_theory,1,201,2);
  f_rho00_Laxis->SetParameter(0,1000.0);
  f_rho00_Laxis->FixParameter(1,450.0);
  f_rho00_Laxis->SetRange(19.0,200.0);
  g_rho_2nd_fit_Laxis->Fit(f_rho00_Laxis,"MNR");

  f_rho00_Laxis->SetLineColor(kRed);
  f_rho00_Laxis->SetLineStyle(1);
  f_rho00_Laxis->SetLineWidth(4);
  f_rho00_Laxis->Draw("l same");

  TGraphAsymmErrors *g_rhoEnergyModel = new TGraphAsymmErrors();
  g_rhoEnergyModel->SetName("g_rhoEnergyModel");
  TGraphAsymmErrors *g_rhoMuBModel = new TGraphAsymmErrors();
  g_rhoMuBModel->SetName("g_rhoMuBModel");
  for(int iEnergy = 0; iEnergy < 5; ++iEnergy)
  {
    g_rhoEnergyModel->SetPoint(iEnergy,beamEnergy[iEnergy],f_rho00_Laxis->Eval(beamEnergy[iEnergy]));
    g_rhoMuBModel->SetPoint(iEnergy,muB[iEnergy],f_rho00_Laxis->Eval(beamEnergy[iEnergy]));
  }
  g_rhoEnergyModel->SetMarkerStyle(24);
  g_rhoEnergyModel->SetMarkerSize(1.0);
  g_rhoEnergyModel->SetMarkerColor(2);
  g_rhoEnergyModel->Draw("p same");

  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/TestMuB/VecMeson/rho00MuB_model.root","RECREATE");
  File_OutPut->cd();
  g_rhoEnergyModel->Write();
  g_rhoMuBModel->Write();
  File_OutPut->Close();
}

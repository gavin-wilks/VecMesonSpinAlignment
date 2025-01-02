#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "Utility/type.h"
#include <string>
//#include "SpecFuncMathMore.h"
//#include "Math.h"
#include <cmath>
#include "Math/SpecFuncMathMore.h"
#include <iostream>
#include <fstream>


using namespace std;

void plotQA_EPDependentEfficiency_Cent(const int energy = 4, const int pid = 0, int etamode = 0, int ypadding = 2, int order = 2, int yspectra = 0) {

  std::string spectra;
  if(yspectra) spectra = "_WithRapiditySpectra";  
  if(!yspectra) spectra = "_NoRapiditySpectra";  

  std::string etastring;
  if(etamode == 0) etastring = "eta1_eta1";
  if(etamode == 1) etastring = "eta1_eta1p5";
  if(etamode == 2) etastring = "eta1p5_eta1p5";
  if(etamode == 3) etastring = "eta0p4";
  if(etamode == 4) etastring = "eta0p6";
  if(etamode == 5) etastring = "eta0p8";

  gROOT->Reset();
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);

  Int_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};

  TH2D *eff[10][10];
  //TFile *eff_file = new TFile(Form("../RcPhiEffCorr/Rebinned/Eff_20GeV_SingleParticle_noToF_Mode2_EtaMode%d.root",etamode),"READ");
  //TFile *eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_NoEPResAndAcceptance/Eff_19GeV_SingleParticle_noToF_Mode2_EtaMode%d.root",etamode),"READ");
  TFile *eff_file;
  if(order == 1) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order1%s_EP/Eff_19GeV_SingleParticle_noToF_Mode1_EtaMode%d.root",spectra.c_str(),etamode),"READ");
  if(order == 2) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2%s_EP/Eff_19GeV_SingleParticle_noToF_Mode1_EtaMode%d.root",spectra.c_str(),etamode),"READ");
  if(energy == 3 && order == 1) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu14GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order1%s_EP/Eff_14GeV_SingleParticle_noToF_Mode1_EtaMode%d.root",spectra.c_str(),etamode),"READ");
  if(energy == 3 && order == 2) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu14GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2%s_EP/Eff_14GeV_SingleParticle_noToF_Mode1_EtaMode%d.root",spectra.c_str(),etamode),"READ");
  eff_file->Print();

  for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
  {
    for(int icent = 0; icent < 9; icent++)
    {
      std::string HistName;
      HistName = Form("h_mEffCosEP_Cent_%d_Pt_%d",icent,ipt);
      cout << HistName << endl;
      eff[ipt][icent] = (TH2D*)eff_file->Get(HistName.c_str());
      eff[ipt][icent]->Print();
    }
  }

  TCanvas *c_fit = new TCanvas("c_fit","c_fit",10,10,900,900);
  c_fit->Divide(3,3);
  for(int i = 0; i < 9; i++)
  {
    c_fit->cd(i+1)->SetLeftMargin(0.15);
    c_fit->cd(i+1)->SetRightMargin(0.2);
    c_fit->cd(i+1)->SetBottomMargin(0.15);
    c_fit->cd(i+1)->SetGrid(0,0);
    c_fit->cd(i+1)->SetTicks(1,1);
  }

  std::string outputname = Form("figures/%s/%s/centstudy_EPefficiency_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order);
  std::string output_start = Form("%s[",outputname.c_str());

  c_fit->Print(output_start.c_str());

  for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
  {
    for(int icent = 0; icent < 9; icent++)
    {
      c_fit->cd(icent+1);
      //eff[ipt][icent]->CenterTitle();
      eff[ipt][icent]->SetTitle(Form("%1.1f<p_{T}<%1.1f, Cent %d-%d",vmsa::pt_low_cent[energy][ipt],vmsa::pt_up_cent[energy][ipt],cent_set[icent+1],cent_set[icent]));
      cout << Form("%1.1f<p_{T}<%1.1f, Cent %d-%d",vmsa::pt_low_cent[energy][ipt],vmsa::pt_up_cent[energy][ipt],cent_set[icent+1],cent_set[icent]) << endl;
      eff[ipt][icent]->GetXaxis()->SetTitle("|cos#theta*|");
      eff[ipt][icent]->GetXaxis()->CenterTitle();
      eff[ipt][icent]->GetYaxis()->SetTitle("#phi-#Psi_{2}");
      eff[ipt][icent]->GetYaxis()->CenterTitle();
      eff[ipt][icent]->Draw("colz");
    }
    c_fit->Update();
    c_fit->Print(outputname.c_str());
  }
  std::string output_stop = Form("%s]",outputname.c_str());
  c_fit->Print(output_stop.c_str()); // close pdf file 
  
}

void Correction(Double_t Res, Double_t Res_error, Double_t obv_rho, Double_t obv_rho_error, Double_t &real_rho, Double_t &real_rho_error) {

  real_rho = (Res-1+4*obv_rho)/(1+3*Res);

  Double_t real_rho_error_1 = obv_rho_error*4/(1+3*Res);
  Double_t real_rho_error_2 = Res_error*(4-12*obv_rho)/(1+3*Res)/(1+3*Res);
  real_rho_error = TMath::Sqrt(real_rho_error_1*real_rho_error_1+real_rho_error_2*real_rho_error_2);

}

Double_t chi(Double_t res) {

  Double_t chi = 2.0;
  Double_t delta = 1.0;

  for(Int_t i=0; i<15; i++) {
    chi = (res1(chi) < res) ? chi + delta : chi - delta;
    delta = delta/2.;
  }

  return chi;
}

Double_t res1(Double_t chi) {

  Double_t con = 0.626657;                   // TMath::Sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;

  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

  return res;
}

Double_t resEventPlane(Double_t chi) {

  Double_t con = 0.626657;                   // TMath::Sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;
  Double_t halfpi = TMath::Pi()/2.;

  Double_t besselOneHalf = TMath::Sqrt(arg/halfpi) * TMath::SinH(arg)/arg;
  Double_t besselThreeHalfs = TMath::Sqrt(arg/halfpi) * (TMath::CosH(arg)/arg - TMath::SinH(arg)/(arg*arg));
  Double_t res = con * chi * exp(-arg) * (besselOneHalf + besselThreeHalfs);

  return res;
}

Double_t resEventPlane1(Double_t chi) {

  Double_t con = 0.626657;                   // TMath::Sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;
  Double_t halfpi = TMath::Pi()/2.;

  TF1* iBesselHalf = new TF1("I_12", "ROOT::Math::cyl_bessel_i([0],x)", 0, 10);
  iBesselHalf->SetParameter(0,1./2.);
  TF1* iBessel3Half = new TF1("I_32", "ROOT::Math::cyl_bessel_i([0],x)", 0, 10);
  iBessel3Half->SetParameter(0,3./2.);
  Double_t besselOneHalf = iBesselHalf->Eval(arg);
  Double_t besselThreeHalfs = iBessel3Half->Eval(arg);
  Double_t res = con * chi * exp(-arg) * (besselOneHalf + besselThreeHalfs);

  return res;
}

double FuncAD(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double D = par[2];
  double R = par[3];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double result = (1.+Bs*D/2.) + (As+D)*CosTheta*CosTheta + (As*D-Bs*D/2.)*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}

double Func4th(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double F = par[1];
  double G = par[2];

  //double result = 1. + (4.*F+3.*G)/8. - (2.*F+3.*G)/4.*CosTheta*CosTheta + 3.*G/8.*CosTheta*CosTheta*CosTheta*CosTheta;
   
  double order0 = 2. + F + 3.*G/4.;
  double order2 = (-1.*F - 3.*G/2.)*CosTheta*CosTheta;
  double order4 = (3.*G/4.)*CosTheta*CosTheta*CosTheta*CosTheta;
 
  double result = order0 + order2 + order4;

  return N*result;

}

double FuncAFG(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double R = par[4];
  //double rho10 = par[5];

  double A = (3.*rho-1.)/(1.-rho);
  //double B = rho10/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));
  //double denom = 4. + A*(1.-R) + B*(-1.+R);
  //double As = (A*(1.+3.*R) + B*(3.-3.*R))/denom;
  //double Bs = (A*(1.-R)    + B*(3.+R)   )/denom;


  double order0 = 2. + F - Bs*F/2. + 3.*G/4. - Bs*G/2.;    
  double order2 = (2.*As - F + As*F + Bs*F - 3.*G/2. + 3.*As*G/4. + 3.*Bs*G/2.)*CosTheta*CosTheta;
  double order4 = (-1.*As*F - Bs*F/2. + 3.*G/4. - 3.*As*G/2. - 3.*Bs*G/2.)*CosTheta*CosTheta*CosTheta*CosTheta;
  double order6 = (3.*As*G/4. + Bs*G/2.)*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;

  double result = order0 + order2 + order4 + order6;

  return N*result;

}

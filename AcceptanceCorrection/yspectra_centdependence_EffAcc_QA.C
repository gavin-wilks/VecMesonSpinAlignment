#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "Utility/type.h"
#include <string>
//#include "SpecFuncMathMore.h"
//#include "Math.h"
#include <cmath>
#include "Math/SpecFuncMathMore.h"

using namespace std;

void yspectra_centdependence_EffAcc_QA(const int energy = 4, const int pid = 0, int mode = 1, int etamode = 0, int ypadding = 2, int order = 2) {

  std::string ordertext[2] = {"","2nd"};

  std::string xtraopt = "_noDelta";

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

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};

  TH1D *h_theta_star_ratio[10][10];
  TH1D *h_theta_ratio[10][10];

  Char_t Char[9][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  Double_t centCent[9] = {75.0,65.0,55.0,45.0,35.0,25.0,15.0,7.5,2.5};

  TF1 *line = new TF1("line","1/3",-0.5,8.5);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
     
  TString *Title;
  {
    //double eff[2][3][vmsa::y_total][7]; // 2 pt bins, 3 centrality bins, 14 |y| bins, 7 cos* bins
    //double eff_error[2][3][vmsa::y_total][7];
    //double eff_y[2][3][vmsa::y_total][7]; // 2 pt bins, 3 centrality bins, 14 |y| bins, 7 cos* bins
    //double eff_error_y[2][3][vmsa::y_total][7];

    TH1D *eff_ratio[vmsa::pt_rebin_cent][9];

    TFile *eff_file;
    if(order == 1) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order1_NoRapiditySpectra/Eff_%s_SingleParticle_noToF_Mode1_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
    if(order == 2) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_NoRapiditySpectra/Eff_%s_SingleParticle_noToF_Mode1_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
    eff_file->Print();
    TFile *eff_file_y;
    if(order == 1) eff_file_y = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order1_WithRapiditySpectra/Eff_%s_SingleParticle_noToF_Mode1_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
    if(order == 2) eff_file_y = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_WithRapiditySpectra/Eff_%s_SingleParticle_noToF_Mode1_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
    eff_file_y->Print();

    for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
    {
      for(int icent = 0; icent < 9; icent++)
      {
        Title = new TString(Form("h_mEffCos_Cent_%d_Pt_%d",icent,ipt));
        cout << Title->Data() << endl;
        TH1D *eff_hist = (TH1D*)eff_file->Get(Title->Data());//->Clone();
        TH1D *eff_hist_y = (TH1D*)eff_file_y->Get(Title->Data());//->Clone();
        //TH1D *eff_hist = (TH1D*)eff_file->Get(Title->Data())->Print();
        //TH1D *eff_hist_y = (TH1D*)eff_file_y->Get(Title->Data())->Print();
        eff_hist->Print();
        eff_hist_y->Print();
        eff_ratio[ipt][icent] = (TH1D*) eff_hist_y->Clone();
        eff_ratio[ipt][icent]->Divide(eff_hist);
        eff_ratio[ipt][icent]->Print();
      }
    }
    //eff_file->Close();
    //eff_file_y->Close();

    TCanvas *c_y = new TCanvas("c_y","c_y",10,10,900,900);
    c_y->Divide(3,3);

    string outputname = Form("output/%s/%s/yspectra_centdpendenceQA_%s_Order%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str(),order);
    string output_start = Form("%s[",outputname.c_str());
    string output_stop = Form("%s]",outputname.c_str());

    c_y->Print(output_start.c_str());

    int centval[10] = {80,70,60,50,40,30,20,10,5,0};
       
    for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
    {
      for(int icent = 0; icent < 9; icent++)
      {
        c_y->cd(icent+1);
        c_y->cd(icent+1)->SetLeftMargin(0.15);
        c_y->cd(icent+1)->SetBottomMargin(0.15);
        c_y->cd(icent+1)->SetTicks(1,1);
        c_y->cd(icent+1)->SetGrid(0,0); 

        cout << centval[icent] << endl;
        cout << centval[icent+1] << endl;
        cout << vmsa::pt_low_cent[energy][ipt] << endl;
	cout << vmsa::pt_up_cent[energy][ipt] << endl;
        string efftitle = Form("Cent %d-%d, %1.1f<p_{T}<%1.1f GeV/c", centval[icent+1], centval[icent], vmsa::pt_low_cent[energy][ipt], vmsa::pt_up_cent[energy][ipt]);
        cout << efftitle << endl;
        eff_ratio[ipt][icent]->Print();
        eff_ratio[ipt][icent]->SetTitle(efftitle.c_str());
        eff_ratio[ipt][icent]->GetXaxis()->SetTitle("|cos#theta*|");
        eff_ratio[ipt][icent]->GetYaxis()->SetTitle("(w/ y-spectra) / (w/o)");
        eff_ratio[ipt][icent]->Draw("pE");
        for(int i = 0; i < 7; i++)
        {
          cout << eff_ratio[ipt][icent]->GetBinContent(i+1) << endl;
        }
      }
      c_y->Update();
      c_y->Print(outputname.c_str());
    }
    c_y->Print(output_stop.c_str());
  }
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

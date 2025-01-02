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

void rapidityEffAccRatio(const int energy = 4, const int pid = 0, bool doall = false, bool isBesI = false, bool random3D = false, int mode = 1, int etamode = 0, int ypadding = 2, int ipt = 1, int icent = 1, std::string ep = "sub", int order = 2) {

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

  TH1D *h_theta_star_ratio[10][10][15];
  TH1D *h_theta_ratio[10][10][15];

  Char_t Char[9][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  Double_t centCent[9] = {75.0,65.0,55.0,45.0,35.0,25.0,15.0,7.5,2.5};

  TF1 *line = new TF1("line","1/3",-0.5,8.5);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
     
  TString *Title;

  double eff[vmsa::y_total][7];
  double eff_error[vmsa::y_total][7];
  double eff_yspec[vmsa::y_total][7];
  double eff_yspec_error[vmsa::y_total][7];
  //TFile *eff_file = new TFile(Form("../RcPhiEffCorr/Rebinned/Eff_20GeV_SingleParticle_noToF_Mode2_EtaMode%d.root",etamode),"READ");
  //TFile *eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_NoEPResAndAcceptance/Eff_19GeV_SingleParticle_noToF_Mode2_EtaMode%d.root",etamode),"READ");
  TFile *eff_file;
  TFile *eff_file_yspec;
  eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_FixedRes/Eff_19GeV_SingleParticle_noToF_Mode2_EtaMode%d.root",etamode),"READ");
  eff_file_yspec = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2_FixedRes_yspec/Eff_19GeV_SingleParticle_noToF_Mode2_EtaMode%d.root",etamode),"READ");
  eff_file->Print();
  eff_file_yspec->Print();

  for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++)
  { 
    Title = new TString(Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy));
    cout << Title->Data() << endl;
    TH1D *eff_hist = (TH1D*)eff_file->Get(Title->Data());
    eff_hist->Print();
    for(int itheta = 0; itheta < 7; itheta++)
    {
      eff[iy][itheta] = eff_hist->GetBinContent(itheta+1);
      eff_error[iy][itheta] = eff_hist->GetBinError(itheta+1)/eff[iy][itheta];
      cout << "icent = " << icent << " ipt = " << ipt << " iy = " << iy << " itheta = " << itheta << " eff = " << eff[iy][itheta] << " +/- " << eff_error[iy][itheta] << endl;
    }
    delete Title;
    delete eff_hist;
  }
  eff_file->Close();

  for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++)
  { 
    Title = new TString(Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy));
    cout << Title->Data() << endl;
    TH1D *eff_hist = (TH1D*)eff_file->Get(Title->Data());
    eff_hist->Print();
    for(int itheta = 0; itheta < 7; itheta++)
    {
      eff_yspec[iy][itheta] = eff_hist->GetBinContent(itheta+1);
      eff_yspec_error[iy][itheta] = eff_hist->GetBinError(itheta+1)/eff_yspec[iy][itheta];
      cout << "icent = " << icent << " ipt = " << ipt << " iy = " << iy << " itheta = " << itheta << " eff = " << eff_yspec[iy][itheta] << " +/- " << eff_yspec_error[iy][itheta] << endl;
    }
    delete Title;
    delete eff_hist;
  }
  eff_file->Close();

  int centVal[4]  = {80,40,10,0};

  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  {
    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
    {
      if( i_dca != 0 && i_sig != 0 ) continue;
      for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
      {
        for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
        {
          for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
          {
            TH1D *rho00_hist;
            TGraphAsymmErrors *g_rho00 = new TGraphAsymmErrors();
  
            double weight_rho00 = 0.0;
            double weight_error_rho00 = 0.0;
            double weight_all = 0.0;

            double y_rho00[vmsa::y_total] = {0.0};
            double y_error_rho00[vmsa::y_total] = {0.0};
            double y_all[vmsa::y_total] = {0.0};
          
            TCanvas *c_fit = new TCanvas("c_fit","c_fit",10,10,1500,600);
            c_fit->Divide(5,2);
            for(int i = 0; i < 10; i++)
            {
              c_fit->cd(i+1)->SetLeftMargin(0.15);
              c_fit->cd(i+1)->SetBottomMargin(0.15);
              c_fit->cd(i+1)->SetGrid(0,0);
              c_fit->cd(i+1)->SetTicks(1,1);
            }
            Title = new TString(Form("fit/Phi/%s/Order%d/yield_pt_%d_cent_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_EtaMode%d_Divided_ratio.pdf",vmsa::mBeamEnergy[energy].c_str(),order,ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),etamode));

            string yieldtitle = Form("yield_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());                             
            TH1D *yieldvsy = new TH1D(yieldtitle.c_str(),yieldtitle.c_str(),10,-1.0,1.0);
            
            TH1F *PtCos[vmsa::y_total];
            TF1 *Func_rho[vmsa::y_total];
            TF1 *Func_obs[vmsa::y_total];
            for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
            {
              string key = Form("eta_%d_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",iy,ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
              cout << key << endl;

              TH1F *PtCos_raw = (TH1F*)input->Get(key.c_str())->Clone("PtCos_raw");
              PtCos[iy] = new TH1F(key.c_str(),key.c_str(), 7, 0, 1);
              PtCos[iy]->Sumw2();
              //delete Title;
              for(int itheta = 0; itheta < 7; itheta++) 
              {
                float inte_mean = PtCos_raw->GetBinContent(itheta+1);
                float inte_mean_error = PtCos_raw->GetBinError(itheta+1);
                cout << "inte_mean  = " << inte_mean << " +/- " << inte_mean_error << endl;
                PtCos[iy]->SetBinContent(itheta+1, inte_mean/eff[iy][itheta]);
                PtCos[iy]->SetBinError(itheta+1, inte_mean/eff[iy][itheta]*TMath::Sqrt(inte_mean_error*inte_mean_error/inte_mean/inte_mean+eff_error[iy][itheta]*eff_error[iy][itheta]));
              }
              PtCos[iy]->Write();

              c_fit->cd(iy+1-ypadding);
              PtCos[iy]->GetXaxis()->SetTitleOffset(1.2);
              PtCos[iy]->SetTitle(Form("%1.1f<y<%1.1f,%1.1f<p_{T}<%1.1fGeV/c,Cent %d-%d",vmsa::y_bin[iy],vmsa::y_bin[iy+1],vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt],centVal[icent+1],centVal[icent]));
              std::string title = Form("%1.1f<y<%1.1f,%1.1f<p_{T}<%1.1fGeV/c,Cent %d-%d",vmsa::y_bin[iy],vmsa::y_bin[iy+1],vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt],centVal[icent+1],centVal[icent]);
              std::cout << title.c_str() << std::endl;
              PtCos[iy]->GetXaxis()->SetTitle("cos#theta*");
              PtCos[iy]->GetYaxis()->SetTitle("EffxAcc Corrected Yield");
              PtCos[iy]->GetYaxis()->SetTitleOffset(1.0);
              PtCos[iy]->SetMarkerColor(kBlack);
              PtCos[iy]->SetMarkerSize(1.8);
              PtCos[iy]->SetMarkerStyle(20);
              PtCos[iy]->Draw("pE");
            
              //Func_rho[iy]->SetParameter(0, PtCos[iy]->GetBinContent(5));
              //Func_rho[iy]->SetParameter(1, 0.333);
              //Func_rho[iy]->SetParLimits(1, 0.0, 1.0);
              //Func_rho[iy]->FixParameter(2, Fval[iy]);
              //Func_rho[iy]->FixParameter(3, Gval[iy]);
              //Func_rho[iy]->FixParameter(4, Res_12[icent]);
              //PtCos[iy]->Fit(Func_rho[iy], "NMRI");
              //Func_rho[iy]->SetLineColor(kRed);
              //Func_rho[iy]->Draw("same");

              Func_rho[iy] = new TF1(Form("Func_rho_%d",iy),FuncAFG,0,1,5);
              Func_obs[iy] = new TF1("Func_obs","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

              PtCos[iy]->Fit(Func_obs[iy],"QNMI"); // fit corrected distribution for observerved rho00
              double rho_obs = Func_obs[iy]->GetParameter(1);
              double rho_obs_err = Func_obs[iy]->GetParError(1);
              Func_obs[iy]->SetLineColor(kRed);
              Func_obs[iy]->Draw("same");
     
              cout << "iy = " << iy << ", rho_obs = " << rho_obs << endl;       

              // Error calculation
              double drhodobs = 4./(1.+3.*Res_12[icent]);  // calculate d(rho)/d(rho_obs)
              double drhodR = -12.*(rho_obs - 1./3.)/(1.+3.*Res_12[icent])/(1.+3.*Res_12[icent]); // calculate d(rho)/d(R)

              double real_rho = 1./3. + 4./(1.+3.*Res_12[icent])*(rho_obs - 1./3.);
              double real_rho_error = TMath::Sqrt((rho_obs_err*rho_obs_err)*(drhodobs*drhodobs) + (Res_12_err[icent]*Res_12_err[icent])*(drhodR*drhodR));
              cout << "iy = " << iy << ", real_rho = " << real_rho << endl;       
         
              //float real_rho = Func_rho[iy]->GetParameter(1);
              //float real_rho_error = Func_rho[iy]->GetParError(1);
              Double_t interr = 0.0;
              float weight = PtCos[iy]->IntegralAndError(1,7,interr);
              double ymean = (vmsa::y_bin[iy] + vmsa::y_bin[iy+1])/2.0;
              int ybin = yieldvsy->FindBin(ymean);
              yieldvsy->SetBinContent(ybin,weight);
              yieldvsy->SetBinError(ybin,interr);
               
              cout << "bin = " << ybin << ", y = " << ymean << ", yield = " << weight << " +/- " << interr << endl;

              weight_rho00 += real_rho*weight;
              weight_error_rho00 += real_rho_error*real_rho_error*weight*weight;
              weight_all += weight;

              y_rho00[iy] += real_rho*weight;
              y_error_rho00[iy] += real_rho_error*real_rho_error*weight*weight;
              y_all[iy] += weight;

              delete PtCos_raw;
            }
            yieldvsy->Write();
            c_fit->SaveAs(Title->Data());

            TCanvas *c_yield = new TCanvas("c_yield","c_yield",10,10,600,600);
            c_yield->SetLeftMargin(0.15);
            c_yield->SetBottomMargin(0.15);
            c_yield->SetGrid(0,0);
            c_yield->SetTicks(1,1);
            yieldvsy->SetMaximum(yieldvsy->GetMaximum()*1.1);
            yieldvsy->SetMinimum(yieldvsy->GetMinimum()*0.9);

            yieldvsy->GetXaxis()->SetTitle("y");
            yieldvsy->GetYaxis()->SetTitle("dN/dy (EffXAcc Corrected)");
            yieldvsy->SetMarkerStyle(20);
            yieldvsy->Draw("pE");
            c_yield->SaveAs(Form("figures/%s/yield_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s.pdf",vmsa::mBeamEnergy[energy].c_str(),ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str()));                             

            //delete c_fit;
            weight_rho00 = weight_rho00/weight_all;
            weight_error_rho00 = TMath::Sqrt(weight_error_rho00)/weight_all;

            for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
            {
              y_rho00[iy] = y_rho00[iy]/y_all[iy];
              y_error_rho00[iy] = TMath::Sqrt(y_error_rho00[iy])/y_all[iy];
            }

            cout << "rho00 = " << weight_rho00<< " +/- " << weight_error_rho00 << endl;

            for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
            {
              cout << y_rho00[iy];
              if(iy == vmsa::y_total-ypadding-1) cout << " " << weight_rho00 << endl;
              else cout << " ";
            }

 
            Title = new TString(Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str()));                             

            g_rho00 = new TGraphAsymmErrors();
            g_rho00->SetName(Title->Data());
            for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++)
            {
              double ymean = (vmsa::y_bin[iy] + vmsa::y_bin[iy+1])/2.0;
              g_rho00->SetPoint(iy-ypadding,ymean,y_rho00[iy]);
              g_rho00->SetPointError(iy-ypadding,0.0,0.0,y_error_rho00[iy],y_error_rho00[iy]);
              cout << "iy = " << iy << " ymean = " << ymean << " y_rho00[iy] = " << y_rho00[iy] << " +/- " << y_error_rho00[iy] << endl;
            }
            delete Title;

            Title = new TString(Form("rhoFinalWeighted_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str()));                             

            rho00_hist = new TH1D(Title->Data(), Title->Data(), vmsa::y_total+1, 0.5, (double)vmsa::y_total+1.5);
            rho00_hist->SetBinContent(vmsa::y_total+1,weight_rho00);
            rho00_hist->SetBinError(vmsa::y_total+1,weight_error_rho00);
            for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
            {
              rho00_hist->SetBinContent(iy+1,y_rho00[iy]);
              rho00_hist->SetBinError(iy+1,y_error_rho00[iy]);
            }

            g_rho00->Write();
            rho00_hist->Write();

            delete Title;
            delete rho00_hist;
            delete g_rho00;
          }
        }
      }
    }
  }
  
  output->Close();
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

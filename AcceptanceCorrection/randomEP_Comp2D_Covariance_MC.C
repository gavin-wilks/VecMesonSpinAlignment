#include "TFile.h"
#include <iostream>
#include "TH3D.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF2.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "../Utility/functions.h"
#include <vector>
#include "TMatrixDSym.h"

double rho00_FromRandom(double R, double rho00_EP, double rho00_Rand)
{
  return rho00_EP/R+(1.-1./R)*rho00_Rand;
}

double drho00_FromRandom(double R, double rho00_EP, double rho00_Rand, double drho00_EP, double drho00_Rand, double cov)
{
  double d_drho00_EP = 1./R;
  double d_drho00_Rand = (1.-1./R); 
 
  return sqrt(d_drho00_EP*d_drho00_EP*drho00_EP*drho00_EP + d_drho00_Rand*d_drho00_Rand*drho00_Rand*drho00_Rand + 2.*d_drho00_EP*d_drho00_Rand*cov);
}

double rerho1n1_FromRandom(double R, double rho00_EP, double rho00_Rand)
{
  return -rho00_EP/(2.*R)+rho00_Rand*(1.+3.*R)/(2.*R)-1./2.;
}

double drerho1n1_FromRandom(double R, double rho00_EP, double rho00_Rand, double drho00_EP, double drho00_Rand, double cov)
{
  double d_drho00_EP = -1./(2.*R);
  double d_drho00_Rand = (1.+3.*R)/(2.*R); 
 
  return sqrt(d_drho00_EP*d_drho00_EP*drho00_EP*drho00_EP + d_drho00_Rand*d_drho00_Rand*drho00_Rand*drho00_Rand + 2.*d_drho00_EP*d_drho00_Rand*cov);
}




double rho00_RP_1D(double R, double rho00_EP)
{
  return 4./(1.+3.*R)*(rho00_EP-1./3.)+1./3.;
}
double drho00_RP_1D(double R, double rho00_EP, double dR, double drho00_EP)
{
  double ddR = -12.*(rho00_EP-1./3.)/(1.+3.*R)/(1.+3.*R);
  double ddrho00_EP = 4./(1.+3.*R);

  return TMath::Sqrt(ddR*ddR*dR*dR + ddrho00_EP*ddrho00_EP*drho00_EP*drho00_EP);
}


double rho00_RP(double R, double rho00_EP, double rerho1n1_EP)
{
  return ((-1.+R)*(1.+2.*rerho1n1_EP)+(3.+R)*rho00_EP)/(4.*R);
}

double drho00_RP(double R, double rho00_EP, double rerho1n1_EP, double dR, double drho00_EP, double drerho1n1_EP, double cov)
{
  double ddR = (1./(4.*R*R))*(1.+2.*rerho1n1_EP)-3./(4.*R*R)*rho00_EP;
  double ddrho00_EP = (3.+R)/(4.*R);
  double ddrerho1n1_EP = (-1.+R)/(2.*R);

  return TMath::Sqrt(ddR*ddR*dR*dR + ddrho00_EP*ddrho00_EP*drho00_EP*drho00_EP + ddrerho1n1_EP*ddrerho1n1_EP*drerho1n1_EP*drerho1n1_EP + 2.*ddrho00_EP*ddrerho1n1_EP*cov);
}


double rerho1n1_RP(double R, double rho00_EP, double rerho1n1_EP)
{
  return (1.+2.*rerho1n1_EP-3.*rho00_EP+R*(-1.+6.*rerho1n1_EP+3.*rho00_EP))/(8.*R);
}

double drerho1n1_RP(double R, double rho00_EP, double rerho1n1_EP, double dR, double drho00_EP, double drerho1n1_EP, double cov)
{
  double ddR = (-1./(8.*R*R))*(1.+2.*rerho1n1_EP-3.*rho00_EP);
  double ddrho00_EP = (-3.+3.*R)/(8.*R);
  double ddrerho1n1_EP = (1.+3.*R)/(4.*R);

  return TMath::Sqrt(ddR*ddR*dR*dR + ddrho00_EP*ddrho00_EP*drho00_EP*drho00_EP + ddrerho1n1_EP*ddrerho1n1_EP*drerho1n1_EP*drerho1n1_EP + 2.*ddrho00_EP*ddrerho1n1_EP*cov);
}

double rhofit(double *x, double *par)
{
  double ep   = (1.-par[0])+(3.*par[0]-1.)*x[0]*x[0];
  double rand = (1.-par[1])+(3.*par[1]-1.)*x[1]*x[1];

  return par[2]*ep*rand;
}


void randomEP_Comp2D_Covariance_MC(int study = 0, double ep_res = 0.4, int order = 1)
{
  int i = 1; 
  string pdftag = "both0p015_withweight_removecut";
  
  double inputrho = 1./3.+0.015;
  double inputre  = 0.015;


  double ep_res1 = 0.0; 
  double ep_res2 = 0.0;

  if(order == 1) 
  {
    TF1* f_res1 = new TF1(Form("resolution1"),EventPlaneResolution,0,80,0);
    TF1* f_res1_k2 = new TF1(Form("resolution1_k2"),EventPlaneResolutionK2,0,80,0);
    double mChi1 = f_res1->GetX(ep_res); // This is for sub  event plane resolution
    cout << "mChi = " << mChi1 << endl;
    mChi1 *= TMath::Sqrt(2.);
    cout << "mChi*sqrt(2) = " << mChi1 << endl;
    ep_res1 = f_res1->Eval(mChi1);
    ep_res2 = f_res1_k2->Eval(mChi1);
    cout << "Approximate method for ep_res1 = " << 0.626657*mChi1-0.09694*pow(mChi1,3)+0.02754*pow(mChi1,4)-0.002283*pow(mChi1,5) << endl;;
    cout << "Approximate method for ep_res2 = " << 0.25*pow(mChi1,2)-0.011414*pow(mChi1,3)-0.034726*pow(mChi1,4)+0.006815*pow(mChi1,5) << endl;;
  }
  if(order == 2)
  {
    TF1* f_res2 = new TF1(Form("resolution2"),EventPlaneResolution,0,80,0);
    double mChi = f_res2->GetX(ep_res);
    //ep_res2 = ep_res;
    ep_res2 = f_res2->Eval(mChi);
  }
  cout << "ep_res1 = " << ep_res1 << endl; 
  cout << "ep_res2 = " << ep_res2 << endl; 

  TFile *file[500];
  //TFile *file[0] = new TFile::Open("randomEPCovTest/random_ep_cent5_psi1_method2_study0.root");
  
  TH3D *h3[2][3][500];
  TH2D *psipsirand[500];
  TProfile *flow[2][2][500];

  std::vector<int> files;
  
  for(int ifile = 0; ifile < 500; ifile++)
  { 
    file[ifile] = new TFile(Form("randomEPCovTest/random_ep_cent5_psi%d_method2_study%d_epres0p4_deltarho0p015_rerho1n10p015_20260309_noweights_actual/Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_pt100.0_y100.0_ptbin5_cent10_%d.root",order,study,ifile),"READ");

    if(!file[ifile] || file[ifile]->IsZombie())
    {
      continue;
    }
    else 
    {
      files.push_back(ifile);
    }

    file[ifile]->Print();
    for(int icut = 0; icut < 3; icut++)
    {
      h3[0][icut][ifile] = (TH3D*)((TH3D*) file[ifile]->Get(Form("random_cut_%d",icut)))->Clone(Form("random_cut_%d_%d",icut,ifile)); 
    }
    
    psipsirand[ifile] = (TH2D*)((TH2D*) file[ifile]->Get("PsiPsiRandom"))->Clone(Form("PsiPsiRandom_%d",ifile));
    
    flow[0][0][ifile] = (TProfile*)((TProfile*) file[ifile]->Get("v1_MC"))->Clone(Form("v1_MC_%d",ifile));
    flow[0][1][ifile] = (TProfile*)((TProfile*) file[ifile]->Get("v1_RC"))->Clone(Form("v1_RC_%d",ifile));
    flow[1][0][ifile] = (TProfile*)((TProfile*) file[ifile]->Get("v2_MC"))->Clone(Form("v2_MC_%d",ifile));
    flow[1][1][ifile] = (TProfile*)((TProfile*) file[ifile]->Get("v2_RC"))->Clone(Form("v2_RC_%d",ifile));

    //h3[1][0][ifile] = (TH3D*)((TH3D*) file->Get("h3_mMcEffCosPhiPrimeYSmearMC_epsmear_0"))->Clone(Form("h3_mMcEffCosPhiPrimeYSmearMC_epsmear_0_%d",ifile));
    //h3[1][1][ifile] = (TH3D*)((TH3D*) file->Get("h3_mMcEffCosPhiPrimeYSmearMC_epsmear_1"))->Clone(Form("h3_mMcEffCosPhiPrimeYSmearMC_epsmear_1_%d",ifile));
    //h3[1][2][ifile] = (TH3D*)((TH3D*) file->Get("h3_mMcEffCosPhiPrimeYSmear_cut1_epsmear_1_ptsmear_41"))->Clone(Form("h3_mMcEffCosPhiPrimeYSmear_cut1_epsmear_1_ptsmear_41_%d",ifile));
  }

  const int nbins = h3[0][0][0]->GetNbinsX(); 
  TH2D *h2_rand[3][500]; // Cut Level, study bin 
  TH2D *h2_2D[3][500]; // Cut Level, study bin 
  TH1D *h1_rand[2][3][500];
  TF2 *fit2_rand[2][500];
  TF2 *fit2_2D[2][500]; 
  TF1 *fit1_rand[2][3][500];

  TH2D *psipsirandAll;

  // Now go to random EP stuff
  for(int ifile : files)
  {
    if(ifile == files[0]) psipsirandAll = (TH2D*) psipsirand[ifile]->Clone("psipsirandAll");
    else                 psipsirandAll->Add(psipsirand[ifile],1.0);
    for(int icut = 0; icut < 3; icut++)
    {
      {
        h2_rand[icut][ifile] = (TH2D*) ((TH2D*) h3[0][icut][ifile]->Project3D("zy"))->Clone(Form("rand2D_%d_%d_%d",icut,i-1,ifile));

        h1_rand[0][icut][ifile] = (TH1D*) h3[0][icut][ifile]->ProjectionY(Form("rand1D_0_%d_%d_%d",icut,i-1,ifile),0,-1,0,-1);
        h1_rand[1][icut][ifile] = (TH1D*) h3[0][icut][ifile]->ProjectionZ(Form("rand1D_1_%d_%d_%d",icut,i-1,ifile),0,-1,0,-1);
      }
    }
  }  

  double rhoEP[3][500];
  double rhoEPerr[3][500];
  double rhoRand[3][500];
  double rhoRanderr[3][500];

  for(int ifile : files)
  {
    for(int ic = 0; ic < 3; ic++)
    {
      {
        for(int icos = 0; icos < 2; icos++)
        {
          fit1_rand[icos][ic][ifile] = new TF1(Form("fit1_rand_%d_%d_%d_%d",icos,ic,i,ifile),SpinDensity,-1,1,2);
          fit1_rand[icos][ic][ifile]->SetParameter(0,1./3.);
          fit1_rand[icos][ic][ifile]->SetParameter(1,h1_rand[icos][ic][ifile]->GetMaximum());
 
          h1_rand[icos][ic][ifile]->Fit(fit1_rand[icos][ic][ifile],"NMRI");
        }

        double rho00EP = fit1_rand[0][ic][ifile]->GetParameter(0);
        double rho00EPerr = fit1_rand[0][ic][ifile]->GetParError(0);
        double rho00rand = fit1_rand[1][ic][ifile]->GetParameter(0);
        double rho00randerr = fit1_rand[1][ic][ifile]->GetParError(0);

        rhoEP[ic][ifile] = rho00EP;
        rhoRand[ic][ifile] = rho00rand;
        rhoEPerr[ic][ifile] = rho00EPerr;
        rhoRanderr[ic][ifile] = rho00randerr;
     
        cout << "ifile = " << ifile << ", ic = " << ic << ", rhoEP = " << rhoEP[ic][ifile] << " +/- " << rhoEPerr[ic][ifile] << endl;
        cout << "ifile = " << ifile << ", ic = " << ic << ", rhoRand = " << rhoRand[ic][ifile] << " +/- " << rhoRanderr[ic][ifile] << endl << endl;
      }
    }
  }

  TCanvas *c = new TCanvas("c","c",10,10,1200,2400);
  c->Divide(3,6);
  for(int i = 0; i < 18; i++)
  {
    c->cd(i+1);
    c->cd(i+1)->SetLeftMargin(0.15);  
    c->cd(i+1)->SetBottomMargin(0.15);
    c->cd(i+1)->SetTicks(1,1);
    c->cd(i+1)->SetGrid(0,0);
  }


  double low  = 0.341;
  double high = 0.349;
  //double low  = 1./3.-0.0025;
  //double high = 1./3.+0.0025;

  TH1D *h = new TH1D("h","h",100,low,high);
  for(int ix = 1; ix <= h->GetNbinsX(); ix++)
  {
    h->SetBinContent(ix,-1); 
  }

  string cuttag[2] = {"BeforeCuts","AfterCuts"};

  for(int ic = 1; ic < 3; ic++)
  {
    TGraphAsymmErrors *gcorr = new TGraphAsymmErrors();
    TH2D *h2 = new TH2D("h2","h2",20,low,high,20,low,high);

    TGraph *g_v1_ep = new TGraph();
    TGraph *g_v2_ep = new TGraph();
    TGraph *g_v1_rand = new TGraph();
    TGraph *g_v2_rand = new TGraph();

    TH2D *h2_v1_ep = new TH2D("h2_v1_ep","h2_v1_ep",20,low,high,20,-0.003,0.003);
    TH2D *h2_v2_ep = new TH2D("h2_v2_ep","h2_v2_ep",20,low,high,20,-0.003,0.003);
    TH2D *h2_v1_rand = new TH2D("h2_v1_rand","h2_v1_rand",20,low,high,20,-0.003,0.003);
    TH2D *h2_v2_rand = new TH2D("h2_v2_rand","h2_v2_rand",20,low,high,20,-0.003,0.003);


    int idx = 0;
    for(int ifile : files) 
    {
      gcorr->SetPoint(idx,rhoEP[ic][ifile],rhoRand[ic][ifile]);
      gcorr->SetPointError(idx,rhoEPerr[ic][ifile],rhoEPerr[ic][ifile],rhoRanderr[ic][ifile],rhoRanderr[ic][ifile]);

      g_v1_ep->SetPoint(idx,rhoEP[ic][ifile],flow[0][ic-1][ifile]->GetBinContent(1));    
      g_v2_ep->SetPoint(idx,rhoEP[ic][ifile],flow[1][ic-1][ifile]->GetBinContent(1));    
      g_v1_rand->SetPoint(idx,rhoRand[ic][ifile],flow[0][ic-1][ifile]->GetBinContent(1));    
      g_v2_rand->SetPoint(idx,rhoRand[ic][ifile],flow[1][ic-1][ifile]->GetBinContent(1));    

      idx++;
      h2->Fill(rhoEP[ic][ifile],rhoRand[ic][ifile]);
      h2_v1_ep->Fill(rhoEP[ic][ifile],flow[0][ic-1][ifile]->GetBinContent(1));    
      h2_v2_ep->Fill(rhoEP[ic][ifile],flow[1][ic-1][ifile]->GetBinContent(1));    
      h2_v1_rand->Fill(rhoRand[ic][ifile],flow[0][ic-1][ifile]->GetBinContent(1));    
      h2_v2_rand->Fill(rhoRand[ic][ifile],flow[1][ic-1][ifile]->GetBinContent(1));    

    }
    ///////////////////////////////////////////////////////////////////// 
    int ipad = 1;
    c->cd(ipad); 
    h->SetStats(0);
    h->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{1}}");
    h->GetYaxis()->SetTitle("#rho_{00}^{#Psi_{Rand}}");
    h->GetYaxis()->SetRangeUser(low,high);

    gcorr->SetMarkerStyle(24);
    gcorr->SetMarkerColor(kBlack);

    TFitResultPtr r = gcorr->Fit("pol1","SW");
    TMatrixDSym cov = r->GetCovarianceMatrix();
    cov.Print();

    TF1 *f = gcorr->GetFunction("pol1");
    h->SetTitle(Form("Corr Factor = %1.3f, slope=%1.3f#pm%1.3f",gcorr->GetCorrelationFactor(),f->GetParameter(1),f->GetParError(1)));
 
    h->DrawCopy();
    gcorr->Draw("PX same");
    ipad++;

    c->cd(ipad);
    h2->SetStats(0);
    h2->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{1}}");
    h2->GetYaxis()->SetTitle("#rho_{00}^{#Psi_{Rand}}");
    h2->Draw("COLZ");
    ipad++;
 
    c->cd(ipad);
    TProfile *p2 = h2->ProfileX();
    p2->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{1}}");
    p2->GetYaxis()->SetTitle("#rho_{00}^{#Psi_{Rand}}");
    p2->GetYaxis()->SetRangeUser(low,high);
    p2->SetStats(0);
    TF1 *f2 = new TF1("f2","[0]+x*[1]",0.3,0.37);
    p2->Approximate();
    p2->Fit(f2,"M");
    p2->SetTitle(Form("slope=%1.3f#pm%1.3f",f2->GetParameter(1),f2->GetParError(1)));
    p2->Draw("PE");
    ipad++;
    ////////////////////////////////////////////////////////////////////////// 
    c->cd(ipad); 
    h->SetStats(0);
    h->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{1}}");
    h->GetYaxis()->SetTitle("v_{1}");
    h->GetYaxis()->SetRangeUser(-0.003,0.003);

    g_v1_ep->SetMarkerStyle(24);
    g_v1_ep->SetMarkerColor(kBlack);
    g_v1_ep->Fit("pol1","SW");

    TF1 *f_v1_ep = g_v1_ep->GetFunction("pol1");
    h->SetTitle(Form("Corr Factor = %1.3f, slope=%1.3f#pm%1.3f",g_v1_ep->GetCorrelationFactor(),f_v1_ep->GetParameter(1),f_v1_ep->GetParError(1)));
    //h->SetTitle(Form("Corr Factor = %1.3f",g_v1_ep->GetCorrelationFactor())); 

    h->DrawCopy();
    g_v1_ep->Draw("PX same");
    ipad++;

    c->cd(ipad);
    h2_v1_ep->SetStats(0);
    h2_v1_ep->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{1}}");
    h2_v1_ep->GetYaxis()->SetTitle("v_{1}");
    h2_v1_ep->Draw("COLZ");
    ipad++;  
 
    c->cd(ipad);
    TProfile *p2_v1_ep = h2_v1_ep->ProfileX();
    p2_v1_ep->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{1}}");
    p2_v1_ep->GetYaxis()->SetTitle("v_{1}");
    p2_v1_ep->GetYaxis()->SetRangeUser(-0.003,0.003);
    TF1 *f2_v1_ep = new TF1("f2_v1_ep","[0]+x*[1]",0.3,0.37);
    p2_v1_ep->Approximate();
    p2_v1_ep->Fit(f2_v1_ep,"M");
    p2_v1_ep->SetTitle(Form("slope=%1.3f#pm%1.3f",f2_v1_ep->GetParameter(1),f2_v1_ep->GetParError(1)));
    p2_v1_ep->Draw("PE");
    p2_v1_ep->SetStats(0);
    ipad++;
    ////////////////////////////////////////////////////////////////////////// 
    c->cd(ipad); 
    h->SetStats(0);
    h->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{Rand}}");
    h->GetYaxis()->SetTitle("v_{1}");
    h->GetYaxis()->SetRangeUser(-0.003,0.003);

    g_v1_rand->SetMarkerStyle(24);
    g_v1_rand->SetMarkerColor(kBlack);
    g_v1_rand->Fit("pol1","SW");

    TF1 *f_v1_rand = g_v1_rand->GetFunction("pol1");
    h->SetTitle(Form("Corr Factor = %1.3f, slope=%1.3f#pm%1.3f",g_v1_rand->GetCorrelationFactor(),f_v1_rand->GetParameter(1),f_v1_rand->GetParError(1)));
    //h->SetTitle(Form("Corr Factor = %1.3f",g_v1_rand->GetCorrelationFactor())); 

    h->DrawCopy();
    g_v1_rand->Draw("PX same");
    ipad++;

    c->cd(ipad);
    h2_v1_rand->SetStats(0);
    h2_v1_rand->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{Rand}}");
    h2_v1_rand->GetYaxis()->SetTitle("v_{1}");
    h2_v1_rand->Draw("COLZ");
    ipad++;  
 
    c->cd(ipad);
    TProfile *p2_v1_rand = h2_v1_rand->ProfileX();
    p2_v1_rand->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{Rand}}");
    p2_v1_rand->GetYaxis()->SetTitle("v_{1}");
    p2_v1_rand->GetYaxis()->SetRangeUser(-0.003,0.003);
    p2_v1_rand->SetStats(0);
    TF1 *f2_v1_rand = new TF1("f2_v1_rand","[0]+x*[1]",0.3,0.37);
    p2_v1_rand->Approximate();
    p2_v1_rand->Fit(f2_v1_rand,"M");
    p2_v1_rand->SetTitle(Form("slope=%1.3f#pm%1.3f",f2_v1_rand->GetParameter(1),f2_v1_rand->GetParError(1)));
    p2_v1_rand->Draw("PE");
    ipad++;
    ////////////////////////////////////////////////////////////////////////// 
    c->cd(ipad); 
    h->SetStats(0);
    h->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{1}}");
    h->GetYaxis()->SetTitle("v_{2}");
    h->GetYaxis()->SetRangeUser(-0.003,0.003);

    g_v2_ep->SetMarkerStyle(24);
    g_v2_ep->SetMarkerColor(kBlack);
    g_v2_ep->Fit("pol1","SW");

    TF1 *f_v2_ep = g_v2_ep->GetFunction("pol1");
    h->SetTitle(Form("Corr Factor = %1.3f, slope=%1.3f#pm%1.3f",g_v2_ep->GetCorrelationFactor(),f_v2_ep->GetParameter(1),f_v2_ep->GetParError(1)));
    //h->SetTitle(Form("Corr Factor = %1.3f",g_v2_ep->GetCorrelationFactor())); 

    h->DrawCopy();
    g_v2_ep->Draw("PX same");
    ipad++;

    c->cd(ipad);
    h2_v2_ep->SetStats(0);
    h2_v2_ep->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{1}}");
    h2_v2_ep->GetYaxis()->SetTitle("v_{2}");
    h2_v2_ep->Draw("COLZ");
    ipad++;  
 
    c->cd(ipad);
    TProfile *p2_v2_ep = h2_v2_ep->ProfileX();
    p2_v2_ep->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{1}}");
    p2_v2_ep->GetYaxis()->SetTitle("v_{2}");
    p2_v2_ep->GetYaxis()->SetRangeUser(-0.003,0.003);
    p2_v2_ep->SetStats(0);
    TF1 *f2_v2_ep = new TF1("f2_v2_ep","[0]+x*[1]",0.3,0.37);
    p2_v2_ep->Approximate();
    p2_v2_ep->Fit(f2_v2_ep,"M");
    p2_v2_ep->SetTitle(Form("slope=%1.3f#pm%1.3f",f2_v2_ep->GetParameter(1),f2_v2_ep->GetParError(1)));
    p2_v2_ep->Draw("PE");
    ipad++;
    ////////////////////////////////////////////////////////////////////////// 
    c->cd(ipad); 
    h->SetStats(0);
    h->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{Rand}}");
    h->GetYaxis()->SetTitle("v_{2}");
    h->GetYaxis()->SetRangeUser(-0.003,0.003);

    g_v2_rand->SetMarkerStyle(24);
    g_v2_rand->SetMarkerColor(kBlack);
    g_v2_rand->Fit("pol1","SW");

    TF1 *f_v2_rand = g_v2_rand->GetFunction("pol1");
    h->SetTitle(Form("Corr Factor = %1.3f, slope=%1.3f#pm%1.3f",g_v2_rand->GetCorrelationFactor(),f_v2_rand->GetParameter(1),f_v2_rand->GetParError(1)));
    //h->SetTitle(Form("Corr Factor = %1.3f",g_v2_rand->GetCorrelationFactor())); 

    h->DrawCopy();
    g_v2_rand->Draw("PX same");
    ipad++;

    c->cd(ipad);
    h2_v2_rand->SetStats(0);
    h2_v2_rand->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{Rand}}");
    h2_v2_rand->GetYaxis()->SetTitle("v_{2}");
    h2_v2_rand->Draw("COLZ");
    ipad++;  
 
    c->cd(ipad);
    TProfile *p2_v2_rand = h2_v2_rand->ProfileX();
    p2_v2_rand->GetXaxis()->SetTitle("#rho_{00}^{#Psi_{Rand}}");
    p2_v2_rand->GetYaxis()->SetTitle("v_{2}");
    p2_v2_rand->GetYaxis()->SetRangeUser(-0.003,0.003);
    p2_v2_rand->SetStats(0);
    TF1 *f2_v2_rand = new TF1("f2_v2_rand","[0]+x*[1]",0.3,0.37);
    p2_v2_rand->Approximate();
    p2_v2_rand->Fit(f2_v2_rand,"M");
    p2_v2_rand->SetTitle(Form("slope=%1.3f#pm%1.3f",f2_v2_rand->GetParameter(1),f2_v2_rand->GetParError(1)));
    p2_v2_rand->Draw("PE");
    ipad++;
    ////////////////////////////////////////////////////////////////////////// 
    ipad++;

    c->cd(ipad);
    psipsirandAll->SetStats(0);
    psipsirandAll->GetXaxis()->SetTitle("#Psi_{1}");
    psipsirandAll->GetYaxis()->SetTitle("#Psi_{rand}");
    psipsirandAll->Draw("COLZ");
    ipad++;  
 
    c->cd(ipad);
    TProfile *p_psipsirandAll = psipsirandAll->ProfileX();
    p_psipsirandAll->GetXaxis()->SetTitle("#Psi_{1}");
    p_psipsirandAll->GetYaxis()->SetTitle("#Psi_{rand}");
    p_psipsirandAll->SetStats(0);
    TF1 *f2_psipsirand = new TF1("f2_psipsirand","[0]+x*[1]");
    p_psipsirandAll->Approximate();
    p_psipsirandAll->Fit(f2_psipsirand,"M");
    p_psipsirandAll->SetTitle(Form("slope=%1.3f#pm%1.3f",f2_psipsirand->GetParameter(1),f2_psipsirand->GetParError(1)));
    p_psipsirandAll->Draw("PE");
    ipad++;
    

 
    c->SaveAs(Form("figures/RandomEP/RandomEP_rhocorr_Order%d_%s_%s_inputrhorerho1n1.pdf",order,pdftag.c_str(),cuttag[ic-1].c_str())); 
  } 
}

#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
//#include "../Utility/draw.h"
#include "../../../draw.h"
#include "../Utility/functions.h"
#include "resolution_pt.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void compareDataRcPhiStarHReal(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 23, std::string simmode = "Rc")//defaultF = 0 is BESII, defaultF = 1 is BESI
{
  const int nvar = 7;
  std::string sigmay_text[nvar] = {"-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"};
  float sigmay[nvar] = {-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3};

  const int defaultfile = 3;

  std::string spectra[nvar];
  for(int i = 0; i < nvar; i++)
  {
    spectra[i] = Form("CEPT_O2_NoYSpec_prelimv2_phieff_phiTPC_20240808_2_ysig1000_rho0.3333_rerho0.0_imrho0.0_r0.0_i0.0_hrho0.3333_hrerho0.0_himrho0.0_hr%s_hi0.0",sigmay_text[i].c_str());
  }
  
  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  
  string inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_%sPhiPtSys_%s_PolySys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1D *h_mCounts[4];  

  TFile *File_InPutRC[nvar]; 
  TH1D *h_mCounts_Cos_RC[4][nvar];  
  TH1D *h_mCounts_RC[4][nvar];  
  TH2D *h_m2D_RC[4][nvar];  
  TH1D *h_mCounts_Ratio[4][nvar];  
  for(int i_pt = 2; i_pt < 6; ++i_pt) // pt loop
  {
    string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    h_mCounts[i_pt-2] = (TH1D*) File_InPut->Get(KEY_counts.c_str());

    float data = h_mCounts[i_pt-2]->Integral(1,10);

    for(int i = 0; i < nvar; i++)
    {
      string inputfileRC = Form("effaccfiles/%s/%s/%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[i].c_str(),vmsa::mBeamEnergy[energy].c_str());
      File_InPutRC[i] = TFile::Open(inputfileRC.c_str());
      string KEY_counts_RC = Form("h_m%sEffPhiS_Cent_9_Pt_%d",simmode.c_str(),i_pt);
      string KEY_2D_RC = Form("h_m%sEffCosHPhiPrime_Cent_9_Pt_%d",simmode.c_str(),i_pt);
      string KEY_counts_Cos_RC = Form("h_m%sEffCos_Cent_9_Pt_%d",simmode.c_str(),i_pt);
      h_mCounts_RC[i_pt-2][i] = (TH1D*) ((TH1D*) File_InPutRC[i]->Get(KEY_counts_RC.c_str()))->Clone();
      h_m2D_RC[i_pt-2][i] = (TH2D*) ((TH2D*) File_InPutRC[i]->Get(KEY_2D_RC.c_str()))->Clone();
      h_mCounts_Cos_RC[i_pt-2][i] = (TH1D*) ((TH1D*) File_InPutRC[i]->Get(KEY_counts_Cos_RC.c_str()))->Clone();
      float rc   = h_mCounts_RC[i_pt-2][i]->Integral(1,10);

      h_mCounts_RC[i_pt-2][i]->Scale(data/rc);

      h_mCounts_Ratio[i_pt-2][i] = (TH1D*) h_mCounts[i_pt-2]->Clone();
      h_mCounts_Ratio[i_pt-2][i]->Divide(h_mCounts_RC[i_pt-2][i]);
    }
  }

  TGraphAsymmErrors *gslope[4];
  TGraphAsymmErrors *gintercept[4]; 
  TGraphAsymmErrors *gmean[4]; 
  TGraphAsymmErrors *gglobal[4]; 
  TGraphAsymmErrors *g2D[5][4]; 
  TF1 *fits[4][nvar];
  TF1 *fits2D[4][nvar];
  TF1 *fits1Dcos[4][nvar];
  float pt[4] = {1.5,2.1,2.7,3.6};

  TCanvas *c = new TCanvas("c","c",10,10,1200,1200);
  c->Divide(3,3);

  std::string outputname = Form("figures/%s/%s/pTstudy/%s_realcomparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[defaultfile].c_str(),simmode.c_str());  
  std::string outputstart = Form("%s[",outputname.c_str()); 
  std::string outputstop = Form("%s]",outputname.c_str()); 
 
  c->Print(outputstart.c_str());

  for(int ipt = 0; ipt < 4; ipt++)
  {
    gslope[ipt] = new TGraphAsymmErrors();
    gintercept[ipt] = new TGraphAsymmErrors();
    gmean[ipt] = new TGraphAsymmErrors();
    gglobal[ipt] = new TGraphAsymmErrors();
    for(int par = 0; par < 5; par++) g2D[par][ipt] = new TGraphAsymmErrors();
    for(int i = 0; i < nvar; i++)
    {
      c->cd(i+1);
      c->cd(i+1)->SetLeftMargin(0.15);
      c->cd(i+1)->SetBottomMargin(0.15);
      c->cd(i+1)->SetTicks(1,1);
      c->cd(i+1)->SetGrid(0,0);
        
      double mean = h_mCounts_RC[ipt][i]->GetMean();
      double mean_err = h_mCounts_RC[ipt][i]->GetMeanError();

      h_mCounts_Ratio[ipt][i]->SetTitle(Form("Re(#rho_{10})-Re(#rho_{0-1})=%s, %1.1f<p_{T}<%1.1f, 20-60 Centrality",sigmay_text[i].c_str(),vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mCounts_Ratio[ipt][i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
      h_mCounts_Ratio[ipt][i]->GetYaxis()->SetTitle(Form("(Data)/(%s from Simulation)",simmode.c_str()));
      h_mCounts_Ratio[ipt][i]->Draw("pE");
      fits1Dcos[ipt][i] = new TF1(Form("fit1Dcos_%d_%d",ipt,i),SpinDensity,-1.0,1.0,2);
      h_mCounts_Cos_RC[ipt][i]->Fit(fits1Dcos[ipt][i],"NMI");
     
      fits[ipt][i] = new TF1(Form("fit_%d_%d",ipt,i),"[0]+[1]*x",-1.0,1.0);
      h_mCounts_Ratio[ipt][i]->Fit(fits[ipt][i],"NMI");
      fits[ipt][i]->SetLineColor(kRed);
      fits[ipt][i]->Draw("l same");

      fits2D[ipt][i] = new TF2(Form("fit2D_%d_%d",ipt,i),SpinDensity2Dcos,-1.0,1.0,-TMath::Pi(),TMath::Pi(),6);
      h_m2D_RC[ipt][i]->Fit(fits2D[ipt][i],"NMI");
      for(int par = 0; par < 5; par++)
      {
        if(par == 0) g2D[par][ipt]->SetPoint(i,sigmay[i],fits2D[ipt][i]->GetParameter(0)-1./3.);
        else g2D[par][ipt]->SetPoint(i,sigmay[i],fits2D[ipt][i]->GetParameter(par));
        g2D[par][ipt]->SetPointError(i,0.0,0.0,fits2D[ipt][i]->GetParError(par),fits2D[ipt][i]->GetParError(par));
      }
      
      gslope[ipt]->SetPoint(i,sigmay[i],fits[ipt][i]->GetParameter(1)); 
      gslope[ipt]->SetPointError(i,0.0,0.0,fits[ipt][i]->GetParError(1),fits[ipt][i]->GetParError(1));
      gintercept[ipt]->SetPoint(i,sigmay[i],fits[ipt][i]->GetParameter(0)); 
      gintercept[ipt]->SetPointError(i,0.0,0.0,fits[ipt][i]->GetParError(0),fits[ipt][i]->GetParError(0));
      gmean[ipt]->SetPoint(i,sigmay[i],mean); 
      gmean[ipt]->SetPointError(i,0.0,0.0,mean_err,mean_err);
      gglobal[ipt]->SetPoint(i,sigmay[i],fits1Dcos[ipt][i]->GetParameter(0)-1./3.);
      gglobal[ipt]->SetPointError(i,0.0,0.0,fits1Dcos[ipt][i]->GetParError(0),fits1Dcos[ipt][i]->GetParError(0));
    }
    c->Update();
    c->Print(outputname.c_str());
  }
  c->Print(outputstop.c_str());


  TCanvas *cfit2D = new TCanvas("cfit2D","cfit2D",10,10,2000,1600);
  cfit2D->Divide(5,4);
  
  for(int i = 0; i < 20; i++)
  {
    cfit2D->cd(i+1);
    cfit2D->cd(i+1)->SetLeftMargin(0.20);  
    cfit2D->cd(i+1)->SetBottomMargin(0.15);
    cfit2D->cd(i+1)->SetTicks(1,1);
    cfit2D->cd(i+1)->SetGrid(0,0);
  }
 
  float ptedges[5] = {1.2,1.8,2.4,3.0,4.2};
  string param[5] = {"#Delta#rho_{00}","Re(#rho_{10})-Re(#rho_{0-1})","Im(#rho_{10})-Im(#rho_{0-1})","Re(#rho_{1-1})","Im(#rho_{1-1})"};

  for(int i = 0; i < 4; i++)
  {
    for(int par = 0; par < 5; par++)
    {
      cfit2D->cd(1+i*5+par);
      //cfit->cd(i+5)->SetLogx();
      g2D[par][i]->SetMarkerStyle(20);
      g2D[par][i]->SetMarkerSize(1.2);
      g2D[par][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
      g2D[par][i]->GetXaxis()->SetTitle("Re(#rho_{10})-Re(#rho_{0-1})");
      g2D[par][i]->GetYaxis()->SetTitle(param[par].c_str());
      g2D[par][i]->Draw("APE"); 
    }
  }
  cfit2D->SaveAs(Form("figures/%s/%s/pTstudy/%s_parametercomparison_realcomparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[defaultfile].c_str(),simmode.c_str()));  

  TCanvas *cfit = new TCanvas("cfit","cfit",10,10,1600,1600);
  cfit->Divide(4,4);
  
  for(int i = 0; i < 16; i++)
  {
    cfit->cd(i+1);
    cfit->cd(i+1)->SetLeftMargin(0.15);  
    cfit->cd(i+1)->SetBottomMargin(0.15);
    cfit->cd(i+1)->SetTicks(1,1);
    cfit->cd(i+1)->SetGrid(0,0);
  }
 

  for(int i = 0; i < 4; i++)
  {
    cfit->cd(i+1);
    //cfit->cd(i+1)->SetLogx();
    gslope[i]->SetMarkerStyle(20);
    gslope[i]->SetMarkerSize(1.2);
    gslope[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    gslope[i]->GetXaxis()->SetTitle("Re(#rho_{10})-Re(#rho_{0-1})");
    gslope[i]->GetYaxis()->SetTitle(Form("d(Data/%s)/dp_{T}",simmode.c_str()));
    gslope[i]->SetMarkerSize(1.2);
    gslope[i]->Draw("APE"); 
  }

  for(int i = 0; i < 4; i++)
  {
    cfit->cd(i+5);
    //cfit->cd(i+5)->SetLogx();
    gintercept[i]->SetMarkerStyle(20);
    gintercept[i]->SetMarkerSize(1.2);
    gintercept[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    gintercept[i]->GetXaxis()->SetTitle("Re(#rho_{10})-Re(#rho_{0-1})");
    gintercept[i]->GetYaxis()->SetTitle(Form("(Data/%s) intercept",simmode.c_str()));
    gintercept[i]->Draw("APE"); 
  }

  for(int i = 0; i < 4; i++)
  {
    cfit->cd(i+9);
    //cfit->cd(i+5)->SetLogx();
    gmean[i]->SetMarkerStyle(20);
    gmean[i]->SetMarkerSize(1.2);
    gmean[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    gmean[i]->GetXaxis()->SetTitle("Re(#rho_{10})-Re(#rho_{0-1})");
    gmean[i]->GetYaxis()->SetTitle("<cos(2#phi*-2#phi)>");
    gmean[i]->Draw("APE"); 
  }
  for(int i = 0; i < 4; i++)
  {
    cfit->cd(i+13);
    //cfit->cd(i+1)->SetLogx();
    gglobal[i]->SetMarkerStyle(20);
    gglobal[i]->SetMarkerSize(1.2);
    gglobal[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    gglobal[i]->GetXaxis()->SetTitle("Re(#rho_{10})-Re(#rho_{0-1})^{helicity} input");
    gglobal[i]->GetYaxis()->SetTitle("#Delta#rho_{00}^{global}");
    gglobal[i]->SetMarkerSize(1.2);
    gglobal[i]->Draw("APE"); 
  }
  cfit->SaveAs(Form("figures/%s/%s/pTstudy/FitParams_%s_realcomparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[defaultfile].c_str(),simmode.c_str()));  


  outputname = Form("figures/%s/%s/pTstudy/Comparison_%s_realcomparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[defaultfile].c_str(),simmode.c_str());  
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 
 
  c->Print(outputstart.c_str());
  TLegend *leg = new TLegend(0.4,0.7,0.6,0.9);
  for(int ipt = 0; ipt < 4; ipt++)
  {
    for(int i = 0; i < nvar; i++)
    { 
      c->cd(i+1); 
      h_mCounts[ipt]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
      h_mCounts[ipt]->GetYaxis()->SetTitle("Counts");
      h_mCounts[ipt]->SetMarkerStyle(20);
      h_mCounts[ipt]->SetMarkerColor(kOrange+7);
      h_mCounts[ipt]->SetLineColor(kOrange+7);

      int min = h_mCounts[ipt]->GetMinimum();
      int max = h_mCounts[ipt]->GetMaximum();
      if(h_mCounts_RC[ipt][i]->GetMinimum() < min) min = h_mCounts_RC[ipt][i]->GetMinimum();
      if(h_mCounts_RC[ipt][i]->GetMaximum() > max) max = h_mCounts_RC[ipt][i]->GetMaximum();
      h_mCounts_RC[ipt][i]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);


      h_mCounts_RC[ipt][i]->SetTitle(Form("Re(#rho_{10})-Re(#rho_{0-1})=%s, %1.1f<p_{T}<%1.1f, 20-60 Centrality",sigmay_text[i].c_str(),vmsa::pt_low[energy][ipt+2],vmsa::pt_up[energy][ipt+2]));
      h_mCounts_RC[ipt][i]->SetMarkerStyle(24);
      h_mCounts_RC[ipt][i]->SetMarkerColor(kBlack);
      h_mCounts_RC[ipt][i]->SetLineColor(kBlack);
      h_mCounts_RC[ipt][i]->Draw("pE");
      h_mCounts[ipt]->Draw("pE same");

      if(i == 0 && ipt == 0)
      {
        leg->AddEntry(h_mCounts[ipt],"Data","p");
        leg->AddEntry(h_mCounts_RC[ipt][i],simmode.c_str(),"p");
      }
      leg->Draw("same");
    }
    c->Update();
    c->Print(outputname.c_str());
  }
  c->Print(outputstop.c_str());


}

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nPt = g_rho->GetN();
  TBox *bSys[nPt];
  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot sys errors
  {
    double pt, rho;
    g_rho->GetPoint(i_pt,pt,rho);
    double err = g_rho->GetErrorYhigh(i_pt);

    bSys[i_pt] = new TBox(pt-0.08,rho-err,pt+0.08,rho+err);
    bSys[i_pt]->SetFillColor(0);
    bSys[i_pt]->SetFillStyle(0);
    bSys[i_pt]->SetLineStyle(1);
    bSys[i_pt]->SetLineWidth(1);
    bSys[i_pt]->SetLineColor(plot_color);
    bSys[i_pt]->Draw("l Same");
  }
}


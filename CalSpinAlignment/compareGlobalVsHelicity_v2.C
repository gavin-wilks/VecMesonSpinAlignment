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

void compareGlobalVsHelicity_v2(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 23, std::string simmode = "Mc")//defaultF = 0 is BESII, defaultF = 1 is BESI
{

  const int start = 0;
  const int stop  = 7;

  const int nvar = 7;
  std::string sigmay_text[nvar] = {"0.2333","0.2667","0.3000","0.3333","0.3667","0.4000","0.4333"};
  float sigmay[nvar] = {0.2333,0.2667,0.3000,0.3333,0.3667,0.4000,0.4333};
  std::string sigmay_text_offdiag[nvar] = {"-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3"};
  float sigmay_offdiag[nvar] = {-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3};

  const int defaultfile = 3;

  const int nopts = 2;

  std::string date[nopts] = {"20240809","20240808"}; // nov2, prelimv2
  std::string v2[nopts] = {"nov2","prelimv2"};

  std::string spectra[nopts][5][nvar];
  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int i = start; i < stop; i++)
    {
      spectra[iopt][0][i] = Form("CEPT_O2_NoYSpec_%s_phieff_phiTPC_%s_2_ysig1000_rho0.3333_rerho0.0_imrho0.0_r0.0_i0.0_hrho%s_hrerho0.0_himrho0.0_hr0.0_hi0.0"   ,v2[iopt].c_str(),date[iopt].c_str(),sigmay_text[i].c_str());
      spectra[iopt][1][i] = Form("CEPT_O2_NoYSpec_%s_phieff_phiTPC_%s_2_ysig1000_rho0.3333_rerho0.0_imrho0.0_r0.0_i0.0_hrho0.3333_hrerho0.0_himrho0.0_hr%s_hi0.0",v2[iopt].c_str(),date[iopt].c_str(),sigmay_text_offdiag[i].c_str());
      spectra[iopt][2][i] = Form("CEPT_O2_NoYSpec_%s_phieff_phiTPC_%s_2_ysig1000_rho0.3333_rerho0.0_imrho0.0_r0.0_i0.0_hrho0.3333_hrerho0.0_himrho0.0_hr0.0_hi%s",v2[iopt].c_str(),date[iopt].c_str(),sigmay_text_offdiag[i].c_str());
      spectra[iopt][3][i] = Form("CEPT_O2_NoYSpec_%s_phieff_phiTPC_%s_2_ysig1000_rho0.3333_rerho0.0_imrho0.0_r0.0_i0.0_hrho0.3333_hrerho%s_himrho0.0_hr0.0_hi0.0",v2[iopt].c_str(),date[iopt].c_str(),sigmay_text_offdiag[i].c_str());
      spectra[iopt][4][i] = Form("CEPT_O2_NoYSpec_%s_phieff_phiTPC_%s_2_ysig1000_rho0.3333_rerho0.0_imrho0.0_r0.0_i0.0_hrho0.3333_hrerho0.0_himrho%s_hr0.0_hi0.0",v2[iopt].c_str(),date[iopt].c_str(),sigmay_text_offdiag[i].c_str());
    }
  }
  
  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  
  string inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_%sPhiPtSys_%s_PolySys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  //TFile *File_InPut = TFile::Open(inputfile.c_str());
  //TH1D *h_mCounts[4];  

  TFile *File_InPutRC[nopts][5][nvar]; 
  //TH1D *h_mCounts_RC[4][nvar];  
  //TH1D *h_mCounts_Cos_RC[4][nvar];  
  //TH1D *h_mCounts_CosH_RC[4][nvar];  
  TH2D *h_m2D_RC[nopts][5][4][nvar];  
  TH2D *h_m2D_RC_Global[nopts][5][4][nvar];  
  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int irho = 0; irho < 5; irho++)
    {
      for(int i_pt = 2; i_pt < 6; ++i_pt) // pt loop
      {
        //string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
        //h_mCounts[i_pt-2] = (TH1D*) File_InPut->Get(KEY_counts.c_str());

        //float data = h_mCounts[i_pt-2]->Integral(1,10);

        for(int i = start; i < stop; i++)
        {
          string inputfileRC = Form("effaccfiles/%s/%s/%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[iopt][irho][i].c_str(),vmsa::mBeamEnergy[energy].c_str());
          File_InPutRC[iopt][irho][i] = TFile::Open(inputfileRC.c_str());
          //string KEY_counts_RC = Form("h_m%sEffPhiS_Cent_9_Pt_%d",simmode.c_str(),i_pt);
          //string KEY_counts_Cos_RC = Form("h_m%sEffCos_Cent_9_Pt_%d",simmode.c_str(),i_pt);
          //string KEY_counts_CosH_RC = Form("h_m%sEffCosH_Cent_9_Pt_%d",simmode.c_str(),i_pt);
          string KEY_2D_RC = Form("h_m%sEffCosHPhiPrime_Cent_9_Pt_%d",simmode.c_str(),i_pt);
          string KEY_2D_RC_Global = Form("h_m%sEffCosPhiPrime_Cent_9_Pt_%d",simmode.c_str(),i_pt);
          string KEY_2D = Form("h_m%sEffCosHPhiPrime_Cent_9_Pt_%d_%d_%d",simmode.c_str(),i_pt,irho,iopt);
          string KEY_2D_Global = Form("h_m%sEffCosPhiPrime_Cent_9_Pt_%d_%d_%d",simmode.c_str(),i_pt,irho,iopt);
          h_m2D_RC[iopt][irho][i_pt-2][i] = (TH2D*) ((TH2D*) File_InPutRC[iopt][irho][i]->Get(KEY_2D_RC.c_str()))->Clone(KEY_2D.c_str());
          h_m2D_RC_Global[iopt][irho][i_pt-2][i] = (TH2D*) ((TH2D*) File_InPutRC[iopt][irho][i]->Get(KEY_2D_RC_Global.c_str()))->Clone(KEY_2D_Global.c_str());
          //h_mCounts_CosH_RC[i_pt-2][i] = (TH1D*) ((TH1D*) File_InPutRC[i]->Get(KEY_counts_CosH_RC.c_str()))->Clone();
          //h_mCounts_RC[i_pt-2][i] = (TH1D*) ((TH1D*) File_InPutRC[i]->Get(KEY_counts_RC.c_str()))->Clone();
          //h_mCounts_Cos_RC[i_pt-2][i] = (TH1D*) ((TH1D*) File_InPutRC[i]->Get(KEY_counts_Cos_RC.c_str()))->Clone();
          //float rc   = h_mCounts_RC[i_pt-2][i]->Integral(1,10);

          //h_mCounts_RC[i_pt-2][i]->Scale(data/rc);

          //h_mCounts_Ratio[i_pt-2][i] = (TH1D*) h_mCounts[i_pt-2]->Clone();
          //h_mCounts_Ratio[i_pt-2][i]->Divide(h_mCounts_RC[i_pt-2][i]);

          h_m2D_RC[iopt][irho][i_pt-2][i]->SetDirectory(0);
          h_m2D_RC_Global[iopt][irho][i_pt-2][i]->SetDirectory(0);

          File_InPutRC[iopt][irho][i]->Close();
          cout << inputfileRC << endl;
        }
      }
    }
  }
  //TGraphAsymmErrors *gslope[4];
  //TGraphAsymmErrors *gintercept[4]; 
  //TGraphAsymmErrors *gmean[4]; 
  //TGraphAsymmErrors *gglobal[4]; 
  TGraphAsymmErrors *g2D[nopts][5][5][4]; 
  TGraphAsymmErrors *g2D_Global[nopts][5][5][4]; 
  //TGraphAsymmErrors *g2D_GlobalFromHelicity[5][4]; 
  //TGraphAsymmErrors *g2D_GlobalFromHelicityVsGlobal[5][4]; 
  //TF1 *fits[4][nvar];
  //TF1 *fits1D[4][nvar];
  //TF1 *fits1Dcos[4][nvar];
  TF1 *fits2D[nopts][5][4][nvar];
  TF1 *fits2D_Global[nopts][5][4][nvar];
  //TF1 *fits2D_GlobalFromHelicity[4][nvar];
  float pt[4] = {1.5,2.1,2.7,3.6};

  //TCanvas *c = new TCanvas("c","c",10,10,1200,1200);
  //c->Divide(3,3);

 
  //c->Print(outputstart.c_str());

  for(int iopt = 0; iopt < nopts; iopt++)
  {
    for(int irho = 0; irho < 5; irho++)
    {
      for(int ipt = 0; ipt < 4; ipt++)
      {
        //gslope[ipt] = new TGraphAsymmErrors();
        //gintercept[ipt] = new TGraphAsymmErrors();
        //gmean[ipt] = new TGraphAsymmErrors();
        //gglobal[ipt] = new TGraphAsymmErrors();
        for(int par = 0; par < 5; par++) 
        {
          g2D[iopt][irho][par][ipt] = new TGraphAsymmErrors();
          g2D_Global[iopt][irho][par][ipt] = new TGraphAsymmErrors();
          //g2D_GlobalFromHelicity[par][ipt] = new TGraphAsymmErrors();
          //g2D_GlobalFromHelicityVsGlobal[par][ipt] = new TGraphAsymmErrors();
        }   
        for(int i = start; i < stop; i++)
        {
          //c->cd(i+1);
          //c->cd(i+1)->SetLeftMargin(0.15);
          //c->cd(i+1)->SetBottomMargin(0.15);
          //c->cd(i+1)->SetTicks(1,1);
          //c->cd(i+1)->SetGrid(0,0);
            
          //double mean = h_mCounts_RC[ipt][i]->GetMean();
          //double mean_err = h_mCounts_RC[ipt][i]->GetMeanError();

          //h_mCounts_Ratio[ipt][i]->SetTitle(Form("#rho_{00}=%s, %1.1f<p_{T}<%1.1f, 20-60 Centrality",sigmay_text[i].c_str(),vmsa::pt_low[energy][ipt+2],vmsa::pt_up[energy][ipt+2]));
          //h_mCounts_Ratio[ipt][i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
          //h_mCounts_Ratio[ipt][i]->GetYaxis()->SetTitle(Form("(Data)/(%s from Simulation)",simmode.c_str()));
          //h_mCounts_Ratio[ipt][i]->Draw("pE");
         
          //fits[ipt][i] = new TF1(Form("fit_%d_%d",ipt,i),"[0]+[1]*x",-1.0,1.0);
          //h_mCounts_Ratio[ipt][i]->Fit(fits[ipt][i],"NMI");
          //fits[ipt][i]->SetLineColor(kRed);
          //fits[ipt][i]->Draw("l same");

          //cout << "About to Fit for 1D Spin Density Matrix" << endl;
          //fits1D[ipt][i] = new TF1(Form("fit1D_%d_%d",ipt,i),SpinDensity,-1.0,1.0,2);
          //fits1Dcos[ipt][i] = new TF1(Form("fit1Dcos_%d_%d",ipt,i),SpinDensity,-1.0,1.0,2);
          //h_mCounts_CosH_RC[ipt][i]->Fit(fits1D[ipt][i],"NMI");
          //h_mCounts_Cos_RC[ipt][i]->Fit(fits1Dcos[ipt][i],"NMI");
          //cout << "Finished Fit for 1D Spin Density Matrix" << endl;

          //cout << "iopt = " << iopt << ", irho = " << irho << ", ipt = " << ipt << ", i = " << i << endl;
          fits2D[iopt][irho][ipt][i] = new TF2(Form("fit2D_%d_%d_%d_%d",ipt,i,irho,iopt),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
          fits2D[iopt][irho][ipt][i]->SetParameter(0,1./3.);
          fits2D[iopt][irho][ipt][i]->SetParLimits(0,0.0,1.0);
          fits2D[iopt][irho][ipt][i]->SetParameter(5,h_m2D_RC[iopt][irho][ipt][i]->GetMaximum());
          h_m2D_RC[iopt][irho][ipt][i]->Fit(fits2D[iopt][irho][ipt][i],"NMI");
          for(int par = 0; par < 5; par++)
          {
            if(par == 0 && irho == 0) g2D[iopt][irho][par][ipt]->SetPoint(i-start,sigmay[i]-1./3.+0.02*1./3.*(0.5-float(iopt)),    fits2D[iopt][irho][ipt][i]->GetParameter(0)-1./3.);
            if(par != 0 && irho == 0) g2D[iopt][irho][par][ipt]->SetPoint(i-start,sigmay[i]-1./3.+0.02*(0.5-float(iopt)),    fits2D[iopt][irho][ipt][i]->GetParameter(par));
            if(par == 0 && irho != 0) g2D[iopt][irho][par][ipt]->SetPoint(i-start,sigmay_offdiag[i]+0.02*(0.5-float(iopt)),  fits2D[iopt][irho][ipt][i]->GetParameter(0)-1./3.);
            if(par != 0 && irho != 0) g2D[iopt][irho][par][ipt]->SetPoint(i-start,sigmay_offdiag[i]+0.02*(0.5-float(iopt)),  fits2D[iopt][irho][ipt][i]->GetParameter(par));
            g2D[iopt][irho][par][ipt]->SetPointError(i-start,0.0,0.0,fits2D[iopt][irho][ipt][i]->GetParError(par),fits2D[iopt][irho][ipt][i]->GetParError(par));
          }

          fits2D_Global[iopt][irho][ipt][i] = new TF2(Form("fit2D_Global_%d_%d_%d_%d",ipt,i,irho,iopt),SpinDensity2Dcos,-1.0,1.0,0,2.0*TMath::Pi(),6);
          fits2D_Global[iopt][irho][ipt][i]->SetParameter(0,1./3.);
          fits2D_Global[iopt][irho][ipt][i]->SetParLimits(0,0.0,1.0);
          fits2D_Global[iopt][irho][ipt][i]->SetParameter(5,h_m2D_RC_Global[iopt][irho][ipt][i]->GetMaximum());
          h_m2D_RC_Global[iopt][irho][ipt][i]->Fit(fits2D_Global[iopt][irho][ipt][i],"NMI");
          for(int par = 0; par < 5; par++)
          {
            if(par == 0 && irho == 0) g2D_Global[iopt][irho][par][ipt]->SetPoint(i-start,sigmay[i]-1./3.+0.02*1./3.*(0.5-float(iopt)),                       fits2D_Global[iopt][irho][ipt][i]->GetParameter(0)-1./3.);
            if(par != 0 && irho == 0) g2D_Global[iopt][irho][par][ipt]->SetPoint(i-start,sigmay[i]-1./3.+0.02*(0.5-float(iopt)),                       fits2D_Global[iopt][irho][ipt][i]->GetParameter(par));
            if(par == 0 && irho != 0) g2D_Global[iopt][irho][par][ipt]->SetPoint(i-start,sigmay_offdiag[i]+0.02*(0.5-float(iopt)),                     fits2D_Global[iopt][irho][ipt][i]->GetParameter(0)-1./3.);
            if(par != 0 && irho != 0) g2D_Global[iopt][irho][par][ipt]->SetPoint(i-start,sigmay_offdiag[i]+0.02*(0.5-float(iopt)),                     fits2D_Global[iopt][irho][ipt][i]->GetParameter(par));
            g2D_Global[iopt][irho][par][ipt]->SetPointError(i-start,0.0,0.0,fits2D_Global[iopt][irho][ipt][i]->GetParError(par),fits2D_Global[iopt][irho][ipt][i]->GetParError(par));
          }

          //fits2D_GlobalFromHelicity[ipt][i] = new TF2(Form("fit2D_%d_%d",ipt,i),SpinDensity2Dcos_GlobalFromHelicity,-1.0,1.0,-TMath::Pi(),TMath::Pi(),6);
          //fits2D_GlobalFromHelicity[ipt][i]->SetParameter(0,1./3.);
          //fits2D_GlobalFromHelicity[ipt][i]->SetParLimits(0,0.0,1.0);
          //fits2D_GlobalFromHelicity[ipt][i]->SetParameter(5,h_m2D_RC[ipt][i]->GetMaximum());
          //h_m2D_RC[ipt][i]->Fit(fits2D_GlobalFromHelicity[ipt][i],"NMI");
          //for(int par = 0; par < 5; par++)
          //{
          //  if(par == 0) g2D_GlobalFromHelicity[par][ipt]->SetPoint(i-start,sigmay[i]-1./3.,fits2D_GlobalFromHelicity[ipt][i]->GetParameter(0)-1./3.);
          //  else g2D_GlobalFromHelicity[par][ipt]->SetPoint(i-start,sigmay[i]-1./3.,fits2D_GlobalFromHelicity[ipt][i]->GetParameter(par));
          //  g2D_GlobalFromHelicity[par][ipt]->SetPointError(i-start,0.0,0.0,fits2D_GlobalFromHelicity[ipt][i]->GetParError(par),fits2D_GlobalFromHelicity[ipt][i]->GetParError(par));
          //}
          //for(int par = 0; par < 5; par++)
          //{
          //  if(par == 0) g2D_GlobalFromHelicityVsGlobal[par][ipt]->SetPoint(i-start,fits2D_Global[ipt][i]->GetParameter(0)-1./3.,fits2D_GlobalFromHelicity[ipt][i]->GetParameter(0)-1./3.);
          //  else g2D_GlobalFromHelicityVsGlobal[par][ipt]->SetPoint(i-start,fits2D_Global[ipt][i]->GetParameter(par),fits2D_GlobalFromHelicity[ipt][i]->GetParameter(par));
          //  g2D_GlobalFromHelicityVsGlobal[par][ipt]->SetPointError(i-start,fits2D_Global[ipt][i]->GetParError(par),fits2D_Global[ipt][i]->GetParError(par),fits2D_GlobalFromHelicity[ipt][i]->GetParError(par),fits2D_GlobalFromHelicity[ipt][i]->GetParError(par));
          //}
          
          //gslope[ipt]->SetPoint(i-start,sigmay[i],fits[ipt][i]->GetParameter(1)); 
          //gslope[ipt]->SetPointError(i-start,0.0,0.0,fits[ipt][i]->GetParError(1),fits[ipt][i]->GetParError(1));
          //gintercept[ipt]->SetPoint(i-start,sigmay[i],fits[ipt][i]->GetParameter(0)); 
          //gintercept[ipt]->SetPointError(i-start,0.0,0.0,fits[ipt][i]->GetParError(0),fits[ipt][i]->GetParError(0));
          //gmean[ipt]->SetPoint(i-start,sigmay[i],mean); 
          //gmean[ipt]->SetPointError(i-start,0.0,0.0,mean_err,mean_err);
        
          //gglobal[ipt]->SetPoint(i-start,sigmay[i]-1./3.,fits1Dcos[ipt][i]->GetParameter(0)-1./3.);
          //gglobal[ipt]->SetPointError(i-start,0.0,0.0,fits1Dcos[ipt][i]->GetParError(0),fits1Dcos[ipt][i]->GetParError(0));
        }
        //c->Update();
        //c->Print(outputname.c_str());
      }
      //c->Print(outputstop.c_str());
    }
  }
  TCanvas *cfit2D = new TCanvas("cfit2D","cfit2D",10,10,2000,2000);
  cfit2D->Divide(5,5);
  
  for(int i = 0; i < 25; i++)
  {
    cfit2D->cd(i+1);
    cfit2D->cd(i+1)->SetLeftMargin(0.20);  
    cfit2D->cd(i+1)->SetBottomMargin(0.15);
    cfit2D->cd(i+1)->SetTicks(1,1);
    cfit2D->cd(i+1)->SetGrid(0,0);
  }
 
  int color[nopts] = {kBlue, kOrange+7}; 
  std::string label[nopts] = {"v_{2} OFF","v_{2} ON (Prelim)"};

  float ptedges[5] = {1.2,1.8,2.4,3.0,4.2};
  string param[5] = {"#Delta#rho_{00}","Re(#rho_{10})-Re(#rho_{0-1})","Im(#rho_{10})-Im(#rho_{0-1})","Re(#rho_{1-1})","Im(#rho_{1-1})"};

  std::string outputname = Form("figures/%s/%s/pTstudy/GlobalvsHelicity_parametercomparison_%s_v2.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str());
  std::string outputstart = Form("%s[",outputname.c_str()); 
  std::string outputstop = Form("%s]",outputname.c_str()); 

  cfit2D->Print(outputstart.c_str());

  TLegend *leg = new TLegend(0.4,0.75,0.7,0.9);
  for(int i = 0; i < 4; i++)
  {
    for(int irho = 0; irho < 5; irho++)
    {
      for(int par = 0; par < 5; par++)
      {
        cfit2D->cd(1+par*5+irho);
        //TLegend *leg = new TLegend(0.4,0.6,0.6,0.8);
        double max = -999999;
        double min = 999999;
        for(int iopt = 0; iopt < nopts; iopt++)
        {
          double tmax = g2D_Global[iopt][irho][par][i]->GetHistogram()->GetMaximum();   
          double tmin = g2D_Global[iopt][irho][par][i]->GetHistogram()->GetMinimum(); 

          if(tmax > max) max = tmax;
          if(tmin < min) min = tmin;

          cout << "max = " << max << ", min = " << min << endl;
        }
        for(int iopt = 0; iopt < nopts; iopt++)
        {

          //cfit->cd(i+5)->SetLogx();
          //g2D[irho][par][i]->Print();
          g2D_Global[iopt][irho][par][i]->SetMarkerStyle(20);
          g2D_Global[iopt][irho][par][i]->SetMarkerSize(1.2);
          g2D_Global[iopt][irho][par][i]->SetMarkerColor(color[iopt]);
          g2D_Global[iopt][irho][par][i]->SetLineColor(color[iopt]);
          g2D_Global[iopt][irho][par][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
          g2D_Global[iopt][irho][par][i]->GetXaxis()->SetTitle(Form("%s^{Helicity} Input",param[irho].c_str()));
          g2D_Global[iopt][irho][par][i]->GetYaxis()->SetTitle(Form("%s^{Global} from 2D Fit",param[par].c_str()));
          if(min < 0.0 && max < 0.0) g2D_Global[iopt][irho][par][i]->GetYaxis()->SetRangeUser(1.5*min,0.5*max);
          if(min < 0.0 && max > 0.0) g2D_Global[iopt][irho][par][i]->GetYaxis()->SetRangeUser(1.5*min,1.5*max);
          if(min > 0.0 && max > 0.0) g2D_Global[iopt][irho][par][i]->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

          if(iopt == 0) g2D_Global[iopt][irho][par][i]->Draw("APE"); 
          else          g2D_Global[iopt][irho][par][i]->Draw("PE same"); 
            
          if(i == 0 && irho == 0 && par == 0) leg->AddEntry(g2D_Global[iopt][irho][par][i],label[iopt].c_str(),"p");
        }
        leg->Draw("same");
      }
    }
    cfit2D->Print(outputname.c_str());
    cfit2D->Update();
  }
  cfit2D->Print(outputstop.c_str());

  outputname = Form("figures/%s/%s/pTstudy/HelicityvsHelicity_parametercomparison_%s_v2.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),simmode.c_str());
  outputstart = Form("%s[",outputname.c_str()); 
  outputstop = Form("%s]",outputname.c_str()); 

  cfit2D->Print(outputstart.c_str());

  for(int i = 0; i < 4; i++)
  {
    for(int irho = 0; irho < 5; irho++)
    {
      for(int par = 0; par < 5; par++)
      {
        cfit2D->cd(1+par*5+irho);
        double max = -999999;
        double min = 999999;
        for(int iopt = 0; iopt < nopts; iopt++)
        {
          double tmax = g2D[iopt][irho][par][i]->GetHistogram()->GetMaximum();   
          double tmin = g2D[iopt][irho][par][i]->GetHistogram()->GetMinimum(); 

          if(tmax > max) max = tmax;
          if(tmin < min) min = tmin;

          cout << "max = " << max << ", min = " << min << endl;
        }
        for(int iopt = 0; iopt < nopts; iopt++)
        {
          //cfit->cd(i+5)->SetLogx();
          //g2D[irho][par][i]->Print();
          g2D[iopt][irho][par][i]->SetMarkerStyle(20);
          g2D[iopt][irho][par][i]->SetMarkerSize(1.2);
          g2D[iopt][irho][par][i]->SetMarkerColor(color[iopt]);
          g2D[iopt][irho][par][i]->SetLineColor(color[iopt]);
          g2D[iopt][irho][par][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
          g2D[iopt][irho][par][i]->GetXaxis()->SetTitle(Form("%s^{Helicity} Input",param[irho].c_str()));
          g2D[iopt][irho][par][i]->GetYaxis()->SetTitle(Form("%s^{Helicity} from 2D Fit",param[par].c_str()));
          if(min < 0.0 && max < 0.0) g2D[iopt][irho][par][i]->GetYaxis()->SetRangeUser(1.5*min,0.5*max);
          if(min < 0.0 && max > 0.0) g2D[iopt][irho][par][i]->GetYaxis()->SetRangeUser(1.5*min,1.5*max);
          if(min > 0.0 && max > 0.0) g2D[iopt][irho][par][i]->GetYaxis()->SetRangeUser(0.5*min,1.5*max);

          if(iopt == 0) g2D[iopt][irho][par][i]->Draw("APE"); 
          else          g2D[iopt][irho][par][i]->Draw("PE same"); 
            
        }
        leg->Draw("same");
      }
    }
    cfit2D->Print(outputname.c_str());
    cfit2D->Update();
  }
  cfit2D->Print(outputstop.c_str());

  //for(int i = 0; i < 4; i++)
  //{
  //  for(int par = 0; par < 5; par++)
  //  {
  //    cfit2D->cd(1+i*5+par);
  //    //cfit->cd(i+5)->SetLogx();
  //    g2D_Global[par][i]->Print();
  //    g2D_Global[par][i]->SetMarkerStyle(20);
  //    g2D_Global[par][i]->SetMarkerSize(1.2);
  //    g2D_Global[par][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
  //    g2D_Global[par][i]->GetXaxis()->SetTitle("#Delta#rho_{00}");
  //    g2D_Global[par][i]->GetYaxis()->SetTitle(param[par].c_str());
  //    g2D_Global[par][i]->Draw("APE"); 
  //  }
  //}
  //cfit2D->SaveAs(Form("figures/%s/%s/pTstudy/%s_parametercomparison_rho00comparison_Global_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[defaultfile].c_str(),simmode.c_str()));  

  //for(int i = 0; i < 4; i++)
  //{
  //  for(int par = 0; par < 5; par++)
  //  {
  //    cfit2D->cd(1+i*5+par);
  //    //cfit->cd(i+5)->SetLogx();
  //    g2D_GlobalFromHelicity[par][i]->Print();
  //    g2D_GlobalFromHelicity[par][i]->SetMarkerStyle(20);
  //    g2D_GlobalFromHelicity[par][i]->SetMarkerSize(1.2);
  //    g2D_GlobalFromHelicity[par][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
  //    g2D_GlobalFromHelicity[par][i]->GetXaxis()->SetTitle("#Delta#rho_{00}");
  //    g2D_GlobalFromHelicity[par][i]->GetYaxis()->SetTitle(param[par].c_str());
  //    g2D_GlobalFromHelicity[par][i]->Draw("APE"); 
  //  }
  //}
  //cfit2D->SaveAs(Form("figures/%s/%s/pTstudy/%s_parametercomparison_rho00comparison_GlobalFromHelicity_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[defaultfile].c_str(),simmode.c_str()));  

  //for(int i = 0; i < 4; i++)
  //{
  //  for(int par = 0; par < 5; par++)
  //  {
  //    cfit2D->cd(1+i*5+par);
  //    //cfit->cd(i+5)->SetLogx();
  //    g2D_GlobalFromHelicityVsGlobal[par][i]->Print();
  //    g2D_GlobalFromHelicityVsGlobal[par][i]->SetMarkerStyle(20);
  //    g2D_GlobalFromHelicityVsGlobal[par][i]->SetMarkerSize(1.2);
  //    g2D_GlobalFromHelicityVsGlobal[par][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
  //    g2D_GlobalFromHelicityVsGlobal[par][i]->GetXaxis()->SetTitle(Form("Global %s",param[par].c_str()));
  //    g2D_GlobalFromHelicityVsGlobal[par][i]->GetYaxis()->SetTitle(Form("Global from Helicity %s",param[par].c_str()));
  //    g2D_GlobalFromHelicityVsGlobal[par][i]->Draw("APE"); 
  //  }
  //}
  //cfit2D->SaveAs(Form("figures/%s/%s/pTstudy/%s_parametercomparison_rho00comparison_GlobalFromHelicityVsGlobal_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[defaultfile].c_str(),simmode.c_str()));  

  //TCanvas *cfit = new TCanvas("cfit","cfit",10,10,1600,1600);
  //cfit->Divide(4,4);
  //
  //for(int i = 0; i < 16; i++)
  //{
  //  cfit->cd(i+1);
  //  cfit->cd(i+1)->SetLeftMargin(0.15);  
  //  cfit->cd(i+1)->SetBottomMargin(0.15);
  //  cfit->cd(i+1)->SetTicks(1,1);
  //  cfit->cd(i+1)->SetGrid(0,0);
  //}
 

  //for(int i = 0; i < 4; i++)
  //{
  //  cfit->cd(i+1);
  //  //cfit->cd(i+1)->SetLogx();
  //  gslope[i]->SetMarkerStyle(20);
  //  gslope[i]->SetMarkerSize(1.2);
  //  gslope[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
  //  gslope[i]->GetXaxis()->SetTitle("#rho_{00}");
  //  gslope[i]->GetYaxis()->SetTitle(Form("d(Data/%s)/dp_{T}",simmode.c_str()));
  //  gslope[i]->SetMarkerSize(1.2);
  //  gslope[i]->Draw("APE"); 
  //}

  //for(int i = 0; i < 4; i++)
  //{
  //  cfit->cd(i+5);
  //  //cfit->cd(i+5)->SetLogx();
  //  gintercept[i]->SetMarkerStyle(20);
  //  gintercept[i]->SetMarkerSize(1.2);
  //  gintercept[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
  //  gintercept[i]->GetXaxis()->SetTitle("#rho_{00}");
  //  gintercept[i]->GetYaxis()->SetTitle(Form("(Data/%s) intercept",simmode.c_str()));
  //  gintercept[i]->Draw("APE"); 
  //}

  //for(int i = 0; i < 4; i++)
  //{
  //  cfit->cd(i+9);
  //  //cfit->cd(i+5)->SetLogx();
  //  gmean[i]->SetMarkerStyle(20);
  //  gmean[i]->SetMarkerSize(1.2);
  //  gmean[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
  //  gmean[i]->GetXaxis()->SetTitle("#rho_{00}");
  //  gmean[i]->GetYaxis()->SetTitle("<cos(2#phi*-2#phi)>");
  //  gmean[i]->Draw("APE"); 
  //}

  //for(int i = 0; i < 4; i++)
  //{
  //  cfit->cd(i+13);
  //  //cfit->cd(i+1)->SetLogx();
  //  gglobal[i]->SetMarkerStyle(20);
  //  gglobal[i]->SetMarkerSize(1.2);
  //  gglobal[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
  //  gglobal[i]->GetXaxis()->SetTitle("#Delta#rho_{00}^{helicity} input");
  //  gglobal[i]->GetYaxis()->SetTitle("#Delta#rho_{00}^{global}");
  //  gglobal[i]->SetMarkerSize(1.2);
  //  gglobal[i]->Draw("APE"); 
  //}

  //cfit->SaveAs(Form("figures/%s/%s/pTstudy/FitParams_%s_rho00comparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[defaultfile].c_str(),simmode.c_str()));  


  //outputname = Form("figures/%s/%s/pTstudy/Comparison_%s_rho00comparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[defaultfile].c_str(),simmode.c_str());  
  //outputstart = Form("%s[",outputname.c_str()); 
  //outputstop = Form("%s]",outputname.c_str()); 
 
  //c->Print(outputstart.c_str());
  //TLegend *leg = new TLegend(0.4,0.7,0.6,0.9);
  //for(int ipt = 0; ipt < 4; ipt++)
  //{
  //  for(int i = start; i < stop; i++)
  //  { 
  //    c->cd(i+1); 
  //    h_mCounts[ipt]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
  //    h_mCounts[ipt]->GetYaxis()->SetTitle("Counts");
  //    h_mCounts[ipt]->SetMarkerStyle(20);
  //    h_mCounts[ipt]->SetMarkerColor(kOrange+7);
  //    h_mCounts[ipt]->SetLineColor(kOrange+7);

  //    int min = h_mCounts[ipt]->GetMinimum();
  //    int max = h_mCounts[ipt]->GetMaximum();
  //    if(h_mCounts_RC[ipt][i]->GetMinimum() < min) min = h_mCounts_RC[ipt][i]->GetMinimum();
  //    if(h_mCounts_RC[ipt][i]->GetMaximum() > max) max = h_mCounts_RC[ipt][i]->GetMaximum();
  //    h_mCounts_RC[ipt][i]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);


  //    h_mCounts_RC[ipt][i]->SetTitle(Form("#rho_{00}=%s, %1.1f<p_{T}<%1.1f, 20-60 Centrality",sigmay_text[i].c_str(),vmsa::pt_low[energy][ipt+2],vmsa::pt_up[energy][ipt+2]));
  //    h_mCounts_RC[ipt][i]->SetMarkerStyle(24);
  //    h_mCounts_RC[ipt][i]->SetMarkerColor(kBlack);
  //    h_mCounts_RC[ipt][i]->SetLineColor(kBlack);
  //    h_mCounts_RC[ipt][i]->Draw("pE");
  //    h_mCounts[ipt]->Draw("pE same");

  //    if(i == 0 && ipt == 0)
  //    {
  //      leg->AddEntry(h_mCounts[ipt],"Data","p");
  //      leg->AddEntry(h_mCounts_RC[ipt][i],simmode.c_str(),"p");
  //    }
  //    leg->Draw("same");
  //  }
  //  c->Update();
  //  c->Print(outputname.c_str());
  //}
  //c->Print(outputstop.c_str());


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


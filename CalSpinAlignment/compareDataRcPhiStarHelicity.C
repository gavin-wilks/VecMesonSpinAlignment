#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
//#include "../Utilitpt/draw.h"
#include "../../../draw.h"
#include "../Utility/functions.h"
//#include "resolution_pt.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

double finter(double *x, double *par) {
   return TMath::Abs(0 - (par[0] + par[1]*x[0]));
}

void compareDataRcPhiStarHelicity(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 23, std::string simmode = "Rc")//defaultF = 0 is BESII, defaultF = 1 is BESI
{
  const int nvar = 21;
  float inputglobalrho00 = 1./3.;
  //std::string sigmay_text[nvar] = {"0.00","0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95","1.00"};
  //float sigmay[nvar] = {0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00};
  std::string sigmay_text[nvar] = {"0.20","0.205","0.21","0.215","0.22","0.225","0.23","0.235","0.24","0.245","0.25","0.255","0.26","0.265","0.27","0.275","0.28","0.285","0.29","0.295","0.30"};
  float sigmay[nvar] = {0.20,0.205,0.21,0.215,0.22,0.225,0.23,0.235,0.24,0.245,0.25,0.255,0.26,0.265,0.27,0.275,0.28,0.285,0.29,0.295,0.30};

  std::string spectrafile = Form("_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240522_fixed2Dfunc_Spectra_rerho1n10.0_ysigma1000");
 
  std::string spectra[nvar];
  for(int i = 0; i < nvar; i++)
  {
    spectra[i] = Form("%s_rho000.3333_r0.0_i0.0_imrho1n10.0_helicityrho%s",spectrafile.c_str(),sigmay_text[i].c_str());
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
  TH1D *h_mCounts_RC[4][nvar];  
  TH2D *h_m2D_RC[4][nvar];  
  TH1D *h_mCosH_RC[4][nvar];  
  TH1D *h_mCos_RC[4][nvar];  
  TH2D *h_mCosCosH_RC[4][nvar];  
  TH1D *h_mCounts_Ratio[4][nvar];  
  for(int i_pt = 2; i_pt < 6; ++i_pt) // pt loop
  {
    string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    h_mCounts[i_pt-2] = (TH1D*) File_InPut->Get(KEY_counts.c_str());

    float data = h_mCounts[i_pt-2]->Integral(1,10);

    for(int i = 0; i < nvar; i++)
    {
      string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra[i].c_str(),vmsa::mBeamEnergy[energy].c_str());
      File_InPutRC[i] = TFile::Open(inputfileRC.c_str());
      string KEY_counts_RC = Form("h_m%sEffPhiS_Cent_9_Pt_%d",simmode.c_str(),i_pt);
      string KEY_2D_RC = Form("h_m%sEffCosPhiPrime_Cent_9_Pt_%d",simmode.c_str(),i_pt);
      string KEY_CosH_RC = Form("h_m%sEffCosH_Cent_9_Pt_%d",simmode.c_str(),i_pt);
      string KEY_CosCosH_RC = Form("h_m%sEffCosCosH_Cent_9_Pt_%d",simmode.c_str(),i_pt);
      string KEY_Cos_RC = Form("h_m%sEffCos_Cent_9_Pt_%d",simmode.c_str(),i_pt);
      h_mCounts_RC[i_pt-2][i] = (TH1D*) ((TH1D*) File_InPutRC[i]->Get(KEY_counts_RC.c_str()))->Clone();
      h_mCosH_RC[i_pt-2][i] = (TH1D*) ((TH1D*) File_InPutRC[i]->Get(KEY_CosH_RC.c_str()))->Clone();
      h_mCosCosH_RC[i_pt-2][i] = (TH2D*) ((TH1D*) File_InPutRC[i]->Get(KEY_CosCosH_RC.c_str()))->Clone();
      h_mCos_RC[i_pt-2][i] = (TH1D*) ((TH1D*) File_InPutRC[i]->Get(KEY_Cos_RC.c_str()))->Clone();
      h_m2D_RC[i_pt-2][i] = (TH2D*) ((TH2D*) File_InPutRC[i]->Get(KEY_2D_RC.c_str()))->Clone();
      float rc   = h_mCounts_RC[i_pt-2][i]->Integral(1,10);

      h_mCounts_RC[i_pt-2][i]->Scale(data/rc);

      h_mCounts_Ratio[i_pt-2][i] = (TH1D*) h_mCounts[i_pt-2]->Clone();
      h_mCounts_Ratio[i_pt-2][i]->Divide(h_mCounts_RC[i_pt-2][i]);
    }
  }

  TGraphAsymmErrors *gslope[4];
  TGraphAsymmErrors *gintercept[4]; 
  TGraphAsymmErrors *gmean[4]; 
  TGraphAsymmErrors *gmeandiff[4]; 
  TGraphAsymmErrors *gCosH[4]; 
  TGraphAsymmErrors *g1D[4]; 
  TGraphAsymmErrors *ginput[4]; 
  TGraphAsymmErrors *g2D[5][4]; 
  TGraphAsymmErrors *gCosCosH[5][4]; 
  TGraphAsymmErrors *gCos2D[4]; 
  TGraphAsymmErrors *gCosH2D[4]; 
  TF1 *fits[4][nvar];
  TF1 *fits1D[4][nvar];
  TF1 *fitsCosH[4][nvar];
  TF1 *fitsCosCosH[4][nvar];
  TF1 *fits2D[4][nvar];
  float pt[4] = {1.5,2.1,2.7,3.6};

  TCanvas *c = new TCanvas("c","c",10,10,2000,2000);
  c->Divide(5,5);

  std::string outputname = Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_Ratio_%s_Order%d_%s%s_helicityrho00comparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectrafile.c_str(),simmode.c_str());  
  std::string outputstart = Form("%s[",outputname.c_str()); 
  std::string outputstop = Form("%s]",outputname.c_str()); 
 
  c->Print(outputstart.c_str());

  for(int ipt = 0; ipt < 4; ipt++)
  {
    gslope[ipt] = new TGraphAsymmErrors();
    gintercept[ipt] = new TGraphAsymmErrors();
    gmean[ipt] = new TGraphAsymmErrors();
    gmeandiff[ipt] = new TGraphAsymmErrors();
    gCosH[ipt] = new TGraphAsymmErrors();
    g1D[ipt] = new TGraphAsymmErrors();
    ginput[ipt] = new TGraphAsymmErrors();
    gCos2D[ipt] = new TGraphAsymmErrors();
    gCosH2D[ipt] = new TGraphAsymmErrors();
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
      double meandata = h_mCounts[ipt]->GetMean();
      double meandata_err = h_mCounts[ipt]->GetMeanError();

      h_mCounts_Ratio[ipt][i]->SetTitle(Form("#rho_{00}=%s, %1.1f<p_{T}<%1.1f, 20-60 Centrality",sigmay_text[i].c_str(),vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mCounts_Ratio[ipt][i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
      h_mCounts_Ratio[ipt][i]->GetYaxis()->SetTitle(Form("(Data)/(%s from Simulation)",simmode.c_str()));
      h_mCounts_Ratio[ipt][i]->Draw("pE");
     
      fits[ipt][i] = new TF1(Form("fit_%d_%d",ipt,i),"[0]+[1]*x",-1.0,1.0);
      fits1D[ipt][i] = new TF1(Form("fitCos_%d_%d",ipt,i),SpinDensity,-1.0,1.0,2);
      fitsCosH[ipt][i] = new TF1(Form("fitCosH_%d_%d",ipt,i),SpinDensity,-1.0,1.0,2);
      fitsCosCosH[ipt][i] = new TF2(Form("fitCosCosH_%d_%d",ipt,i),SpinDensityCosCosH,-1.0,1.0,-1.0,1.0,3);
      fitsCosH[ipt][i]->SetLineColor(kRed);
      h_mCounts_Ratio[ipt][i]->Fit(fits[ipt][i],"NMI");
      h_mCosH_RC[ipt][i]->Fit(fitsCosH[ipt][i],"NMI");
      h_mCosCosH_RC[ipt][i]->Fit(fitsCosCosH[ipt][i],"NMI");
      h_mCos_RC[ipt][i]->Fit(fits1D[ipt][i],"NMI");
      fits[ipt][i]->SetLineColor(kRed);
      fits[ipt][i]->Draw("l same");

      fits2D[ipt][i] = new TF2(Form("fit2D_%d_%d",ipt,i),SpinDensity2Dcos,-1.0,1.0,0,2.0*TMath::Pi(),6);
      h_m2D_RC[ipt][i]->Fit(fits2D[ipt][i],"NMI");
      for(int par = 0; par < 5; par++)
      {
        if(par == 0) g2D[par][ipt]->SetPoint(i,sigmay[i],fits2D[ipt][i]->GetParameter(0));
        else g2D[par][ipt]->SetPoint(i,sigmay[i],fits2D[ipt][i]->GetParameter(par));
        g2D[par][ipt]->SetPointError(i,0.0,0.0,fits2D[ipt][i]->GetParError(par),fits2D[ipt][i]->GetParError(par));
      }
      
      gslope[ipt]->SetPoint(i,sigmay[i],fits[ipt][i]->GetParameter(1)); 
      gslope[ipt]->SetPointError(i,0.0,0.0,fits[ipt][i]->GetParError(1),fits[ipt][i]->GetParError(1));
      gintercept[ipt]->SetPoint(i,sigmay[i],fits[ipt][i]->GetParameter(0)); 
      gintercept[ipt]->SetPointError(i,0.0,0.0,fits[ipt][i]->GetParError(0),fits[ipt][i]->GetParError(0));
      gmean[ipt]->SetPoint(i,sigmay[i],mean); 
      gmean[ipt]->SetPointError(i,0.0,0.0,mean_err,mean_err);
      gmeandiff[ipt]->SetPoint(i,sigmay[i],meandata-mean); 
      double differr = TMath::Sqrt(meandata_err*meandata_err+mean_err*mean_err);
      gmeandiff[ipt]->SetPointError(i,0.0,0.0,differr,differr);
      gCosH[ipt]->SetPoint(i,sigmay[i],fitsCosH[ipt][i]->GetParameter(0)); 
      gCosH[ipt]->SetPointError(i,0.0,0.0,fitsCosH[ipt][i]->GetParError(0),fitsCosH[ipt][i]->GetParError(0));
      g1D[ipt]->SetPoint(i,sigmay[i],fits1D[ipt][i]->GetParameter(0)-1./3.); 
      g1D[ipt]->SetPointError(i,0.0,0.0,fits1D[ipt][i]->GetParError(0),fits1D[ipt][i]->GetParError(0));
      ginput[ipt]->SetPoint(i,sigmay[i],fitsCosH[ipt][i]->GetParameter(0)-sigmay[i]); 
      ginput[ipt]->SetPointError(i,0.0,0.0,fitsCosH[ipt][i]->GetParError(0),fitsCosH[ipt][i]->GetParError(0));
      // 2D fits to Cos(theta*) global and Cos(theta*) helicity 
      cout << "Filling 2D fits" << endl;
      gCos2D[ipt]->SetPoint(i,sigmay[i],fitsCosCosH[ipt][i]->GetParameter(0)-inputglobalrho00); 
      gCos2D[ipt]->SetPointError(i,0.0,0.0,fitsCosCosH[ipt][i]->GetParError(0),fitsCosCosH[ipt][i]->GetParError(0));
      gCosH2D[ipt]->SetPoint(i,sigmay[i],fitsCosCosH[ipt][i]->GetParameter(1)-sigmay[i]); 
      gCosH2D[ipt]->SetPointError(i,0.0,0.0,fitsCosCosH[ipt][i]->GetParError(1),fitsCosCosH[ipt][i]->GetParError(1));
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
  string param[5] = {"#rho_{00}","Re(#rho_{10})-Re(#rho_{0-1})","Im(#rho_{10})-Im(#rho_{0-1})","Re(#rho_{1-1})","Im(#rho_{1-1})"};

  for(int i = 0; i < 4; i++)
  {
    for(int par = 0; par < 5; par++)
    {
      cfit2D->cd(1+i*5+par);
      //cfit->cd(i+5)->SetLogx();
      g2D[par][i]->SetMarkerStyle(20);
      g2D[par][i]->SetMarkerSize(1.2);
      g2D[par][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
      g2D[par][i]->GetXaxis()->SetTitle("local #rho_{00}");
      g2D[par][i]->GetYaxis()->SetTitle(param[par].c_str());
      g2D[par][i]->Draw("APE"); 
    }
  }
  cfit2D->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_parametercomparison_%s_Order%d_%s%s_helicityrho00comparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectrafile.c_str(),simmode.c_str()));  

  TCanvas *cfit = new TCanvas("cfit","cfit",10,10,1600,2800);
  cfit->Divide(4,7);
  
  for(int i = 0; i < 28; i++)
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
    gslope[i]->GetXaxis()->SetTitle("local #rho_{00}");
    gslope[i]->GetYaxis()->SetTitle(Form("d(Data/%s)/d(cos(2#phi*-2#phi))",simmode.c_str()));
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
    gintercept[i]->GetXaxis()->SetTitle("local #rho_{00}");
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
    gmean[i]->GetXaxis()->SetTitle("local #rho_{00}");
    gmean[i]->GetYaxis()->SetTitle("<cos(2#phi*-2#phi)>");
    gmean[i]->Draw("APE"); 
  }

  double datarho[4] = {0.307032,0.295729,0.288274,0.246577}; 
  double intersection[4];

  for(int i = 0; i < 4; i++)
  {
    cfit->cd(i+13);
    //cfit->cd(i+5)->SetLogx();
    gmeandiff[i]->SetMarkerStyle(20);
    gmeandiff[i]->SetMarkerSize(1.2);
    gmeandiff[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    gmeandiff[i]->GetXaxis()->SetTitle("local #rho_{00}");
    gmeandiff[i]->GetYaxis()->SetTitle(Form("<cos(2#phi*-2#phi)>_{Data}-<cos(2#phi*-2#phi)>_{%s}",simmode.c_str()));
    gmeandiff[i]->Draw("APE"); 
    TF1 *zero = new TF1("zero","0",-1.0,1.0);
    zero->SetLineColor(kRed);
    zero->Draw("same");
    TF1 *poly1 = new TF1("poly1","[0]+[1]*x",sigmay[0],sigmay[nvar-1]);
    gmeandiff[i]->Fit(poly1,"NMI");
    poly1->SetLineColor(kRed);
    poly1->Draw("same");

    cout << "pt = " << i << ", delta = " << poly1->Eval(datarho[i]) << endl;   

    TF1 *fint = new TF1("fint",finter,sigmay[0],sigmay[nvar-1],2);
    fint->SetParameter(0,poly1->GetParameter(0));
    fint->SetParError(0,poly1->GetParError(0));
    fint->SetParameter(1,poly1->GetParameter(1));
    fint->SetParError(1,poly1->GetParError(1));
    double xint = fint->GetMinimumX();
    TMarker *m = new TMarker(xint, poly1->Eval(xint),24);
    m->SetMarkerColor(kRed);
    m->SetMarkerSize(3);
    m->Draw();
    printf("xint=%g\n",xint);
    intersection[i] = xint;    
  }
  double intersectionrc[4] = {0.292105,0.280381,0.262049,0.232921};
  TLegend *leg1 = new TLegend(0.2,0.6,0.4,0.8);
  for(int i = 0; i < 4; i++)
  {
    cfit->cd(i+17);
    //cfit->cd(i+5)->SetLogx();
    gCosH[i]->SetMarkerStyle(20);
    gCosH[i]->SetMarkerSize(1.2);
    gCosH[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    gCosH[i]->GetXaxis()->SetTitle("local MC #rho_{00}");
    gCosH[i]->GetYaxis()->SetTitle("local RC #rho_{00}");
    gCosH[i]->Draw("APE"); 
    TF1 *line = new TF1("line","x",sigmay[0],sigmay[nvar-1]);
    TF1 *poly1 = new TF1("poly1","[0]+x",sigmay[0],sigmay[nvar-1]);
    gCosH[i]->Fit(poly1,"NMI");
    //line->SetLineColor(kRed);

    //line->Draw("same");
    //if(i == 0) leg1->AddEntry(line,"RC #rho_{00} = MC #rho_{00}","l");
    poly1->SetLineColor(kRed);
    poly1->Draw("same");
    if(i == 0) leg1->AddEntry(poly1,"poly1","l");
    leg1->Draw("same");
  }
  for(int i = 0; i < 4; i++)
  {
    cfit->cd(i+21);
    //cfit->cd(i+5)->SetLogx();
    g1D[i]->SetMarkerStyle(20);
    g1D[i]->SetMarkerSize(1.2);
    g1D[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    g1D[i]->GetXaxis()->SetTitle("local #rho_{00}");
    g1D[i]->GetYaxis()->SetTitle("global #Delta#rho_{00}");
    g1D[i]->Draw("APE"); 
    
    PlotLine(intersectionrc[i],intersectionrc[i],g1D[i]->GetHistogram()->GetMinimum(),g1D[i]->GetHistogram()->GetMaximum(),kRed,2,2);
  }
  for(int i = 0; i < 4; i++)
  {
    cfit->cd(i+25);
    //cfit->cd(i+5)->SetLogx();
    ginput[i]->SetMarkerStyle(20);
    ginput[i]->SetMarkerSize(1.2);
    ginput[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    ginput[i]->GetXaxis()->SetTitle("local input #rho_{00}");
    ginput[i]->GetYaxis()->SetTitle(Form("local %s #rho_{00} - input #rho_{00}",simmode.c_str()));
    ginput[i]->Draw("APE"); 
  }
  cfit->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_FitParams_%s_Order%d_%s%s_helicityrho00comparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectrafile.c_str(),simmode.c_str()));  

  TCanvas *cCosCosH = new TCanvas("cCosCosH","cCosCosH",10,10,1600,800);
  cCosCosH->Divide(4,2);
  
  for(int i = 0; i < 8; i++)
  {
    cCosCosH->cd(i+1);
    cCosCosH->cd(i+1)->SetLeftMargin(0.15);  
    cCosCosH->cd(i+1)->SetBottomMargin(0.15);
    cCosCosH->cd(i+1)->SetTicks(1,1);
    cCosCosH->cd(i+1)->SetGrid(0,0);
  }

  for(int i = 0; i < 4; i++)
  {
    cCosCosH->cd(i+1);
    //cfit->cd(i+5)->SetLogx();
    gCos2D[i]->SetMarkerStyle(20);
    gCos2D[i]->SetMarkerSize(1.2);
    gCos2D[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    gCos2D[i]->GetXaxis()->SetTitle("input local (helicity) #rho_{00}");
    gCos2D[i]->GetYaxis()->SetTitle(Form("%s global #rho_{00} - input global #rho_{00}",simmode.c_str()));
    gCos2D[i]->Draw("APE"); 
  }
  for(int i = 0; i < 4; i++)
  {
    cCosCosH->cd(i+5);
    //cfit->cd(i+5)->SetLogx();
    gCosH2D[i]->SetMarkerStyle(20);
    gCosH2D[i]->SetMarkerSize(1.2);
    gCosH2D[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    gCosH2D[i]->GetXaxis()->SetTitle("input local (helicity) #rho_{00}");
    gCosH2D[i]->GetYaxis()->SetTitle(Form("%s local #rho_{00} - input local #rho_{00}",simmode.c_str()));
    gCosH2D[i]->Draw("APE"); 
  }
  cCosCosH->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_CosCosH_%s_Order%d_%s%s_helicityrho00comparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectrafile.c_str(),simmode.c_str()));  


  outputname = Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_Comparison_%s_Order%d_%s%s_helicityrho00comparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectrafile.c_str(),simmode.c_str());  
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


      h_mCounts_RC[ipt][i]->SetTitle(Form("local #rho_{00}=%s, %1.1f<p_{T}<%1.1f, 20-60 Centrality",sigmay_text[i].c_str(),vmsa::pt_low[energy][ipt+2],vmsa::pt_up[energy][ipt+2]));
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
  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot spts errors
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


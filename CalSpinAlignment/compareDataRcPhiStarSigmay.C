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

void compareDataRcPhiStarSigmay(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 23, std::string simmode = "Mc")//defaultF = 0 is BESII, defaultF = 1 is BESI
{
  const int nvar = 8;
  std::string sigmay_text[nvar] = {"1000","100","10","1","0.5","0.25","0.125","0.0625"};
  float sigmay[nvar] = {1000,100,10,1,0.5,0.25,0.125,0.0625};

  std::string spectrafile = Form("_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240505_fixed2Dfunc_rerho1n10.00");
 
  std::string spectra[nvar];
  for(int i = 0; i < nvar; i++)
  {
    spectra[i] = Form("%s_ysigma%s_rho000.3333",spectrafile.c_str(),sigmay_text[i].c_str());
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
  TH1D *h_mCounts_Ratio[4][nvar];  
  for(int i_pt = 2; i_pt < 6; ++i_pt) // pt loop
  {
    string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    h_mCounts[i_pt-2] = (TH1D*) File_InPut->Get(KEY_counts.c_str());

    float data = h_mCounts[i_pt-2]->Integral(1,10);

    for(int i = 0; i < nvar; i++)
    {
      string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra[i].c_str(),vmsa::mBeamEnergy[energy].c_str());
      File_InPutRC[i] = TFile::Open(inputfileRC.c_str());
      string KEY_counts_RC = Form("h_m%sEffPhiS_Cent_9_Pt_%d",simmode.c_str(),i_pt);
      h_mCounts_RC[i_pt-2][i] = (TH1D*) ((TH1D*) File_InPutRC[i]->Get(KEY_counts_RC.c_str()))->Clone();
      float rc   = h_mCounts_RC[i_pt-2][i]->Integral(1,10);

      h_mCounts_RC[i_pt-2][i]->Scale(data/rc);

      h_mCounts_Ratio[i_pt-2][i] = (TH1D*) h_mCounts[i_pt-2]->Clone();
      h_mCounts_Ratio[i_pt-2][i]->Divide(h_mCounts_RC[i_pt-2][i]);
    }
  }

  TGraphAsymmErrors *gslope[4];
  TGraphAsymmErrors *gintercept[4]; 
  TF1 *fits[4][nvar];
  float pt[4] = {1.5,2.1,2.7,3.6};

  TCanvas *c = new TCanvas("c","c",10,10,1200,1200);
  c->Divide(3,3);

  std::string outputname = Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_Ratio_%s_Order%d_%s%s_sigmaycomparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectrafile.c_str(),simmode.c_str());  
  std::string outputstart = Form("%s[",outputname.c_str()); 
  std::string outputstop = Form("%s]",outputname.c_str()); 
 
  c->Print(outputstart.c_str());

  for(int ipt = 0; ipt < 4; ipt++)
  {
    gslope[ipt] = new TGraphAsymmErrors();
    gintercept[ipt] = new TGraphAsymmErrors();
    for(int i = 0; i < nvar; i++)
    {
      c->cd(i+1);
      c->cd(i+1)->SetLeftMargin(0.15);
      c->cd(i+1)->SetBottomMargin(0.15);
      c->cd(i+1)->SetTicks(1,1);
      c->cd(i+1)->SetGrid(0,0);
      
      h_mCounts_Ratio[ipt][i]->SetTitle(Form("#sigma_{y}=%s, %1.1f<p_{T}<%1.1f, 20-60 Centrality",sigmay_text[i].c_str(),vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mCounts_Ratio[ipt][i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
      h_mCounts_Ratio[ipt][i]->GetYaxis()->SetTitle(Form("(Data)/(%s from Simulation)",simmode.c_str()));
      h_mCounts_Ratio[ipt][i]->Draw("pE");
     
      fits[ipt][i] = new TF1(Form("fit_%d_%d",ipt,i),"[0]+[1]*x",-1.0,1.0);
      h_mCounts_Ratio[ipt][i]->Fit(fits[ipt][i],"NMI");
      fits[ipt][i]->SetLineColor(kRed);
      fits[ipt][i]->Draw("l same");
      
      gslope[ipt]->SetPoint(i,sigmay[i],fits[ipt][i]->GetParameter(1)); 
      gslope[ipt]->SetPointError(i,0.0,0.0,fits[ipt][i]->GetParError(1),fits[ipt][i]->GetParError(1));
      gintercept[ipt]->SetPoint(i,sigmay[i],fits[ipt][i]->GetParameter(0)); 
      gintercept[ipt]->SetPointError(i,0.0,0.0,fits[ipt][i]->GetParError(0),fits[ipt][i]->GetParError(0));
    }
    c->Update();
    c->Print(outputname.c_str());
  }
  c->Print(outputstop.c_str());



  TCanvas *cfit = new TCanvas("cfit","cfit",10,10,1600,800);
  cfit->Divide(4,2);
  
  for(int i = 0; i < 8; i++)
  {
    cfit->cd(i+1);
    cfit->cd(i+1)->SetLeftMargin(0.15);  
    cfit->cd(i+1)->SetBottomMargin(0.15);
    cfit->cd(i+1)->SetTicks(1,1);
    cfit->cd(i+1)->SetGrid(0,0);
  }
 
  float ptedges[5] = {1.2,1.8,2.4,3.0,4.2};

  for(int i = 0; i < 4; i++)
  {
    cfit->cd(i+1);
    cfit->cd(i+1)->SetLogx();
    gslope[i]->SetMarkerStyle(20);
    gslope[i]->SetMarkerSize(1.2);
    gslope[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    gslope[i]->GetXaxis()->SetTitle("#sigma_{y}");
    gslope[i]->GetYaxis()->SetTitle(Form("d(Data/%s)/dp_{T}",simmode.c_str()));
    gslope[i]->SetMarkerSize(1.2);
    gslope[i]->Draw("APE"); 
  }

  for(int i = 0; i < 4; i++)
  {
    cfit->cd(i+5);
    cfit->cd(i+5)->SetLogx();
    gintercept[i]->SetMarkerStyle(20);
    gintercept[i]->SetMarkerSize(1.2);
    gintercept[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f GeV/c",ptedges[i],ptedges[i+1]));
    gintercept[i]->GetXaxis()->SetTitle("#sigma_{y}");
    gintercept[i]->GetYaxis()->SetTitle(Form("(Data/%s) intercept",simmode.c_str()));
    gintercept[i]->Draw("APE"); 
  }
  cfit->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_FitParams_%s_Order%d_%s%s_sigmaycomparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectrafile.c_str(),simmode.c_str()));  


  outputname = Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_Comparison_%s_Order%d_%s%s_sigmaycomparison_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectrafile.c_str(),simmode.c_str());  
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


      h_mCounts_RC[ipt][i]->SetTitle(Form("#sigma_{y}=%s, %1.1f<p_{T}<%1.1f, 20-60 Centrality",sigmay_text[i].c_str(),vmsa::pt_low[energy][ipt+2],vmsa::pt_up[energy][ipt+2]));
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


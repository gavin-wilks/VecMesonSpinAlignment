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

void compareDataRcThetaStar_Spectra(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 25)//defaultF = 0 is BESII, defaultF = 1 is BESI
{

  std::string comp = "GlobalRhoPt";
  const int nopts = 2; 
  //std::string spectra[4] = 
  //{
  //  "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240523_fixed2Dfunc_FlatSpectra_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0_helicityrho0.3333",
  //  "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240507_fixed2Dfunc_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0", 
  //  "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240515_fixed2Dfunc_FinerPhiBins_Spectra_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0",
  //  "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240515_fixed2Dfunc_FinerPhiBins_SpectraPhiAndKaon_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0"
  //};
  //std::string spectra[nopts] = 
  //{
  //  "Pythia_NoWeights_20240731",
  //  "Pythia_PhiMesonWeights_20240731_Retry",
  //  "Pythia_PhiMesonWeights_KaonWeights_20240731"
  //};

  //std::string spectra[nopts] = 
  //{
  //  "Pythia_NoWeights_20240731",
  //  "Pythia_Helicity2DPt_20240813",
  //  "Pythia_Helicity2DY_20240813",
  //};
  std::string spectra[nopts] = 
  {
    "Pythia_NoWeights_20240731",
    "Pythia_Global1DPt_20240814",
  };

  //std::string spectra[nopts] = 
  //{
  //  "Pythia_NoWeights_20240731",
  //  "Pythia_Helicity2DPt_20240813",
  //};
  int mode[nopts] = {0,0};
  //int mode[nopts] = {0,0,0};

  //std::string weight_text[nopts] = {"S0) p_{T} spectra", "S1)S0+#phi-meson weight", "S2)S1+kaon weight"};
  //std::string weight_text[nopts] = {"No Helicity #rho", "Helicity #rho(p_{T})", "Helicity #rho(y)"};
  std::string weight_text[nopts] = {"No Global #rho", "Global #rho(p_{T})"};
  //std::string weight_text[nopts] = {"No Helicity #rho", "Helicity #rho(p_{T})"};
  
  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_%s_PolySys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/CosTheta_%sPhiPtSys_%s_PolySys_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1D *h_mCounts[4];  
  TH1D *h_mCounts_RC[4][4];  
  TH1D *h_mCounts_Ratio[4][4];  

  TFile *File_InPutRC[4];
 
  for(int i_pt = 2; i_pt < 6; ++i_pt) // pt loop
  {
    string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    h_mCounts[i_pt-2] = (TH1D*) File_InPut->Get(KEY_counts.c_str());
    for(int ifile = 0; ifile < nopts; ifile++)
    {
      //string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra[ifile].c_str(),vmsa::mBeamEnergy[energy].c_str());
      string inputfileRC = Form("effaccfiles/%s/%s/%s/Eff_%s_SingleParticle_noToF_Mode%d_EtaMode0_PtBins2_5.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[ifile].c_str(),vmsa::mBeamEnergy[energy].c_str(),mode[ifile]);
      cout << inputfileRC << endl;
      File_InPutRC[ifile] = TFile::Open(inputfileRC.c_str());
      //string KEY_counts_RC = Form("h_mRc5EffPhiS_Cent_9_Pt_%d",i_pt);
      string KEY_counts_RC = Form("h_mRc5EffCos_Cent_9_Pt_%d",i_pt);
      TH1D *h_temp = (TH1D*) ((TH1D*)File_InPutRC[ifile]->Get(KEY_counts_RC.c_str()))->Clone();
      //h_mCounts_RC[ifile][i_pt-2] = (TH1D*) ((TH1D*)File_InPutRC[ifile]->Get(KEY_counts_RC.c_str()))->Clone();
      string coshist = Form("name%d%d",i_pt,ifile);
      h_mCounts_RC[ifile][i_pt-2] = new TH1D(coshist.c_str(),coshist.c_str(),7,0.0,1.0);//(TH1D*) ((TH1D*)File_InPutRC[ifile]->Get(KEY_counts_RC.c_str()))->Clone();
      for(int i = 0; i < 7; i++)
      {
        double val1 = h_temp->GetBinContent(i+1);
        double val2 = h_temp->GetBinContent(14-i);
        double err1 = h_temp->GetBinError(i+1);
        double err2 = h_temp->GetBinError(14-i);

        cout << "i = " << i << ", val1 = " << val1 << ", val2 = " << val2 << endl;

        double val = val1 + val2; 
        double err = TMath::Sqrt(err1*err1+err2*err2);

        h_mCounts_RC[ifile][i_pt-2]->SetBinContent(i+1,val);
        h_mCounts_RC[ifile][i_pt-2]->SetBinError(i+1,err);
      }
     

      float data = h_mCounts[i_pt-2]->Integral(1,7);
      float rc   = h_mCounts_RC[ifile][i_pt-2]->Integral(1,7);
      cout << "i_pt = " << i_pt << ", ifile = " << ifile << ", counts = " << rc << endl; 

      h_mCounts_RC[ifile][i_pt-2]->Scale(data/rc);

      h_mCounts_Ratio[ifile][i_pt-2] = (TH1D*) h_mCounts[i_pt-2]->Clone();
      h_mCounts_Ratio[ifile][i_pt-2]->Divide(h_mCounts_RC[ifile][i_pt-2]);
    }
  }

  int color[4] = {kBlack, kBlue, kOrange+7, kGray+2};

  TGraphAsymmErrors *gslope[4];// = new TGraphAsymmErrors();
  TGraphAsymmErrors *gintercept[4];// = new TGraphAsymmErrors();
  TF1 *fits[4][4];
  float pt[4] = {1.5,2.1,2.7,3.6};

  TCanvas *c = new TCanvas("c","c",10,10,800,800);
  c->Divide(2,2);

  TLegend *leg = new TLegend(0.2,0.2,0.6,0.4);
  for(int i = 0; i < 4; i++)
  {
    c->cd(i+1);
    c->cd(i+1)->SetLeftMargin(0.15);
    c->cd(i+1)->SetBottomMargin(0.15);
    c->cd(i+1)->SetTicks(1,1);
    c->cd(i+1)->SetGrid(0,0);
    
    for(int ifile = 0; ifile < nopts; ifile++)
    {   
      if(i == 0) gslope[ifile] = new TGraphAsymmErrors();
      if(i == 0) gintercept[ifile] = new TGraphAsymmErrors();
      h_mCounts_Ratio[ifile][i]->SetStats(0);
      h_mCounts_Ratio[ifile][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mCounts_Ratio[ifile][i]->GetXaxis()->SetTitle("cos(#theta*)");
      h_mCounts_Ratio[ifile][i]->GetYaxis()->SetTitle("(Data)/(RC from Simulation)");
      h_mCounts_Ratio[ifile][i]->SetMarkerStyle(20);
      h_mCounts_Ratio[ifile][i]->SetMarkerColor(color[ifile]);
      h_mCounts_Ratio[ifile][i]->SetLineColor(color[ifile]);
      if(ifile == 0) h_mCounts_Ratio[ifile][i]->Draw("pE");
      if(ifile != 0) h_mCounts_Ratio[ifile][i]->Draw("pE same");
      if(i == 0) leg->AddEntry(h_mCounts_Ratio[ifile][i],weight_text[ifile].c_str(),"p");
   
      fits[ifile][i] = new TF1(Form("fit_%d_%d",ifile,i),"[0]+[1]*x",-1.0,1.0);
      h_mCounts_Ratio[ifile][i]->Fit(fits[ifile][i],"NMI");
      fits[ifile][i]->SetLineColor(kRed);

      gslope[ifile]->SetPoint(i,pt[i]+ifile*0.05,fits[ifile][i]->GetParameter(1)); 
      gslope[ifile]->SetPointError(i,0.0,0.0,fits[ifile][i]->GetParError(1),fits[ifile][i]->GetParError(1));
      gintercept[ifile]->SetPoint(i,pt[i]+ifile*0.05,fits[ifile][i]->GetParameter(0)); 
      gintercept[ifile]->SetPointError(i,0.0,0.0,fits[ifile][i]->GetParError(0),fits[ifile][i]->GetParError(0));
    }
    leg->Draw("same");
  }
  c->SaveAs(Form("figures/%s/%s/pTstudy/CosTheta_Ratio%sComparison_%s_Order%d_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),comp.c_str(),correction.c_str(),order,etamode.c_str()));  

  TCanvas *cfit = new TCanvas("cfit","cfit",10,10,800,400);
  cfit->Divide(2,1);
  
  cfit->cd(1);
  cfit->cd(1)->SetLeftMargin(0.15);  
  cfit->cd(1)->SetBottomMargin(0.15);
  cfit->cd(1)->SetTicks(1,1);
  cfit->cd(1)->SetGrid(0,0);

  double max = -999999999;
  double min =  999999999;
  for(int ifile = 0; ifile < nopts; ifile++)
  {
    double tempmax = gslope[ifile]->GetHistogram()->GetMaximum();
    double tempmin = gslope[ifile]->GetHistogram()->GetMinimum();

    if(tempmax > max) max = tempmax;
    if(tempmin < min) min = tempmin;
  }

  for(int ifile = 0; ifile < nopts; ifile++)
  {
    gslope[ifile]->SetMarkerStyle(20);
    gslope[ifile]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gslope[ifile]->GetYaxis()->SetTitle("d(Data/RC)/d(cos(theta*))");
    gslope[ifile]->SetMarkerSize(1.2);
    gslope[ifile]->GetHistogram()->SetMaximum(max);
    gslope[ifile]->GetHistogram()->SetMinimum(min);
    gslope[ifile]->SetMarkerColor(color[ifile]);
    gslope[ifile]->SetLineColor(color[ifile]);
    if(ifile == 0) gslope[ifile]->Draw("APE"); 
    else gslope[ifile]->Draw("pE same"); 
  }
  leg->Draw("same");

  cfit->cd(2);
  cfit->cd(2)->SetLeftMargin(0.15);  
  cfit->cd(2)->SetBottomMargin(0.15);
  cfit->cd(2)->SetTicks(1,1);
  cfit->cd(2)->SetGrid(0,0);

  for(int ifile = 0; ifile < nopts; ifile++)
  {
    gintercept[ifile]->SetMarkerStyle(20);
    gintercept[ifile]->SetMarkerSize(1.2);
    gintercept[ifile]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    gintercept[ifile]->GetYaxis()->SetTitle("(Data/RC) intercept");
    gintercept[ifile]->SetMarkerColor(color[ifile]);
    gintercept[ifile]->SetLineColor(color[ifile]);
    if(ifile == 0) gintercept[ifile]->Draw("APE"); 
    else gintercept[ifile]->Draw("pE same"); 
  }
  leg->Draw("same");

  cfit->SaveAs(Form("figures/%s/%s/pTstudy/CosTheta_FitParams%sComparison_%s_Order%d_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),comp.c_str(),correction.c_str(),order,etamode.c_str()));  

  //TLegend *leg = new TLegend(0.4,0.7,0.6,0.9);
  //for(int i = 0; i < 4; i++)
  //{
  //  c->cd(i+1);
  //  c->cd(i+1)->SetLeftMargin(0.15);
  //  c->cd(i+1)->SetBottomMargin(0.15);
  //  c->cd(i+1)->SetTicks(1,1);
  //  c->cd(i+1)->SetGrid(0,0);
  //  
  //  h_mCounts[i]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
  //  h_mCounts[i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
  //  h_mCounts[i]->GetYaxis()->SetTitle("Counts");
  //  h_mCounts[i]->SetMarkerStyle(20);
  //  h_mCounts[i]->SetMarkerColor(kOrange+7);
  //  h_mCounts[i]->SetLineColor(kOrange+7);

  //  int min = h_mCounts[i]->GetMinimum();
  //  int max = h_mCounts[i]->GetMaximum();
  //  if(h_mCounts_RC[i]->GetMinimum() < min) min = h_mCounts_RC[i]->GetMinimum();
  //  if(h_mCounts_RC[i]->GetMaximum() > max) max = h_mCounts_RC[i]->GetMaximum();
  //  h_mCounts[i]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

  //  h_mCounts[i]->Draw("pE");

  //  h_mCounts_RC[i]->SetMarkerStyle(24);
  //  h_mCounts_RC[i]->SetMarkerColor(kBlack);
  //  h_mCounts_RC[i]->SetLineColor(kBlack);
  //  h_mCounts_RC[i]->Draw("pE same");

  //  if(i == 0)
  //  {
  //    leg->AddEntry(h_mCounts[i],"Data","p");
  //    leg->AddEntry(h_mCounts_RC[i],"RC","p");
  //  }
  //  leg->Draw("same");
  //}

  //c->SaveAs(Form("figures/%s/%s/pTstudy/CosTheta_Comparison_%s_Order%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),spectra.c_str()));  


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


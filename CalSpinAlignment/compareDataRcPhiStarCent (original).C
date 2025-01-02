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

void compareDataRcPhiStarCent(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 5)//defaultF = 0 is BESII, defaultF = 1 is BESI
{
 
  std::string spectra = "";
  if(yspectra == 0) spectra = "_NoRapiditySpectra";
  if(yspectra == 1) spectra = "_NoRapiditySpectra_EP";
  if(yspectra == 2) spectra = "_NoRapiditySpectra_EP_noV2";
  if(yspectra == 3) spectra = "_WithRapiditySpectra";
  if(yspectra == 4) spectra = "_WithRapiditySpectra_HalfSigma";
  if(yspectra == 5) spectra = "_NoRapiditySpectra_FixedFirstEP";
  //if(yspectra == 5) spectra = "_NoRapiditySpectra_FixedFirstEP";
  
  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  string inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_%sRhoCentSys_%s_Poly.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_%sRhoCentSys_%s_Poly_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1D *h_mCounts[vmsa::pt_rebin_cent][9];  

  //string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra.c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode1_EtaMode0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra.c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPutRC = TFile::Open(inputfileRC.c_str());
  TH1D *h_mCounts_RC[vmsa::pt_rebin_cent][9];  

  TH1D *h_mCounts_Ratio[vmsa::pt_rebin_cent][9];  
  
  cout << inputfile << endl;
  cout << inputfileRC << endl;
 
  for(int i_pt = vmsa::pt_rebin_first_cent[energy]; i_pt <= vmsa::pt_rebin_last_cent[energy]; ++i_pt) // pt loop
  {
    for(int i_cent = 0; i_cent < 9; i_cent++)
    {
      string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
      cout << KEY_counts << endl;
      h_mCounts[i_pt][i_cent] = (TH1D*) File_InPut->Get(KEY_counts.c_str());
      string KEY_counts_RC = Form("h_mRcEffPhiS_Cent_%d_Pt_%d",i_cent,i_pt);
      cout << KEY_counts_RC << endl;
      h_mCounts_RC[i_pt][i_cent] = (TH1D*) File_InPutRC->Get(KEY_counts_RC.c_str());
      h_mCounts_RC[i_pt][i_cent]->Rebin(2);
  
      float data = h_mCounts[i_pt][i_cent]->Integral(1,10);
      float rc   = h_mCounts_RC[i_pt][i_cent]->Integral(1,10);
      cout << data << "    " << rc << endl;  

      h_mCounts_RC[i_pt][i_cent]->Scale(data/rc);
  
      h_mCounts_Ratio[i_pt][i_cent] = (TH1D*) h_mCounts[i_pt][i_cent]->Clone();
      h_mCounts_Ratio[i_pt][i_cent]->Divide(h_mCounts_RC[i_pt][i_cent]);


    }
  }
  cout << "Loaded histograms" << endl;

  TCanvas *c = new TCanvas("c","c",600,10,1200,1200);
  c->Divide(3,3);

  for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
  {
    for(int i = 0; i < 9; i++)
    {
      c->cd(i+1);
      c->cd(i+1)->SetLeftMargin(0.15);
      c->cd(i+1)->SetBottomMargin(0.15);
      c->cd(i+1)->SetTicks(1,1);
      c->cd(i+1)->SetGrid(0,0);
      
      h_mCounts_Ratio[ipt][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %2.0f-%2.0f Centrality",vmsa::pt_low_cent[energy][ipt],vmsa::pt_up_cent[energy][ipt],vmsa::centVal[i+1],vmsa::centVal[i]));
      h_mCounts_Ratio[ipt][i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
      h_mCounts_Ratio[ipt][i]->GetYaxis()->SetTitle("(Data)/(RC from Simulation)");
      h_mCounts_Ratio[ipt][i]->Draw("pE");
    }
    c->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_CentRatio_%s_Order%d_%s_pt%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),ipt));  
  }

  TLegend *leg = new TLegend(0.4,0.7,0.6,0.9);
  for(int ipt = 0; ipt < vmsa::pt_rebin_cent; ipt++)
  {
    for(int i = 0; i < 9; i++)
    {
      c->cd(i+1);
      c->cd(i+1)->SetLeftMargin(0.15);
      c->cd(i+1)->SetBottomMargin(0.15);
      c->cd(i+1)->SetTicks(1,1);
      c->cd(i+1)->SetGrid(0,0);
      
      h_mCounts[ipt][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %2.0f-%2.0f Centrality",vmsa::pt_low_cent[energy][ipt],vmsa::pt_up_cent[energy][ipt],vmsa::centVal[i+1],vmsa::centVal[i]));
      h_mCounts[ipt][i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
      h_mCounts[ipt][i]->GetYaxis()->SetTitle("Counts");
      h_mCounts[ipt][i]->SetMarkerStyle(20);
      h_mCounts[ipt][i]->SetMarkerColor(kOrange+7);
      h_mCounts[ipt][i]->SetLineColor(kOrange+7);

      int min = h_mCounts[ipt][i]->GetMinimum();
      int max = h_mCounts[ipt][i]->GetMaximum();
      if(h_mCounts_RC[ipt][i]->GetMinimum() < min) min = h_mCounts_RC[ipt][i]->GetMinimum();
      if(h_mCounts_RC[ipt][i]->GetMaximum() > max) max = h_mCounts_RC[ipt][i]->GetMaximum();
      h_mCounts[ipt][i]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

      h_mCounts[ipt][i]->Draw("pE");

      h_mCounts_RC[ipt][i]->SetMarkerStyle(24);
      h_mCounts_RC[ipt][i]->SetMarkerColor(kBlack);
      h_mCounts_RC[ipt][i]->SetLineColor(kBlack);
      h_mCounts_RC[ipt][i]->Draw("pE same");

      if(i == 0)
      {
        leg->AddEntry(h_mCounts[ipt][i],"Data","p");
        leg->AddEntry(h_mCounts_RC[ipt][i],"MC","p");
      }
      leg->Draw("same");
    }

    c->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_CentComparison_%s_Order%d_%s_pt%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),ipt));  
  }

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


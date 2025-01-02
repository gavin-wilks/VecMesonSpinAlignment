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

void compareDataRcPhiStarY(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 5)//defaultF = 0 is BESII, defaultF = 1 is BESI
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

  string inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_%sRhoEtaSys_%s_Poly.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/Cos2PhiStarPhi_%sRhoEtaSys_%s_Poly_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1D *h_mCounts[vmsa::pt_rebin_y][vmsa::cent_rebin_total][vmsa::eta_total];  

  //string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra.c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode2_EtaMode0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra.c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPutRC = TFile::Open(inputfileRC.c_str());
  TH1D *h_mCounts_RC[vmsa::pt_rebin_y][vmsa::cent_rebin_total][vmsa::eta_total];  

  TH1D *h_mCounts_Ratio[vmsa::pt_rebin_y][vmsa::cent_rebin_total][vmsa::eta_total];  
  
  cout << inputfile << endl;
  cout << inputfileRC << endl;
 
  int centidx[3] = {0,4,7};

  for(int i_pt = 0; i_pt < vmsa::pt_rebin_y; ++i_pt) // pt loop
  {
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total; i_cent++)
    {
      for(int i_y = 0; i_y < vmsa::eta_total; i_y++)
      {
        string KEY_counts = Form("eta_%d_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_y,i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
        h_mCounts[i_pt][i_cent][i_y] = (TH1D*) File_InPut->Get(KEY_counts.c_str());
        string KEY_counts_RC = Form("h_mRcEffPhiS_Cent_%d_Pt_%d_Y_%d",centidx[i_cent],i_pt,i_y);
        h_mCounts_RC[i_pt][i_cent][i_y] = (TH1D*) File_InPutRC->Get(KEY_counts_RC.c_str());
        h_mCounts_RC[i_pt][i_cent][i_y]->Rebin(2);
  
        float data = h_mCounts[i_pt][i_cent][i_y]->Integral(1,10);
        float rc   = h_mCounts_RC[i_pt][i_cent][i_y]->Integral(1,10);
  
        h_mCounts_RC[i_pt][i_cent][i_y]->Scale(data/rc);
  
        h_mCounts_Ratio[i_pt][i_cent][i_y] = (TH1D*) h_mCounts[i_pt][i_cent][i_y]->Clone();
        h_mCounts_Ratio[i_pt][i_cent][i_y]->Divide(h_mCounts_RC[i_pt][i_cent][i_y]);
     }
    }
  }

  cout << "Loaded Histograms" << endl;

  TCanvas *c = new TCanvas("c","c",600,10,1500,900);
  c->Divide(5,3);

  for(int ipt = 0; ipt < vmsa::pt_rebin_y; ipt++)
  {
    for(int icent = 0; icent < vmsa::cent_rebin_total; icent++)
    {
      for(int i = 0; i < vmsa::eta_total; i++)
      {
        c->cd(i+1);
        c->cd(i+1)->SetLeftMargin(0.15);
        c->cd(i+1)->SetBottomMargin(0.15);
        c->cd(i+1)->SetTicks(1,1);
        c->cd(i+1)->SetGrid(0,0);
        
        h_mCounts_Ratio[ipt][icent][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %2.0f-%2.0f Centrality",vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt],vmsa::centVal[vmsa::cent_rebin[icent+1]],vmsa::centVal[vmsa::cent_rebin[icent]]));
        h_mCounts_Ratio[ipt][icent][i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
        h_mCounts_Ratio[ipt][icent][i]->GetYaxis()->SetTitle("(Data)/(RC from Simulation)");
        h_mCounts_Ratio[ipt][icent][i]->Draw("pE");
      }
      c->SaveAs(Form("figures/%s/%s/rapiditystudy/Cos2PhiStarPhi_YRatio_%s_Order%d_%s_pt%d_cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),ipt,icent));  
    }
  }

  for(int ipt = 0; ipt < vmsa::pt_rebin_y; ipt++)
  {
    for(int icent = 0; icent < vmsa::cent_rebin_total; icent++)
    {
      TLegend *leg = new TLegend(0.4,0.7,0.6,0.9);
      for(int i = 0; i < vmsa::eta_total; i++)
      {
        c->cd(i+1);
        c->cd(i+1)->SetLeftMargin(0.15);
        c->cd(i+1)->SetBottomMargin(0.15);
        c->cd(i+1)->SetTicks(1,1);
        c->cd(i+1)->SetGrid(0,0);
        
        h_mCounts[ipt][icent][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f, %2.0f-%2.0f Centrality",vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt],vmsa::centVal[vmsa::cent_rebin[icent+1]],vmsa::centVal[vmsa::cent_rebin[icent]]));
        h_mCounts[ipt][icent][i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
        h_mCounts[ipt][icent][i]->GetYaxis()->SetTitle("Counts");
        h_mCounts[ipt][icent][i]->SetMarkerStyle(20);
        h_mCounts[ipt][icent][i]->SetMarkerColor(kOrange+7);
        h_mCounts[ipt][icent][i]->SetLineColor(kOrange+7);

        int min = h_mCounts[ipt][icent][i]->GetMinimum();
        int max = h_mCounts[ipt][icent][i]->GetMaximum();
        if(h_mCounts_RC[ipt][icent][i]->GetMinimum() < min) min = h_mCounts_RC[ipt][icent][i]->GetMinimum();
        if(h_mCounts_RC[ipt][icent][i]->GetMaximum() > max) max = h_mCounts_RC[ipt][icent][i]->GetMaximum();
        h_mCounts[ipt][icent][i]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

        h_mCounts[ipt][icent][i]->Draw("pE");

        h_mCounts_RC[ipt][icent][i]->SetMarkerStyle(24);
        h_mCounts_RC[ipt][icent][i]->SetMarkerColor(kBlack);
        h_mCounts_RC[ipt][icent][i]->SetLineColor(kBlack);
        h_mCounts_RC[ipt][icent][i]->Draw("pE same");

        if(i == 0)
        {
          leg->AddEntry(h_mCounts[ipt][icent][i],"Data","p");
          leg->AddEntry(h_mCounts_RC[ipt][icent][i],"MC","p");
        }
        leg->Draw("same");
      }

      c->SaveAs(Form("figures/%s/%s/rapiditystudy/Cos2PhiStarPhi_YComparison_%s_Order%d_%s_pt%d_cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),ipt,icent));  
    }
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


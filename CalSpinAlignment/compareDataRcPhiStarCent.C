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

void compareRcPhiStarPt(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 16)//defaultF = 0 is BESII, defaultF = 1 is BESI
{
 
  std::string spectra[4] = {"_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc",
                            "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc_pT0p2",
                            "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc_pT0p2_TPC",
                            "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc_pT0p2_TPC_TOF",};
  
  std::string legtitle[4] = {"RC0: |#eta| < 0.9", 
                             "RC1: RC0 + p_{T}<0.2", 
                             "RC2: RC1 + TPC Eff",
                             "RC3: RC2 + TOF Eff"};
  
  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  int marker[4] = {24,20,21,22};
  int color[4] = {1,2,1,2};

  TFile *File_InPutRC[4];
  TH1D *h_mEff[4][4];  

  for(int i_s = 0; i_s < 4; i_s++)
  { 
    string inputfileRC = Form("effaccfiles/%s/%s/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order%d%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,spectra[i_s].c_str(),vmsa::mBeamEnergy[energy].c_str());
    File_InPutRC[i_s] = TFile::Open(inputfileRC.c_str());
    for(int i_pt = 2; i_pt < 6; ++i_pt) // pt loop
    {
      string KEY_counts_RC = Form("h_mEffPhiS_Cent_9_Pt_%d",i_pt);
      h_mEff[i_s][i_pt-2] = (TH1D*) File_InPutRC[i_s]->Get(KEY_counts_RC.c_str());
    }
  }

  TGraphAsymmErrors *g_a2[4];
  float a2[16] = {0.0};

  TCanvas *c = new TCanvas("c","c",600,10,800,800);
  c->Divide(2,2);

  std::string outputname = Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_EfficiencyFits_%s_Order%d_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str());  
  std::string output_start = Form("%s[",outputname.c_str());
  std::string output_stop = Form("%s]",outputname.c_str());

  c->Print(output_start.c_str());

  for(int i_s = 0; i_s < 4; i_s++)
  {
    g_a2[i_s] = new TGraphAsymmErrors();
    TF1* f_eff[4];    
    for(int i = 0; i < 4; i++)
    {
      f_eff[i] = new TF1(Form("f_eff_%d",i),"[0]*(1.+2.*[1]*x)",-1.0,1.0);
      c->cd(i+1);
      c->cd(i+1)->SetLeftMargin(0.15);
      c->cd(i+1)->SetBottomMargin(0.15);
      c->cd(i+1)->SetTicks(1,1);
      c->cd(i+1)->SetGrid(0,0);
      
      h_mEff[i_s][i]->SetTitle(Form("%1.1f<p_{T}<%1.1f, 20-60 Centrality",vmsa::pt_low[energy][i+2],vmsa::pt_up[energy][i+2]));
      h_mEff[i_s][i]->GetXaxis()->SetTitle("cos(2#phi*-2#phi)");
      h_mEff[i_s][i]->GetYaxis()->SetTitle("Efficiency");
      h_mEff[i_s][i]->Draw("pE");
    
      h_mEff[i_s][i]->Fit(f_eff[i],"QNMRI");
      f_eff[i]->SetLineColor(kRed);
      f_eff[i]->Draw("l same");

      float ptmean = (vmsa::pt_low[energy][i+2]+vmsa::pt_up[energy][i+2])/2.0;
      g_a2[i_s]->SetPoint(i,ptmean+i_s*0.05,f_eff[i]->GetParameter(1));
      g_a2[i_s]->SetPointError(i,0.0,0.0,f_eff[i]->GetParError(1),f_eff[i]->GetParError(1));

      a2[i_s*4+i] = f_eff[i]->GetParameter(1);
      cout << "method = " << i_s << ", pT = " << i+1 << ", a2 = " << f_eff[i]->GetParameter(1) << endl;
    }
    c->Update();
    c->Print(outputname.c_str());
  }
  
  c->Print(output_stop.c_str());


  float min = TMath::MinElement(16,a2);
  float max = TMath::MaxElement(16,a2);

  TCanvas *c_a2 = new TCanvas("c_a2","c_a2",600,10,400,400);
  c_a2->cd();
  c_a2->cd()->SetLeftMargin(0.15);
  c_a2->cd()->SetBottomMargin(0.15);
  c_a2->cd()->SetTicks(1,1);
  c_a2->cd()->SetGrid(0,0);

  TLegend *leg = new TLegend(0.5,0.2,0.8,0.5);
  for(int i = 0; i < 4; i++)
  {
    g_a2[i]->SetTitle(Form("%s 20-60 Centrality",vmsa::mBeamEnergy[energy].c_str()));
    g_a2[i]->GetHistogram()->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    g_a2[i]->GetHistogram()->GetYaxis()->SetTitle("a_{2}");
    if(max > 0) g_a2[i]->GetHistogram()->GetYaxis()->SetRangeUser(min*1.2, max*1.2);
    else        g_a2[i]->GetHistogram()->GetYaxis()->SetRangeUser(min*1.2, max*0.8);
    if(max > 0) g_a2[i]->SetMaximum(max*1.2);
    else        g_a2[i]->SetMaximum(max*0.8);
    g_a2[i]->SetMinimum(min*1.2);
    g_a2[i]->SetMarkerStyle(marker[i]);
    g_a2[i]->SetMarkerColor(color[i]);
    g_a2[i]->SetLineColor(color[i]);

    if(i == 0) g_a2[i]->Draw("APE");
    else       g_a2[i]->Draw("PE same");

    leg->AddEntry(g_a2[i],legtitle[i].c_str(),"p");
  }
  leg->Draw("same");

  c_a2->SaveAs(Form("figures/%s/%s/pTstudy/Cos2PhiStarPhi_a2Comparison_%s_Order%d_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str()));  


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


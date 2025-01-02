#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "../Utility/draw.h"
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void ptEtaDependence(Int_t energy = 4, Int_t pid = 0, std::string setting = "AccRes")
{
  string inputfileF1 = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_0_eta0p6.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),setting.c_str());

  TFile *File_InPutF1 = TFile::Open(inputfileF1.c_str());

  TH1F *h_frame = (TH1F*)File_InPutF1->Get("h_frame");

  TGraphAsymmErrors *g_before = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_after = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_before1 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_after1 = new TGraphAsymmErrors();

  float eta[4]       = {1.0,0.8,0.6,0.4};
 
  // 2nd Order
  float before[4]    = {0.3423,0.3454,0.3546,0.3750};
  float beforeErr[4] = {0.0012,0.0012,0.0015,0.0020};
  float after[4]     = {0.3520,0.3461,0.3417,0.3366}; // v2 on and pt on, Eff then Acc+Res correction
  float afterErr[4]  = {0.0024,0.0025,0.0029,0.0039}; // v2 on and pt on, Eff then Acc+Res correction
 
  // 1st Order
  float before1[4]    = {0.3379,0.3405,0.3504,0.3716};
  float beforeErr1[4] = {0.0012,0.0013,0.0015,0.0019};
  float after1[4]     = {0.3458,0.3387,0.3341,0.3312}; // v2 on and pt on, Eff then Acc+Res correction
  float afterErr1[4]  = {0.0029,0.0031,0.0036,0.0046}; // v2 on and pt on, Eff then Acc+Res correction

  for( int i = 0; i < 4; i++ )
  {
    g_after->SetPoint(i, eta[i], after[i]);
    g_after->SetPointError(i, 0.0, 0.0, afterErr[i], afterErr[i]);
    g_before->SetPoint(i, eta[i], before[i]);
    g_before->SetPointError(i, 0.0, 0.0, beforeErr[i], beforeErr[i]);
    g_after1->SetPoint(i, eta[i], after1[i]);
    g_after1->SetPointError(i, 0.0, 0.0, afterErr1[i], afterErr1[i]);
    g_before1->SetPoint(i, eta[i], before1[i]);
    g_before1->SetPointError(i, 0.0, 0.0, beforeErr1[i], beforeErr1[i]);
  }


  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,800,800);
  c_rho_SysError->cd();
  c_rho_SysError->cd()->SetLeftMargin(0.15);
  c_rho_SysError->cd()->SetBottomMargin(0.15);
  c_rho_SysError->cd()->SetTicks(1,1);
  c_rho_SysError->cd()->SetGrid(0,0);
  h_frame->GetYaxis()->SetRangeUser(0.32,0.40);
  h_frame->GetXaxis()->SetRangeUser(0.2,1.2);
  h_frame->GetXaxis()->SetTitle("|#eta| cut" );
  h_frame->DrawCopy("pE");

  PlotLine(0.2,1.2,1.0/3.0,1.0/3.0,1,2,2);

  gStyle->SetErrorX(0);

  g_before->SetMarkerStyle(20);
  g_before->SetMarkerColor(kOrange+7);
  g_before->SetLineColor(kOrange+7);
  g_before->SetMarkerSize(1.3);
  g_before->Draw("pE same");
  
  g_after->SetMarkerStyle(20);
  g_after->SetMarkerColor(kBlue);
  g_after->SetLineColor(kBlue);
  g_after->SetMarkerSize(1.3);
  g_after->Draw("pE same");

  TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);
  leg->AddEntry(g_before,"Before Corrections","p");
  leg->AddEntry(g_after,"After Corrections","p");
  leg->Draw("same");
  
  c_rho_SysError->SaveAs("./figures/Phi/19GeV/ptDependence_EtaCutPlot.pdf");


  TCanvas *c_b = new TCanvas("c_b","c_b",600,10,800,800);
  c_b->cd();
  c_b->cd()->SetLeftMargin(0.15);
  c_b->cd()->SetBottomMargin(0.15);
  c_b->cd()->SetTicks(1,1);
  c_b->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");

  PlotLine(0.2,1.2,1.0/3.0,1.0/3.0,1,2,2);

  gStyle->SetErrorX(0);

  g_before1->SetMarkerStyle(20);
  g_before1->SetMarkerColor(kOrange+7);
  g_before1->SetLineColor(kOrange+7);
  g_before1->SetMarkerSize(1.3);
  g_before1->Draw("pE same");
  
  g_before->SetMarkerStyle(20);
  g_before->SetMarkerColor(kBlue);
  g_before->SetLineColor(kBlue);
  g_before->SetMarkerSize(1.3);
  g_before->Draw("pE same");

  TLegend *leg2 = new TLegend(0.5,0.6,0.8,0.8);
  leg2->AddEntry(g_before1,"1^{st} Order Before Corrections","p");
  leg2->AddEntry(g_before, "2^{nd} Order Before Corrections","p");
  leg2->Draw("same");
  
  c_b->SaveAs("./figures/Phi/19GeV/beforeComparison_ptDependence_EtaCutPlot.pdf");



  TCanvas *c_a = new TCanvas("c_a","c_a",600,10,800,800);
  c_a->cd();
  c_a->cd()->SetLeftMargin(0.15);
  c_a->cd()->SetBottomMargin(0.15);
  c_a->cd()->SetTicks(1,1);
  c_a->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");

  PlotLine(0.2,1.2,1.0/3.0,1.0/3.0,1,2,2);

  gStyle->SetErrorX(0);

  g_after1->SetMarkerStyle(20);
  g_after1->SetMarkerColor(kOrange+7);
  g_after1->SetLineColor(kOrange+7);
  g_after1->SetMarkerSize(1.3);
  g_after1->Draw("pE same");
  
  g_after->SetMarkerStyle(20);
  g_after->SetMarkerColor(kBlue);
  g_after->SetLineColor(kBlue);
  g_after->SetMarkerSize(1.3);
  g_after->Draw("pE same");

  TLegend *leg3 = new TLegend(0.5,0.6,0.8,0.8);
  leg3->AddEntry(g_after1,"1^{st} Order After Corrections","p");
  leg3->AddEntry(g_after, "2^{nd} Order After Corrections","p");
  leg3->Draw("same");
  
  c_rho_SysError->SaveAs("./figures/Phi/19GeV/afterComparison_ptDependence_EtaCutPlot.pdf");

}


void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE)
{
  const int nEnergy = g_rho->GetN();
  TBox *bSys[nEnergy];
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);
    
    //bSys[i_energy] = new TBox(energy/1.08,rho-err,energy*1.08,rho+err);
    bSys[i_energy] = new TBox(vmsa::pt_low[beamE][i_energy],rho-err,vmsa::pt_up[beamE][i_energy],rho+err);
    bSys[i_energy]->SetFillColor(0);
    bSys[i_energy]->SetFillStyle(0);
    bSys[i_energy]->SetLineStyle(1);
    bSys[i_energy]->SetLineWidth(1);
    bSys[i_energy]->SetLineColor(plot_color);
    bSys[i_energy]->Draw("l Same");
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

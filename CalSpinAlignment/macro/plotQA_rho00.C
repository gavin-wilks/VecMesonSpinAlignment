#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "../../Utility/draw.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

void plotQA_rho00()
{
  TString InPut_new = "/global/homes/x/xusun/AuAu200GeV/SpinAlignment/Phi/rho00/RawRhoPtSys.root";
  TFile *File_new = TFile::Open(InPut_new.Data());
  TH1F *h_frame = (TH1F*)File_new->Get("h_frame");
  TGraphAsymmErrors *g_new = (TGraphAsymmErrors*)File_new->Get("rhoRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Count")->Clone();

  TString InPut_old = "/global/homes/x/xusun/AuAu200GeV/SpinAlignment/Phi/rho00/backup/RawRhoPtSys.root";
  TFile *File_old = TFile::Open(InPut_old.Data());
  TGraphAsymmErrors *g_old = (TGraphAsymmErrors*)File_old->Get("rhoRaw_Centrality_9_2nd_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Count")->Clone();

  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  h_frame->GetYaxis()->SetRangeUser(0.3,0.4);
  h_frame->Draw("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  g_new->SetMarkerStyle(24);
  g_new->SetMarkerColor(kGray+2);
  g_new->SetMarkerSize(1.2);
  g_new->Draw("pE same");

  g_old->SetMarkerStyle(24);
  g_old->SetMarkerColor(2);
  g_old->SetMarkerSize(1.2);
  g_old->Draw("pE same");

  TLegend *leg = new TLegend(0.2,0.7,0.5,0.85);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(g_new,"ToF always","p");
  leg->AddEntry(g_old,"old cuts","p");
  leg->Draw("same");

  c_rho->SaveAs("../../figures/c_rho_ToF.eps");
}

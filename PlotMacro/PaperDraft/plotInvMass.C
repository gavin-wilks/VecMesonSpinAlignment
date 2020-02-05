#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/functions.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"

using namespace std;

void plotInvMass()
{
  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/BESII/invmass.root");
  TH1D *h_mMass_SE = (TH1D*)File_Input->Get("imass_all_clone")->Clone("h_mMass_SE");
  TH1D *h_mMass_ME = (TH1D*)File_Input->Get("imass_all_bg")->Clone("h_mMass_ME");
  TH1D *h_mMass_SM = (TH1D*)File_Input->Get("imass_all")->Clone("h_mMass_SM");

  TCanvas *c_peak = new TCanvas("c_peak","c_peak",10,10,800,800);
  c_peak->cd();
  c_peak->cd()->SetLeftMargin(0.15);
  c_peak->cd()->SetBottomMargin(0.15);
  c_peak->cd()->SetTicks(1,1);
  c_peak->cd()->SetGrid(0,0);
  h_mMass_SE->SetTitle("");
  h_mMass_SE->SetStats(0);
  h_mMass_SE->GetXaxis()->SetTitle("M(K^{+},K^{-})(GeV/c^{2})");
  h_mMass_SE->GetXaxis()->CenterTitle();
  h_mMass_SE->GetXaxis()->SetTitleSize(0.06);
  h_mMass_SE->GetXaxis()->SetTitleOffset(0.9);
  h_mMass_SE->GetXaxis()->SetLabelSize(0.04);
  h_mMass_SE->SetNdivisions(505,"X");

  h_mMass_SE->GetYaxis()->SetTitle("Yields");
  h_mMass_SE->GetYaxis()->CenterTitle();
  h_mMass_SE->GetYaxis()->SetTitleOffset(1.14);
  h_mMass_SE->GetYaxis()->SetTitleSize(0.06);
  h_mMass_SE->GetYaxis()->SetLabelSize(0.04);
  h_mMass_SE->GetYaxis()->SetRangeUser(-0.1*h_mMass_SE->GetMaximum(),1.15*h_mMass_SE->GetMaximum());
  h_mMass_SE->SetNdivisions(505,"Y");
  h_mMass_SE->SetMarkerStyle(24);
  h_mMass_SE->SetMarkerColor(1);
  h_mMass_SE->SetMarkerSize(1.4);
  h_mMass_SE->DrawCopy("PE");
  // PlotLine(vmsa::InvMass_low[0],vmsa::InvMass_high[0],0.0,0.0,1,5,2);

  h_mMass_ME->SetMarkerStyle(1);
  // h_mMass_ME->SetLineColor(kAzure+2);
  // h_mMass_ME->SetFillColor(kAzure+2);
  h_mMass_ME->SetLineColor(kGray+2);
  h_mMass_ME->SetFillColor(kGray);
  h_mMass_ME->SetFillStyle(3003);
  h_mMass_ME->DrawCopy("h same");

  h_mMass_SM->SetMarkerStyle(20);
  h_mMass_SM->SetMarkerColor(2);
  h_mMass_SM->SetMarkerSize(1.2);
  h_mMass_SM->DrawCopy("pE same");

  TF1 *f_sig = new TF1("f_sig",PolyBreitWigner,vmsa::BW_Start[0],vmsa::BW_Stop[0],5);
  f_sig->SetParameter(0,vmsa::InvMass[0]);
  f_sig->SetParLimits(0,vmsa::InvMass[0]-1.5*vmsa::Width[0],vmsa::InvMass[0]+1.5*vmsa::Width[0]);
  f_sig->SetParameter(1,vmsa::Width[0]);
  f_sig->SetParLimits(1,0.004,0.070);
  f_sig->SetParameter(2,1.0);
  f_sig->SetParameter(3,-1.0);
  f_sig->SetParameter(4,1.0);
  f_sig->SetParameter(2,h_mMass_SM->GetMaximum()/f_sig->GetMaximum());
  f_sig->SetRange(vmsa::BW_Start[0],vmsa::BW_Stop[0]);
  h_mMass_SM->Fit(f_sig,"MNR");
  f_sig->SetLineColor(4);
  f_sig->SetLineStyle(1);
  f_sig->SetLineWidth(2);
  f_sig->DrawCopy("l same");

  TF1 *f_bg = new TF1("f_bg",Poly,vmsa::BW_Start[0],vmsa::BW_Stop[0],2);
  for(int i_par = 0; i_par < 2;++i_par)
  {
    f_bg->SetParameter(i_par,f_sig->GetParameter(i_par+3));
  }
  f_bg->SetLineColor(4);
  f_bg->SetLineStyle(2);
  f_bg->SetLineWidth(4);
  f_bg->DrawCopy("l same");

  TLegend *leg = new TLegend(0.5,0.4,0.8,0.55);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(h_mMass_SE,"Same Event (Sig+Bg)","p");
  leg->AddEntry(h_mMass_ME,"Mixed Event (Bg)","f");
  leg->AddEntry(h_mMass_SM,"Signal","p");
  leg->Draw("same");

  string leg_energy = Form("Au+Au %s & 20-60%%",vmsa::mBeamEnergy[6].c_str());
  plotTopLegend((char*)leg_energy.c_str(),0.45,0.84,0.04,1,0.0,42,1);
  plotTopLegend((char*)"1.2 < p_{T} < 1.8 GeV/c",0.50,0.79,0.04,1,0.0,42,1);

  c_peak->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/c_InvMass.eps");
  c_peak->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_InvMass.png");
}

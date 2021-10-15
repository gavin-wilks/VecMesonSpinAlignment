#include <string>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGaxis.h>

#include "../../../Utility/draw.h"
#include "../../../Utility/StSpinAlignmentCons.h"
#include "../../../Utility/functions.h"

using namespace std;

void plotFig4_SignalExtraction()
{
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;

  const int style_Kstr = 20;
  const int color_Kstr = kAzure-9;

  const float size_marker = 1.4;

  TFile *File_InputPhi_BeforeSub = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/fig_2ndEP_27GeV_Run18_BeforeSub.root");
  TH1D *h_mMassPhi = (TH1D*)File_InputPhi_BeforeSub->Get("imass_all_pt2")->Clone("h_mMassPhi");

  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/fig_2ndEP_27GeV_Run18.root");
  TH1D *h_mYieldsPhi = (TH1D*)File_InputPhi->Get("fit_cos_pt2")->Clone("h_mYieldsPhi");
  float cos[7] = {1.0/14.0,3.0/14.0,5.0/14.0,7.0/14.0,9.0/14.0,11.0/14.0,13.0/14.0};
  // float yields[7] = {598.598, 596.721, 606.366, 606.203, 615.265, 626.465,  644.074};
  // float errors[7] = {4.06987, 4.06908, 4.12435, 4.1327, 4.1449, 4.2176, 4.26549};

  const long numOfEvent = 176995424;
  const double binWidth = h_mMassPhi->GetBinWidth(1);
  // cout << "phi-meson: numOfEvent = " << numOfEvent << ", binWidth = " << binWidth << endl;

  TGraphAsymmErrors *g_yieldsPhi = new TGraphAsymmErrors();
  for(int i_cos = 0; i_cos < 7; ++i_cos)
  {
    g_yieldsPhi->SetPoint(i_cos,cos[i_cos],h_mYieldsPhi->GetBinContent(i_cos+1)/(numOfEvent*binWidth));
    g_yieldsPhi->SetPointError(i_cos,1.0/14.0,1.0/14.0,h_mYieldsPhi->GetBinError(i_cos+1)/(numOfEvent*binWidth),h_mYieldsPhi->GetBinError(i_cos+1)/(numOfEvent*binWidth));
  }

  // TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/Kstar_Fig_signal_54_rap0p1.root");
  TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/Kstar_Fig_signal_54_rap0p1_Oct2021.root");
  TH1D *h_mMassKstar = (TH1D*)File_InputKstar->Get("ThetaIntegratedSignal_54GeV_PtBin2")->Clone("h_mMassKstar");
  TF1 *f_SigKstar = (TF1*)File_InputKstar->Get("BWFunc_54GeV_PtBin2")->Clone("f_SigKstar");
  TF1 *f_BkgKstar = (TF1*)File_InputKstar->Get("Res_Func54GeV_PtBin2")->Clone("f_BkgKstar");
  // TH1F *h_mYieldsKstar = (TH1F*)File_InputKstar->Get("Yield_vs_Theta_54GeV_PtBin2")->Clone("h_mYieldsKstar");
  // TF1 *f_rhoKstar = (TF1*)File_InputKstar->Get("Rho00_Fit_Func_54GeV_PtBin2")->Clone("f_rhoKstar");
  TH1F *h_mYieldsKstar = (TH1F*)File_InputKstar->Get("hist_Raw_Yield_theta_54GeV_Pt_2.00_2.50")->Clone("h_mYieldsKstar");
  TF1 *f_rhoKstar = (TF1*)File_InputKstar->Get("fit_Raw_Yield_theta_54GeV_Pt_2.00_2.50")->Clone("f_rhoKstar");

  TGraphAsymmErrors *g_yieldsKstar = new TGraphAsymmErrors();
  for(int i_cos = 0; i_cos < 5; ++i_cos)
  {
    g_yieldsKstar->SetPoint(i_cos,0.1+0.2*i_cos,h_mYieldsKstar->GetBinContent(i_cos+1));
    g_yieldsKstar->SetPointError(i_cos,0.1,0.1,h_mYieldsKstar->GetBinError(i_cos+1),h_mYieldsKstar->GetBinError(i_cos+1));
  }

  TGaxis::SetMaxDigits(3);

  TCanvas *c_Signal = new TCanvas("c_Signal","c_Signal",10,10,800,800);
  c_Signal->Divide(2,2);
  for(int i_pad = 0; i_pad < 4; ++i_pad)
  {
    c_Signal->cd(i_pad+1);
    c_Signal->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Signal->cd(i_pad+1)->SetTopMargin(0.05);
    c_Signal->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Signal->cd(i_pad+1)->SetTicks(1,1);
    c_Signal->cd(i_pad+1)->SetGrid(0,0);
  }

  c_Signal->cd(1);
  {
    h_mMassPhi->Scale(1.0/numOfEvent);
    h_mMassPhi->SetTitle("");
    h_mMassPhi->SetStats(0);
    h_mMassPhi->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
    h_mMassPhi->GetXaxis()->CenterTitle();
    h_mMassPhi->GetXaxis()->SetTitleSize(0.06);
    h_mMassPhi->GetXaxis()->SetTitleOffset(0.9);
    h_mMassPhi->GetXaxis()->SetTitleFont(42);
    h_mMassPhi->GetXaxis()->SetLabelSize(0.04);
    h_mMassPhi->SetNdivisions(505,"X");

    h_mMassPhi->GetYaxis()->SetTitle("Raw Yields");
    h_mMassPhi->GetYaxis()->CenterTitle();
    h_mMassPhi->GetYaxis()->SetTitleOffset(1.14);
    h_mMassPhi->GetYaxis()->SetTitleSize(0.06);
    h_mMassPhi->GetYaxis()->SetTitleFont(42);
    h_mMassPhi->GetYaxis()->SetLabelSize(0.04);
    h_mMassPhi->GetYaxis()->SetRangeUser(-0.1*h_mMassPhi->GetMaximum(),1.6*h_mMassPhi->GetMaximum());
    h_mMassPhi->SetNdivisions(505,"Y");
    h_mMassPhi->SetMarkerStyle(24);
    h_mMassPhi->SetMarkerColor(1);
    h_mMassPhi->SetMarkerSize(1.2);
    h_mMassPhi->DrawCopy("PE");
    PlotLine(vmsa::InvMass_low[0],vmsa::InvMass_high[0],0.0,0.0,1,1,2);

    TF1 *f_SigPhi = new TF1("f_SigPhi",PolyBreitWigner,vmsa::BW_Start[0],vmsa::BW_Stop[0],5);
    f_SigPhi->SetParameter(0,vmsa::InvMass[0]);
    f_SigPhi->SetParLimits(0,vmsa::InvMass[0]-1.5*vmsa::Width[0],vmsa::InvMass[0]+1.5*vmsa::Width[0]);
    f_SigPhi->SetParameter(1,vmsa::Width[0]);
    f_SigPhi->SetParLimits(1,0.004,0.070);
    f_SigPhi->SetParameter(2,1.0/numOfEvent);
    f_SigPhi->SetParameter(3,-1.0);
    f_SigPhi->SetParameter(4,1.0);
    f_SigPhi->SetParameter(2,h_mMassPhi->GetMaximum()/f_SigPhi->GetMaximum());
    f_SigPhi->SetRange(vmsa::BW_Start[0],vmsa::BW_Stop[0]);
    h_mMassPhi->Fit(f_SigPhi,"MNQR");
    f_SigPhi->SetLineColor(kGray+2);
    f_SigPhi->SetLineStyle(1);
    f_SigPhi->SetLineWidth(2);
    f_SigPhi->DrawCopy("l same");

    TF1 *f_BkgPhi = new TF1("f_BkgPhi",Poly,vmsa::BW_Start[0],vmsa::BW_Stop[0],2);
    for(int i_par = 0; i_par < 2;++i_par)
    {
      f_BkgPhi->SetParameter(i_par,f_SigPhi->GetParameter(i_par+3));
    }
    f_BkgPhi->SetLineColor(kGray+4);
    f_BkgPhi->SetLineStyle(2);
    f_BkgPhi->SetLineWidth(2);
    f_BkgPhi->DrawCopy("l same");

    plotTopLegend((char*)"a) #phi",0.20,0.87,0.05,1,0.0,42,1);
    plotTopLegend((char*)"Au+Au 27 GeV & 20-60%",0.20,0.80,0.05,1,0.0,42,1);
    plotTopLegend((char*)"|y| < 1.0 & 1.2 < p_{T} < 1.8 GeV/c",0.20,0.73,0.05,1,0.0,42,1);

    TLegend *leg = new TLegend(0.18,0.55,0.55,0.70);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->AddEntry(f_SigPhi,"BW+res. bkg.","l");
    leg->AddEntry(f_BkgPhi,"res. bkg.","l");
    leg->Draw("same");
  }

  c_Signal->cd(2);
  {
    h_mMassKstar->SetTitle("");
    h_mMassKstar->SetStats(0);
    h_mMassKstar->GetXaxis()->SetTitle("M(#pi^{#pm},K^{#mp}) (GeV/c^{2})");
    h_mMassKstar->GetXaxis()->CenterTitle();
    h_mMassKstar->GetXaxis()->SetTitleSize(0.06);
    h_mMassKstar->GetXaxis()->SetTitleOffset(0.9);
    h_mMassKstar->GetXaxis()->SetTitleFont(42);
    h_mMassKstar->GetXaxis()->SetLabelSize(0.04);
    h_mMassKstar->SetNdivisions(505,"X");

    h_mMassKstar->GetYaxis()->SetTitle("Raw Yields");
    h_mMassKstar->GetYaxis()->CenterTitle();
    h_mMassKstar->GetYaxis()->SetTitleOffset(1.14);
    h_mMassKstar->GetYaxis()->SetTitleSize(0.06);
    h_mMassKstar->GetYaxis()->SetTitleFont(42);
    h_mMassKstar->GetYaxis()->SetLabelSize(0.04);
    h_mMassKstar->GetYaxis()->SetRangeUser(-0.1*h_mMassKstar->GetMaximum(),1.6*h_mMassKstar->GetMaximum());
    h_mMassKstar->SetNdivisions(505,"Y");
    h_mMassKstar->SetMarkerStyle(24);
    h_mMassKstar->SetMarkerColor(1);
    h_mMassKstar->SetMarkerSize(1.2);
    h_mMassKstar->DrawCopy("PE");
    h_mMassKstar->Draw("pE");
    PlotLine(0.74,1.06,0.0,0.0,1,1,2);

    f_SigKstar->SetLineColor(kGray+2);
    f_SigKstar->SetLineStyle(1);
    f_SigKstar->SetLineWidth(2);
    f_SigKstar->Draw("l same");

    f_BkgKstar->SetLineColor(kGray+4);
    f_BkgKstar->SetLineStyle(2);
    f_BkgKstar->SetLineWidth(2);
    f_BkgKstar->Draw("l same");

    plotTopLegend((char*)"b) K^{*0}",0.20,0.87,0.05,1,0.0,42,1);
    plotTopLegend((char*)"Au+Au 54.4 GeV & 20-60%",0.20,0.80,0.05,1,0.0,42,1);
    plotTopLegend((char*)"|y| < 1.0 & 2.0 < p_{T} < 2.5 GeV/c",0.20,0.73,0.05,1,0.0,42,1);

    TLegend *leg = new TLegend(0.18,0.20,0.55,0.35);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    // leg->AddEntry(h_mMassKstar,"K^{*0}(54.4 GeV)","p");
    leg->AddEntry(f_SigKstar,"BW+res. bkg.","l");
    leg->AddEntry(f_BkgKstar,"res. bkg.","l");
    leg->Draw("same");
  }

  c_Signal->cd(3);
  {
    TH1F *h_framePhi = new TH1F("h_framePhi","h_framePhi",100,0.0,1.0);
    for(int i_bin = 0; i_bin < 100; ++i_bin)
    {
      h_framePhi->SetBinContent(i_bin+1,-10.0);
      h_framePhi->SetBinError(i_bin+1,1.0);
    }
    h_framePhi->SetTitle("");
    h_framePhi->SetStats(0);
    h_framePhi->GetXaxis()->SetTitle("cos(#theta*)");
    h_framePhi->GetXaxis()->CenterTitle();
    h_framePhi->GetXaxis()->SetTitleSize(0.06);
    h_framePhi->GetXaxis()->SetTitleOffset(0.9);
    h_framePhi->GetXaxis()->SetTitleFont(42);
    h_framePhi->GetXaxis()->SetLabelSize(0.04);
    h_framePhi->SetNdivisions(505,"X");

    h_framePhi->GetYaxis()->SetTitle("Raw Yields");
    h_framePhi->GetYaxis()->CenterTitle();
    h_framePhi->GetYaxis()->SetTitleSize(0.06);
    h_framePhi->GetYaxis()->SetTitleOffset(1.14);
    h_framePhi->GetYaxis()->SetTitleFont(42);
    h_framePhi->GetYaxis()->SetLabelSize(0.04);
    h_framePhi->GetYaxis()->SetRangeUser(575.0/(numOfEvent*binWidth),675.0/(numOfEvent*binWidth));
    h_framePhi->SetNdivisions(505,"Y");
    h_framePhi->SetMarkerStyle(24);
    h_framePhi->SetMarkerColor(1);
    h_framePhi->SetMarkerSize(1.8);
    h_framePhi->DrawCopy("PE");

    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_yieldsPhi,style_phi_2nd,color_phi_2nd,size_marker+0.2);

    TF1 *f_rhoPhi = new TF1("f_rhoPhi",SpinDensity,0.0,1.0,2);
    f_rhoPhi->ReleaseParameter(0);
    f_rhoPhi->ReleaseParameter(1);
    f_rhoPhi->SetParameter(0,0.33);
    f_rhoPhi->SetParameter(1,1000/(numOfEvent*binWidth));
    f_rhoPhi->SetRange(0.0,1.0);
    g_yieldsPhi->Fit(f_rhoPhi,"NQRI");
    f_rhoPhi->SetLineColor(color_phi_2nd);
    f_rhoPhi->SetLineStyle(2);
    f_rhoPhi->SetLineWidth(4);
    // f_rhoPhi->DrawCopy("l same");

    plotTopLegend((char*)"c) #phi",0.20,0.87,0.05,1,0.0,42,1);
    plotTopLegend((char*)"Au+Au 27 GeV & 20-60%",0.20,0.80,0.05,1,0.0,42,1);
    plotTopLegend((char*)"|y| < 1.0 & 1.2 < p_{T} < 1.8 GeV/c",0.20,0.73,0.05,1,0.0,42,1);

    TLegend *leg = new TLegend(0.18,0.18,0.85,0.28);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->AddEntry(f_rhoPhi,"N#times[(1-#rho_{00})+(3#rho_{00}-1)cos^{2}#theta*]","l");
    // leg->Draw("same");
  }

  c_Signal->cd(4);
  {
    TH1F *h_frameKstar = new TH1F("h_frameKstar","h_frameKstar",100,0.0,1.0);
    for(int i_bin = 0; i_bin < 100; ++i_bin)
    {
      h_frameKstar->SetBinContent(i_bin+1,-10.0);
      h_frameKstar->SetBinError(i_bin+1,1.0);
    }
    h_frameKstar->SetTitle("");
    h_frameKstar->SetStats(0);
    h_frameKstar->GetXaxis()->SetTitle("cos(#theta*)");
    h_frameKstar->GetXaxis()->CenterTitle();
    h_frameKstar->GetXaxis()->SetTitleSize(0.06);
    h_frameKstar->GetXaxis()->SetTitleOffset(0.9);
    h_frameKstar->GetXaxis()->SetTitleFont(42);
    h_frameKstar->GetXaxis()->SetLabelSize(0.04);
    h_frameKstar->SetNdivisions(505,"X");

    h_frameKstar->GetYaxis()->SetTitle("Raw Yields");
    h_frameKstar->GetYaxis()->CenterTitle();
    h_frameKstar->GetYaxis()->SetTitleOffset(1.14);
    h_frameKstar->GetYaxis()->SetTitleSize(0.06);
    h_frameKstar->GetYaxis()->SetTitleFont(42);
    h_frameKstar->GetYaxis()->SetLabelSize(0.04);
    h_frameKstar->GetYaxis()->SetRangeUser(4.2e-3,5.0e-3);
    h_frameKstar->SetNdivisions(505,"Y");
    h_frameKstar->SetMarkerStyle(24);
    h_frameKstar->SetMarkerColor(1);
    h_frameKstar->SetMarkerSize(1.8);
    h_frameKstar->DrawCopy("PE");

    // g_yieldsKstar->Draw("pE same");
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_yieldsKstar,style_Kstr,color_Kstr,size_marker-0.4);
    f_rhoKstar->SetLineColor(color_Kstr);
    f_rhoKstar->SetLineStyle(2);
    f_rhoKstar->SetLineWidth(4);
    // f_rhoKstar->Draw("l same");

    plotTopLegend((char*)"d) K^{*0}",0.20,0.87,0.05,1,0.0,42,1);
    plotTopLegend((char*)"Au+Au 54.4 GeV & 20-60%",0.20,0.80,0.05,1,0.0,42,1);
    plotTopLegend((char*)"|y| < 1.0 & 2.0 < p_{T} < 2.5 GeV/c",0.20,0.73,0.05,1,0.0,42,1);

    TLegend *leg = new TLegend(0.18,0.18,0.85,0.28);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->AddEntry(f_rhoKstar,"N#times[(1-#rho_{00})+(3#rho_{00}-1)cos^{2}#theta*]","l");
    // leg->Draw("same");
  }

  c_Signal->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig4_SignalExtraction.eps");
  c_Signal->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig4_SignalExtraction.png");
}

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

double FuncAD(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double D = par[2];
  double R = par[3];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double result = (1.+Bs*D/2.) + (As+D)*CosTheta*CosTheta + (As*D-Bs*D/2.)*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}

void plotFig4_SignalCorrection()
{
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;

  const int style_Kstr = 20;
  const int color_Kstr = kAzure-9;

  const float size_marker = 1.4;
  TGaxis::SetMaxDigits(3);

  // plot phi signal
  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/Second_ActualFit_27GeV_pT_2.root");
  TH1D *h_mYieldsPhi = (TH1D*)File_InputPhi->Get("PtCos")->Clone("h_mYieldsPhi");
  TF1 *f_rhoPhiInput = (TF1*)File_InputPhi->Get("Func_rho")->Clone("f_rhoPhiInput");
  float cos[7] = {1.0/14.0,3.0/14.0,5.0/14.0,7.0/14.0,9.0/14.0,11.0/14.0,13.0/14.0};
  // float yields[7] = {598.598, 596.721, 606.366, 606.203, 615.265, 626.465,  644.074};
  // float errors[7] = {4.06987, 4.06908, 4.12435, 4.1327, 4.1449, 4.2176, 4.26549};

  const long numOfEvent = 229186355;
  // const long numOfEvent = 176995424;
  const double binWidth = 0.00035;
  // cout << "phi-meson: numOfEvent = " << numOfEvent << ", binWidth = " << binWidth << endl;

  TGraphAsymmErrors *g_yieldsPhi = new TGraphAsymmErrors();
  for(int i_cos = 0; i_cos < 7; ++i_cos)
  {
    g_yieldsPhi->SetPoint(i_cos,cos[i_cos],h_mYieldsPhi->GetBinContent(i_cos+1)/(numOfEvent*binWidth));
    g_yieldsPhi->SetPointError(i_cos,1.0/14.0,1.0/14.0,h_mYieldsPhi->GetBinError(i_cos+1)/(numOfEvent*binWidth),h_mYieldsPhi->GetBinError(i_cos+1)/(numOfEvent*binWidth));
  }

  {
    TCanvas *c_SignalPhi = new TCanvas("c_SignalPhi","c_SignalPhi",10,10,800,800);
    c_SignalPhi->cd();
    c_SignalPhi->cd()->SetLeftMargin(0.15);
    c_SignalPhi->cd()->SetTopMargin(0.05);
    c_SignalPhi->cd()->SetBottomMargin(0.15);
    c_SignalPhi->cd()->SetTicks(1,1);
    c_SignalPhi->cd()->SetGrid(0,0);
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

    h_framePhi->GetYaxis()->SetTitle("Yields");
    // h_framePhi->GetYaxis()->SetTitle("Yields (Efficiency Corr.)");
    h_framePhi->GetYaxis()->CenterTitle();
    h_framePhi->GetYaxis()->SetTitleSize(0.05);
    h_framePhi->GetYaxis()->SetTitleOffset(1.14);
    h_framePhi->GetYaxis()->SetTitleFont(42);
    h_framePhi->GetYaxis()->SetLabelSize(0.04);
    h_framePhi->GetYaxis()->SetRangeUser(575.0/(numOfEvent*binWidth),675.0/(numOfEvent*binWidth));
    h_framePhi->SetNdivisions(505,"Y");
    h_framePhi->SetMarkerStyle(24);
    h_framePhi->SetMarkerColor(1);
    h_framePhi->SetMarkerSize(1.8);
    h_framePhi->DrawCopy("PE");

    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_yieldsPhi,style_phi_2nd,color_phi_2nd,size_marker+0.6);

    float scale = f_rhoPhiInput->GetParameter(0)/(numOfEvent*binWidth);
    TF1 *f_rhoPhi = new TF1("f_rhoPhi",FuncAD,0,1.0,4);
    for(int i_par = 0; i_par < 4; ++i_par)
    {
      f_rhoPhi->ReleaseParameter(i_par);
    }
    f_rhoPhi->SetParameter(0,scale);
    f_rhoPhi->SetParameter(1,f_rhoPhiInput->GetParameter(1));
    f_rhoPhi->SetParameter(2,f_rhoPhiInput->GetParameter(2));
    f_rhoPhi->SetParameter(3,f_rhoPhiInput->GetParameter(3));
    f_rhoPhi->SetLineColor(color_phi_2nd);
    f_rhoPhi->SetLineStyle(2);
    f_rhoPhi->SetLineWidth(4);
    f_rhoPhi->Draw("l same");

    plotTopLegend((char*)"#phi",0.20,0.87,0.05,1,0.0,42,1);
    plotTopLegend((char*)"Au+Au 27 GeV & 20-60%",0.20,0.80,0.05,1,0.0,42,1);
    plotTopLegend((char*)"|y| < 1.0 & 1.2 < p_{T} < 1.8 GeV/c",0.20,0.73,0.05,1,0.0,42,1);

    TLegend *leg = new TLegend(0.15,0.15,0.85,0.30);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->AddEntry(f_rhoPhi,"N#times[(1+#frac{B'F}{2})+(A'+F)cos^{2}#theta*+(A'F-#frac{B'F}{2})cos^{4}#theta*]","l");
    // leg->Draw("same");

    c_SignalPhi->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig4_SignalCorrectionPhi.eps");
    c_SignalPhi->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig4_SignalCorrectionPhi.png");
  }

  // TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/Kstar_Fig_signal_54_rap0p1.root");
  TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/Kstar_Fig_signal_54_rap0p1_Oct2021.root");
  TH1F *h_mYieldsKstar = (TH1F*)File_InputKstar->Get("hist_EffAcc_Corrected_Yield_theta_54GeV_Pt_2.00_2.50")->Clone("h_mYieldsKstar");
  TF1 *f_rhoKstar = (TF1*)File_InputKstar->Get("fit_EffAcc_Corrected_Yield_theta_54GeV_Pt_2.00_2.50")->Clone("f_rhoKstar");

  TGraphAsymmErrors *g_yieldsKstar = new TGraphAsymmErrors();
  for(int i_cos = 0; i_cos < 5; ++i_cos)
  {
    g_yieldsKstar->SetPoint(i_cos,0.1+0.2*i_cos,h_mYieldsKstar->GetBinContent(i_cos+1));
    g_yieldsKstar->SetPointError(i_cos,0.1,0.1,h_mYieldsKstar->GetBinError(i_cos+1),h_mYieldsKstar->GetBinError(i_cos+1));
  }

  {
    TCanvas *c_SignalKstar = new TCanvas("c_SignalKstar","c_SignalKstar",10,10,800,800);
    c_SignalKstar->cd();
    c_SignalKstar->cd()->SetLeftMargin(0.15);
    c_SignalKstar->cd()->SetTopMargin(0.05);
    c_SignalKstar->cd()->SetBottomMargin(0.15);
    c_SignalKstar->cd()->SetTicks(1,1);
    c_SignalKstar->cd()->SetGrid(0,0);

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

    // h_frameKstar->GetYaxis()->SetTitle("Yields (Efficiency*Acceptance Corr.)");
    h_frameKstar->GetYaxis()->SetTitle("Yields");
    h_frameKstar->GetYaxis()->CenterTitle();
    h_frameKstar->GetYaxis()->SetTitleOffset(1.14);
    h_frameKstar->GetYaxis()->SetTitleSize(0.05);
    h_frameKstar->GetYaxis()->SetTitleFont(42);
    h_frameKstar->GetYaxis()->SetLabelSize(0.04);
    h_frameKstar->GetYaxis()->SetRangeUser(8.6e-3,10.6e-3);
    h_frameKstar->SetNdivisions(505,"Y");
    h_frameKstar->SetMarkerStyle(24);
    h_frameKstar->SetMarkerColor(1);
    h_frameKstar->SetMarkerSize(1.8);
    h_frameKstar->DrawCopy("PE");

    // g_yieldsKstar->Draw("pE same");
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_yieldsKstar,style_Kstr,color_Kstr,size_marker);
    f_rhoKstar->SetLineColor(color_Kstr);
    f_rhoKstar->SetLineStyle(2);
    f_rhoKstar->SetLineWidth(4);
    f_rhoKstar->Draw("l same");

    plotTopLegend((char*)"K^{*0}",0.20,0.87,0.05,1,0.0,42,1);
    plotTopLegend((char*)"Au+Au 54.4 GeV & 20-60%",0.20,0.80,0.05,1,0.0,42,1);
    plotTopLegend((char*)"|y| < 1.0 & 2.0 < p_{T} < 2.5 GeV/c",0.20,0.73,0.05,1,0.0,42,1);

    TLegend *leg = new TLegend(0.18,0.18,0.85,0.28);
    leg->SetBorderSize(0);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->AddEntry(f_rhoKstar,"N#times[(1-#rho_{00})+(3#rho_{00}-1)cos^{2}#theta*]","l");
    // leg->Draw("same");

    c_SignalKstar->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig4_SignalCorrectionKstar.eps");
    c_SignalKstar->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/Nature/fig4_SignalCorrectionKstar.png");
  }
}

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

#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/functions.h"

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

void plotQA_SignalCorrection()
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
  TH1D *h_mYieldsRaw = (TH1D*)File_InputPhi->Get("raw_yield")->Clone("h_mYieldsPhi");
  TH1D *h_mYieldsEff = (TH1D*)File_InputPhi->Get("PtCos")->Clone("h_mYieldsPhi");
  TF1 *f_rhoPhiInput = (TF1*)File_InputPhi->Get("Func_rho")->Clone("f_rhoPhiInput");
  float cos[7] = {1.0/14.0,3.0/14.0,5.0/14.0,7.0/14.0,9.0/14.0,11.0/14.0,13.0/14.0};
  // float yields[7] = {598.598, 596.721, 606.366, 606.203, 615.265, 626.465,  644.074};
  // float errors[7] = {4.06987, 4.06908, 4.12435, 4.1327, 4.1449, 4.2176, 4.26549};

  const long numOfEvent = 229186355;
  // const long numOfEvent = 176995424;
  const double binWidth = 0.00035;
  // cout << "phi-meson: numOfEvent = " << numOfEvent << ", binWidth = " << binWidth << endl;

  TGraphAsymmErrors *g_yieldsRaw = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_yieldsEff = new TGraphAsymmErrors();
  for(int i_cos = 0; i_cos < 7; ++i_cos)
  {
    g_yieldsRaw->SetPoint(i_cos,cos[i_cos],h_mYieldsRaw->GetBinContent(i_cos+1)/(numOfEvent*binWidth));
    g_yieldsRaw->SetPointError(i_cos,1.0/14.0,1.0/14.0,h_mYieldsRaw->GetBinError(i_cos+1)/(numOfEvent*binWidth),h_mYieldsRaw->GetBinError(i_cos+1)/(numOfEvent*binWidth));

    g_yieldsEff->SetPoint(i_cos,cos[i_cos],h_mYieldsEff->GetBinContent(i_cos+1)/(numOfEvent*binWidth));
    g_yieldsEff->SetPointError(i_cos,1.0/14.0,1.0/14.0,h_mYieldsEff->GetBinError(i_cos+1)/(numOfEvent*binWidth),h_mYieldsEff->GetBinError(i_cos+1)/(numOfEvent*binWidth));
  }

  TCanvas *c_SignalPhi = new TCanvas("c_SignalPhi","c_SignalPhi",10,10,1000,500);
  c_SignalPhi->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_SignalPhi->cd(i_pad+1);
    c_SignalPhi->cd(i_pad+1)->SetLeftMargin(0.15);
    c_SignalPhi->cd(i_pad+1)->SetTopMargin(0.05);
    c_SignalPhi->cd(i_pad+1)->SetBottomMargin(0.15);
    c_SignalPhi->cd(i_pad+1)->SetTicks(1,1);
    c_SignalPhi->cd(i_pad+1)->SetGrid(0,0);
  }
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
  h_framePhi->GetYaxis()->CenterTitle();
  h_framePhi->GetYaxis()->SetTitleSize(0.05);
  h_framePhi->GetYaxis()->SetTitleOffset(1.14);
  h_framePhi->GetYaxis()->SetTitleFont(42);
  h_framePhi->GetYaxis()->SetLabelSize(0.04);
  h_framePhi->GetYaxis()->SetRangeUser(110.0/(numOfEvent*binWidth),140.0/(numOfEvent*binWidth));
  h_framePhi->SetNdivisions(505,"Y");
  h_framePhi->SetMarkerStyle(24);
  h_framePhi->SetMarkerColor(1);
  h_framePhi->SetMarkerSize(1.8);

  c_SignalPhi->cd(1);
  h_framePhi->DrawCopy("pE");

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_yieldsRaw,style_phi_2nd,color_phi_2nd,size_marker+0.6);

  plotTopLegend((char*)"#phi",0.20,0.87,0.05,1,0.0,42,1);
  plotTopLegend((char*)"Au+Au 27 GeV & 20-60%",0.20,0.80,0.05,1,0.0,42,1);
  plotTopLegend((char*)"|y| < 1.0 & 1.2 < p_{T} < 1.8 GeV/c",0.20,0.73,0.05,1,0.0,42,1);

  c_SignalPhi->cd(2);
  h_framePhi->GetYaxis()->SetTitle("Yields (Efficiency Corr.)");
  h_framePhi->GetYaxis()->SetRangeUser(575.0/(numOfEvent*binWidth),675.0/(numOfEvent*binWidth));
  h_framePhi->DrawCopy("pE");

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_yieldsEff,style_phi_2nd,color_phi_2nd,size_marker+0.6);

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
  leg->Draw("same");

  c_SignalPhi->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/phiMeson/c_SigCorrectionSteps.eps");
}

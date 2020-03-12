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

void plotYields_27GeV()
{
  float cos[7] = {1.0/14.0,3.0/14.0,5.0/14.0,7.0/14.0,9.0/14.0,11.0/14.0,13.0/14.0};
  float yields[7] = {598.598, 596.721, 606.366, 606.203, 615.265, 626.465,  644.074};
  float errors[7] = {4.06987, 4.06908, 4.12435, 4.1327, 4.1449, 4.2176, 4.26549};

  TGraphAsymmErrors *g_yields = new TGraphAsymmErrors();
  for(int i_cos = 0; i_cos < 7; ++i_cos)
  {
    g_yields->SetPoint(i_cos,cos[i_cos],yields[i_cos]);
    g_yields->SetPointError(i_cos,1.0/14.0,1.0/14.0,errors[i_cos],errors[i_cos]);
  }

  TCanvas *c_yields = new TCanvas("c_yields","c_yields",10,10,800,800);
  c_yields->cd();
  c_yields->cd()->SetLeftMargin(0.15);
  c_yields->cd()->SetBottomMargin(0.15);
  c_yields->cd()->SetTicks(1,1);
  c_yields->cd()->SetGrid(0,0);

  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,1.0);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("cos(#theta*)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetTitleOffset(0.9);
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->SetNdivisions(505,"X");

  h_play->GetYaxis()->SetTitle("Yields");
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetTitleOffset(1.14);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetRangeUser(580.0,660.0);
  h_play->SetNdivisions(505,"Y");
  h_play->SetMarkerStyle(24);
  h_play->SetMarkerColor(1);
  h_play->SetMarkerSize(1.8);
  h_play->DrawCopy("PE");

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_yields,30,kRed,2.4);

  TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
  f_rho->SetParameter(0,0.33);
  f_rho->SetParameter(1,1000);
  g_yields->Fit(f_rho,"NQMRI");
  f_rho->SetLineColor(2);
  f_rho->SetLineStyle(2);
  f_rho->SetLineWidth(4);
  f_rho->DrawCopy("l same");

  string leg_energy = Form("Au+Au %s & 20-60%% & 2^{nd}-order EP",vmsa::mBeamEnergy[3].c_str());
  plotTopLegend((char*)leg_energy.c_str(),0.2,0.83,0.04,1,0.0,42,1);
  plotTopLegend((char*)"1.2 < p_{T} < 1.8 GeV/c",0.30,0.77,0.04,1,0.0,42,1);

  c_yields->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/BESII/c_yields_27GeV.eps");
  c_yields->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_yields_27GeV.png");
}

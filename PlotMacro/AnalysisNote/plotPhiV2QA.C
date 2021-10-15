#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "../../Utility/functions.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"

using namespace std;

std::pair<double, double> const momentumRange(0.2,5.4);
int const pT_low = 1;
int const pT_high = 14;
int const MarkerColorRP = kGray+2;
int const MarkerStyleRP = 24;

void plotPhiV2QA(int energy = 6)
{
  string InPutFile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Phi/MonteCarlo/Data/v2_200.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  TGraphAsymmErrors *g_v2 = (TGraphAsymmErrors*)File_InPut->Get("Graph");
  TF1 *f_v2 = (TF1*)File_InPut->Get("f_v2");
  // TF1 *f_v2 = new TF1("f_v2",v2_pT_FitFunc,momentumRange.first,momentumRange.second,5);
  // f_v2->FixParameter(0,2);
  // f_v2->SetParameter(1,0.1);
  // f_v2->SetParameter(2,0.1);
  // f_v2->SetParameter(3,0.1);
  // f_v2->SetParameter(4,0.1);
  // f_v2->SetLineColor(kGray+2);
  // f_v2->SetLineWidth(2);
  // f_v2->SetLineStyle(2);
  // g_v2->Fit(f_v2,"N");

  TCanvas *c_v2 = new TCanvas("c_v2","c_v2",10,10,800,800);
  c_v2->cd()->SetLeftMargin(0.15);
  c_v2->cd()->SetBottomMargin(0.15);
  c_v2->cd()->SetTicks(1,1);
  c_v2->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_play->SetBinContent(i_bin,-10.0);
    h_play->SetBinError(i_bin,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetXaxis()->SetRangeUser(momentumRange.first,momentumRange.second);
  h_play->GetXaxis()->SetNdivisions(505);

  h_play->GetYaxis()->SetTitle("v_{2}");
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetRangeUser(0.0,0.2);
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->Draw("pE");
  f_v2->Draw("l same");

  g_v2->SetMarkerStyle(30);
  g_v2->SetMarkerColor(2);
  g_v2->SetLineColor(2);
  g_v2->SetMarkerSize(2.4);
  g_v2->Draw("pE same");

  plotTopLegend((char*)"AuAu 200 GeV 0%-80%",0.4,0.35,0.04,1,0.0,42,1);
  TLegend *legv2 = new TLegend(0.2,0.7,0.65,0.85);
  legv2->SetBorderSize(0.0);
  legv2->SetFillColor(10);
  legv2->AddEntry(g_v2,"#phi STAR PRL 116 062301","p");
  legv2->AddEntry(f_v2,"fit","l");
  legv2->Draw("same");

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/phiMeson/c_V2QA_%s.eps",vmsa::mBeamEnergy[6].c_str());
  c_v2->SaveAs(FigName.c_str());
}

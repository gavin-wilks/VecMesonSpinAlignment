#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"

using namespace std;

void plotRho00EnergySideBand()
{
  gStyle->SetOptDate(0);
  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/SideBand/rhoEnergySideBand.root");
  TGraphAsymmErrors *g_rhoSideBand_2nd = (TGraphAsymmErrors*)File_Input->Get("g_rhoSideBand_2nd");

  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,0,1000);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  c_rho00->cd()->SetLogx();
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(9.0,240.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.25,0.38);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(9.0,240.0,1.0/3.0,1.0/3.0,1,3,2);

  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoSideBand_2nd,30,kRed,1.8);
  Draw_TGAE_Point_new_Symbol(80,0.36,0.0,0.0,0.0,0.0,30,kRed,1.8);
  plotTopLegend((char*)"2^{nd}-order EP",90,0.3585,0.04,1,0.0,42,0);
  plotTopLegend((char*)"Side Band",90,0.35,0.04,1,0.0,42,0);

  // plotTopLegend((char*)"#rho_{00} = 1/3",100,0.328,0.04,1,0.0,42,0);

  plotTopLegend((char*)"Au+Au (20-60\% & |#eta| < 1)",0.45,0.30,0.04,1,0.0,42,1);
  plotTopLegend((char*)"#phi-meson (1.2 < p_{T}< 5.4 GeV/c)",0.4,0.25,0.04,1,0.0,42,1);

  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_rhoEnergySideBand.eps");
  c_rho00->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_rhoEnergySideBand.png");
}


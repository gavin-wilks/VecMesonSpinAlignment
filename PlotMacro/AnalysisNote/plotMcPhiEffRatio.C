#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/draw.h"
#include "../../Utility/functions.h"

using namespace std;

void plotMcPhiEffRatio(int cent = 9, int ptBin = 1)
{
  string inputfile_flow = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu200GeV/Phi/Efficiency/Eff_200GeV_SingleKaon_2060_withFlowSpecPtCut.root");
  TFile *File_InPutFlow = TFile::Open(inputfile_flow.c_str());
  TH1D *h_mEffCosFlow[10][vmsa::BinPt];
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_pt = vmsa::pt_rebin_first[6]; i_pt < vmsa::pt_rebin_last[6]; ++i_pt) // use rebinned pt
    {
      string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mEffCosFlow[i_cent][i_pt] = (TH1D*)File_InPutFlow->Get(HistName.c_str());
    }
  }

  string inputfile_noflow = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu200GeV/Phi/Efficiency/Eff_200GeV_SingleKaon_2060_woFlowSpecPtCut.root");
  TFile *File_InPutNoFlow = TFile::Open(inputfile_noflow.c_str());
  TH1D *h_mEffCosNoFlow[10][vmsa::BinPt];
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_pt = vmsa::pt_rebin_first[6]; i_pt < vmsa::pt_rebin_last[6]; ++i_pt) // use rebinned pt
    {
      string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mEffCosNoFlow[i_cent][i_pt] = (TH1D*)File_InPutNoFlow->Get(HistName.c_str());
    }
  }

  TH1D *h_mEffCosRatio[10][vmsa::BinPt];
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_pt = vmsa::pt_rebin_first[6]; i_pt < vmsa::pt_rebin_last[6]; ++i_pt) // use rebinned pt
    {
      string HistName = Form("h_mEffCosRatio_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mEffCosRatio[i_cent][i_pt] = (TH1D*)h_mEffCosFlow[i_cent][i_pt]->Clone(HistName.c_str());
      h_mEffCosRatio[i_cent][i_pt]->Divide(h_mEffCosNoFlow[i_cent][i_pt]);
    }
  }

  TCanvas *c_eff = new TCanvas("c_eff","c_eff",10,10,800,800);
  c_eff->cd()->SetLeftMargin(0.15);
  c_eff->cd()->SetBottomMargin(0.15);
  c_eff->cd()->SetGrid(0,0);
  c_eff->cd()->SetTicks(1,1);

  c_eff->cd(1);
  string title = Form("#phi-meson: %1.1f < p_{T} < %1.1f GeV/c",vmsa::pt_low[6][ptBin],vmsa::pt_up[6][ptBin]);
  h_mEffCosRatio[cent][ptBin]->SetTitle(title.c_str());
  h_mEffCosRatio[cent][ptBin]->SetStats(0);
  h_mEffCosRatio[cent][ptBin]->GetXaxis()->SetTitle("cos(#theta*)");
  h_mEffCosRatio[cent][ptBin]->GetXaxis()->SetTitleSize(0.06);
  h_mEffCosRatio[cent][ptBin]->GetXaxis()->CenterTitle();
  h_mEffCosRatio[cent][ptBin]->GetXaxis()->SetLabelSize(0.04);
  h_mEffCosRatio[cent][ptBin]->GetXaxis()->SetNdivisions(505);
  h_mEffCosRatio[cent][ptBin]->GetXaxis()->SetRangeUser(0.0,1.0);

  h_mEffCosRatio[cent][ptBin]->GetYaxis()->SetTitle("eff^{w. v_{2}}/eff^{w/o v_{2}}");
  h_mEffCosRatio[cent][ptBin]->GetYaxis()->SetTitleSize(0.06);
  h_mEffCosRatio[cent][ptBin]->GetYaxis()->CenterTitle();
  h_mEffCosRatio[cent][ptBin]->GetYaxis()->SetLabelSize(0.04);
  h_mEffCosRatio[cent][ptBin]->GetYaxis()->SetNdivisions(505);
  h_mEffCosRatio[cent][ptBin]->GetYaxis()->SetTitleOffset(1.2);
  h_mEffCosRatio[cent][ptBin]->GetYaxis()->SetRangeUser(0.8,1.2);
  h_mEffCosRatio[cent][ptBin]->SetMarkerStyle(20);
  h_mEffCosRatio[cent][ptBin]->SetMarkerColor(kGray+2);
  h_mEffCosRatio[cent][ptBin]->SetMarkerSize(1.4);
  h_mEffCosRatio[cent][ptBin]->SetLineColor(kGray+2);
  h_mEffCosRatio[cent][ptBin]->Draw("pE");
  PlotLine(0.0,1.0,1.0,1.0,1,2,2);

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/phiMeson/c_effMcPhiRatioAuAu%sPt%d.eps",vmsa::mBeamEnergy[6].c_str(),ptBin);
  c_eff->SaveAs(FigName.c_str());
}

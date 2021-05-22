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

void plotMcPhiEffEP(int cent = 9, int ptBin = 2, bool isFlow = false)
{
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu200GeV/Phi/Efficiency/Eff_200GeV_SingleKaon_2060_withFlowSpecPtCut.root");
  if(!isFlow) inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu200GeV/Phi/Efficiency/Eff_200GeV_SingleKaon_2060_woFlowSpecPtCut.root");
  // string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu200GeV/Phi/Efficiency/Eff_200GeV_SingleKaon_2060_withFlow.root");
  // if(!isFlow) inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu200GeV/Phi/Efficiency/Eff_200GeV_SingleKaon_2060_woFlow.root");
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH2D *h_mMcCosEP[10][vmsa::BinPt]; // efficiency vs CosThetaStar & EP as a function of centrality and pt
  TH2D *h_mRcCosEP[10][vmsa::BinPt];
  TH2D *h_mEffCosEP[10][vmsa::BinPt];
  TH1D *h_mEffCos[10][vmsa::BinPt];

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    for(int i_pt = vmsa::pt_rebin_first[6]; i_pt < vmsa::pt_rebin_last[6]; ++i_pt) // use rebinned pt
    {
      string HistName;
      HistName = Form("h_mMcCosEP_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mMcCosEP[i_cent][i_pt] = (TH2D*)File_InPut->Get(HistName.c_str());
      h_mMcCosEP[i_cent][i_pt]->Sumw2();

      HistName = Form("h_mRcCosEP_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mRcCosEP[i_cent][i_pt] = (TH2D*)File_InPut->Get(HistName.c_str());
      h_mRcCosEP[i_cent][i_pt]->Sumw2();
      
      HistName = Form("h_mEffCosEP_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mEffCosEP[i_cent][i_pt] = (TH2D*)h_mRcCosEP[i_cent][i_pt]->Clone(HistName.c_str());
      h_mEffCosEP[i_cent][i_pt]->Divide(h_mMcCosEP[i_cent][i_pt]);

      HistName = Form("h_mEffCos_Cent_%d_Pt_%d",i_cent,i_pt);
      h_mEffCos[i_cent][i_pt] = (TH1D*)File_InPut->Get(HistName.c_str());
    }
  }

  TCanvas *c_eff = new TCanvas("c_eff","c_eff",10,10,800,400);
  c_eff->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_eff->cd(i_pad+1)->SetLeftMargin(0.15);
    c_eff->cd(i_pad+1)->SetBottomMargin(0.15);
    c_eff->cd(i_pad+1)->SetGrid(0,0);
    c_eff->cd(i_pad+1)->SetTicks(1,1);
  }

  c_eff->cd(1);
  c_eff->cd(1)->SetRightMargin(0.15);
  h_mEffCosEP[cent][ptBin]->SetStats(0);
  string title = Form("%1.1f < p_{T} < %1.1f GeV/c",vmsa::pt_low[6][ptBin],vmsa::pt_up[6][ptBin]);
  h_mEffCosEP[cent][ptBin]->SetTitle(title.c_str());
  h_mEffCosEP[cent][ptBin]->GetXaxis()->SetTitle("cos(#theta*)");
  h_mEffCosEP[cent][ptBin]->GetXaxis()->SetTitleSize(0.06);
  h_mEffCosEP[cent][ptBin]->GetXaxis()->CenterTitle();
  h_mEffCosEP[cent][ptBin]->GetYaxis()->SetTitle("#phi-#Psi_{2}");
  h_mEffCosEP[cent][ptBin]->GetYaxis()->SetTitleSize(0.06);
  h_mEffCosEP[cent][ptBin]->GetYaxis()->CenterTitle();
  h_mEffCosEP[cent][ptBin]->Draw("colz");

  c_eff->cd(2);
  string legEnergy = Form("AuAu %s 20%%-60%%",vmsa::mBeamEnergy[6].c_str());
  h_mEffCos[cent][ptBin]->SetTitle("#phi-meson");
  h_mEffCos[cent][ptBin]->SetStats(0);
  h_mEffCos[cent][ptBin]->GetXaxis()->SetTitle("cos(#theta*)");
  h_mEffCos[cent][ptBin]->GetXaxis()->SetTitleSize(0.06);
  h_mEffCos[cent][ptBin]->GetXaxis()->CenterTitle();
  h_mEffCos[cent][ptBin]->GetXaxis()->SetLabelSize(0.04);
  h_mEffCos[cent][ptBin]->GetXaxis()->SetNdivisions(505);
  h_mEffCos[cent][ptBin]->GetXaxis()->SetRangeUser(0.0,1.0);

  h_mEffCos[cent][ptBin]->GetYaxis()->SetTitle("efficiency");
  h_mEffCos[cent][ptBin]->GetYaxis()->SetTitleSize(0.06);
  h_mEffCos[cent][ptBin]->GetYaxis()->CenterTitle();
  h_mEffCos[cent][ptBin]->GetYaxis()->SetLabelSize(0.04);
  h_mEffCos[cent][ptBin]->GetYaxis()->SetNdivisions(505);
  h_mEffCos[cent][ptBin]->GetYaxis()->SetTitleOffset(1.2);
  h_mEffCos[cent][ptBin]->GetYaxis()->SetRangeUser(0.0,0.3);
  h_mEffCos[cent][ptBin]->SetMarkerStyle(20);
  h_mEffCos[cent][ptBin]->SetMarkerColor(kGray+2);
  h_mEffCos[cent][ptBin]->SetMarkerSize(1.4);
  h_mEffCos[cent][ptBin]->SetLineColor(kGray+2);
  h_mEffCos[cent][ptBin]->Draw("pE");

  TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
  f_rho->SetParameter(0,0.33);
  f_rho->SetParameter(1,h_mEffCos[cent][ptBin]->GetMaximum());
  h_mEffCos[cent][ptBin]->Fit(f_rho,"NMRI");
  f_rho->SetLineColor(2);
  f_rho->SetLineStyle(2);
  f_rho->Draw("l same");

  string leg = Form("#rho_{00} = %1.3f #pm %1.3f",f_rho->GetParameter(0),f_rho->GetParError(0));
  plotTopLegend((char*)leg.c_str(),0.25,0.75,0.05,1,0.0,42,1);

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/phiMeson/c_effMcPhiEpAuAu%sPt%d_withFlowSpecPtCut.eps",vmsa::mBeamEnergy[6].c_str(),ptBin);
  if(!isFlow) FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/phiMeson/c_effMcPhiEpAuAu%sPt%d_woFlowSpecPtCut.eps",vmsa::mBeamEnergy[6].c_str(),ptBin);
  // string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/phiMeson/c_effMcPhiEpAuAu%sPt%d_withFlow.eps",vmsa::mBeamEnergy[6].c_str(),ptBin);
  // if(!isFlow) FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/phiMeson/c_effMcPhiEpAuAu%sPt%d_woFlow.eps",vmsa::mBeamEnergy[6].c_str(),ptBin);
  c_eff->SaveAs(FigName.c_str());
}

#include <string>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLegend.h>
#include <map>
#include "../../Utility/draw.h"

typedef std::map<std::string,TH1D*> TH1DMap;

using namespace std;

double effFit(double *x_val, double *par)
{
  double x = x_val[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];

  double y = p0 * TMath::Exp(-1.0*TMath::Power(p1/x,p2));

  return y;
}

void plotQA_TpcEffRatio(int mEnergy = 1, int mPID = 0, int year = 1)
{
  string const mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  string const mParType[2] = {"Kplus","Kminus"};
  string const mParTex[2] = {"K^{+}","K^{-}"};
  string const Centrality[10] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%","20%-60%"}; // Centrality bin
  string const mYear[2] = {"run11","run10"};
  int const mPlotStyle[10] = {0,0,1,1,1,1,0,0,0,1};
  int const mPlotColor[10] = {0,0,kGray+3,kAzure+4,kOrange+1,kCyan+1,0,0,0,kRed};

  const double ptStart = 0.0;
  const double ptStop  = 5.0;

  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr_2060.root",mBeamEnergy[mEnergy].c_str(),mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str(),mYear[year].c_str());
  cout << "open input file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH1DMap h_mEfficiency;
  TH1DMap h_mEffRatio;
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistName;

    HistName = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
    cout << "read in => " << HistName.c_str() << endl;
    h_mEffRatio[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1600,800);
  c_play->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_play->SetLeftMargin(0.15);
    c_play->SetBottomMargin(0.15);
    c_play->SetGrid(0,0);
    c_play->SetTicks(1,1);
  }

  c_play->cd(1);
  string HistNameEff = "h_mEffPt_Cent_9";
  string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  h_mEfficiency[HistNameEff]->SetTitle(title.c_str());
  h_mEfficiency[HistNameEff]->SetStats(0);

  h_mEfficiency[HistNameEff]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mEfficiency[HistNameEff]->GetXaxis()->CenterTitle();
  h_mEfficiency[HistNameEff]->GetXaxis()->SetRangeUser(ptStart,ptStop);

  h_mEfficiency[HistNameEff]->GetYaxis()->SetTitle("Efficiency");
  h_mEfficiency[HistNameEff]->GetYaxis()->CenterTitle();
  h_mEfficiency[HistNameEff]->GetYaxis()->SetRangeUser(0.0,1.05);
  h_mEfficiency[HistNameEff]->SetLineColor(mPlotColor[9]);
  h_mEfficiency[HistNameEff]->SetLineWidth(4);
  h_mEfficiency[HistNameEff]->DrawCopy("HIST");

  for(int i_cent = 2; i_cent <= 5; ++i_cent)
  {
    HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEfficiency[HistNameEff]->SetLineColor(mPlotColor[i_cent]);
    h_mEfficiency[HistNameEff]->SetLineWidth(2);
    h_mEfficiency[HistNameEff]->DrawCopy("HIST same");
  }
  HistNameEff = "h_mEffPt_Cent_9";
  h_mEfficiency[HistNameEff]->DrawCopy("HIST same");

  TF1 *f_eff = new TF1("f_eff",effFit,0,10,3);
  f_eff->SetParameter(0,1.0);
  f_eff->SetParameter(1,1.0);
  f_eff->SetParameter(2,1.0);
  f_eff->SetRange(ptStart,ptStop);
  h_mEfficiency[HistNameEff]->Fit(f_eff,"NR");
  f_eff->SetLineColor(2);
  f_eff->SetLineStyle(2);
  f_eff->SetLineWidth(4);
  f_eff->Draw("l same");

  TLegend *leg = new TLegend(0.5,0.2,0.8,0.5);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    if(mPlotStyle[i_cent] > 0)
    {
      HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
      leg->AddEntry(h_mEfficiency[HistNameEff],Centrality[i_cent].c_str());
    }
  }
  leg->AddEntry(f_eff,"fit to 20-60%","l");
  leg->Draw("same");


  // get the ratio to the fit
  c_play->cd(2);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetXaxis()->SetRangeUser(ptStart,ptStop);
  h_play->GetYaxis()->SetTitle("Ratio");
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetRangeUser(0.8,1.2);
  h_play->Draw("pE");

  for(int i_cent = 2; i_cent <= 5; ++i_cent)
  {
    HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEffRatio[HistNameEff]->Divide(f_eff);
    h_mEffRatio[HistNameEff]->SetLineColor(mPlotColor[i_cent]);
    h_mEffRatio[HistNameEff]->SetLineWidth(2);
    h_mEffRatio[HistNameEff]->DrawCopy("HIST same");
  }
  HistNameEff = "h_mEffPt_Cent_9";
  h_mEffRatio[HistNameEff]->Divide(f_eff);
  h_mEffRatio[HistNameEff]->SetLineColor(mPlotColor[9]);
  h_mEffRatio[HistNameEff]->SetLineWidth(4);
  h_mEffRatio[HistNameEff]->DrawCopy("HIST same");

  PlotLine(ptStart, ptStop, 1.0, 1.0, 2, 4, 2);

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/TPC/c_TpcEffRatio%s%s.eps",mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_play->SaveAs(FigName.c_str());
}

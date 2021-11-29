#include <string>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <map>
#include "../../Utility/draw.h"

typedef std::map<std::string,TH1D*> TH1DMap;

using namespace std;

// tof matching efficiency
double tof_Kaon(double* x, double* par)
{
  return par[0]*(1.0 / (pow(x[0] - par[1], 2) + par[2]) - par[4] / (exp(x[0] - par[3]) + par[5]) + par[6]);
}

void plotQA_TofEffRatio(int mEnergy = 6, int mPID = 0)
{
  string const mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  string const mParType[2] = {"Kplus","Kminus"};
  string const mParTex[2] = {"K^{+}","K^{-}"};
  string const Centrality[10] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%","20%-60%"}; // Centrality bin
  int const mPlotStyle[10] = {0,0,1,1,1,1,0,0,0,1};
  int const mPlotColor[10] = {0,0,kGray+3,kAzure+4,kOrange+1,kCyan+1,0,0,0,kRed};

  const double ptStart = 0.0;
  const double ptStop  = 8.0;

  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/Eff_%s_ToFMatch_2060.root",mBeamEnergy[mEnergy].c_str(),mBeamEnergy[mEnergy].c_str());
  cout << "open input file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH1DMap h_mEfficiency;
  TH1DMap h_mEffRatio;
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistName;
    HistName = Form("h_mEfficiency_%s_Cent_%d",mParType[mPID].c_str(),i_cent);
    h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str());
    cout << "read in => " << HistName.c_str() << endl;
    h_mEffRatio[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
  }

  string inputEff = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/FitPar_AuAu%s_%s_first_2060.root",mBeamEnergy[mEnergy].c_str(),mBeamEnergy[mEnergy].c_str(),mParType[mPID].c_str());
  TFile *File_Eff = TFile::Open(inputEff.c_str());
  cout << "OPEN ToF Matching Efficiency Fit File : " << inputEff.c_str() << endl;
  string KEY = Form("h_mFitParameters_%s_Cent_9",mParType[mPID].c_str());
  TH1D *h_TofFit = (TH1D*)File_Eff->Get(KEY.c_str());
  TF1 *f_TofFit = new TF1("f_TofFit",tof_Kaon,0.2,10,7);
  for(int i_par = 0; i_par < 7; ++i_par)
  {
    f_TofFit->FixParameter(i_par,h_TofFit->GetBinContent(i_par+1));
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1600,800);
  c_play->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_play->cd(i_pad+1)->SetLeftMargin(0.15);
    c_play->cd(i_pad+1)->SetBottomMargin(0.15);
    c_play->cd(i_pad+1)->SetGrid(0,0);
    c_play->cd(i_pad+1)->SetTicks(1,1);
  }

  {
    c_play->cd(1);
    string HistName = Form("h_mEfficiency_%s_Cent_9",mParType[mPID].c_str());

    string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
    h_mEfficiency[HistName]->SetTitle(title.c_str());
    h_mEfficiency[HistName]->SetStats(0);

    h_mEfficiency[HistName]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_mEfficiency[HistName]->GetXaxis()->CenterTitle();

    h_mEfficiency[HistName]->GetYaxis()->SetTitle("Efficiency");
    h_mEfficiency[HistName]->GetYaxis()->CenterTitle();
    h_mEfficiency[HistName]->GetYaxis()->SetRangeUser(0.0,1.05);
    h_mEfficiency[HistName]->SetLineColor(mPlotColor[9]);
    h_mEfficiency[HistName]->SetLineWidth(4);
    h_mEfficiency[HistName]->DrawCopy("HIST");

    for(int i_cent = 2; i_cent <= 5; ++i_cent)
    {
      HistName = Form("h_mEfficiency_%s_Cent_%d",mParType[mPID].c_str(),i_cent);
      h_mEfficiency[HistName]->SetLineColor(mPlotColor[i_cent]);
      h_mEfficiency[HistName]->SetLineWidth(2);
      h_mEfficiency[HistName]->DrawCopy("HIST same");
    }
    HistName = Form("h_mEfficiency_%s_Cent_9",mParType[mPID].c_str());
    h_mEfficiency[HistName]->DrawCopy("HIST same");

    f_TofFit->SetLineColor(2);
    f_TofFit->SetLineWidth(3);
    f_TofFit->SetLineStyle(2);
    f_TofFit->Draw("l same");

    TLegend *leg = new TLegend(0.5,0.2,0.8,0.5);
    leg->SetFillColor(10);
    leg->SetBorderSize(0);
    for(int i_cent = 0; i_cent < 10; ++i_cent)
    {
      if(mPlotStyle[i_cent] > 0)
      {
	HistName = Form("h_mEfficiency_%s_Cent_%d",mParType[mPID].c_str(),i_cent);
	leg->AddEntry(h_mEfficiency[HistName],Centrality[i_cent].c_str());
      }
    }
    leg->AddEntry(f_TofFit,"fit to 20-60%","l");
    leg->Draw("same");
  }

  {
    // get the ratio to the fit
    c_play->cd(2);
    TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
    for(int i_bin = 0; i_bin < 100; ++i_bin)
    {
      h_play->SetBinContent(i_bin+1,-10.0);
      h_play->SetBinError(i_bin+1,1.0);
    }
    h_play->SetTitle("Ratio to Fit");
    h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_play->GetXaxis()->CenterTitle();
    h_play->GetXaxis()->SetRangeUser(ptStart,ptStop);
    h_play->GetYaxis()->SetTitle("Ratio");
    h_play->GetYaxis()->CenterTitle();
    h_play->GetYaxis()->SetRangeUser(0.8,1.2);
    h_play->Draw("pE");

    string HistNameEff;
    for(int i_cent = 2; i_cent <= 5; ++i_cent)
    {
      HistNameEff = Form("h_mEfficiency_%s_Cent_%d",mParType[mPID].c_str(),i_cent);
      h_mEffRatio[HistNameEff]->Divide(f_TofFit);
      h_mEffRatio[HistNameEff]->SetLineColor(mPlotColor[i_cent]);
      h_mEffRatio[HistNameEff]->SetLineWidth(2);
      h_mEffRatio[HistNameEff]->DrawCopy("HIST same");
    }
    HistNameEff = Form("h_mEfficiency_%s_Cent_9",mParType[mPID].c_str());
    h_mEffRatio[HistNameEff]->Divide(f_TofFit);
    h_mEffRatio[HistNameEff]->SetLineColor(mPlotColor[9]);
    h_mEffRatio[HistNameEff]->SetLineWidth(4);
    h_mEffRatio[HistNameEff]->DrawCopy("HIST same");

    PlotLine(ptStart, ptStop, 1.0, 1.0, 4, 4, 2);

    PlotLine(0.4, 0.4, 0.8, 1.2, kGray+2, 4, 2);
    PlotLine(0.8, 0.8, 0.8, 1.2, kGray+2, 4, 2);
    PlotLine(1.0, 1.0, 0.8, 1.2, kGray, 4, 2);
    PlotLine(2.5, 2.5, 0.8, 1.2, kGray, 4, 2);
  }


  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/ToF/c_TofEffRatio%s%s.eps",mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_play->SaveAs(FigName.c_str());
}

#include <string>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <map>

typedef std::map<std::string,TH1D*> TH1DMap;

using namespace std;

void plotQA_TpcTrackingEfficiency(int mEnergy = 4, int mPID = 0, int year = 1)
{
  string const mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  string const mParType[2] = {"Kplus","Kminus"};
  string const mParTex[2] = {"K^{+}","K^{-}"};
  string const Centrality[10] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%","20%-60%"}; // Centrality bin
  string const mYear[2] = {"run11","run10"};
  int const mPlotStyle[10] = {0,0,1,1,1,1,0,0,0,1};
  int const mPlotColor[10] = {0,0,kGray+3,kAzure+4,kOrange+1,kCyan+1,0,0,0,kRed};

  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent_%s_pr_2060.root",mBeamEnergy[mEnergy].c_str(),mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str(),mYear[year].c_str());
  cout << "open input file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH1DMap h_mEfficiency;
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistName;

    HistName = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str());
    cout << "read in => " << HistName.c_str() << endl;

    HistName = Form("h_mEffEta_Cent_%d",i_cent);
    h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str());
    cout << "read in => " << HistName.c_str() << endl;

    HistName = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str());
    cout << "read in => " << HistName.c_str() << endl;
    for(int i_eta = 0; i_eta < 10; ++i_eta)
    {
      for(int i_phi = 0; i_phi < 12; ++i_phi)
      {
	HistName = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str());
	cout << "read in => " << HistName.c_str() << endl;
      }
    }
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetGrid(0,0);
  c_play->SetTicks(1,1);

  string HistNameEff = "h_mEffPt_Cent_9";
  string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  h_mEfficiency[HistNameEff]->SetTitle(title.c_str());
  h_mEfficiency[HistNameEff]->SetStats(0);

  h_mEfficiency[HistNameEff]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mEfficiency[HistNameEff]->GetXaxis()->CenterTitle();
  h_mEfficiency[HistNameEff]->GetXaxis()->SetRangeUser(0.0,5.00);

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
  leg->Draw("same");

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/TPC/c_TpcEffCentCom_%s%s.eps",mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_play->SaveAs(FigName.c_str());

  int const eta_bin = 4;

  TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1600,1200);
  c_diff->Divide(4,3);
  for(int i_pad = 0; i_pad < 12; ++i_pad)
  {
    c_diff->cd(i_pad+1)->SetLeftMargin(0.15);
    c_diff->cd(i_pad+1)->SetBottomMargin(0.15);
    c_diff->cd(i_pad+1)->SetGrid(0,0);
    c_diff->cd(i_pad+1)->SetTicks(1,1);
    string HistName = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",eta_bin,i_pad);

    string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
    h_mEfficiency[HistName]->SetTitle(title.c_str());
    h_mEfficiency[HistName]->SetStats(0);

    h_mEfficiency[HistName]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_mEfficiency[HistName]->GetXaxis()->CenterTitle();
    h_mEfficiency[HistName]->GetXaxis()->SetRangeUser(0.0,5.00);

    h_mEfficiency[HistName]->GetYaxis()->SetTitle("Efficiency");
    h_mEfficiency[HistName]->GetYaxis()->CenterTitle();
    h_mEfficiency[HistName]->GetYaxis()->SetRangeUser(0.0,1.05);
    h_mEfficiency[HistName]->SetLineColor(mPlotColor[9]);
    h_mEfficiency[HistName]->SetLineWidth(2);
    h_mEfficiency[HistName]->DrawCopy("HIST");

    for(int i_cent = 2; i_cent <= 5; ++i_cent)
    {
      HistName = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,eta_bin,i_pad);
      h_mEfficiency[HistName]->SetLineColor(mPlotColor[i_cent]);
      h_mEfficiency[HistName]->SetLineWidth(1);
      h_mEfficiency[HistName]->DrawCopy("HIST same");
    }
    HistName = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",eta_bin,i_pad);
    h_mEfficiency[HistName]->DrawCopy("HIST same");

    string phi_leg = Form("phi bin: %d",i_pad);
    TLegend *leg_phi = new TLegend(0.3,0.2,0.8,0.3);
    leg_phi->SetFillColor(10);
    leg_phi->SetBorderSize(0);
    leg_phi->AddEntry(h_mEfficiency[HistName],phi_leg.c_str());
    leg_phi->Draw("same");
  }
  FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/TPC/c_TpcEffDiffComEta_%d_%s%s.eps",eta_bin,mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_diff->SaveAs(FigName.c_str());

  int const phi_bin = 4;
  TCanvas *c_diff_phi = new TCanvas("c_diff_phi","c_diff_phi",10,10,1000,500);
  c_diff_phi->Divide(5,2);
  for(int i_pad = 0; i_pad < 10; ++i_pad)
  {
    c_diff_phi->cd(i_pad+1)->SetLeftMargin(0.15);
    c_diff_phi->cd(i_pad+1)->SetBottomMargin(0.15);
    c_diff_phi->cd(i_pad+1)->SetGrid(0,0);
    c_diff_phi->cd(i_pad+1)->SetTicks(1,1);
    string HistName = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",i_pad,phi_bin);

    string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
    h_mEfficiency[HistName]->SetTitle(title.c_str());
    h_mEfficiency[HistName]->SetStats(0);

    h_mEfficiency[HistName]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_mEfficiency[HistName]->GetXaxis()->CenterTitle();
    h_mEfficiency[HistName]->GetXaxis()->SetRangeUser(0.0,5.00);

    h_mEfficiency[HistName]->GetYaxis()->SetTitle("Efficiency");
    h_mEfficiency[HistName]->GetYaxis()->CenterTitle();
    h_mEfficiency[HistName]->GetYaxis()->SetRangeUser(0.0,1.05);
    h_mEfficiency[HistName]->SetLineColor(mPlotColor[9]);
    h_mEfficiency[HistName]->SetLineWidth(2);
    h_mEfficiency[HistName]->DrawCopy("HIST");

    for(int i_cent = 2; i_cent <= 5; ++i_cent)
    {
      HistName = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_pad,phi_bin);
      h_mEfficiency[HistName]->SetLineColor(mPlotColor[i_cent]);
      h_mEfficiency[HistName]->SetLineWidth(1);
      h_mEfficiency[HistName]->DrawCopy("HIST same");
    }
    HistName = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",i_pad,phi_bin);
    h_mEfficiency[HistName]->DrawCopy("HIST same");

    string eta_leg = Form("#eta bin: %d",i_pad);
    TLegend *leg_eta = new TLegend(0.3,0.2,0.8,0.3);
    leg_eta->SetFillColor(10);
    leg_eta->SetBorderSize(0);
    leg_eta->AddEntry(h_mEfficiency[HistName],eta_leg.c_str());
    leg_eta->Draw("same");
  }
  FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/Efficiency/TPC/c_TpcEffDiffComPhi_%d_%s%s.eps",phi_bin,mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_diff_phi->SaveAs(FigName.c_str());
}

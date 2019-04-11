#include <string>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <map>

typedef std::map<std::string,TH1D*> TH1DMap;

using namespace std;

void plotQA_TpcTrackingEfficiency(int mEnergy = 6, int mPID = 0)
{
  string const mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  string const mParType[2] = {"Kplus","Kminus"};
  string const mParTex[2] = {"K^{+}","K^{-}"};
  string const Centrality[10] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%","20%-60%"}; // Centrality bin
  int const mPlotStyle[10] = {0,0,1,1,1,1,0,0,0,1};

  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Embedding/%s/Efficiency/Eff_%s_StMcEvent_run11_pr_2060.root",mBeamEnergy[mEnergy].c_str(),mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  cout << "open input file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH1DMap h_mEfficiency;
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistNameEff;

    HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEfficiency[HistNameEff] = (TH1D*)File_InPut->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;

    HistNameEff = Form("h_mEffEta_Cent_%d",i_cent);
    h_mEfficiency[HistNameEff] = (TH1D*)File_InPut->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;

    HistNameEff = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEfficiency[HistNameEff] = (TH1D*)File_InPut->Get(HistNameEff.c_str());
    cout << "read in => " << HistNameEff.c_str() << endl;
    for(int i_eta = 0; i_eta < 10; ++i_eta)
    {
      for(int i_phi = 0; i_phi < 12; ++i_phi)
      {
	HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	h_mEfficiency[HistNameEff] = (TH1D*)File_InPut->Get(HistNameEff.c_str());
	cout << "read in => " << HistNameEff.c_str() << endl;
      }
    }
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1500,500);
  c_play->Divide(3,1);
  for(int i_pad = 0; i_pad < 3; ++i_pad)
  {
    c_play->cd(i_pad+1)->SetLeftMargin(0.15);
    c_play->cd(i_pad+1)->SetBottomMargin(0.15);
    c_play->cd(i_pad+1)->SetGrid(0,0);
    c_play->cd(i_pad+1)->SetTicks(1,1);

    string HistNameEff;

    for(int i_cent = 2; i_cent <= 5; ++i_cent)
    {
      if(i_pad == 0) HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
      if(i_pad == 1) HistNameEff = Form("h_mEffEta_Cent_%d",i_cent);
      if(i_pad == 2) HistNameEff = Form("h_mEffPhi_Cent_%d",i_cent);

      h_mEfficiency[HistNameEff]->SetLineColor(i_cent+1);
      if(i_cent == 2)
      {
	string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
	h_mEfficiency[HistNameEff]->SetTitle(title.c_str());
	h_mEfficiency[HistNameEff]->SetStats(0);

	if(i_pad == 0) h_mEfficiency[HistNameEff]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	if(i_pad == 1) h_mEfficiency[HistNameEff]->GetXaxis()->SetTitle("#eta");
	if(i_pad == 2) h_mEfficiency[HistNameEff]->GetXaxis()->SetTitle("#phi");
	h_mEfficiency[HistNameEff]->GetXaxis()->CenterTitle();

	h_mEfficiency[HistNameEff]->GetYaxis()->SetTitle("Efficiency");
	h_mEfficiency[HistNameEff]->GetYaxis()->CenterTitle();
	h_mEfficiency[HistNameEff]->GetYaxis()->SetRangeUser(0.0,1.05);
	h_mEfficiency[HistNameEff]->DrawCopy("HIST");
      }
      else
      {
	h_mEfficiency[HistNameEff]->DrawCopy("HIST same");
      }
    }
    if(i_pad == 0) HistNameEff = "h_mEffPt_Cent_9";
    if(i_pad == 1) HistNameEff = "h_mEffEta_Cent_9";
    if(i_pad == 2) HistNameEff = "h_mEffPhi_Cent_9";
    h_mEfficiency[HistNameEff]->SetLineColor(2);
    h_mEfficiency[HistNameEff]->SetLineWidth(2);
    h_mEfficiency[HistNameEff]->DrawCopy("HIST same");

    if(i_pad == 0)
    {
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
    }
  }

  string FigName = Form("/Users/xusun/WorkSpace/Papers/VecMesonSpinAlignment/figures/Efficiency/TPC/c_TpcEffCentCom_%s%s.eps",mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
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
    string HistNameEff = Form("h_mEff_Cent_9_Eta_%d_Phi_%d",eta_bin,i_pad);

    string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
    h_mEfficiency[HistNameEff]->SetTitle(title.c_str());
    h_mEfficiency[HistNameEff]->SetStats(0);

    h_mEfficiency[HistNameEff]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_mEfficiency[HistNameEff]->GetXaxis()->CenterTitle();

    h_mEfficiency[HistNameEff]->GetYaxis()->SetTitle("Efficiency");
    h_mEfficiency[HistNameEff]->GetYaxis()->CenterTitle();
    h_mEfficiency[HistNameEff]->GetYaxis()->SetRangeUser(0.0,1.05);
    h_mEfficiency[HistNameEff]->SetLineColor(2);
    h_mEfficiency[HistNameEff]->SetLineWidth(2);
    // h_mEfficiency[HistNameEff]->SetMarkerStyle(20);
    // h_mEfficiency[HistNameEff]->SetMarkerSize(1.4);
    // h_mEfficiency[HistNameEff]->SetMarkerColor(2);
    h_mEfficiency[HistNameEff]->DrawCopy("HIST");

    for(int i_cent = 2; i_cent <= 5; ++i_cent)
    {
      HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,eta_bin,i_pad);
      h_mEfficiency[HistNameEff]->SetLineColor(i_cent+1);
      h_mEfficiency[HistNameEff]->DrawCopy("HIST same");
    }
    h_mEfficiency[HistNameEff]->DrawCopy("HIST same");
  }
  FigName = Form("/Users/xusun/WorkSpace/Papers/VecMesonSpinAlignment/figures/Efficiency/TPC/c_TpcEffDiffCom_%s%s.eps",mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_diff->SaveAs(FigName.c_str());
}

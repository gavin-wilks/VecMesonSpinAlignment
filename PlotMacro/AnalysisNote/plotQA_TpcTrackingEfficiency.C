#include <string>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLegend.h>
#include <map>
#include "../../anapicodst/StRoot/Utility/StSpinAlignmentCons.h"

typedef std::map<std::string,TH1D*> TH1DMap;

using namespace std;

TH1D* CalEffError(TH1D *h_Mc, TH1D *h_Rc, std::string HistName)
{
  cout << "Before clone" << endl;
  TH1D* h_ratio = (TH1D*)h_Rc->Clone();
  cout << "Cloned RC hist" << endl;
  //h_ratio->Divide(h_Mc);
  h_ratio->Divide(h_Rc,h_Mc,1,1,"B");
  cout << "Divided hists" << endl;
  //for(int i_bin = 1; i_bin < h_ratio->GetNbinsX()+1; ++i_bin)
  //{
  //  double n = h_Mc->GetBinContent(i_bin);
  //  double k = h_Rc->GetBinContent(i_bin);
  //  double variance = (k+1.0)*(k+2.0)/((n+2.0)*(n+3.0))-(k+1.0)*(k+1.0)/((n+2.0)*(n+2.0));
  //  double sigma = TMath::Sqrt(variance);
  //  if(n > 0.0 && k > 0.0) h_ratio->SetBinError(i_bin,sigma);
  //}
  h_ratio->SetName(HistName.c_str());

  return h_ratio;
}



void plotQA_TpcTrackingEfficiency(int mEnergy = 4, int mPID = 0, int year = 1)
{
  string const mBeamEnergy[7] = {"7GeV","9GeV","11GeV","14GeV","19GeV"};
  string const mParType[2] = {"Kplus","Kminus"};
  string const mParTex[2] = {"K^{+}","K^{-}"};
  string const Centrality[10] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%","20%-60%"}; // Centrality bin
  string const mYear[2] = {"run19","run19"};
  int const mPlotStyle[10] = {0,0,1,1,1,1,0,0,0,1};
  int const mPlotColor[10] = {0,0,kGray+3,kAzure+4,kOrange+1,kCyan+1,0,0,0,kRed};

  string inputfile = Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/Embedding/Phi/%s_embedding_%s.root",mBeamEnergy[mEnergy].c_str(),mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  cout << "open input file: " << inputfile.c_str() << endl;
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH3D* h_mMcTracks[10];
  TH3D* h_mRcTracks[10];
  TH1D* h_mMcEffPt[10];
  TH1D* h_mRcEffPt[10];
  TH1D* h_mEffPt[10];
  TH1D* h_mMcEffEta[10];
  TH1D* h_mRcEffEta[10];
  TH1D* h_mEffEta[10];
  TH1D* h_mMcEffPhi[10];
  TH1D* h_mRcEffPhi[10];
  TH1D* h_mEffPhi[10];

  TH1DMap h_mMcEffPEP;
  TH1DMap h_mRcEffPEP;
  TH1DMap h_mEfficiency;

  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    string HistName;
    if(i_cent != 9)
    {
      HistName = Form("h_mMcTracks_%d",i_cent);
      h_mMcTracks[i_cent] = (TH3D*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mRcTracks_%d",i_cent);
      h_mRcTracks[i_cent] = (TH3D*)File_InPut->Get(HistName.c_str());
      cout << "Read in track counts" << endl;
    }

    if(i_cent == 2)
    {
      HistName = "h_mMcTracks_9"; 
      h_mMcTracks[9] = (TH3D*)h_mMcTracks[2]->Clone(HistName.c_str());
      HistName = "h_mRcTracks_9"; 
      h_mRcTracks[9] = (TH3D*)h_mRcTracks[2]->Clone(HistName.c_str());
      cout << "Created 20-60% track counts" << endl;
    }
    if(i_cent > 2 && i_cent <= 5)
    {
      h_mMcTracks[9]->Add(h_mMcTracks[i_cent]);
      h_mRcTracks[9]->Add(h_mRcTracks[i_cent]);
      cout << "Add to 20-60% track counts" << endl;
    }
    //HistName = Form("h_mEffPt_Cent_%d",i_cent);
    //h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str());
    //cout << "read in => " << HistName.c_str() << endl;

    //HistName = Form("h_mEffEta_Cent_%d",i_cent);
    //h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str());
    //cout << "read in => " << HistName.c_str() << endl;

    //HistName = Form("h_mEffPhi_Cent_%d",i_cent);
    //h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str());
    //cout << "read in => " << HistName.c_str() << endl;
    //for(int i_eta = 0; i_eta < 15; ++i_eta)
    //{
    //for(int i_phi = 0; i_phi < 12; ++i_phi)
    //  {
	//HistName = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
	//h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str());
	//cout << "read in => " << HistName.c_str() << endl;
     // }
   // }
    
    //std::string HistName;

    HistName = Form("h_mMcEffPt_Cent_%d",i_cent);
    h_mMcEffPt[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffPt_Cent_%d",i_cent);
    h_mRcEffPt[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("x")->Clone(HistName.c_str());
    HistName = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEffPt[i_cent] = CalEffError(h_mMcEffPt[i_cent],h_mRcEffPt[i_cent],HistName.c_str());

    HistName = Form("h_mMcEffEta_Cent_%d",i_cent);
    h_mMcEffEta[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffEta_Cent_%d",i_cent);
    h_mRcEffEta[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("y")->Clone(HistName.c_str());
    HistName = Form("h_mEffEta_Cent_%d",i_cent);
    h_mEffEta[i_cent] = CalEffError(h_mMcEffEta[i_cent],h_mRcEffEta[i_cent],HistName.c_str());

    HistName = Form("h_mMcEffPhi_Cent_%d",i_cent);
    h_mMcEffPhi[i_cent] = (TH1D*)h_mMcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
    HistName = Form("h_mRcEffPhi_Cent_%d",i_cent);
    h_mRcEffPhi[i_cent] = (TH1D*)h_mRcTracks[i_cent]->Project3D("z")->Clone(HistName.c_str());
    HistName = Form("h_mEffPhi_Cent_%d",i_cent);
    h_mEffPhi[i_cent] = CalEffError(h_mMcEffPhi[i_cent],h_mRcEffPhi[i_cent],HistName.c_str());

    cout << "Calculate Efficiency Errors for integrated bins" << endl;

    for(int i_eta = 0; i_eta < vmsa::BinEta; ++i_eta)
    {
      for(int i_phi = 0; i_phi < vmsa::BinPhi; ++i_phi)
      {
        string HistNameMc = Form("h_mMcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
        h_mMcEffPEP[HistNameMc] = (TH1D*)h_mMcTracks[i_cent]->ProjectionX(HistNameMc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);
        cout << "Project Mc Tracks" << endl;
        string HistNameRc = Form("h_mRcEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
        h_mRcEffPEP[HistNameRc] = (TH1D*)h_mRcTracks[i_cent]->ProjectionX(HistNameRc.c_str(),i_eta+1,i_eta+1,i_phi+1,i_phi+1);
        cout << "Project Rc Tracks" << endl;
        string HistNameEff = Form("h_mEff_Cent_%d_Eta_%d_Phi_%d",i_cent,i_eta,i_phi);
        cout << HistNameEff << endl;
        h_mEfficiency[HistNameEff] = (TH1D*) CalEffError(h_mMcEffPEP[HistNameMc],h_mRcEffPEP[HistNameRc],HistNameEff.c_str());
        cout << "Calculate Efficiency errors for individual phi and eta bins" << endl;
      }
    }    
  }
  {
  TCanvas *c_raw = new TCanvas("c_raw","c_raw",10,10,800,800);
  c_raw->SetLeftMargin(0.15);
  c_raw->SetBottomMargin(0.15);
  c_raw->SetGrid(0,0);
  c_raw->SetTicks(1,1);

  string HistNameEff = "h_mRcEffPt_Cent_9";
  string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  h_mRcEffPt[9]->SetTitle(title.c_str());
  h_mRcEffPt[9]->SetStats(0);

  h_mRcEffPt[9]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mRcEffPt[9]->GetXaxis()->CenterTitle();
  h_mRcEffPt[9]->GetXaxis()->SetRangeUser(0.0,5.00);

  h_mRcEffPt[9]->GetYaxis()->SetTitle("Efficiency");
  h_mRcEffPt[9]->GetYaxis()->CenterTitle();
  //h_mRcEffPt[9]->GetYaxis()->SetRangeUser(0.0,1.05);
  h_mRcEffPt[9]->SetLineColor(mPlotColor[9]);
  h_mRcEffPt[9]->SetLineWidth(4);
  h_mRcEffPt[9]->DrawCopy("HIST");

  for(int i_cent = 2; i_cent <= 5; ++i_cent)
  {
    HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
    h_mRcEffPt[i_cent]->SetLineColor(mPlotColor[i_cent]);
    h_mRcEffPt[i_cent]->SetLineWidth(2);
    h_mRcEffPt[i_cent]->DrawCopy("HIST same");
  }
  HistNameEff = "h_mRcEffPt_Cent_9";
  h_mRcEffPt[9]->DrawCopy("HIST same");

  TLegend *leg = new TLegend(0.5,0.2,0.8,0.5);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    if(mPlotStyle[i_cent] > 0)
    {
      HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
      leg->AddEntry(h_mRcEffPt[i_cent],Centrality[i_cent].c_str());
    }
  }
  leg->Draw("same");

  string FigName = Form("figures/AnalysisNote/Efficiency/TPC/c_RcTpcEffRawHits_%s%s.pdf",mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_raw->SaveAs(FigName.c_str());
  }
  {
  TCanvas *c_raw = new TCanvas("c_raw","c_raw",10,10,800,800);
  c_raw->SetLeftMargin(0.15);
  c_raw->SetBottomMargin(0.15);
  c_raw->SetGrid(0,0);
  c_raw->SetTicks(1,1);

  string HistNameEff = "h_mMcEffPt_Cent_9";
  string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  h_mMcEffPt[9]->SetTitle(title.c_str());
  h_mMcEffPt[9]->SetStats(0);

  h_mMcEffPt[9]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mMcEffPt[9]->GetXaxis()->CenterTitle();
  h_mMcEffPt[9]->GetXaxis()->SetRangeUser(0.0,5.00);

  h_mMcEffPt[9]->GetYaxis()->SetTitle("Efficiency");
  h_mMcEffPt[9]->GetYaxis()->CenterTitle();
  //h_mMcEffPt[9]->GetYaxis()->SetRangeUser(0.0,1.05);
  h_mMcEffPt[9]->SetLineColor(mPlotColor[9]);
  h_mMcEffPt[9]->SetLineWidth(4);
  h_mMcEffPt[9]->DrawCopy("HIST");

  for(int i_cent = 2; i_cent <= 5; ++i_cent)
  {
    HistNameEff = Form("h_mMcEffPt_Cent_%d",i_cent);
    h_mMcEffPt[i_cent]->SetLineColor(mPlotColor[i_cent]);
    h_mMcEffPt[i_cent]->SetLineWidth(2);
    h_mMcEffPt[i_cent]->DrawCopy("HIST same");
  }
  HistNameEff = "h_mMcEffPt_Cent_9";
  h_mMcEffPt[9]->DrawCopy("HIST same");

  TLegend *leg = new TLegend(0.5,0.2,0.8,0.5);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    if(mPlotStyle[i_cent] > 0)
    {
      HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
      leg->AddEntry(h_mMcEffPt[i_cent],Centrality[i_cent].c_str());
    }
  }
  leg->Draw("same");

  string FigName = Form("figures/AnalysisNote/Efficiency/TPC/c_McTpcEffRawHits_%s%s.pdf",mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_raw->SaveAs(FigName.c_str());
  }

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetGrid(0,0);
  c_play->SetTicks(1,1);

  string HistNameEff = "h_mEffPt_Cent_9";
  string title = Form("%s @ Au+Au %s",mParTex[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  h_mEffPt[9]->SetTitle(title.c_str());
  h_mEffPt[9]->SetStats(0);

  h_mEffPt[9]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mEffPt[9]->GetXaxis()->CenterTitle();
  h_mEffPt[9]->GetXaxis()->SetRangeUser(0.0,5.00);

  h_mEffPt[9]->GetYaxis()->SetTitle("Efficiency");
  h_mEffPt[9]->GetYaxis()->CenterTitle();
  h_mEffPt[9]->GetYaxis()->SetRangeUser(0.0,1.05);
  h_mEffPt[9]->SetLineColor(mPlotColor[9]);
  h_mEffPt[9]->SetLineWidth(4);
  h_mEffPt[9]->DrawCopy("HIST");

  for(int i_cent = 2; i_cent <= 5; ++i_cent)
  {
    HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
    h_mEffPt[i_cent]->SetLineColor(mPlotColor[i_cent]);
    h_mEffPt[i_cent]->SetLineWidth(2);
    h_mEffPt[i_cent]->DrawCopy("HIST same");
  }
  HistNameEff = "h_mEffPt_Cent_9";
  h_mEffPt[9]->DrawCopy("HIST same");

  TLegend *leg = new TLegend(0.5,0.2,0.8,0.5);
  leg->SetFillColor(10);
  leg->SetBorderSize(0);
  for(int i_cent = 0; i_cent < 10; ++i_cent)
  {
    if(mPlotStyle[i_cent] > 0)
    {
      HistNameEff = Form("h_mEffPt_Cent_%d",i_cent);
      leg->AddEntry(h_mEffPt[i_cent],Centrality[i_cent].c_str());
    }
  }
  leg->Draw("same");

  string FigName = Form("figures/AnalysisNote/Efficiency/TPC/c_TpcEffCentCom_%s%s.pdf",mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_play->SaveAs(FigName.c_str());

  int const eta_bin = 7;

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
  FigName = Form("figures/AnalysisNote/Efficiency/TPC/c_TpcEffDiffComEta_%d_%s%s.pdf",eta_bin,mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
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
  FigName = Form("figures/AnalysisNote/Efficiency/TPC/c_TpcEffDiffComPhi_%d_%s%s.pdf",phi_bin,mParType[mPID].c_str(),mBeamEnergy[mEnergy].c_str());
  c_diff_phi->SaveAs(FigName.c_str());
}

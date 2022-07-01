#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/draw.h"

using namespace std;

float pos_x[8] = {0.05,0.54,0.05,0.54,0.05,0.54,0.05,0.54};
float pos_y[8] = {0.95,0.95,0.90,0.90,0.85,0.85,0.80,0.80};

void comEffMcPhi(int energy = 4, int cent = 9)
{
  gStyle->SetOptDate(0);

  // string input_first = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/Phi/Efficiency/Eff_%s_SingleKaon_first.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string input_first = Form("/star/u/gwilks3/Workspace/VectorMesonSpinAlignment/VecMesonSpinAlignment/RcPhiEffCorr/Eff_%s_SingleKaon_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_fisrt = TFile::Open(input_first.c_str());
  TH1D *h_mEff_first[vmsa::pt_rebin];
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",cent,i_pt);
    h_mEff_first[i_pt] = (TH1D*)File_fisrt->Get(HistName.c_str());
  }

  // string input_second = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/Phi/Efficiency/Eff_%s_SingleKaon_second.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //string input_second = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Phi/Efficiency/Eff_%s_SingleKaon_second.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //TFile *File_second = TFile::Open(input_second.c_str());
  //TH1D *h_mEff_second[vmsa::pt_rebin];
  //for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  //{
  //  string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",cent,i_pt);
  //  h_mEff_second[i_pt] = (TH1D*)File_second->Get(HistName.c_str());
  //}

  // string input_eta = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/Phi/Efficiency/Eff_%s_SingleKaon_eta.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //string input_eta = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Phi/Efficiency/Eff_%s_SingleKaon_eta.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  //TFile *File_eta = TFile::Open(input_eta.c_str());
  //TH1D *h_mEff_eta[vmsa::pt_rebin];
  //for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  //{
  //  string HistName = Form("h_mEffCos_Cent_%d_Pt_%d",cent,i_pt);
  //  h_mEff_eta[i_pt] = (TH1D*)File_eta->Get(HistName.c_str());
  //}

  TCanvas *c_eff = new TCanvas("c_eff","c_eff",10,10,800,800);
  c_eff->SetLeftMargin(0.15);
  c_eff->SetBottomMargin(0.15);
  c_eff->SetGrid(0,0);
  c_eff->SetTicks(1,1);
  TH1D *h_frame = new TH1D("h_frame","h_frame",100,0.0,1.0);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10);
    h_frame->SetBinError(i_bin+1,1);
  }
  string legEnergy = Form("AuAu %s 20%%-60%%",vmsa::mBeamEnergy[energy].c_str());
  h_frame->SetTitle(legEnergy.c_str());
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetTitle("cos(#theta*)");
  h_frame->GetXaxis()->CenterTitle();
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetNdivisions(505);

  h_frame->GetYaxis()->SetTitle("efficiency");
  h_frame->GetYaxis()->SetTitleSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->SetNdivisions(505);
  h_frame->GetYaxis()->SetTitleOffset(1.2);
  h_frame->GetYaxis()->SetRangeUser(0.0,0.5);
  h_frame->Draw("pE");
  for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt)
  {
    //h_mEff_eta[i_pt]->SetMarkerStyle(29);
    //h_mEff_eta[i_pt]->SetMarkerColor(2);
    //h_mEff_eta[i_pt]->SetMarkerSize(1.4);
    //h_mEff_eta[i_pt]->Draw("pEX0 same");

    h_mEff_first[i_pt]->SetMarkerStyle(24);
    h_mEff_first[i_pt]->SetMarkerColor(6);
    h_mEff_first[i_pt]->SetMarkerSize(1.4);
    h_mEff_first[i_pt]->Draw("pEX0 same");

    //h_mEff_second[i_pt]->SetMarkerStyle(25);
    //h_mEff_second[i_pt]->SetMarkerColor(4);
    //h_mEff_second[i_pt]->SetMarkerSize(1.4);
    //h_mEff_second[i_pt]->Draw("pEX0 same");
  }

  //TLegend *leg = new TLegend(0.6,0.7,0.8,0.9);
  //leg->SetBorderSize(0);
  //leg->SetFillColor(10);
  //leg->AddEntry(h_mEff_eta[0],"eta","p");
  //leg->AddEntry(h_mEff_first[0],"","p");
  //leg->AddEntry(h_mEff_second[0],"phi & plateau","p");
  //leg->Draw("same");

  // c_eff->SaveAs("/star/data01/pwg/sunxuhit/AuAu39GeV/SpinAlignment/Phi/Efficiency/c_eff_com.eps");
  string FigName = Form("figures/AnalysisNote/Efficiency/c_effAuAu%s_com.pdf",vmsa::mBeamEnergy[energy].c_str());
  c_eff->SaveAs(FigName.c_str());
}

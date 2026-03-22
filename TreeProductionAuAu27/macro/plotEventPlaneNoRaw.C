#include "../StRoot/Utility/StSpinAlignmentCons.h"
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"

using namespace std;

void plotEventPlaneNoRaw(int beamEnergy = 5)
{

  //--------------TPC EP---------------
  // x-axis: runIndex | y-axis: EP
  TH1F *h_mTpcRawEpEast; // raw EP
  TH1F *h_mTpcRawEpWest;
  TH1F *h_mTpcRawEpFull;

  TH1F *h_mTpcReCenterEpEast; // recenter EP
  TH1F *h_mTpcReCenterEpWest;
  TH1F *h_mTpcReCenterEpFull;

  TH1F *h_mTpcShiftEpEast; // shift EP
  TH1F *h_mTpcShiftEpWest;
  TH1F *h_mTpcShiftEpRanA;
  TH1F *h_mTpcShiftEpRanB;
  TH1F *h_mTpcShiftEpFull;
  //--------------TPC EP---------------

  string inputRaw = Form("../StRoot/Utility/ReCenterParameter/file_%s_ReCenterPar.root",vmsa::mBeamEnergy[beamEnergy].c_str());
  TFile *File_InPutRaw = TFile::Open(inputRaw.c_str());
  
  string HistName;

  
 cout << "Are these loaded in" << endl;

  string inputShift = Form("../StRoot/Utility/Resolution/file_%s_Resolution.root",vmsa::mBeamEnergy[beamEnergy].c_str());
  TFile *File_InPutShift = TFile::Open(inputShift.c_str());

  HistName = "h_mEastRaw";
  h_mTpcRawEpEast = (TH1F*)File_InPutShift->Get(HistName.c_str());
  HistName = "h_mWestRaw";
  h_mTpcRawEpWest = (TH1F*)File_InPutShift->Get(HistName.c_str());
  HistName = "h_mFullRaw";
  h_mTpcRawEpFull = (TH1F*)File_InPutShift->Get(HistName.c_str());  
  HistName = "h_mEastReCenter";
  h_mTpcReCenterEpEast = (TH1F*)File_InPutShift->Get(HistName.c_str());
  HistName = "h_mWestReCenter";
  h_mTpcReCenterEpWest = (TH1F*)File_InPutShift->Get(HistName.c_str());
  HistName = "h_mFullReCenter";
  h_mTpcReCenterEpFull = (TH1F*)File_InPutShift->Get(HistName.c_str());
  HistName = "h_mEastShift";
  h_mTpcShiftEpEast = (TH1F*)File_InPutShift->Get(HistName.c_str());
  HistName = "h_mWestShift";
  h_mTpcShiftEpWest = (TH1F*)File_InPutShift->Get(HistName.c_str());
  HistName = "h_mFullShift";
  h_mTpcShiftEpFull = (TH1F*)File_InPutShift->Get(HistName.c_str());

  cout << "Are these loaded in 2" << endl;

  TCanvas *c_Tpc = new TCanvas("c_Tpc","c_Tpc",10,10,800,400);
  c_Tpc->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_Tpc->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Tpc->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Tpc->cd(i_pad+1)->SetGrid(0,0);
    c_Tpc->cd(i_pad+1)->SetTicks(1,1);
  }

  TLegend *legE;
  TLegend *legW;
 
  gStyle->SetOptStat("n");

  c_Tpc->cd(1);
  //h_mTpcRawEpEast->SetTitle("2^{nd} Order East EP TPC");
  //h_mTpcRawEpEast->GetXaxis()->SetTitle("#Psi_{2,E}");
  //h_mTpcRawEpEast->GetYaxis()->SetTitle("Count");
  //h_mTpcRawEpEast->GetXaxis()->SetRangeUser(-TMath::Pi()/2.,TMath::Pi()/2.);
  //h_mTpcRawEpEast->GetYaxis()->SetRangeUser(950*1e3,1500*1e3);
  //h_mTpcRawEpEast->SetLineColor(kBlack);
  //h_mTpcRawEpEast->SetMarkerColor(kBlack);
  h_mTpcReCenterEpEast->SetTitle("2^{nd} Order East EP TPC");
  h_mTpcReCenterEpEast->GetXaxis()->SetTitle("#Psi_{2,E}");
  h_mTpcReCenterEpEast->GetYaxis()->SetTitle("Count");
  h_mTpcReCenterEpEast->GetXaxis()->SetRangeUser(-TMath::Pi()/2.,TMath::Pi()/2.);
  h_mTpcReCenterEpEast->SetLineColor(kBlack);
  h_mTpcReCenterEpEast->SetMarkerColor(kBlack);
  h_mTpcReCenterEpEast->SetLineColor(kBlue);
  h_mTpcReCenterEpEast->SetMarkerColor(kBlue);
  h_mTpcShiftEpEast->SetLineColor(kRed);
  h_mTpcShiftEpEast->SetMarkerColor(kRed);

  //h_mTpcRawEpEast->Draw("pE");
  //h_mTpcReCenterEpEast->Draw("pE same");
  h_mTpcReCenterEpEast->Draw("pE");
  h_mTpcShiftEpEast->Draw("pE same");
  h_mTpcRawEpEast->Print();
  h_mTpcReCenterEpEast->Print();
  //h_mTpcReCenterEpEast->Draw("pE");
  h_mTpcShiftEpEast->Print();

  legE = new TLegend(0.65,0.25,0.85,0.45);
  legE->SetFillColor(10);
  legE->SetBorderSize(0);
  //legE->AddEntry(h_mTpcRawEpEast,"Raw EP","l");
  legE->AddEntry(h_mTpcReCenterEpEast,"ReCenter EP","l");
  legE->AddEntry(h_mTpcShiftEpEast,"Shift EP","l"); 
  legE->Draw("same");


  c_Tpc->cd(2);
  h_mTpcReCenterEpWest->SetTitle("2^{nd} Order West EP TPC");
  h_mTpcReCenterEpWest->GetXaxis()->SetTitle("#Psi_{2,W}");
  h_mTpcReCenterEpWest->GetYaxis()->SetTitle("Count");
  h_mTpcReCenterEpWest->GetXaxis()->SetRangeUser(-TMath::Pi()/2.,TMath::Pi()/2.);
  h_mTpcReCenterEpWest->SetLineColor(kBlack);
  h_mTpcReCenterEpWest->SetMarkerColor(kBlack);
  //h_mTpcRawEpWest->SetTitle("2^{nd} Order West EP TPC");
  //h_mTpcRawEpWest->GetXaxis()->SetTitle("#Psi_{2,W}");
  //h_mTpcRawEpWest->GetYaxis()->SetTitle("Count");
  //h_mTpcRawEpWest->GetXaxis()->SetRangeUser(-TMath::Pi()/2.,TMath::Pi()/2.);
  //h_mTpcRawEpWest->GetYaxis()->SetRangeUser(950*1e3,1500*1e3);
  //h_mTpcRawEpWest->SetLineColor(kBlack);
  //h_mTpcRawEpWest->SetMarkerColor(kBlack);
  h_mTpcReCenterEpWest->SetLineColor(kBlue);
  h_mTpcReCenterEpWest->SetMarkerColor(kBlue);
  h_mTpcShiftEpWest->SetLineColor(kRed);
  h_mTpcShiftEpWest->SetMarkerColor(kRed);

  //h_mTpcRawEpWest->Draw("pE");
  h_mTpcReCenterEpWest->Draw("pE");
  //h_mTpcReCenterEpWest->Draw("pE same");
  h_mTpcShiftEpWest->Draw("pE same");
  h_mTpcRawEpWest->Print();
  h_mTpcReCenterEpWest->Print();
  h_mTpcShiftEpWest->Print();

  legW = new TLegend(0.65,0.25,0.85,0.45);
  legW->SetFillColor(10);
  legW->SetBorderSize(0);
  //legW->AddEntry(h_mTpcRawEpWest,"Raw EP","l");
  legW->AddEntry(h_mTpcReCenterEpWest,"ReCenter EP","l");
  legW->AddEntry(h_mTpcShiftEpWest,"Shift EP","l"); 
  legW->Draw("same");

  string FigureName = Form("./figures/c_mTpcEventPlane_%s.pdf",vmsa::mBeamEnergy[beamEnergy].c_str());
  c_Tpc->SaveAs(FigureName.c_str());
}

#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include "../../Utility/draw.h"

void plotQA_rawYields(int energy = 3)
{
  string mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/Phi/rho00/eta_1/RawYieldSys.root",mBeamEnergy[energy].c_str());
  cout << "inputfile is set to " << inputfile.c_str() << endl;

  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_rawYields = (TH1F*)File_InPut->Get("Yield_Centrality_9_Dca_0_Sig_0_Phi_Norm_0_Sigma_0_Count");

  TCanvas *c_rawYields = new TCanvas("c_rawYields","c_rawYields",10,10,800,800);
  c_rawYields->cd();
  c_rawYields->cd()->SetLeftMargin(0.15);
  c_rawYields->cd()->SetBottomMargin(0.15);
  c_rawYields->cd()->SetTicks(1,1);
  c_rawYields->cd()->SetGrid(0,0);
  c_rawYields->cd()->SetLogy(1);

  h_rawYields->SetStats(0);
  h_rawYields->SetTitle("raw Yields at 27GeV 20-60% & |v_{z}| < 70");
  h_rawYields->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_rawYields->GetXaxis()->CenterTitle();
  h_rawYields->GetYaxis()->SetTitle("Yields");
  h_rawYields->GetYaxis()->CenterTitle();
  h_rawYields->SetMarkerSize(2.0);
  h_rawYields->SetMarkerColor(kGray+2);
  h_rawYields->Draw("pE");

  string FigName = Form("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/AnalysisNote/RhoSig/c_rawYields_AuAu%s.eps",mBeamEnergy[energy].c_str());
  c_rawYields->SaveAs(FigName.c_str());
}

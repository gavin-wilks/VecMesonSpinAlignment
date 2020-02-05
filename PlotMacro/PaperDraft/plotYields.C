#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/functions.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"

using namespace std;

void plotYields()
{
  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/BESII/counts.root");
  TH1D *h_mYields = (TH1D*)File_Input->Get("PtCos")->Clone("h_mYields");

  TCanvas *c_yields = new TCanvas("c_yields","c_yields",10,10,800,800);
  c_yields->cd();
  c_yields->cd()->SetLeftMargin(0.15);
  c_yields->cd()->SetBottomMargin(0.15);
  c_yields->cd()->SetTicks(1,1);
  c_yields->cd()->SetGrid(0,0);

  h_mYields->SetTitle("");
  h_mYields->SetStats(0);
  h_mYields->GetXaxis()->SetTitle("cos(#theta*)");
  h_mYields->GetXaxis()->CenterTitle();
  h_mYields->GetXaxis()->SetTitleSize(0.06);
  h_mYields->GetXaxis()->SetTitleOffset(0.9);
  h_mYields->GetXaxis()->SetLabelSize(0.04);
  h_mYields->SetNdivisions(505,"X");

  h_mYields->GetYaxis()->SetTitle("Yields");
  h_mYields->GetYaxis()->CenterTitle();
  h_mYields->GetYaxis()->SetTitleOffset(1.14);
  h_mYields->GetYaxis()->SetTitleSize(0.06);
  h_mYields->GetYaxis()->SetLabelSize(0.04);
  h_mYields->GetYaxis()->SetRangeUser(0.99*h_mYields->GetMinimum(),1.01*h_mYields->GetMaximum());
  h_mYields->SetNdivisions(505,"Y");
  h_mYields->SetMarkerStyle(24);
  h_mYields->SetMarkerColor(1);
  h_mYields->SetMarkerSize(1.8);
  h_mYields->DrawCopy("PE");

  TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
  f_rho->SetParameter(0,0.33);
  f_rho->SetParameter(1,1000);
  h_mYields->Fit(f_rho,"NQMRI");
  f_rho->SetLineColor(4);
  f_rho->SetLineStyle(2);
  f_rho->SetLineWidth(4);
  f_rho->DrawCopy("l same");

  string leg_energy = Form("Au+Au %s & 20-60%%",vmsa::mBeamEnergy[6].c_str());
  plotTopLegend((char*)leg_energy.c_str(),0.2,0.75,0.04,1,0.0,42,1);
  plotTopLegend((char*)"1.2 < p_{T} < 1.8 GeV/c",0.25,0.70,0.04,1,0.0,42,1);

  c_yields->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperDraft/c_yields.eps");
  c_yields->SaveAs("/Users/xusun/WorkSpace/STAR/figures/SpinAlignment/PaperProposal/c_yields.png");
}

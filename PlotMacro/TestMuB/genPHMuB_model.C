#include <string>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TStyle.h>
#include <TF1.h>

void genPHMuB_model()
{
  const double amptEnergy[6] = {3.0, 7.7,11.5,14.5,19.6,27.0};
  const double amptMuB[6]    = {750, 420, 315, 260, 205, 155};
  const double amptUpper[6]  = {1.304,1.499,1.374,1.262,1.140,0.985};
  const double amptLower[6]  = {1.150,1.486,1.364,1.245,1.126,0.971};

  const double urqmdEnergy[5] = {7.7,19.6,39.0,62.4,200.0};
  const double urqmdMuB[5]    = {420, 205, 115,  70,   20};
  const double urqmdMean[5]   = {0.01752,0.00987,0.00584,0.00400,0.00172};

  const double ckEnergy[9] = {7.7,11.5,14.5,19.6,27.0,39.0,62.4,200.0};
  const double ckMuB[9]    = {420, 315, 260, 205, 155, 115,  70,   20};
  const double ckMean[9]   = {2.43925,2.30841,2.23364,2.03738,1.91589,1.71028,1.47664,0.86916};

  const double fdEnergy[12] = {2.4, 2.7, 3.3, 3.8, 4.8, 5.6, 7.7, 11.5, 14.5, 19.6, 27.0, 39.0};
  const double fdMuB[12]    = {790, 753, 688, 642, 566, 517, 420,  315,  260,  205,  155,  115};
  const double fdUpper[12]  = {5.728,7.006,10.201,7.077,4.965,4.503,2.994,1.822,1.592,1.343,1.041,0.669};
  const double fdLower[12]  = {5.142,5.515, 6.509,5.408,3.935,3.207,2.320,1.450,1.166,0.882,0.651,0.278};

  TGraphAsymmErrors *g_ampt = new TGraphAsymmErrors();
  g_ampt->SetName("g_ampt");
  for(int iEnergy = 0; iEnergy < 6; ++iEnergy)
  {
    double amptMean = 0.5*(amptUpper[iEnergy]+amptLower[iEnergy]);
    double amptDiff = 0.5*(amptUpper[iEnergy]-amptLower[iEnergy]);
    g_ampt->SetPoint(iEnergy,amptMuB[iEnergy],amptMean);
    g_ampt->SetPointError(iEnergy,0.0,0.0,amptDiff,amptDiff);
  }

  TGraphAsymmErrors *g_urqmd = new TGraphAsymmErrors();
  g_urqmd->SetName("g_urqmd");
  for(int iEnergy = 0; iEnergy < 5; ++iEnergy)
  {
    g_urqmd->SetPoint(iEnergy,urqmdMuB[iEnergy],urqmdMean[iEnergy]*100.0);
  }

  TGraphAsymmErrors *g_ck = new TGraphAsymmErrors();
  g_ck->SetName("g_ck");
  for(int iEnergy = 0; iEnergy < 9; ++iEnergy)
  {
    g_ck->SetPoint(iEnergy,ckMuB[iEnergy],ckMean[iEnergy]);
  }

  TGraphAsymmErrors *g_fd = new TGraphAsymmErrors();
  g_fd->SetName("g_fd");
  for(int iEnergy = 0; iEnergy < 12; ++iEnergy)
  {
    double fdMean = 0.5*(fdUpper[iEnergy]+fdLower[iEnergy]);
    double fdDiff = 0.5*(fdUpper[iEnergy]-fdLower[iEnergy]);
    g_fd->SetPoint(iEnergy,fdMuB[iEnergy],fdMean);
    g_fd->SetPointError(iEnergy,0.0,0.0,fdDiff,fdDiff);
  }

  TH1F *h_frame = new TH1F("h_frame","h_frame",1000,0,1000);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,850.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetTitle("Baryon Chemical Potential #mu_{B} (MeV)");
  h_frame->GetXaxis()->SetTitleSize(0.08);
  h_frame->GetXaxis()->SetTitleOffset(0.9);
  h_frame->GetXaxis()->CenterTitle();
  h_frame->GetXaxis()->SetLabelSize(0.06);

  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitleSize(0.08);
  h_frame->GetYaxis()->SetTitleOffset(0.8);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->GetYaxis()->SetLabelSize(0.06);
  h_frame->GetYaxis()->SetRangeUser(-0.5,12.1);
  h_frame->GetYaxis()->SetTitle("Spin Polarization P_{H} (%)");
  h_frame->DrawCopy("PE");

  // g_ampt->SetMarkerColor(kBlue);
  // g_ampt->SetLineColor(kBlue);
  // g_ampt->SetFillColor(kBlue);
  // g_ampt->SetFillStyle(3001);

  g_fd->SetMarkerColor(kYellow);
  g_fd->SetLineColor(kYellow);
  g_fd->SetFillColor(kYellow);
  g_fd->SetFillStyle(3001);
  g_fd->Draw("pE3 same");

  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/TestMuB/Lambda/PHMuB_model.root","RECREATE");
  File_OutPut->cd();
  g_ampt->Write();
  g_urqmd->Write();
  g_ck->Write();
  g_fd->Write();
  File_OutPut->Close();
}

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
// #include "../../Utility/functions.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/type.h"
#include "../StRoot/StToFMatchMaker/StToFMatchCons.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

TH1D* hTmp;

double tof_Kaon(double* x, double* par)
{
   if ((x[0] > 0.4 && x[0] < 0.9) || (x[0] > 1.1 && x[0] < 2.5))
   {
      int bin = hTmp->FindBin(x[0]);
      return hTmp->GetBinContent(bin);
   }
   else
   {
      return par[0] / (pow(x[0] - par[1], 2) + par[2]) - par[4] / (exp(x[0] - par[3]) + par[5]) + par[6];
   }

}


void fitToFMatchEfficiency(int energy = 6)
{
  // string inputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ToFMatch/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  string inputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/ToFMatch/Eff_%s_ToFMatch.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  const int Cent_Start = 9;
  const int Cent_Stop = 10;

  TH1DMap h_mEfficiency_cent; // ToF matching efficiency
  TH1DMap h_mEfficiency; // ToF matching efficiency
  for(int i_pid = tof::mPID_Start; i_pid < tof::mPID_Stop; ++i_pid)
  {
    for(int i_charge = 0; i_charge < 2; ++i_charge)
    {
      for(int i_cent = Cent_Start; i_cent < Cent_Stop; ++i_cent)
      {
	string HistName = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent);
	h_mEfficiency_cent[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
	for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
	{
	  for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
	  {
	    string HistName = Form("h_mEfficiency_%s%s_Cent_%d_Eta_%d_Phi_%d",tof::mPID_ToF[i_pid].c_str(),tof::mCharge[i_charge].c_str(),i_cent,i_eta,i_phi);
	    h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
	  }
	}
      }
    }
  }

#if _PlotQA_
  int pidQA = tof::pidQA;
  int chargeQA = tof::chargeQA;
  int centQA = 9;
  int etaQA = tof::etaQA;
  int phiQA = tof::phiQA;

  TCanvas *c_Efficiency = new TCanvas("c_Efficiency","c_Efficiency",10,10,800,800);
  c_Efficiency->cd();
  c_Efficiency->cd()->SetLeftMargin(0.15);
  c_Efficiency->cd()->SetBottomMargin(0.15);
  c_Efficiency->cd()->SetTicks(1,1);
  c_Efficiency->cd()->SetGrid(0,0);

  string HistEff = Form("h_mEfficiency_%s%s_Cent_%d",tof::mPID_ToF[pidQA].c_str(),tof::mCharge[chargeQA].c_str(),centQA);
  // h_mEfficiency[HistEff]->SetTitle("TH1::Divide(\"B\")");
  h_mEfficiency_cent[HistEff]->SetMarkerStyle(24);
  h_mEfficiency_cent[HistEff]->SetMarkerSize(1.2);
  h_mEfficiency_cent[HistEff]->SetMarkerColor(2);
  h_mEfficiency_cent[HistEff]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mEfficiency_cent[HistEff]->GetXaxis()->CenterTitle();
  h_mEfficiency_cent[HistEff]->GetYaxis()->SetTitle("Efficiency");
  h_mEfficiency_cent[HistEff]->GetYaxis()->SetRangeUser(0.0,1.2);
  h_mEfficiency_cent[HistEff]->Draw("pE");

  // fit parameters

  double par0[9] = {-0.018979 , -0.0322611 , -0.0680754 , -0.0698575 , -0.0315267 , -0.00589929 , -0.00226724 , -0.00212137 , -0.00389514};
  double par1[9] = {0.0308943 , -0.0939411 , -0.14377 , -0.19003 , -0.116323 , 0.180593 , 0.207874 , 0.208863 , 0.194876};
  double par2[9] = {-0.00592033 , -0.0600635 , -0.0515391 , -0.0708703 , -0.0756912 , 0.00912449 , 0.00500487 , 0.00497987 , 0.00824164 };
  double par3[9] = {1.28883, 1.53952 , 1.52213 , 1.01707 , 1.5415 , 2.75657 , -0.326349 , 1.11207 , 1.37717 };
  double par4[9] = {1.58923e-06 , -0.00130657 , -0.00973403 , -0.0163526 , -0.00162583 , -2.20034e-05 , 0.773984 , 0.119933 , -1.14531e-06 };
  double par5[9] = {-0.340008, -0.261115 , -0.236246 , -0.345869 , -0.260416 , -0.0777638 , 6.00519 , 1.05048 , -0.311133 };
  double par6[9] = {0.685642, 0.687831 , 0.682733 , 0.674683 , 0.659559 , 0.639223 , 0.658294 , 0.626411 , 0.578902};

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(30000);
  hTmp = h_mEfficiency_cent[HistEff];
  TF1 *f_km = new TF1("f_km",tof_Kaon,0.2,10,7);
  f_km->SetParameters(par0[3], par1[3], par2[3], par3[3], par4[3], par5[3], par6[3]);
  h_mEfficiency_cent[HistEff]->Draw();
  h_mEfficiency_cent[HistEff]->Fit(f_km, "NR", "", 0.2, 3);
  f_km->SetLineColor(2);
  f_km->SetLineWidth(2);
  f_km->SetLineStyle(1);
  f_km->SetNpx(1000);
  f_km->Draw("l same");

  TF1 *funkm = new TF1("funkm", "[0]/(pow(x-[1],2)+[2])-[4]/(exp(x-[3])+[5])+[6]", 0.2, 10);
  for (int ii = 0; ii < 7; ii++)
  {
    funkm->SetParameter(ii, f_km->GetParameter(ii));
  }
  funkm->SetNpx(1000);
  funkm->SetLineColor(2);
  funkm->SetLineWidth(3);
  funkm->SetLineStyle(2);
  funkm->Draw("same");

#endif 

}

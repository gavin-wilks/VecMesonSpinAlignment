#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "../../Utility/type.h"
#include "../StRoot/StToFMatchMaker/StToFMatchCons.h"

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

// TH1D* hTmp;

/*
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
*/
double tof_Kaon(double* x, double* par)
{
  return par[0]*(1.0 / (pow(x[0] - par[1], 2) + par[2]) - par[4] / (exp(x[0] - par[3]) + par[5]) + par[6]);
}

void plotQA_ToFMatchEff_2060(const int energy = 6, const int charge = 0)
{
  int const mColor[12] = {11,11,2,2,4,4,7,7,8,8,38,38};
  int const mStyle[12] = {30,24,30,24,30,24,30,24,30,24,30,24};
  int const mLineStyle[12] = {3,2,3,2,3,2,3,2,3,2,3,2};
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);

  std::string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/Eff_%s_ToFMatch_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH1DMap h_mEffCent;
  TH1DMap h_mEffCent_plot;
  TH1DMap h_mEffEta;
  TH1DMap h_mEffEta_plot;
  TH1DMap h_mEfficiency;
  TH1DMap h_mEfficiency_plot;

  string HistName = Form("h_mEfficiency_K%s_Cent_9",tof::mCharge[charge].c_str());
  h_mEffCent[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
  h_mEffCent[HistName]->SetMarkerStyle(24);
  h_mEffCent[HistName]->SetMarkerSize(1.2);
  h_mEffCent[HistName]->SetMarkerColor(2);
  h_mEffCent[HistName]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mEffCent[HistName]->GetXaxis()->CenterTitle();
  h_mEffCent[HistName]->GetYaxis()->SetTitle("Efficiency");
  h_mEffCent[HistName]->GetYaxis()->SetRangeUser(0.0,1.2);
  h_mEffCent_plot[HistName] = (TH1D*)h_mEffCent[HistName]->Clone();
  for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
  {
    HistName = Form("h_mEfficiency_K%s_Cent_9_Eta_%d",tof::mCharge[charge].c_str(),i_eta);
    h_mEffEta[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
    h_mEffEta[HistName]->SetMarkerStyle(29);
    h_mEffEta[HistName]->SetMarkerSize(1.6);
    h_mEffEta[HistName]->SetMarkerColor(2);
    h_mEffEta[HistName]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_mEffEta[HistName]->GetXaxis()->CenterTitle();
    h_mEffEta[HistName]->GetYaxis()->SetTitle("Efficiency");
    h_mEffEta[HistName]->GetYaxis()->SetRangeUser(0.0,1.2);
    h_mEffEta_plot[HistName] = (TH1D*)h_mEffEta[HistName]->Clone();
    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
    {
      HistName = Form("h_mEfficiency_K%s_Cent_9_Eta_%d_Phi_%d",tof::mCharge[charge].c_str(),i_eta,i_phi);
      h_mEfficiency[HistName] = (TH1D*)File_InPut->Get(HistName.c_str())->Clone();
      h_mEfficiency[HistName]->SetMarkerStyle(mStyle[i_phi]);
      h_mEfficiency[HistName]->SetMarkerSize(1.2);
      h_mEfficiency[HistName]->SetMarkerColor(mColor[i_phi]);
      h_mEfficiency[HistName]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      h_mEfficiency[HistName]->GetXaxis()->CenterTitle();
      h_mEfficiency[HistName]->GetYaxis()->SetTitle("Efficiency");
      h_mEfficiency[HistName]->GetYaxis()->SetRangeUser(0.0,1.2);
      h_mEfficiency_plot[HistName] = (TH1D*)h_mEfficiency[HistName]->Clone();
    }
  }

  TH1DMap h_mFitPar; // initialize histogram for fit parameters
  std::string HistName_FitPar;
  HistName_FitPar = Form("h_mFitParameters_K%s_Cent_9",tof::mCharge[charge].c_str());
  h_mFitPar[HistName_FitPar] = new TH1D(HistName_FitPar.c_str(),HistName_FitPar.c_str(),7,-0.5,6.5);
  for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
  {
    HistName_FitPar = Form("h_mFitParameters_K%s_Cent_9_Eta_%d",tof::mCharge[charge].c_str(),i_eta);
    h_mFitPar[HistName_FitPar] = new TH1D(HistName_FitPar.c_str(),HistName_FitPar.c_str(),7,-0.5,6.5);
    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi)
    {
      HistName_FitPar = Form("h_mFitParameters_K%s_Cent_9_Eta_%d_Phi_%d",tof::mCharge[charge].c_str(),i_eta,i_phi);
      h_mFitPar[HistName_FitPar] = new TH1D(HistName_FitPar.c_str(),HistName_FitPar.c_str(),7,-0.5,6.5);
    }
  }

  float dip_1st[7] = {0.4,0.4,0.4,0.4,0.4,0.4,0.4};
  float dip_2nd[7] = {0.8,0.8,0.8,0.8,0.8,0.8,0.8};
  float dip_3rd[7] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  float dip_4th[7] = {2.5,2.5,2.5,2.5,2.5,2.5,2.5}; // Kplus
  // float dip_4th[7] = {2.5,2.5,2.5,2.3,2.5,2.5,2.5}; // Kminus

  double range[7] = {4.0,4.0,5.0,4.0,4.7,4.8,5.0}; // Kplus
  // double range[7] = {4.0,4.0,4.5,4.0,4.7,4.8,5.0}; // Kminus

  int cent_kaon = -1;
  if(charge == 0) cent_kaon = 2;
  if(charge == 1) cent_kaon = 2;
  double par0[9] = {-0.018979, -0.0322611, -0.0680754, -0.0698575, -0.0315267, -0.00589929, -0.00226724, -0.00212137, -0.00389514 };
  double par1[9] = {0.0308943, -0.0939411, -0.14377, -0.19003, -0.116323, 0.180593, 0.207874, 0.208863, 0.194876};
  double par2[9] = {-0.00592033, -0.0600635, -0.0515391, -0.0708703, -0.0756912, 0.00912449, 0.00500487, 0.00497987, 0.00824164};
  double par3[9] = {1.28883, 1.53952, 1.52213, 1.01707, 1.5415, 2.75657, -0.326349, 1.11207, 1.37717};
  double par4[9] = {1.58923e-06, -0.00130657, -0.00973403, -0.0163526, -0.00162583, -2.20034e-05, 0.773984, 0.119933, -1.14531e-06};
  double par5[9] = {-0.340008, -0.261115, -0.236246, -0.345869, -0.260416, -0.0777638, 6.00519, 1.05048, -0.311133};
  double par6[9] = {0.685642, 0.687831, 0.682733, 0.674683, 0.659559, 0.639223, 0.658294, 0.626411, 0.578902};

  HistName = Form("h_mEfficiency_K%s_Cent_9",tof::mCharge[charge].c_str()); // ToF efficiency at cent 9
  for(int i_bin = 0; i_bin < h_mEffCent[HistName]->GetNbinsX(); ++i_bin) 
  {
    float pt = h_mEffCent[HistName]->GetBinCenter(i_bin+1); 
    if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy])) h_mEffCent[HistName]->SetBinError(i_bin+1,0.0); 
    if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy])) h_mEffCent[HistName]->SetBinContent(i_bin+1,0.0); 
  }

  TF1 *f_kaon_cent = new TF1("f_kaon_cent",tof_Kaon,0.2,10,7);
  f_kaon_cent->SetParameter(0,par0[cent_kaon]);
  f_kaon_cent->SetParameter(1,par1[cent_kaon]);
  f_kaon_cent->SetParameter(2,par2[cent_kaon]);
  f_kaon_cent->SetParameter(3,par3[cent_kaon]);
  f_kaon_cent->SetParameter(4,par4[cent_kaon]);
  f_kaon_cent->SetParameter(5,par5[cent_kaon]);
  f_kaon_cent->SetParameter(6,par6[cent_kaon]);
  for(int i_par = 0; i_par < 7; ++i_par)
  {
    f_kaon_cent->SetParLimits(i_par,-10.0,10.0);
  }
  f_kaon_cent->SetRange(0.2,range[energy]);

  h_mEffCent[HistName]->Fit(f_kaon_cent, "MWNR");

  double par_fit[7];
  HistName_FitPar = Form("h_mFitParameters_K%s_Cent_9",tof::mCharge[charge].c_str());
  for(int i_par = 0; i_par < 7; ++i_par)
  {
    par_fit[i_par] = f_kaon_cent->GetParameter(i_par);
    h_mFitPar[HistName_FitPar]->SetBinContent(i_par+1,f_kaon_cent->GetParameter(i_par));
  }

  // fit eta distribution
  TF1 *f_kaon_eta[tof::BinEta];
  for(int i_eta = 1; i_eta < tof::BinEta-1; ++i_eta) // ToF efficiency at different eta and phi bin at cent 9
  {
    HistName = Form("h_mEfficiency_K%s_Cent_9_Eta_%d",tof::mCharge[charge].c_str(),i_eta);
    for(int i_bin = 0; i_bin < h_mEffEta[HistName]->GetNbinsX(); ++i_bin) 
    {
      float pt = h_mEffEta[HistName]->GetBinCenter(i_bin+1); 
      if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]+0.5)) h_mEffEta[HistName]->SetBinError(i_bin+1,0.0); 
      if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]+0.5)) h_mEffEta[HistName]->SetBinContent(i_bin+1,0.0); 
    }

    std::string FuncName = Form("f_kaon_cent_9_eta_%d",i_eta);
    f_kaon_eta[i_eta] = new TF1(FuncName.c_str(),tof_Kaon,0.2,10,7);
    f_kaon_eta[i_eta]->SetParameter(0,par_fit[0]);
    for(int i_par = 1; i_par < 7; ++i_par)
    {
      f_kaon_eta[i_eta]->FixParameter(i_par,par_fit[i_par]);
    }
    f_kaon_eta[i_eta]->SetRange(0.2,range[energy]);

    h_mEffEta[HistName]->Fit(f_kaon_eta[i_eta], "NR");

    HistName_FitPar = Form("h_mFitParameters_K%s_Cent_9_Eta_%d",tof::mCharge[charge].c_str(),i_eta);
    for(int i_bin = 0; i_bin < 7; ++i_bin)
    {
      h_mFitPar[HistName_FitPar]->SetBinContent(i_bin+1,f_kaon_eta[i_eta]->GetParameter(i_bin));
    }
  }

  // fit eta & phi distribution
  TF1 *f_kaon_phi[tof::BinEta][tof::BinPhi];
  for(int i_eta = 1; i_eta < tof::BinEta-1; ++i_eta) // ToF efficiency at different eta and phi bin at cent 9
  {
    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi) // eta & phi dependence
    {
      HistName = Form("h_mEfficiency_K%s_Cent_9_Eta_%d_Phi_%d",tof::mCharge[charge].c_str(),i_eta,i_phi);
      for(int i_bin = 0; i_bin < h_mEfficiency[HistName]->GetNbinsX(); ++i_bin) 
      {
	float pt = h_mEfficiency[HistName]->GetBinCenter(i_bin+1); 
	if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]+0.5)) h_mEfficiency[HistName]->SetBinError(i_bin+1,0.0); 
	if ((pt > dip_1st[energy] && pt < dip_2nd[energy]) || (pt > dip_3rd[energy] && pt < dip_4th[energy]+0.5)) h_mEfficiency[HistName]->SetBinContent(i_bin+1,0.0); 
      }

      std::string FuncName = Form("f_kaon_cent_9_eta_%d_phi_%d",i_eta,i_phi);
      f_kaon_phi[i_eta][i_phi] = new TF1(FuncName.c_str(),tof_Kaon,0.2,10,7);
      f_kaon_phi[i_eta][i_phi]->SetParameter(0,par_fit[0]);
      for(int i_par = 1; i_par < 7; ++i_par)
      {
	f_kaon_phi[i_eta][i_phi]->FixParameter(i_par,par_fit[i_par]);
      }
      f_kaon_phi[i_eta][i_phi]->SetRange(0.2,range[energy]);
      /* if(i_phi == 4 || i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(0.2,dip_3rd[energy]); // first */
      /* if(i_phi == 4 || i_phi == 5 || i_phi == 6) f_kaon_phi[i_eta][i_phi]->SetRange(dip_4th[energy]+0.5,range[energy]); // second */

      h_mEfficiency[HistName]->Fit(f_kaon_phi[i_eta][i_phi], "NR");

      HistName_FitPar = Form("h_mFitParameters_K%s_Cent_9_Eta_%d_Phi_%d",tof::mCharge[charge].c_str(),i_eta,i_phi);
      for(int i_bin = 0; i_bin < 7; ++i_bin)
      {
	h_mFitPar[HistName_FitPar]->SetBinContent(i_bin+1,f_kaon_phi[i_eta][i_phi]->GetParameter(i_bin));
      }
    }
  }

#if _PlotQA_

  TCanvas *c_Efficiency = new TCanvas("c_Efficiency","c_Efficiency",10,10,800,800);
  c_Efficiency->cd()->SetLeftMargin(0.1);
  c_Efficiency->cd()->SetBottomMargin(0.1);
  c_Efficiency->cd()->SetGrid(0,0);
  c_Efficiency->cd()->SetTicks(1,1);

  string output_start = Form("/Users/xusun/WorkSpace/Papers/VecMesonSpinAlignment/figures/Efficiency/ToF/QAs/EffFit_AuAu%s_K%s_QA_2060.pdf[",vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str());
  c_Efficiency->Print(output_start.c_str());

  string outputname = Form("/Users/xusun/WorkSpace/Papers/VecMesonSpinAlignment/figures/Efficiency/ToF/QAs/EffFit_AuAu%s_K%s_QA_2060.pdf",vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str());

  HistName = Form("h_mEfficiency_K%s_Cent_9",tof::mCharge[charge].c_str()); // ToF efficiency at cent 9
  h_mEffCent[HistName]->Draw("pE");
  h_mEffCent_plot[HistName]->Draw("pE same");
  f_kaon_cent->SetLineColor(4);
  f_kaon_cent->SetLineWidth(2);
  f_kaon_cent->SetLineStyle(2);
  f_kaon_cent->SetNpx(1000);
  f_kaon_cent->SetRange(0.2,8.0);
  f_kaon_cent->Draw("l same");
  c_Efficiency->Update();
  c_Efficiency->Print(outputname.c_str());

  for(int i_eta = 1; i_eta < tof::BinEta-1; ++i_eta) // ToF efficiency at different eta and phi bin at cent 9
  {
    TLegend *leg = new TLegend(0.5,0.7,0.8,0.9);
    leg->SetFillColor(10);
    leg->SetBorderSize(0);
    HistName = Form("h_mEfficiency_K%s_Cent_9_Eta_%d",tof::mCharge[charge].c_str(),i_eta);
    // h_mEffEta[HistName]->Draw("pE");
    h_mEffEta_plot[HistName]->Draw("pE");
    f_kaon_eta[i_eta]->SetLineColor(2);
    f_kaon_eta[i_eta]->SetLineWidth(4);
    f_kaon_eta[i_eta]->SetLineStyle(2);
    f_kaon_eta[i_eta]->SetNpx(1000);
    f_kaon_eta[i_eta]->SetRange(0.2,8.0);
    f_kaon_eta[i_eta]->Draw("l same");
    leg->AddEntry(h_mEffEta[HistName],"eta & phi integrated","p");
    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi) 
    {
      HistName = Form("h_mEfficiency_K%s_Cent_9_Eta_%d_Phi_%d",tof::mCharge[charge].c_str(),i_eta,i_phi);
      // h_mEfficiency[HistName]->Draw("pE same");
      h_mEfficiency_plot[HistName]->Draw("pE same");
      f_kaon_phi[i_eta][i_phi]->SetLineColor(mColor[i_phi]);
      f_kaon_phi[i_eta][i_phi]->SetLineWidth(2);
      f_kaon_phi[i_eta][i_phi]->SetLineStyle(mLineStyle[i_phi]);
      f_kaon_phi[i_eta][i_phi]->SetNpx(1000);
      // f_kaon_phi[i_eta][i_phi]->SetRange(0.2,range[energy]);
      f_kaon_phi[i_eta][i_phi]->SetRange(0.2,8.0);
      f_kaon_phi[i_eta][i_phi]->Draw("l same");
      string Leg_phi = Form("%d bin",i_phi);
      leg->AddEntry(h_mEfficiency[HistName],Leg_phi.c_str(),"p");
    }
    HistName = Form("h_mEfficiency_K%s_Cent_9_Eta_%d",tof::mCharge[charge].c_str(),i_eta);
    // h_mEffEta[HistName]->Draw("pE same");
    h_mEffEta_plot[HistName]->Draw("pE same");

    leg->Draw("same");
    c_Efficiency->Update();
    c_Efficiency->Print(outputname.c_str());
  }

  string output_stop = Form("/Users/xusun/WorkSpace/Papers/VecMesonSpinAlignment/figures/Efficiency/ToF/QAs/EffFit_AuAu%s_K%s_QA_2060.pdf]",vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str());
  c_Efficiency->Print(output_stop.c_str());

#endif

  string outputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/ToFMatch/FitPar_AuAu%s_K%s_2060.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),tof::mCharge[charge].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  HistName_FitPar = Form("h_mFitParameters_K%s_Cent_9",tof::mCharge[charge].c_str());
  h_mFitPar[HistName_FitPar]->Write();
  for(int i_eta = 0; i_eta < tof::BinEta; ++i_eta)
  {
    HistName_FitPar = Form("h_mFitParameters_K%s_Cent_9_Eta_%d",tof::mCharge[charge].c_str(),i_eta);
    h_mFitPar[HistName_FitPar]->Write();
    for(int i_phi = 0; i_phi < tof::BinPhi; ++i_phi) 
    {
      HistName_FitPar = Form("h_mFitParameters_K%s_Cent_9_Eta_%d_Phi_%d",tof::mCharge[charge].c_str(),i_eta,i_phi);
      h_mFitPar[HistName_FitPar]->Write();
    }
  }
}

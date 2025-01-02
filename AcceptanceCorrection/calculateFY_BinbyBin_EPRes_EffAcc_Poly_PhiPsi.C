#include "StRoot/Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "Utility/type.h"
#include <string>
//#include "SpecFuncMathMore.h"
//#include "Math.h"
#include <cmath>
#include "Math/SpecFuncMathMore.h"
#include <iostream>
#include <fstream>


using namespace std;

void calculateFY_BinbyBin_EPRes_EffAcc_Poly_PhiPsi(const int energy = 4, const int pid = 0, bool doall = false, bool isBesI = false, bool random3D = false, int mode = 1, int etamode = 0, int ypadding = 2, int ipt = 1, int icent = 1, int i_poly = 0, std::string ep = "sub", int order = 1, int yspectra = 0, int EP = 1, int v2 = 0) {

  std::string spectra;
  if(yspectra == 2) spectra = "_WithRapiditySpectra_HalfSigma";  
  if(yspectra == 1) spectra = "_WithRapiditySpectra";  
  if(yspectra == 0) spectra = "_NoRapiditySpectra";  
  if(yspectra == -1) spectra = "_NoRapiditySpectra_FixedFirstEP";  
  if(yspectra == -2) spectra = "_NoRapiditySpectra_FixedFirstEP_PhiPsi";  
  if(yspectra == -3) spectra = "_NoRapiditySpectra_PhiPsi";  

  std::string eptext = "";
  if(EP == 1) eptext = "_EP";
 
  std::string v2text = "";
  if(v2 == 0) v2text = "_noV2";

  std::string ordertext[2] = {"","2nd"};

  std::string xtraopt = "_noDelta";

  std::string etastring;
  if(etamode == 0) etastring = "eta1_eta1";
  if(etamode == 1) etastring = "eta1_eta1p5";
  if(etamode == 2) etastring = "eta1p5_eta1p5";
  if(etamode == 3) etastring = "eta0p4";
  if(etamode == 4) etastring = "eta0p6";
  if(etamode == 5) etastring = "eta0p8";

  gROOT->Reset();
  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};

  TH1D *h_theta_star_ratio[10][10][15];
  TH1D *h_theta_ratio[10][10][15];

  //TFile *MCFile;
  //double Fval[15] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  //double Ferr[15] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  //double Gval[15] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};
  //double Gerr[15] = {0.0};//{-0.0104367, -0.0104367, -0.0112319, -0.00879225, -0.00634735, -0.0053006};//{0.0}; //0.08555382,0.0660837,0.0363711,-0.00175735,0.0126171,-0.0257385};

  //MCFile = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/AcceptanceProccessTupleEPRes_%s%s/Acceptance_%s_19GeV_Mode%d_EtaMode%d.root",ep.c_str(),xtraopt.c_str(),vmsa::mPID[pid].c_str(),mode,etamode),"READ");

  //TCanvas *c_play = new TCanvas("c_play","c_play",10,10,1500,900);
  //c_play->Divide(5,3);
  //for(int i = 0; i < 15; i++)
  //{
  //  //c_play->cd(i+1)->SetTopMargin(0.10);
  //  c_play->cd(i+1)->SetLeftMargin(0.15);
  //  c_play->cd(i+1)->SetBottomMargin(0.15);
  //  c_play->cd(i+1)->SetGrid(0,0);
  //  c_play->cd(i+1)->SetTicks(1,1);
  //}
  ///////////////////////////// COS(THETA) /////////////////////////////////////
  //std::string outputname = Form("figures/Fplots_%s_%s_Mode%d_EtaMode%d_rapidity_theta_cent%d_pt%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),mode,etamode,icent,ipt,ep.c_str(),xtraopt.c_str());
  //std::string output_start = Form("%s[",outputname.c_str());
  //c_play->Print(output_start.c_str());
  //cout << "c_play->Print()" << endl;
  //
  //for(int iy = 0; iy < vmsa::y_total; iy++)
  //{
  //  c_play->cd(iy+1);
  //  h_theta_ratio[ipt][icent][iy] = (TH1D*) MCFile->Get(Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy))->Clone();
  //  cout << Form("Cent %d  %1.1f<pT<%1.1f",icent,vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]) << endl;
  //  h_theta_ratio[ipt][icent][iy]->GetXaxis()->SetTitle(Form("Cent %d  %1.1f <y<%1.1f  %1.1f<pT<%1.1f    cos#theta",icent,vmsa::y_bin[iy],vmsa::y_bin[iy+1],vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]));
  //  TF1 *Func_A = new TF1("Func_A","[0]*(1+[1]*x*x+[2]*x*x*x*x)",0,1);

  //  Func_A->SetParameter(0,h_theta_ratio[ipt][icent][iy]->GetBinContent(1));
  //  Func_A->SetParameter(1,0);
  //  h_theta_ratio[ipt][icent][iy]->Fit(Func_A,"ER");

  //  h_theta_ratio[ipt][icent][iy]->Draw("pE");
  //  Func_A->SetLineColor(kRed);
  //  Func_A->SetLineWidth(1);
  //  Func_A->Draw("same");

  //  Fval[iy] = Func_A->GetParameter(1);
  //  Ferr[iy] = Func_A->GetParError(1);
  //  Gval[iy] = Func_A->GetParameter(2);
  //  Gerr[iy] = Func_A->GetParError(2);

  //}
  //c_play->Update();
  //c_play->Print(outputname.c_str());

  //string output_stop = Form("%s]",outputname.c_str());
  //c_play->Print(output_stop.c_str()); // close pdf file
  ///////////////////////////// COS(THETA) /////////////////////////////////////

  ///////////////////////////// COS(THETA*) /////////////////////////////////////
  //outputname = Form("figures/Fplots_%s_%s_Mode%d_EtaMode%d_rapidity_thetastar_cent%d_pt%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),mode,etamode,icent,ipt,ep.c_str(),xtraopt.c_str());
  //output_start = Form("%s[",outputname.c_str());
  //c_play->Print(output_start.c_str());
  //cout << "c_play->Print()" << endl;
  //
  //for(int iy = 0; iy < vmsa::y_total; iy++)
  //{
  //  c_play->cd(iy+1);
  //  h_theta_star_ratio[ipt][icent][iy] = (TH1D*) MCFile->Get(Form("h_mEffCosS_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy))->Clone();
  //  h_theta_star_ratio[ipt][icent][iy]->Sumw2();
  //  cout << Form("Cent %d  %1.1f<pT<%1.1f",icent,vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]) << endl;
  //  h_theta_star_ratio[ipt][icent][iy]->GetXaxis()->SetTitle(Form("Cent %d  %1.1f <y<%1.1f  %1.1f<pT<%1.1f    cos#theta*",icent,vmsa::y_bin[iy],vmsa::y_bin[iy+1],vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]));
  //  TF1 *Func_A = new TF1("Func_A",Func4th,0,1,3);

  //  Func_A->SetParameter(0,h_theta_star_ratio[ipt][icent][iy]->GetBinContent(1));
  //  Func_A->SetParameter(1,0);
  //  h_theta_star_ratio[ipt][icent][iy]->Fit(Func_A,"ER");

  //  //Fval[iy] = Func_A->GetParameter(1);
  //  //Ferr[iy] = Func_A->GetParError(1);
  //  //Gval[iy] = Func_A->GetParameter(2);
  //  //Gerr[iy] = Func_A->GetParError(2);

  //  h_theta_star_ratio[ipt][icent][iy]->Draw("pE");
  //  Func_A->SetLineColor(kRed);
  //  Func_A->SetLineWidth(1);
  //  Func_A->Draw("same");
  //}
  //c_play->Update();
  //c_play->Print(outputname.c_str());

  //string output_stop = Form("%s]",outputname.c_str());
  //c_play->Print(output_stop.c_str()); // close pdf file
  ///////////////////////////// COS(THETA*) /////////////////////////////////////

  //for(int iy = 0; iy < vmsa::y_total; iy++)
  //{
  //  cout << "F = " << Fval[iy] << " +/- " << Ferr[iy]  << endl;
  //  cout << "G = " << Gval[iy] << " +/- " << Gerr[iy]  << endl;
  //}
  ///////////////////////////// GRAPHING F AND G VALUES /////////////////////////////////////
  //  
  //TGraphAsymmErrors *gFvalues[2];

  //TCanvas *c_play2 = new TCanvas("c_play2","c_play2",10,10,800,800);
  //c_play2->cd()->SetLeftMargin(0.15);
  //c_play2->cd()->SetBottomMargin(0.15);
  //c_play2->cd()->SetGrid(0,0);
  //c_play2->cd()->SetTicks(1,1);

  //string outputname = Form("figures/Fplots_%s_%s_Mode%d_EtaMode%d_rapidity_Fvalues_pt%d_cent%d_%s%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),mode,etamode,ipt,icent,ep.c_str(),xtraopt.c_str());
  //string output_start = Form("%s[",outputname.c_str());
  //c_play2->Print(output_start.c_str());
  //cout << "c_play2->Print()" << endl;

  //TH1F *h_play = new TH1F("h_play","h_play",100,-2.0,2.0);
  //for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  //{
  //  h_play->SetBinContent(i_bin+1,-1000.0);
  //  h_play->SetBinError(i_bin+1,1.0);
  //}
  //h_play->SetTitle("");
  //h_play->SetStats(0);
  //h_play->GetYaxis()->SetTitle("");
  //h_play->GetXaxis()->CenterTitle();
  //h_play->GetYaxis()->CenterTitle();
  //h_play->GetXaxis()->SetTitleSize(0.06);
  //h_play->GetYaxis()->SetTitleSize(0.06);
  //h_play->GetXaxis()->SetRangeUser(-1.5,1.5);
  //h_play->GetXaxis()->SetLabelSize(0.04);
  //h_play->GetYaxis()->SetLabelSize(0.04);
  //h_play->SetNdivisions(505,"X");
  //h_play->SetNdivisions(505,"Y");

  //for(int iF = 0; iF < 2; iF++)
  //{
  //  c_play2->cd();
  //  gFvalues[iF] = new TGraphAsymmErrors();
  //  cout << "icent = " << icent << endl;
  //  for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++)
  //  {
  //    double yval = (vmsa::y_bin[iy]+vmsa::y_bin[iy+1])/2.;
  //    if(iF == 0) 
  //    {
  //      gFvalues[iF]->SetPoint(iy-ypadding,yval,Fval[iy]);
  //      gFvalues[iF]->SetPointError(iy-ypadding,0.0,0.0,Ferr[iy],Ferr[iy]);
  //    }
  //    if(iF == 1) 
  //    {
  //      gFvalues[iF]->SetPoint(iy-ypadding,yval,Gval[iy]);
  //      gFvalues[iF]->SetPointError(iy-ypadding,0.0,0.0,Gerr[iy],Gerr[iy]);
  //    }
  //  }
  //  double min = TMath::MinElement(gFvalues[iF]->GetN(),gFvalues[iF]->GetY());
  //  double max = TMath::MaxElement(gFvalues[iF]->GetN(),gFvalues[iF]->GetY());
  //  if(min < 0.0)  min =  1.2*min;
  //  if(min > 0.0)  min =  0.6*min;
  //  if(min == 0.0) min = -0.2*max;
  //  if(max < 0.0)  max =  0.8*max;
  //  if(max > 0.0)  max =  1.2*max;
  //  if(max == 0.0) max = -0.2*min;
  //  h_play->GetYaxis()->SetRangeUser(min,max);
  //  h_play->GetXaxis()->SetTitle(Form("Cent %d  %1.1f<pT<%1.1f      y",icent,vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt]));
  //  h_play->DrawCopy("pE");
  //  gFvalues[iF]->SetMarkerStyle(20); 
  //  gFvalues[iF]->SetLineColor(kBlue); 
  //  gFvalues[iF]->SetMarkerColor(kBlue); 
  //  gFvalues[iF]->Draw("pE same"); 
  //   
  //  c_play2->Update();
  //  c_play2->Print(outputname.c_str());
  //}
  //string output_stop = Form("%s]",outputname.c_str());
  //c_play2->Print(output_stop.c_str()); // close pdf file
  /////////////////////////// GRAPHING F AND G VALUES /////////////////////////////////////
    
  Char_t Char[9][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  Double_t centCent[9] = {75.0,65.0,55.0,45.0,35.0,25.0,15.0,7.5,2.5};

  TF1 *line = new TF1("line","1/3",-0.5,8.5);
  line->SetLineStyle(2);
  line->SetLineWidth(1);
     
  TFile *fres1 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/Resolution_file_19GeV_EpdCorrections_4.root","READ");
  if(energy == 3) fres1 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/Resolution_file_14GeV_EpdCorrections_6.root","READ");
  TFile *fres2 = new TFile("../TreeProduction/UtilityFilesEta1p5OfficialCent/file_19GeV_Resolution.root","READ");
  if(energy == 3) fres2 = new TFile("../TreeProduction/StRoot/Utility/Resolution/file_14GeV_Resolution.root","READ");
  TFile  *fd12 = new TFile("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/SpinAlignment/Phi/file_19GeV_EpdFlow_4.root","READ");
  if(energy == 3) fd12 = new TFile("../TreeProduction/StRoot/Utility/EpdResolution/file_14GeV_EpdFlow_6.root","READ");
  
  TH1D *DeltaPsi1  = (TH1D*) fres1->Get("AveCosDeltaPsi1");
  TH1D *DeltaPsi2  = (TH1D*) fres2->Get("p_mRes2_Sub");
  TH1D *DeltaPsi12 = (TH1D*)  fd12->Get("p_mD12");

  TH1D *resolution_1 = new TH1D("resolution_1","resolution_1", 9, -0.5, 8.5);
  TH1D *resolution_12 = new TH1D("resolution_12","resolution_12", 9, -0.5, 8.5);
  TH1D *resolution_2 = new TH1D("resolution_2","resolution_2", 9, -0.5, 8.5);
  TH1D *resolution_2sub = new TH1D("resolution_2sub","resolution_2sub", 9, -0.5, 8.5);

  TGraphAsymmErrors *g_res1  = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_res2  = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_res12 = new TGraphAsymmErrors();

  std::ofstream outfile(Form("resolution_%s_%d.h",vmsa::mBeamEnergy[energy].c_str(),order));

  outfile << "const double r" << vmsa::mBeamEnergy[energy] << "ep" << order << "[3][2] = {" << std::endl;

  for(Int_t cent=0; cent<9; cent++) {

    Double_t CosMean = DeltaPsi1->GetBinContent(cent+1); // load from resolution file
    Double_t CosError = DeltaPsi1->GetBinError(cent+1);
    Double_t ZDCSMD_resSub = CosMean>0.? (TMath::Sqrt(CosMean)) : 0.;
    Double_t ZDCSMD_resSubErr = CosMean>0.? CosError/2./CosMean : 0.;
    Double_t ZDCSMD_chiSub = chi(ZDCSMD_resSub);
    cout << "chi value = " << chi(ZDCSMD_resSub) << endl;
    cout << "chi value * sqrt(2)  = " << chi(ZDCSMD_resSub)*TMath::Sqrt(2.) << endl;
    Double_t ZDCSMD_chiSubDelta = chi(ZDCSMD_resSub+0.005);
    Double_t ResZ = resEventPlane(TMath::Sqrt(2.)*ZDCSMD_chiSub);
    Double_t ResZDelta = resEventPlane(TMath::Sqrt(2.)*ZDCSMD_chiSubDelta);
    Double_t ResErrZ = ZDCSMD_resSubErr * TMath::Abs((ResZ-ResZDelta)/0.005);

    resolution_1->SetBinContent(cent+1, ResZ);
    resolution_1->SetBinError(cent+1, ResErrZ);

    cout << "cent = " << cent << ",   <cos(psi_e-psi_w)> = " << CosMean << ",  sub EP Res = " << ZDCSMD_resSub << ",  <cos2(psi-psiRP)> = " << ResZ << endl; 

    Double_t Cos12 = DeltaPsi12->GetBinContent(cent+1)/TMath::Sqrt(2.); // load file for D12 
    Double_t CosError12 = DeltaPsi12->GetBinError(cent+1)/TMath::Sqrt(2.);
    Double_t mean = Cos12/ResZ;
    Double_t error = ResErrZ*ResErrZ/ResZ/ResZ + CosError12*CosError12/Cos12/Cos12;

    resolution_12->SetBinContent(cent+1,mean);
    resolution_12->SetBinError(cent+1, mean*TMath::Sqrt(error));

    Double_t CosMean_2 = DeltaPsi2->GetBinContent(cent+1);
    Double_t CosError_2 = DeltaPsi2->GetBinError(cent+1);
    Double_t ZDCSMD_resSub_2 = CosMean_2>0.? (TMath::Sqrt(CosMean_2)) : 0.;
    Double_t ZDCSMD_resSubErr_2 = CosMean_2>0.? CosError_2/2./CosMean_2 : 0.;
    Double_t ZDCSMD_chiSub_2 = chi(ZDCSMD_resSub_2);
    Double_t ZDCSMD_chiSubDelta_2 = chi(ZDCSMD_resSub_2+0.005);
    Double_t ResZ_2 = res1(TMath::Sqrt(2.)*ZDCSMD_chiSub_2);
    Double_t ResZDelta_2 = res1(TMath::Sqrt(2.)*ZDCSMD_chiSubDelta_2);
    Double_t ResErrZ_2 = ZDCSMD_resSubErr_2 * TMath::Abs((ResZ_2-ResZDelta_2)/0.005);

    resolution_2->SetBinContent(cent+1,ResZ_2);
    resolution_2->SetBinError(cent+1,ResErrZ_2);
    resolution_2sub->SetBinContent(cent+1,ZDCSMD_resSub_2);
    resolution_2sub->SetBinError(cent+1,ZDCSMD_resSubErr_2);

    cout<<"Cent"<<cent_set[cent+1]<<"_"<<cent_set[cent]<<": "<<ResZ<<" +/- "<<ResErrZ<<" "<<ResZ_2<<" +/- "<<ResErrZ_2<<" "<<mean<<" +/- "<<mean*TMath::Sqrt(error)<<endl;

    g_res1->SetPoint(cent,centCent[cent],ResZ);
    g_res1->SetPointError(cent,0.0,0.0,ResErrZ,ResErrZ);
    g_res2->SetPoint(cent,centCent[cent],ResZ_2);
    g_res2->SetPointError(cent,0.0,0.0,ResErrZ_2,ResErrZ_2);
    g_res12->SetPoint(cent,centCent[cent],mean);
    g_res12->SetPointError(cent,0.0,0.0,mean*TMath::Sqrt(error),mean*TMath::Sqrt(error));
  }

  
  double Res_12[vmsa::cent_rebin_total] = {0.0};
  double Res_12_err[vmsa::cent_rebin_total] = {0.0};
  double Res_12_weight[vmsa::cent_rebin_total] = {0.0};
  if(order == 1)
  {
    for(int icent_rebin = 0; icent_rebin < vmsa::cent_rebin_total; icent_rebin++)
    {
      for(int i_cent = vmsa::cent_rebin[icent_rebin]; i_cent < vmsa::cent_rebin[icent_rebin+1]; i_cent++) {
        Res_12[icent_rebin] += resolution_1->GetBinContent(i_cent+1)/resolution_1->GetBinError(i_cent+1)/resolution_1->GetBinError(i_cent+1);
        Res_12_weight[icent_rebin] += 1./resolution_1->GetBinError(i_cent+1)/resolution_1->GetBinError(i_cent+1);
      }
      Res_12[icent_rebin] = Res_12[icent_rebin]/Res_12_weight[icent_rebin];
      Res_12_err[icent_rebin] = TMath::Sqrt(1./Res_12_weight[icent_rebin]); 
      cout << "centbin = " << icent_rebin << "    for cent range = [" << vmsa::cent_rebin[icent_rebin] << "," << vmsa::cent_rebin[icent_rebin+1]-1 << "]       with res12 = " << Res_12[icent_rebin] << " +/- " << Res_12_err[icent_rebin] << endl;
      outfile << "{" << Res_12[icent_rebin] << "," << Res_12_err[icent_rebin] << "}" << std::endl;
    }
  }

  if(order == 2)
  {
    for(int icent_rebin = 0; icent_rebin < vmsa::cent_rebin_total; icent_rebin++)
    {
      for(int i_cent = vmsa::cent_rebin[icent_rebin]; i_cent < vmsa::cent_rebin[icent_rebin+1]; i_cent++) {
        Res_12[icent_rebin] += resolution_12->GetBinContent(i_cent+1)/resolution_12->GetBinError(i_cent+1)/resolution_12->GetBinError(i_cent+1);
        Res_12_weight[icent_rebin] += 1./resolution_12->GetBinError(i_cent+1)/resolution_12->GetBinError(i_cent+1);
      }
      Res_12[icent_rebin] = Res_12[icent_rebin]/Res_12_weight[icent_rebin];
      Res_12_err[icent_rebin] = TMath::Sqrt(1./Res_12_weight[icent_rebin]); 
      cout << "centbin = " << icent_rebin << "    for cent range = [" << vmsa::cent_rebin[icent_rebin] << "," << vmsa::cent_rebin[icent_rebin+1]-1 << "]       with res12 = " << Res_12[icent_rebin] << " +/- " << Res_12_err[icent_rebin] << endl;
      outfile << "{" << Res_12[icent_rebin] << "," << Res_12_err[icent_rebin] << "}," << std::endl;
    }
  }
  if(random3D)
  {
    for(int icent_rebin = 0; icent_rebin < vmsa::cent_rebin_total; icent_rebin++)
    {
      Res_12[icent_rebin] = 0.0;//Res_12[icent_rebin]/Res_12_weight[icent_rebin];
      Res_12_err[icent_rebin] = 0.0;//TMath::Sqrt(1./Res_12_weight[icent_rebin]); 
      cout << "centbin = " << icent_rebin << "    for cent range = [" << vmsa::cent_rebin[icent_rebin] << "," << vmsa::cent_rebin[icent_rebin+1]-1 << "]       with res12 = " << Res_12[icent_rebin] << " +/- " << Res_12_err[icent_rebin] << endl;
    }
  }
  outfile << "};" << std::endl;
  outfile.close();
  //for(int icent_rebin = 0; icent_rebin < vmsa::cent_rebin_total; icent_rebin++)
  //{
  //  for(int i_cent = vmsa::cent_rebin[icent_rebin]; i_cent < vmsa::cent_rebin[icent_rebin+1]; i_cent++) {
  //    Res_12[icent_rebin] += resolution_2sub->GetBinContent(i_cent+1)/resolution_2sub->GetBinError(i_cent+1)/resolution_2sub->GetBinError(i_cent+1);
  //    Res_12_weight[icent_rebin] += 1./resolution_2sub->GetBinError(i_cent+1)/resolution_2sub->GetBinError(i_cent+1);
  //  }
  //  Res_12[icent_rebin] = Res_12[icent_rebin]/Res_12_weight[icent_rebin];
  //  cout << "CentBin = " << icent_rebin << "    for cent range = [" << vmsa::cent_rebin[icent_rebin] << "," << vmsa::cent_rebin[icent_rebin+1]-1 << "]       with Res12 = " << Res_12[icent_rebin] << endl;
  //}
  
  TString *Title;
  if(doall)
  {
    double eff[vmsa::y_total][7][7];
    double eff_error[vmsa::y_total][7][7];
    //TFile *eff_file = new TFile(Form("../RcPhiEffCorr/Rebinned/Eff_20GeV_SingleParticle_noToF_Mode2_EtaMode%d.root",etamode),"READ");
    //TFile *eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_NoEPResAndAcceptance/Eff_19GeV_SingleParticle_noToF_Mode2_EtaMode%d.root",etamode),"READ");
    TFile *eff_file;
    if(order == 1) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order1%s%s%s/Eff_%s_SingleParticle_noToF_Mode2_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
    if(order == 2) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu%s_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order2%s%s%s/Eff_%s_SingleParticle_noToF_Mode2_EtaMode%d.root",vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str(),vmsa::mBeamEnergy[energy].c_str(),etamode),"READ");
    if(random3D) eff_file = new TFile(Form("/gpfs01/star/pwg/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/CosEfficiencyProccessTuple_EPRes_EffAcc_Ycut_Order3/Eff_19GeV_SingleParticle_noToF_Mode2_EtaMode%d.root",etamode),"READ");
    eff_file->Print();

    for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++)
    { 
      Title = new TString(Form("h_mEffCos_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy));
      cout << Title->Data() << endl;
      TH2D *eff_hist = (TH2D*)eff_file->Get(Title->Data());
      eff_hist->Print();
      for(int itheta = 0; itheta < 7; itheta++)
      {
        for(int iphipsi = 0; iphipsi < 7; iphipsi++)
        {
          eff[iy][itheta][iphipsi] = eff_hist->GetBinContent(itheta+1,iphipsi+1);
          eff_error[iy][itheta][iphipsi] = eff_hist->GetBinError(itheta+1,iphipsi+1)/eff[iy][itheta][iphipsi];
          cout << "icent = " << icent << " ipt = " << ipt << " iy = " << iy << " itheta = " << itheta << " eff = " << eff[iy][itheta][iphipsi] << " +/- " << eff_error[iy][itheta][iphipsi] << endl;
        }
      }
      delete Title;
      delete eff_hist;
    }
    eff_file->Close();

    TFile *input;
    if(order == 1) input = new TFile(Form("rho00/%s/%s/Poly/PhiPsi_RawRhoEtaSys_%s_Poly_FirstOrder.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str()),"READ");
    if(order == 2) input = new TFile(Form("rho00/%s/%s/Poly/PhiPsi_RawRhoEtaSys_%s_Poly.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str()),"READ");
    if(random3D) input = new TFile(Form("rho00/%s/%s/3DRandom/RawRhoEtaSys_%s_FirstOrder.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str()),"READ");
    
    //string outputname = Form("output/%s/%s/AccResRhoEtaSys_%s_cent_%d_ipt_%d_%s%s_2ndSubEPRes.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str(),icent,ipt,ep.c_str(),xtraopt.c_str());
    string outputname;
    if(order == 1) outputname = Form("output/%s/%s/PhiPsi_AccResRhoEtaSys_%s_Poly_FirstOrder_cent_%d_ipt_%d_poly%d%s%s%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str(),icent,ipt,i_poly+1,spectra.c_str(),eptext.c_str(),v2text.c_str());
    if(order == 2) outputname = Form("output/%s/%s/PhiPsi_AccResRhoEtaSys_%s_Poly_cent_%d_ipt_%d_poly%d%s%s%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str(),icent,ipt,i_poly+1,spectra.c_str(),eptext.c_str(),v2text.c_str());
    if(random3D) outputname = Form("output/%s/%s/3DRandom/AccResRhoEtaSys_%s_FirstOrder_cent_%d_ipt_%d.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str(),icent,ipt);
    TFile *output = new TFile(outputname.c_str(),"RECREATE");
    output->cd();

    int centVal[4]  = {80,40,10,0};

    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    {
      for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
      {
        if( i_dca != 0 && i_sig != 0 ) continue;
        for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
        {
          for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
          {
            for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
            {

                TH1D *rho00_hist;
                TGraphAsymmErrors *g_rho00 = new TGraphAsymmErrors();
    
                double weight_rho00 = 0.0;
                double weight_error_rho00 = 0.0;
                double weight_all = 0.0;

                double y_rho00[vmsa::y_total] = {0.0};
                double y_error_rho00[vmsa::y_total] = {0.0};
                double y_all[vmsa::y_total] = {0.0};
            
                TCanvas *c_fit = new TCanvas("c_fit","c_fit",10,10,1500,600);
                c_fit->Divide(5,2);
                for(int i = 0; i < 10; i++)
                {
                  c_fit->cd(i+1)->SetLeftMargin(0.15);
                  c_fit->cd(i+1)->SetBottomMargin(0.15);
                  c_fit->cd(i+1)->SetGrid(0,0);
                  c_fit->cd(i+1)->SetTicks(1,1);
                }
                Title = new TString(Form("fit/Phi/%s/Order%d/yield_pt_%d_cent_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_EtaMode%d_Divided.pdf",vmsa::mBeamEnergy[energy].c_str(),order,ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1,etamode));

                string yieldtitle = Form("yield_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);                             
                TH1D *yieldvsy = new TH1D(yieldtitle.c_str(),yieldtitle.c_str(),10,-1.0,1.0);
                
                TH2F *PtCos_2D[vmsa::y_total];
                TH1F *PtCos[vmsa::y_total];
                TF1 *Func_rho[vmsa::y_total];
                TF1 *Func_obs[vmsa::y_total];
                for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
                {
                  string key = Form("eta_%d_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",iy,ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  cout << key << endl;

                  TH2F *PtCos_raw = (TH2F*)input->Get(key.c_str())->Clone("PtCos_raw");
                  PtCos_2D[iy] = new TH2F(key.c_str(),key.c_str(), 7, 0, 1, 7, -TMath::Pi()/2.0, TMath::Pi()/2.0);
                  PtCos_2D[iy]->Sumw2();
                  //delete Title;
                  for(int itheta = 0; itheta < 7; itheta++) 
                  {
                    for(int iphipsi = 0; iphipsi < 7; iphipsi++) 
                    {
                      float inte_mean = PtCos_raw->GetBinContent(itheta+1,iphipsi+1);
                      float inte_mean_error = PtCos_raw->GetBinError(itheta+1,iphipsi+1);
                      //cout << "inte_mean  = " << inte_mean << " +/- " << inte_mean_error << endl;
                      PtCos_2D[iy]->SetBinContent(itheta+1, iphipsi+1, inte_mean/eff[iy][itheta][iphipsi]);
                      PtCos_2D[iy]->SetBinError(itheta+1, iphipsi+1, inte_mean/eff[iy][itheta][iphipsi]*TMath::Sqrt(inte_mean_error*inte_mean_error/inte_mean/inte_mean+eff_error[iy][itheta][iphipsi]*eff_error[iy][itheta][iphipsi]));
                    }
                  }
                  PtCos_2D[iy]->Write();
                 
                  PtCos[iy] = (TH1F*) PtCos_2D[iy]->ProjectionX();

                  c_fit->cd(iy+1-ypadding);
                  PtCos[iy]->GetXaxis()->SetTitleOffset(1.2);
                  PtCos[iy]->SetTitle(Form("%1.1f<y<%1.1f,%1.1f<p_{T}<%1.1fGeV/c,Cent %d-%d",vmsa::y_bin[iy],vmsa::y_bin[iy+1],vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt],centVal[icent+1],centVal[icent]));
                  std::string title = Form("%1.1f<y<%1.1f,%1.1f<p_{T}<%1.1fGeV/c,Cent %d-%d",vmsa::y_bin[iy],vmsa::y_bin[iy+1],vmsa::pt_low_y[energy][ipt],vmsa::pt_up_y[energy][ipt],centVal[icent+1],centVal[icent]);
                  std::cout << title.c_str() << std::endl;
                  PtCos[iy]->GetXaxis()->SetTitle("cos#theta*");
                  PtCos[iy]->GetYaxis()->SetTitle("EffxAcc Corrected Yield");
                  PtCos[iy]->GetYaxis()->SetTitleOffset(1.0);
                  PtCos[iy]->SetMarkerColor(kBlack);
                  PtCos[iy]->SetMarkerSize(1.8);
                  PtCos[iy]->SetMarkerStyle(20);
                  PtCos[iy]->Draw("pE");
                
                  //Func_rho[iy]->SetParameter(0, PtCos[iy]->GetBinContent(5));
                  //Func_rho[iy]->SetParameter(1, 0.333);
                  //Func_rho[iy]->SetParLimits(1, 0.0, 1.0);
                  //Func_rho[iy]->FixParameter(2, Fval[iy]);
                  //Func_rho[iy]->FixParameter(3, Gval[iy]);
                  //Func_rho[iy]->FixParameter(4, Res_12[icent]);
                  //PtCos[iy]->Fit(Func_rho[iy], "NMRI");
                  //Func_rho[iy]->SetLineColor(kRed);
                  //Func_rho[iy]->Draw("same");

                  Func_rho[iy] = new TF1(Form("Func_rho_%d",iy),FuncAFG,0,1,5);
                  Func_obs[iy] = new TF1("Func_obs","[0]*(1.-[1]+(3.*[1]-1)*(x*x))",0,1);

                  PtCos[iy]->Fit(Func_obs[iy],"QNMI"); // fit corrected distribution for observerved rho00
                  double rho_obs = Func_obs[iy]->GetParameter(1);
                  double rho_obs_err = Func_obs[iy]->GetParError(1);
                  Func_obs[iy]->SetLineColor(kRed);
                  Func_obs[iy]->Draw("same");
       
                  cout << "iy = " << iy << ", rho_obs = " << rho_obs << endl;       

                  // Error calculation
                  double drhodobs = 4./(1.+3.*Res_12[icent]);  // calculate d(rho)/d(rho_obs)
                  double drhodR = -12.*(rho_obs - 1./3.)/(1.+3.*Res_12[icent])/(1.+3.*Res_12[icent]); // calculate d(rho)/d(R)

                  double real_rho = 1./3. + 4./(1.+3.*Res_12[icent])*(rho_obs - 1./3.);
                  double real_rho_error = TMath::Sqrt((rho_obs_err*rho_obs_err)*(drhodobs*drhodobs) + (Res_12_err[icent]*Res_12_err[icent])*(drhodR*drhodR));
                  cout << "iy = " << iy << ", real_rho = " << real_rho << endl;       
           
                  //float real_rho = Func_rho[iy]->GetParameter(1);
                  //float real_rho_error = Func_rho[iy]->GetParError(1);
                  Double_t interr = 0.0;
                  float weight = PtCos[iy]->IntegralAndError(1,7,interr);
                  double ymean = (vmsa::y_bin[iy] + vmsa::y_bin[iy+1])/2.0;
                  int ybin = yieldvsy->FindBin(ymean);
                  yieldvsy->SetBinContent(ybin,weight);
                  yieldvsy->SetBinError(ybin,interr);
                   
                  cout << "bin = " << ybin << ", y = " << ymean << ", yield = " << weight << " +/- " << interr << endl;

                  weight_rho00 += real_rho*weight;
                  weight_error_rho00 += real_rho_error*real_rho_error*weight*weight;
                  weight_all += weight;

                  y_rho00[iy] += real_rho*weight;
                  y_error_rho00[iy] += real_rho_error*real_rho_error*weight*weight;
                  y_all[iy] += weight;

                  delete PtCos_raw;
                }
                yieldvsy->Write();
                c_fit->SaveAs(Title->Data());

                TCanvas *c_yield = new TCanvas("c_yield","c_yield",10,10,600,600);
                c_yield->SetLeftMargin(0.15);
                c_yield->SetBottomMargin(0.15);
                c_yield->SetGrid(0,0);
                c_yield->SetTicks(1,1);
                yieldvsy->SetMaximum(yieldvsy->GetMaximum()*1.1);
                yieldvsy->SetMinimum(yieldvsy->GetMinimum()*0.9);

                yieldvsy->GetXaxis()->SetTitle("y");
                yieldvsy->GetYaxis()->SetTitle("dN/dy (EffXAcc Corrected)");
                yieldvsy->SetMarkerStyle(20);
                yieldvsy->Draw("pE");
                c_yield->SaveAs(Form("figures/%s/yield_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d.pdf",vmsa::mBeamEnergy[energy].c_str(),ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             

                //delete c_fit;
                weight_rho00 = weight_rho00/weight_all;
                weight_error_rho00 = TMath::Sqrt(weight_error_rho00)/weight_all;

                for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
                {
                  y_rho00[iy] = y_rho00[iy]/y_all[iy];
                  y_error_rho00[iy] = TMath::Sqrt(y_error_rho00[iy])/y_all[iy];
                }

                cout << "rho00 = " << weight_rho00<< " +/- " << weight_error_rho00 << endl;

                for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
                {
                  cout << y_rho00[iy];
                  if(iy == vmsa::y_total-ypadding-1) cout << " " << weight_rho00 << endl;
                  else cout << " ";
                }

 
                Title = new TString(Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             

                g_rho00 = new TGraphAsymmErrors();
                g_rho00->SetName(Title->Data());
                for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++)
                {
                  double ymean = (vmsa::y_bin[iy] + vmsa::y_bin[iy+1])/2.0;
                  g_rho00->SetPoint(iy-ypadding,ymean,y_rho00[iy]);
                  g_rho00->SetPointError(iy-ypadding,0.0,0.0,y_error_rho00[iy],y_error_rho00[iy]);
                  cout << "iy = " << iy << " ymean = " << ymean << " y_rho00[iy] = " << y_rho00[iy] << " +/- " << y_error_rho00[iy] << endl;
                }
                delete Title;

                Title = new TString(Form("rhoFinalWeighted_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             

                rho00_hist = new TH1D(Title->Data(), Title->Data(), vmsa::y_total+1, 0.5, (double)vmsa::y_total+1.5);
                rho00_hist->SetBinContent(vmsa::y_total+1,weight_rho00);
                rho00_hist->SetBinError(vmsa::y_total+1,weight_error_rho00);
                for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
                {
                  rho00_hist->SetBinContent(iy+1,y_rho00[iy]);
                  rho00_hist->SetBinError(iy+1,y_error_rho00[iy]);
                }

                g_rho00->Write();
                rho00_hist->Write();

                delete Title;
                delete rho00_hist;
                delete g_rho00;
                delete yieldvsy;
                //delete PtCos;
               
                for(int iy = 0+ypadding; iy < vmsa::y_total-ypadding; iy++) 
                {
                  delete PtCos[iy];
                }
             
            }
          }
        }
      }
    }
    
    output->Close();
  }
}

void Correction(Double_t Res, Double_t Res_error, Double_t obv_rho, Double_t obv_rho_error, Double_t &real_rho, Double_t &real_rho_error) {

  real_rho = (Res-1+4*obv_rho)/(1+3*Res);

  Double_t real_rho_error_1 = obv_rho_error*4/(1+3*Res);
  Double_t real_rho_error_2 = Res_error*(4-12*obv_rho)/(1+3*Res)/(1+3*Res);
  real_rho_error = TMath::Sqrt(real_rho_error_1*real_rho_error_1+real_rho_error_2*real_rho_error_2);

}

Double_t chi(Double_t res) {

  Double_t chi = 2.0;
  Double_t delta = 1.0;

  for(Int_t i=0; i<15; i++) {
    chi = (res1(chi) < res) ? chi + delta : chi - delta;
    delta = delta/2.;
  }

  return chi;
}

Double_t res1(Double_t chi) {

  Double_t con = 0.626657;                   // TMath::Sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;

  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

  return res;
}

Double_t resEventPlane(Double_t chi) {

  Double_t con = 0.626657;                   // TMath::Sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;
  Double_t halfpi = TMath::Pi()/2.;

  Double_t besselOneHalf = TMath::Sqrt(arg/halfpi) * TMath::SinH(arg)/arg;
  Double_t besselThreeHalfs = TMath::Sqrt(arg/halfpi) * (TMath::CosH(arg)/arg - TMath::SinH(arg)/(arg*arg));
  Double_t res = con * chi * exp(-arg) * (besselOneHalf + besselThreeHalfs);

  return res;
}

Double_t resEventPlane1(Double_t chi) {

  Double_t con = 0.626657;                   // TMath::Sqrt(pi/2)/2
  Double_t arg = chi * chi / 4.;
  Double_t halfpi = TMath::Pi()/2.;

  TF1* iBesselHalf = new TF1("I_12", "ROOT::Math::cyl_bessel_i([0],x)", 0, 10);
  iBesselHalf->SetParameter(0,1./2.);
  TF1* iBessel3Half = new TF1("I_32", "ROOT::Math::cyl_bessel_i([0],x)", 0, 10);
  iBessel3Half->SetParameter(0,3./2.);
  Double_t besselOneHalf = iBesselHalf->Eval(arg);
  Double_t besselThreeHalfs = iBessel3Half->Eval(arg);
  Double_t res = con * chi * exp(-arg) * (besselOneHalf + besselThreeHalfs);

  return res;
}

double FuncAD(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double D = par[2];
  double R = par[3];

  double A = (3.*rho-1.)/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));

  double result = (1.+Bs*D/2.) + (As+D)*CosTheta*CosTheta + (As*D-Bs*D/2.)*CosTheta*CosTheta*CosTheta*CosTheta;

  return N*result;

}

double Func4th(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double F = par[1];
  double G = par[2];

  //double result = 1. + (4.*F+3.*G)/8. - (2.*F+3.*G)/4.*CosTheta*CosTheta + 3.*G/8.*CosTheta*CosTheta*CosTheta*CosTheta;
   
  double order0 = 2. + F + 3.*G/4.;
  double order2 = (-1.*F - 3.*G/2.)*CosTheta*CosTheta;
  double order4 = (3.*G/4.)*CosTheta*CosTheta*CosTheta*CosTheta;
 
  double result = order0 + order2 + order4;

  return N*result;

}

double FuncAFG(double *x_val, double *par) {

  double CosTheta = x_val[0];
  double N = par[0];
  double rho = par[1];
  double F = par[2];
  double G = par[3];
  double R = par[4];
  //double rho10 = par[5];

  double A = (3.*rho-1.)/(1.-rho);
  //double B = rho10/(1.-rho);
  double As = A*(1.+3.*R)/(4.+A*(1.-R));
  double Bs = A*(1.-R)/(4.+A*(1.-R));
  //double denom = 4. + A*(1.-R) + B*(-1.+R);
  //double As = (A*(1.+3.*R) + B*(3.-3.*R))/denom;
  //double Bs = (A*(1.-R)    + B*(3.+R)   )/denom;


  double order0 = 2. + F - Bs*F/2. + 3.*G/4. - Bs*G/2.;    
  double order2 = (2.*As - F + As*F + Bs*F - 3.*G/2. + 3.*As*G/4. + 3.*Bs*G/2.)*CosTheta*CosTheta;
  double order4 = (-1.*As*F - Bs*F/2. + 3.*G/4. - 3.*As*G/2. - 3.*Bs*G/2.)*CosTheta*CosTheta*CosTheta*CosTheta;
  double order6 = (3.*As*G/4. + Bs*G/2.)*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta*CosTheta;

  double result = order0 + order2 + order4 + order6;

  return N*result;

}

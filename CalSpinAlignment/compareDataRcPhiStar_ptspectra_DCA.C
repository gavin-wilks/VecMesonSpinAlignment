#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
//#include "../Utility/draw.h"
#include "../../../draw.h"
#include "../Utility/functions.h"
#include "resolution_pt.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void compareDataRcPhiStar_ptspectra_DCA(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 40)//defaultF = 0 is BESII, defaultF = 1 is BESI
{
  std::string filename = ""; 
  std::string spectra = "";
  if(yspectra == 0) spectra = "_NoRapiditySpectra";
  if(yspectra == 1) spectra = "_NoRapiditySpectra_EP";
  if(yspectra == 2) spectra = "_NoRapiditySpectra_EP_noV2";
  if(yspectra == 3) spectra = "_WithRapiditySpectra";
  if(yspectra == 4) spectra = "_WithRapiditySpectra_HalfSigma";
  if(yspectra == 5) spectra = "_NoRapiditySpectra_FixedFirstEP";
  if(yspectra == 6) spectra = "_NoRapiditySpectra_FixedFirstEP_InputRho";
  if(yspectra == 7) spectra = "_NoRapiditySpectra_InputRho0p4_Fixed";
  if(yspectra == 8) spectra = "_NoRapiditySpectra_v2times3";
  if(yspectra == 9) spectra = "_NoRapiditySpectra_FixedFirstEP_AcceptanceOnly";
  if(yspectra == 10) spectra = "_NoRapiditySpectra_FixedFirstEP_ptmin0p5";
  if(yspectra == 11) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins";
  if(yspectra == 12) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_InputRho0p5";
  if(yspectra == 13) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_InputRho0p5_3xV2";
  if(yspectra == 14) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_Eta0p5";
  if(yspectra == 15) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_rhoforeachstudy";
  if(yspectra == 16) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP";
  if(yspectra == 17) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc";
  if(yspectra == 18) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc_pT0p2";
  if(yspectra == 19) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc_pT0p2_TPC";
  if(yspectra == 20) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_20240319_Acc_pT0p2_TPC_TOF";
  if(yspectra == 21) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_EPeff_flatRP_Acc1Rc_pT0p1_TPC_TOF_kaons_20240327";
  if(yspectra == 22) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_phis_20240423";
  if(yspectra == 23) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_phis_20240426_withRho_deltar";
  if(yspectra == 24) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240508_fixed2Dfunc_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0";
  if(yspectra == 25) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240508_fixed2Dfunc_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0_SpectraWeight";
  if(yspectra == 26) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240523_fixed2Dfunc_PhiSpectra_rerho1n10.0_ysigma1000_rho000.3333_r0.0_i0.0_imrho1n10.0_helicityrho0.3333";
  if(yspectra == 27) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_NoWeight_rerho1n10.0_ysigma1000_noalignment";
  if(yspectra == 28) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_PhiWeight_rerho1n10.0_ysigma1000_noalignment";
  if(yspectra == 29) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_PhiAndKaonWeight_rerho1n10.0_ysigma1000_noalignment";
  if(yspectra == 30) spectra = "_NoRapiditySpectra_FixedFirstEP_RCBins_prelimv2_phieff_flatRP_Acc1Rc_pT0p1_TPC_TOF_everything_20240529_fixed2Dfunc_FinerPhiBins_PhiAndKaonWeightAgain_rerho1n10.0_ysigma1000_noalignment";
  if(yspectra == 35) spectra = "Pythia_PhiEmbedding_PIDEff"; 
  if(yspectra == 39) spectra = "PhiEmbeddingWithWeights_PIDEff_20240726"; 
  if(yspectra == 40)
  {
    spectra = "Pythia_PIDEff_Helicity2D_20240729";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 41)
  {
    spectra = "Pythia_PIDEff_Helicity2Dplusdelta_20240729";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 42)
  {
    spectra = "Pythia_PIDEff_Helicity2Dminusdelta_20240729";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 43)
  {
    spectra = "Pythia_PIDEff_Helicity2Dminusdeltaplusrhodelta_20240729";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 44)
  {
    spectra = "Pythia_PIDEff_20240729";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 45)
  {
    spectra = "Pythia_PIDEff_Helicity2Dminusdeltaplusrhodelta_20240730";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
  if(yspectra == 46)
  {
    spectra = "Pythia_PIDEff_Sigmay1p5_20240731";
    filename = "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"; 
  }
 
  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  string inputfile = Form("../output/AuAu%s/%s/ptyspectra_%sPhiPtSys_%s_PolySys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  string inputfileRC = Form("effaccfiles/%s/%s/%s/Eff_%s_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPutRC = TFile::Open(inputfileRC.c_str());

  TH2D *h_mCounts;  
  TH2D *h_mCounts_RC;  
  TH2D *h_mCounts_Ratio;  
  TH1D *h_mCounts_pt;  
  TH1D *h_mCounts_RC_pt;  
  TH1D *h_mCounts_Ratio_pt;  
  TH1D *h_mCounts_y;  
  TH1D *h_mCounts_RC_y;  
  TH1D *h_mCounts_Ratio_y;  
  TH1D *h_mCounts_ybins[vmsa::pt_total];  
  TH1D *h_mCounts_RC_ybins[vmsa::pt_total];  
  TH1D *h_mCounts_Ratio_ybins[vmsa::pt_total];  
  
  cout << inputfile << endl;
  cout << inputfileRC << endl;
 
  string KEY_counts = Form("ptyspectra_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  h_mCounts = (TH2D*) File_InPut->Get(KEY_counts.c_str());
  //h_mCounts->RebinY(2);
  string KEY_counts_RC = Form("ptyspectra_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),2,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  h_mCounts_RC = (TH2D*) File_InPut->Get(KEY_counts_RC.c_str());
  //h_mCounts_RC->RebinY(2);


  ///////////// 2D y vs pT /////////////
  float data = h_mCounts->Integral(0,-1,0,-1);
  float rc   = h_mCounts_RC->Integral(0,-1,0,-1);

  h_mCounts_RC->Scale(data/rc);
  cout << "NBinsX = " << h_mCounts_RC->GetNbinsX() << endl;
  cout << "NBinsY = " << h_mCounts_RC->GetNbinsY() << endl;

  h_mCounts_Ratio = (TH2D*) h_mCounts->Clone();
  h_mCounts_Ratio->Divide(h_mCounts_RC);
  ///////////////////////////////////////

  for(int ipt = 0; ipt < vmsa::pt_total; ipt++)
  {
    string KEY_counts = Form("counts_pt_%d",ipt);
    h_mCounts_ybins[ipt] = (TH1D*) h_mCounts->ProjectionY(KEY_counts.c_str(),ipt+1,ipt+1,"e");
    string KEY_counts_RC = Form("counts_rc_pt_%d",ipt);
    h_mCounts_RC_ybins[ipt] = (TH1D*) h_mCounts_RC->ProjectionY(KEY_counts_RC.c_str(),ipt+1,ipt+1,"e");
    
    string KEY_counts_RC_ratio = Form("countsratio_rc_pt_%d",ipt);
    h_mCounts_Ratio_ybins[ipt] = (TH1D*) h_mCounts_ybins[ipt]->Clone(KEY_counts_RC_ratio.c_str());
    h_mCounts_Ratio_ybins[ipt]->Divide(h_mCounts_RC_ybins[ipt]);
  }
 
  ///////////// 1D pT /////////////
  h_mCounts_pt = (TH1D*) h_mCounts->ProjectionX("pt",0,-1,"e"); 
  h_mCounts_RC_pt = (TH1D*) h_mCounts_RC->ProjectionX("ptrc",0,-1,"e"); 
  //float datapt = h_mCounts_pt->Integral(0,-1);
  //float rcpt   = h_mCounts_RC_pt->Integral(0,-1);

  //h_mCounts_RC_pt->Scale(datapt/rcpt);

  h_mCounts_Ratio_pt = (TH1D*) h_mCounts_pt->Clone();
  h_mCounts_Ratio_pt->Divide(h_mCounts_RC_pt);
  ///////////////////////////////////////
  
  ///////////// 1D y /////////////
  h_mCounts_y = (TH1D*) h_mCounts->ProjectionY("y",0,-1,"e"); 
  h_mCounts_RC_y = (TH1D*) h_mCounts_RC->ProjectionY("yrc",0,-1,"e"); 
  //float datay = h_mCounts_y->Integral(0,-1);
  //float rcy   = h_mCounts_RC_y->Integral(0,-1);

  //h_mCounts_RC_y->Scale(datay/rcy);

  h_mCounts_Ratio_y = (TH1D*) h_mCounts_y->Clone();
  h_mCounts_Ratio_y->Divide(h_mCounts_RC_y);
  ///////////////////////////////////////

  TCanvas *c = new TCanvas("c","c",10,10,400,400);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);
  c->SetTicks(1,1);
  c->SetGrid(0,0);

  c->SetLogy();
  h_mCounts_pt->SetTitle(Form("p_{T} 20-60 Centrality"));
  h_mCounts_pt->GetXaxis()->SetTitle("p_{T}");
  h_mCounts_pt->GetYaxis()->SetTitle("Yields");
  h_mCounts_pt->SetMarkerStyle(20);
  h_mCounts_pt->SetMarkerColor(kOrange+7);
  h_mCounts_pt->SetLineColor(kOrange+7);
  h_mCounts_RC_pt->SetMarkerStyle(20);
  h_mCounts_RC_pt->SetMarkerColor(kBlue);
  h_mCounts_RC_pt->SetLineColor(kBlue);

  double max = h_mCounts_pt->GetMaximum();
  double min = h_mCounts_pt->GetMinimum();
  if(h_mCounts_RC_pt->GetMaximum() > max) max = h_mCounts_RC_pt->GetMaximum();
  if(h_mCounts_RC_pt->GetMinimum() < min) min = h_mCounts_RC_pt->GetMinimum();
 
  TLegend *leg = new TLegend(0.6,0.5,0.8,0.7);
  leg->AddEntry(h_mCounts_pt,"|DCA| < 2 cm","p");
  leg->AddEntry(h_mCounts_RC_pt,"|DCA| < 3 cm","p");

  h_mCounts_pt->Draw("pE");
  h_mCounts_RC_pt->Draw("pE same");
  leg->Draw("same");  

  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_ptspectra_%s_Order%d_%s_DATADCA.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),correction.c_str(),order,etamode.c_str()));  

  c->SetLogy(0);
  h_mCounts_y->SetTitle(Form("y 20-60 Centrality"));
  h_mCounts_y->GetXaxis()->SetTitle("y");
  h_mCounts_y->GetYaxis()->SetTitle("Yields");
  h_mCounts_y->SetMarkerStyle(20);
  h_mCounts_y->SetMarkerColor(kOrange+7);
  h_mCounts_y->SetLineColor(kOrange+7);
  h_mCounts_RC_y->SetMarkerStyle(20);
  h_mCounts_RC_y->SetMarkerColor(kBlue);
  h_mCounts_RC_y->SetLineColor(kBlue);

  max = h_mCounts_y->GetMaximum();
  min = h_mCounts_y->GetMinimum();
  if(h_mCounts_RC_y->GetMaximum() > max) max = h_mCounts_RC_y->GetMaximum();
  if(h_mCounts_RC_y->GetMinimum() < min) min = h_mCounts_RC_y->GetMinimum();
 
  TLegend *legy = new TLegend(0.4,0.2,0.6,0.4);
  legy->AddEntry(h_mCounts_y,"|DCA| < 2.0 cm","p");
  legy->AddEntry(h_mCounts_RC_y,"|DCA| < 3.0 cm","p");

  h_mCounts_y->Draw("pE");
  h_mCounts_RC_y->Draw("pE same");
  legy->Draw("same");  

  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_yspectra_%s_Order%d_%s_DATADCA.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),correction.c_str(),order,etamode.c_str()));  
  
  h_mCounts_Ratio->SetTitle(Form("y vs pT ratio 20-60 Centrality"));
  h_mCounts_Ratio->GetXaxis()->SetTitle("p_{T}");
  h_mCounts_Ratio->GetYaxis()->SetTitle("y");
  h_mCounts_Ratio->Draw("colz");
  
  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_ptyspectra2Dratio_%s_Order%d_%s_DATADCA.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),correction.c_str(),order,etamode.c_str()));  

  h_mCounts_Ratio_pt->SetTitle(Form("20-60 Centrality"));
  h_mCounts_Ratio_pt->GetXaxis()->SetTitle("p_{T}");
  h_mCounts_Ratio_pt->GetYaxis()->SetTitle("(|DCA|<2cm)/(|DCA|<3cm)");
  h_mCounts_Ratio_pt->GetYaxis()->SetRangeUser(0.9,1.1);
  h_mCounts_Ratio_pt->Draw("pE");
  
  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_ptspectra1Dratio_%s_Order%d_%s_DATADCA.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),correction.c_str(),order,etamode.c_str()));  

  h_mCounts_Ratio_y->SetTitle(Form("20-60 Centrality"));
  h_mCounts_Ratio_y->GetXaxis()->SetTitle("y");
  h_mCounts_Ratio_y->GetYaxis()->SetTitle("(|DCA|<2cm)/(|DCA|<3cm)");
  h_mCounts_Ratio_y->GetYaxis()->SetRangeUser(0.6,1.4);
  h_mCounts_Ratio_y->Draw("pE");
  
  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_yspectra1Dratio_%s_Order%d_%s_DATADCA.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),correction.c_str(),order,etamode.c_str()));  

  string outputfile = "ptyspectra_datarcratio_19GeV.root";
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_mCounts_Ratio->SetName("pty_datarcratio");
  h_mCounts_Ratio->Write();  
  File_OutPut->Close();  

  TCanvas *cy = new TCanvas("c","c",10,10,2500,2500);
  cy->Divide(5,5);
  
  for(int ipt = 0; ipt < vmsa::pt_total; ipt++)
  {
    cy->cd(ipt+1);
    cy->cd(ipt+1)->SetLeftMargin(0.15);
    cy->cd(ipt+1)->SetBottomMargin(0.15);
    cy->cd(ipt+1)->SetTicks(1,1);
    cy->cd(ipt+1)->SetGrid(0,0);

    h_mCounts_ybins[ipt]->SetTitle(Form("%1.3f<p_{T}<%1.3f, 20-60 Centrality",vmsa::ptRawStart[ipt],vmsa::ptRawStop[ipt]));
    h_mCounts_ybins[ipt]->GetXaxis()->SetTitle("y");
    h_mCounts_ybins[ipt]->GetYaxis()->SetTitle("Yields");
    h_mCounts_ybins[ipt]->SetMarkerStyle(20);
    h_mCounts_ybins[ipt]->SetMarkerColor(kOrange+7);
    h_mCounts_ybins[ipt]->SetLineColor(kOrange+7);
    h_mCounts_RC_ybins[ipt]->SetMarkerStyle(20);
    h_mCounts_RC_ybins[ipt]->SetMarkerColor(kBlue);
    h_mCounts_RC_ybins[ipt]->SetLineColor(kBlue);

    double max = h_mCounts_ybins[ipt]->GetMaximum();
    double min = h_mCounts_ybins[ipt]->GetMinimum();
    if(h_mCounts_RC_ybins[ipt]->GetMaximum() > max) max = h_mCounts_RC_ybins[ipt]->GetMaximum();
    if(h_mCounts_RC_ybins[ipt]->GetMinimum() < min) min = h_mCounts_RC_ybins[ipt]->GetMinimum();
 
    TLegend *leg = new TLegend(0.6,0.5,0.8,0.7);
    leg->AddEntry(h_mCounts_ybins[ipt],"|DCA| < 2 cm","p");
    leg->AddEntry(h_mCounts_RC_ybins[ipt],"|DCA| < 3 cm","p");

    h_mCounts_ybins[ipt]->Draw("pE");
    h_mCounts_RC_ybins[ipt]->Draw("pE same");
    leg->Draw("same");  
  }
  cy->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_ptspectra_ybincomp_%s_Order%d_%s_DATADCA.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),correction.c_str(),order,etamode.c_str()));  

  for(int ipt = 0; ipt < vmsa::pt_total; ipt++)
  {
    cy->cd(ipt+1);
    h_mCounts_Ratio_ybins[ipt]->SetTitle(Form("%1.3f<p_{T}<%1.3f, 20-60 Centrality",vmsa::ptRawStart[ipt],vmsa::ptRawStop[ipt]));
    h_mCounts_Ratio_ybins[ipt]->GetXaxis()->SetTitle("y");
    h_mCounts_Ratio_ybins[ipt]->GetYaxis()->SetTitle("Data/RC");
    h_mCounts_Ratio_ybins[ipt]->GetYaxis()->SetRangeUser(0.6,1.4);
    h_mCounts_Ratio_ybins[ipt]->Draw("pE");  
  }
  cy->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_ptspectraratio_ybincomp_%s_Order%d_%s_DATADCA.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),correction.c_str(),order,etamode.c_str()));  


}


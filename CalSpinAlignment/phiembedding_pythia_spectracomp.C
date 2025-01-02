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

void phiembedding_pythia_spectracomp(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 40, string sim_py = "Rc5", string sim_em = "Rc8")//defaultF = 0 is BESII, defaultF = 1 is BESI
{
 
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
  
  if(yspectra == 35) spectra = "PhiEmbeddingWithWeights";
  if(yspectra == 36) spectra = "PhiEmbeddingWithWeights_PIDEff_20240722";
  if(yspectra == 37) spectra = "PhiEmbeddingNoWeights";
  if(yspectra == 38) spectra = "PhiEmbeddingWithWeights_PIDEff_20240725";
  if(yspectra == 39) spectra = "PhiEmbeddingWithWeights_PIDEff_20240726";
  

  if(yspectra == 40) spectra = "PhiEmbeddingNoWeights_PIDEff_20240730";

  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  string inputfile = Form("effaccfiles/%s/%s/%s/Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_Weights_PIDEff_2.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_WithWeights_PIDEff_20240722_1.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_NoWeights_PIDEff_20240723_1.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_NoWeights_PIDEff_20240724_fourMomentum_1.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  //string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_WithWeights_PIDEff_20240725_FixedMC_0.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
  string inputfileRC = Form("effaccfiles/%s/%s/%s/PhiEmbed_NoWeights_20240731_Full.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str());
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
  
  cout << inputfile << endl;
  cout << inputfileRC << endl;
 
  string KEY_Pythia = Form("h_m%sEffPtY_Cent_%d",sim_py.c_str(),i_cent);
  string KEY_Embedding = Form("h_m%sEffPtY_Cent_%d",sim_em.c_str(),i_cent);
  cout << KEY_Pythia << endl;
  cout << KEY_Embedding << endl;
  string KEY_counts = Form("h_mRcEffPtY_Cent_%d_Pythia",i_cent);
  h_mCounts = (TH2D*)((TH2D*) File_InPut->Get(KEY_Pythia.c_str()))->Clone(KEY_counts.c_str());
  cout << "Loaded Pythia" << endl;
  string KEY_counts_RC = Form("h_mRcEffPtY_Cent_%d_Embedding",i_cent);
  h_mCounts_RC = (TH2D*)((TH2D*) File_InPutRC->Get(KEY_Embedding.c_str()))->Clone(KEY_counts_RC.c_str());
  cout << "Loaded Embedding" << endl;
 
  cout << "Loaded Histograms" << endl;

  h_mCounts->Print();
  h_mCounts_RC->Print();

  ///////////// 2D y vs pT /////////////
  //float data = h_mcounts->integral(0,-1,0,-1);
  //float rc   = h_mcounts_rc->integral(0,-1,0,-1);
  cout << "NBinsX = " << h_mCounts_RC->GetNbinsX() << endl;
  cout << "NBinsY = " << h_mCounts_RC->GetNbinsY() << endl;
  int nbinsx = h_mCounts_RC->GetNbinsX();
  int nbinsy = h_mCounts_RC->GetNbinsY();
  float data = h_mCounts->Integral(1,nbinsx,1,nbinsy);
  float rc   = h_mCounts_RC->Integral(1,nbinsx,1,nbinsy);

  h_mCounts_RC->Scale(data/rc);

  h_mCounts_Ratio = (TH2D*) h_mCounts->Clone();
  h_mCounts_Ratio->Divide(h_mCounts_RC);
  ///////////////////////////////////////

  ///////////// 1D pT /////////////
  h_mCounts_pt = (TH1D*) h_mCounts->ProjectionX("pt",1,nbinsy,"e"); 
  h_mCounts_RC_pt = (TH1D*) h_mCounts_RC->ProjectionX("ptrc",1,nbinsy,"e"); 
  //float datapt = h_mCounts_pt->Integral(0,-1);
  //float rcpt   = h_mCounts_RC_pt->Integral(0,-1);

  //h_mCounts_RC_pt->Scale(datapt/rcpt);

  h_mCounts_Ratio_pt = (TH1D*) h_mCounts_pt->Clone();
  h_mCounts_Ratio_pt->Divide(h_mCounts_RC_pt);
  ///////////////////////////////////////
  
  ///////////// 1D y /////////////
  h_mCounts_y = (TH1D*) h_mCounts->ProjectionY("y",1,nbinsx,"e"); 
  h_mCounts_RC_y = (TH1D*) h_mCounts_RC->ProjectionY("yrc",1,nbinsx,"e"); 
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
  
  h_mCounts_Ratio->SetTitle(Form("y vs pT ratio Centrality %d",i_cent));
  h_mCounts_Ratio->GetXaxis()->SetTitle("p_{T}");
  h_mCounts_Ratio->GetYaxis()->SetTitle("y");
  h_mCounts_Ratio->Draw("colz");
  
  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/PythiaEmbeddingRatio_ptyspectra2D_%s_Order%d_%s_Cent%d_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra.c_str(),correction.c_str(),order,etamode.c_str(),i_cent,sim_py.c_str(),sim_em.c_str()));  

  h_mCounts_Ratio_pt->SetTitle(Form("Centrality %d",i_cent));
  h_mCounts_Ratio_pt->GetXaxis()->SetTitle("p_{T}");
  h_mCounts_Ratio_pt->GetYaxis()->SetTitle("Pythia/Embedding");
  h_mCounts_Ratio_pt->Draw("pE");
  
  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/PythiaEmbeddingRatio_ptspectra1D_%s_Order%d_%s_Cent%d_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(), spectra.c_str(),correction.c_str(),order,etamode.c_str(),i_cent,sim_py.c_str(),sim_em.c_str()));  
  TLegend *leg = new TLegend(0.5,0.2,0.7,0.4);

  h_mCounts_pt->SetTitle(Form("Centrality %d",i_cent));
  h_mCounts_pt->GetXaxis()->SetTitle("p_{T}");
  h_mCounts_pt->GetYaxis()->SetTitle("Counts");
  h_mCounts_pt->SetMarkerStyle(20);
  h_mCounts_pt->SetMarkerColor(kOrange+7);
  h_mCounts_RC_pt->SetMarkerStyle(20);
  h_mCounts_RC_pt->SetMarkerColor(kBlue);
  double max = h_mCounts_pt->GetMaximum();
  double min = h_mCounts_pt->GetMinimum();
  if(h_mCounts_RC_pt->GetMaximum() > max) max = h_mCounts_RC_pt->GetMaximum();
  if(h_mCounts_RC_pt->GetMinimum() < min) min = h_mCounts_RC_pt->GetMinimum();
  cout << "Max = " << max << ", min = " << min << endl;
  h_mCounts_pt->GetYaxis()->SetRangeUser(min*0.9,max*1.1);
  h_mCounts_pt->Draw("pE");
  h_mCounts_RC_pt->Draw("pE same");
  leg->AddEntry(h_mCounts_pt,"Pythia","p");
  leg->AddEntry(h_mCounts_RC_pt,"Embedding","p");
  leg->Draw("same"); 
 
  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/PythiaEmbedding_ptspectra1D_%s_Order%d_%s_Cent%d_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(), spectra.c_str(),correction.c_str(),order,etamode.c_str(),i_cent,sim_py.c_str(),sim_em.c_str()));  

  h_mCounts_Ratio_y->SetTitle(Form("Centrality %d",i_cent));
  h_mCounts_Ratio_y->GetXaxis()->SetTitle("y");
  h_mCounts_Ratio_y->GetYaxis()->SetTitle("Pythia/Embedding");
  h_mCounts_Ratio_y->Draw("pE");
  
  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/PythiaEmbeddingRatio_yspectra1D_%s_Order%d_%s_Cent%d_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),  spectra.c_str(),correction.c_str(),order,etamode.c_str(),i_cent,sim_py.c_str(),sim_em.c_str()));  

  h_mCounts_y->SetTitle(Form("Centrality %d",i_cent));
  h_mCounts_y->GetXaxis()->SetTitle("y");
  h_mCounts_y->GetYaxis()->SetTitle("Counts");
  h_mCounts_y->SetMarkerStyle(20);
  h_mCounts_y->SetMarkerColor(kOrange+7);
  h_mCounts_RC_y->SetMarkerStyle(20);
  h_mCounts_RC_y->SetMarkerColor(kBlue);
  max = h_mCounts_y->GetMaximum();
  min = h_mCounts_y->GetMinimum();
  if(h_mCounts_RC_y->GetMaximum() > max) max = h_mCounts_RC_y->GetMaximum();
  if(h_mCounts_RC_y->GetMinimum() < min) min = h_mCounts_RC_y->GetMinimum();
  cout << "Max = " << max << ", min = " << min << endl;
  h_mCounts_y->GetYaxis()->SetRangeUser(min*0.9,max*1.1);
  h_mCounts_y->Draw("pE");
  h_mCounts_RC_y->Draw("pE same");
  leg->Draw("same"); 
  
  c->SaveAs(Form("figures/%s/%s/pTstudy/%s/PythiaEmbedding_yspectra1D_%s_Order%d_%s_Cent%d_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(), spectra.c_str(),correction.c_str(),order,etamode.c_str(),i_cent,sim_py.c_str(),sim_em.c_str()));  

  string outputfile = "ptyspectra_datarcratio_19GeV.root";
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_mCounts_Ratio->SetName("pty_datarcratio");
  h_mCounts_Ratio->Write();  
  File_OutPut->Close();  

}


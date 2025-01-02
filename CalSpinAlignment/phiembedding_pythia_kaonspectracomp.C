#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TProfile2D.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "phi_data_constants_19GeV.h"
//#ifdef MAKECINT
//#pragma link C++ class std::map<std::string,TH1F*>+;
//#endif

#ifndef _PlotQA_
#define _PlotQA_  1
#endif

#ifndef _SaveQA_
#define _SaveQA_  0
#endif

using namespace std;

void phiembedding_pythia_kaonspectracomp(int energy = 4, int pid = 0, int year = 0, string date = "20240715", bool random3D = false, int order = 2, string etamode = "eta1_eta1", int deltaonly = 0, string sim_py = "rc5", string sim_em = "rc8")
{
  std::string EP[2] = {"","2nd"};
  TGaxis::SetMaxDigits(4);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);


  string folder = "PhiEmbeddingNoWeights_PIDEff_20240730";
  string InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",folder.c_str());  
  TFile *File_RC_Pythia = TFile::Open(InPutFile_RC.c_str());

  InPutFile_RC = Form("effaccfiles/Phi/19GeV/%s/PhiEmbed_NoWeights_PIDEff_20240730_FixedMC_FixedAngle_0.root",folder.c_str());  
  TFile *File_RC_Embed = TFile::Open(InPutFile_RC.c_str());

  cout << "Loaded All Files" << endl;

  TH2FMap h_mMass_RC;

  string KEY_minus_RC_py = Form("%s_kminus_cent9",sim_py.c_str());
  string KEY_minus_RC_em = Form("%s_kminus_cent9",sim_em.c_str());
  string KEY_minus_RC_Pythia = Form("kminus_Pythia");
  string KEY_minus_RC_Embed = Form("kminus_Embed");
  h_mMass_RC[KEY_minus_RC_Pythia] = (TH2F*)((TH2F*)File_RC_Pythia->Get(KEY_minus_RC_py.c_str()))->Clone(KEY_minus_RC_Pythia.c_str());
  h_mMass_RC[KEY_minus_RC_Embed]  = (TH2F*)((TH2F*)File_RC_Embed->Get(KEY_minus_RC_em.c_str()))->Clone(KEY_minus_RC_Embed.c_str());

  string KEY_plus_RC_py = Form("%s_kplus_cent9",sim_py.c_str());
  string KEY_plus_RC_em = Form("%s_kplus_cent9",sim_em.c_str());
  string KEY_plus_RC_Pythia = Form("kplus_Pythia");
  string KEY_plus_RC_Embed = Form("kplus_Embed");
  h_mMass_RC[KEY_plus_RC_Pythia] = (TH2F*)((TH2F*)File_RC_Pythia->Get(KEY_plus_RC_py.c_str()))->Clone(KEY_plus_RC_Pythia.c_str());
  h_mMass_RC[KEY_plus_RC_Embed]  = (TH2F*)((TH2F*)File_RC_Embed->Get(KEY_plus_RC_em.c_str()))->Clone(KEY_plus_RC_Embed.c_str());

  cout << "Loaded All Hists" << endl;

  TH2FMap h_mMass_phiy;   
  TH2FMap h_mRcMass_phiy;   
  TH1FMap h_mMass_y;   
  TH1FMap h_mRcMass_y;

  TH2FMap h_mMass_pty;
  TH2FMap h_mRcMass_pty;
  TH2FMap h_mRatio_pty;
  TH1FMap h_mMass_pt;
  TH1FMap h_mRcMass_pt;
  TH1FMap h_mRatio_pt;
  TH1FMap h_mRatio_y;

  int rebinconstant = 3;

  h_mMass_RC[KEY_plus_RC_Pythia]->RebinY(2);
  h_mMass_RC[KEY_plus_RC_Pythia]->RebinX(rebinconstant);
  h_mMass_RC[KEY_plus_RC_Embed]->RebinY(2);
  h_mMass_RC[KEY_plus_RC_Embed]->RebinX(rebinconstant);
  h_mMass_RC[KEY_minus_RC_Pythia]->RebinY(2);
  h_mMass_RC[KEY_minus_RC_Pythia]->RebinX(rebinconstant);
  h_mMass_RC[KEY_minus_RC_Embed]->RebinY(2);
  h_mMass_RC[KEY_minus_RC_Embed]->RebinX(rebinconstant);

  cout << "Rebinned" << endl;

  int nbx = h_mMass_RC[KEY_plus_RC_Pythia]->GetNbinsX();
  int nby = h_mMass_RC[KEY_plus_RC_Pythia]->GetNbinsY();

  double yields_plus_Pythia = h_mMass_RC[KEY_plus_RC_Pythia]->Integral(1,nbx,1,nby);
  double yields_minus_Pythia = h_mMass_RC[KEY_minus_RC_Pythia]->Integral(1,nbx,1,nby);
  double yields_plus_Embed  = h_mMass_RC[KEY_plus_RC_Embed]->Integral(1,nbx,1,nby);
  double yields_minus_Embed  = h_mMass_RC[KEY_minus_RC_Embed]->Integral(1,nbx,1,nby);
 
  h_mMass_RC[KEY_plus_RC_Embed]->Scale(yields_plus_Pythia/yields_plus_Embed);
  h_mMass_RC[KEY_minus_RC_Embed]->Scale(yields_minus_Pythia/yields_minus_Embed);

  string KEY_plus_ratio_pty = Form("kplus_ratio_pty");   
  h_mRatio_pty[KEY_plus_ratio_pty] = (TH2F*) h_mMass_RC[KEY_plus_RC_Pythia]->Clone(KEY_plus_ratio_pty.c_str());
  h_mRatio_pty[KEY_plus_ratio_pty]->Divide(h_mMass_RC[KEY_plus_RC_Embed]);

  string KEY_minus_ratio_pty = Form("kminus_ratio_pty");   
  h_mRatio_pty[KEY_minus_ratio_pty] = (TH2F*) h_mMass_RC[KEY_minus_RC_Pythia]->Clone(KEY_minus_ratio_pty.c_str());
  h_mRatio_pty[KEY_minus_ratio_pty]->Divide(h_mMass_RC[KEY_minus_RC_Embed]);
  
  string KEY_plus_RC_pt_Pythia = Form("kplus_pt_Pythia");   
  string KEY_plus_RC_pt_Embed  = Form("kplus_pt_Embed");   
  string KEY_plus_ratio_pt = Form("kplus_ratio_pt");   
  h_mRcMass_pt[KEY_plus_RC_pt_Pythia] = (TH1F*) h_mMass_RC[KEY_plus_RC_Pythia]->ProjectionX(KEY_plus_RC_pt_Pythia.c_str(),1,nby,"e");
  h_mRcMass_pt[KEY_plus_RC_pt_Embed]  = (TH1F*) h_mMass_RC[KEY_plus_RC_Embed]->ProjectionX(KEY_plus_RC_pt_Embed.c_str(),1,nby,"e");
  h_mRatio_pt[KEY_plus_ratio_pt] = (TH1F*) h_mRcMass_pt[KEY_plus_RC_pt_Pythia]->Clone(KEY_plus_ratio_pt.c_str());
  h_mRatio_pt[KEY_plus_ratio_pt]->Divide(h_mRcMass_pt[KEY_plus_RC_pt_Embed]);
  
  string KEY_plus_RC_y_Pythia = Form("kplus_y_Pythia");   
  string KEY_plus_RC_y_Embed  = Form("kplus_y_Embed");   
  string KEY_plus_ratio_y = Form("kplus_ratio_y");   
  h_mRcMass_y[KEY_plus_RC_y_Pythia] = (TH1F*) h_mMass_RC[KEY_plus_RC_Pythia]->ProjectionY(KEY_plus_RC_y_Pythia.c_str(),1,nbx,"e");
  h_mRcMass_y[KEY_plus_RC_y_Embed]  = (TH1F*) h_mMass_RC[KEY_plus_RC_Embed]->ProjectionY(KEY_plus_RC_y_Embed.c_str(),1,nbx,"e");
  h_mRatio_y[KEY_plus_ratio_y] = (TH1F*) h_mRcMass_y[KEY_plus_RC_y_Pythia]->Clone(KEY_plus_ratio_y.c_str());
  h_mRatio_y[KEY_plus_ratio_y]->Divide(h_mRcMass_y[KEY_plus_RC_y_Embed]);

  string KEY_minus_RC_pt_Pythia = Form("kminus_pt_Pythia");   
  string KEY_minus_RC_pt_Embed  = Form("kminus_pt_Embed");   
  string KEY_minus_ratio_pt = Form("kminus_ratio_pt");   
  h_mRcMass_pt[KEY_minus_RC_pt_Pythia] = (TH1F*) h_mMass_RC[KEY_minus_RC_Pythia]->ProjectionX(KEY_minus_RC_pt_Pythia.c_str(),1,nby,"e");
  h_mRcMass_pt[KEY_minus_RC_pt_Embed]  = (TH1F*) h_mMass_RC[KEY_minus_RC_Embed]->ProjectionX(KEY_minus_RC_pt_Embed.c_str(),1,nby,"e");
  h_mRatio_pt[KEY_minus_ratio_pt] = (TH1F*) h_mRcMass_pt[KEY_minus_RC_pt_Pythia]->Clone(KEY_minus_ratio_pt.c_str());
  h_mRatio_pt[KEY_minus_ratio_pt]->Divide(h_mRcMass_pt[KEY_minus_RC_pt_Embed]);
  
  string KEY_minus_RC_y_Pythia = Form("kminus_y_Pythia");   
  string KEY_minus_RC_y_Embed  = Form("kminus_y_Embed");   
  string KEY_minus_ratio_y = Form("kminus_ratio_y");   
  h_mRcMass_y[KEY_minus_RC_y_Pythia] = (TH1F*) h_mMass_RC[KEY_minus_RC_Pythia]->ProjectionY(KEY_minus_RC_y_Pythia.c_str(),1,nbx,"e");
  h_mRcMass_y[KEY_minus_RC_y_Embed]  = (TH1F*) h_mMass_RC[KEY_minus_RC_Embed]->ProjectionY(KEY_minus_RC_y_Embed.c_str(),1,nbx,"e");
  h_mRatio_y[KEY_minus_ratio_y] = (TH1F*) h_mRcMass_y[KEY_minus_RC_y_Pythia]->Clone(KEY_minus_ratio_y.c_str());
  h_mRatio_y[KEY_minus_ratio_y]->Divide(h_mRcMass_y[KEY_minus_RC_y_Embed]);

  cout << "Created all Hists, Time to Plot" << endl;

  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,400);     
  c1->Divide(2,1);
  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    c1->cd(i+1)->SetLeftMargin(0.15);
    c1->cd(i+1)->SetBottomMargin(0.15);
    c1->cd(i+1)->SetTicks(1,1);
    c1->cd(i+1)->SetGrid(0,0);  
  } 

  string histname[2] = {"kplus","kminus"};
  string histlabel[2] = {"K^{+}","K^{-}"};
  /////////////////////////////// No Cos2PhiStarPhi separation //////////////////////////////
  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_RC_pty = Form("%s_Pythia",histname[i].c_str());   
    h_mMass_RC[KEY_RC_pty]->SetTitle(Form("%s Pythia %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mMass_RC[KEY_RC_pty]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mMass_RC[KEY_RC_pty]->GetYaxis()->SetTitle("y");

    h_mMass_RC[KEY_RC_pty]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pty_Pythia_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim_py.c_str(),sim_em.c_str()));

  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_RC_pty = Form("%s_Embed",histname[i].c_str());   
    h_mMass_RC[KEY_RC_pty]->SetTitle(Form("%s Embedding  %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mMass_RC[KEY_RC_pty]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mMass_RC[KEY_RC_pty]->GetYaxis()->SetTitle("y");

    h_mMass_RC[KEY_RC_pty]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pty_Embedding_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim_py.c_str(),sim_em.c_str()));

  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_ratio_pty = Form("%s_ratio_pty",histname[i].c_str());   
    h_mRatio_pty[KEY_ratio_pty]->SetTitle(Form("%s Pythia/Embedding %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mRatio_pty[KEY_ratio_pty]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mRatio_pty[KEY_ratio_pty]->GetYaxis()->SetTitle("y");

    h_mRatio_pty[KEY_ratio_pty]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_pty_PythiaEbmeddginRatio_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim_py.c_str(),sim_em.c_str()));

  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_ratio_pt = Form("%s_ratio_pt",histname[i].c_str());   
    h_mRatio_pt[KEY_ratio_pt]->SetTitle(Form("%s Pythia/Embedding %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mRatio_pt[KEY_ratio_pt]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mRatio_pt[KEY_ratio_pt]->GetYaxis()->SetTitle("Pythia/Embedding");

    h_mRatio_pt[KEY_ratio_pt]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_ptonly_PythiaEmbeddingRatio_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim_py.c_str(),sim_em.c_str()));

  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_ratio_y = Form("%s_ratio_y",histname[i].c_str());   
    h_mRatio_y[KEY_ratio_y]->SetTitle(Form("%s Pythia/Embedding %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mRatio_y[KEY_ratio_y]->GetXaxis()->SetTitle("y");
    h_mRatio_y[KEY_ratio_y]->GetYaxis()->SetTitle("Pythia/Embedding");

    h_mRatio_y[KEY_ratio_y]->Draw("colz");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_yonly_PythiaEmbeddingRatio_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim_py.c_str(),sim_em.c_str()));

  TLegend *legy = new TLegend(0.4,0.2,0.6,0.4);
  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    string KEY_y = Form("%s_y_Pythia",histname[i].c_str());   
    string KEY_RC_y = Form("%s_y_Embed",histname[i].c_str());   
    h_mRcMass_y[KEY_y]->SetTitle(Form("%s %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mRcMass_y[KEY_y]->GetXaxis()->SetTitle("y");
    h_mRcMass_y[KEY_y]->GetYaxis()->SetTitle("Counts");
    h_mRcMass_y[KEY_y]->SetMarkerStyle(20);
    h_mRcMass_y[KEY_y]->SetMarkerColor(kOrange+7);
    h_mRcMass_y[KEY_y]->SetLineColor(kOrange+7);

    int min = h_mRcMass_y[KEY_y]->GetMinimum();
    int max = h_mRcMass_y[KEY_y]->GetMaximum();
    if(h_mRcMass_y[KEY_RC_y]->GetMinimum() < min) min = h_mRcMass_y[KEY_RC_y]->GetMinimum();
    if(h_mRcMass_y[KEY_RC_y]->GetMaximum() > max) max = h_mRcMass_y[KEY_RC_y]->GetMaximum();
    h_mRcMass_y[KEY_y]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);

    h_mRcMass_y[KEY_y]->Draw("pE");

    h_mRcMass_y[KEY_RC_y]->SetMarkerStyle(24);
    h_mRcMass_y[KEY_RC_y]->SetMarkerColor(kBlack);
    h_mRcMass_y[KEY_RC_y]->SetLineColor(kBlack);
    h_mRcMass_y[KEY_RC_y]->Draw("pE same");

    if(i == 0)
    {
      legy->AddEntry(h_mRcMass_y[KEY_y],"Pythia","p");
      legy->AddEntry(h_mRcMass_y[KEY_RC_y],"Embedding","p");
    }
    legy->Draw("same");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_yonly_PythiaEmbedding_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim_py.c_str(),sim_em.c_str()));
   
  TLegend *legpt = new TLegend(0.4,0.2,0.6,0.4);
  for(int i = 0; i < 2; i++)
  {
    c1->cd(i+1);
    c1->cd(i+1)->SetLogy();
    string KEY_pt = Form("%s_pt_Pythia",histname[i].c_str());   
    string KEY_RC_pt = Form("%s_pt_Embed",histname[i].c_str());   
    h_mRcMass_pt[KEY_pt]->SetTitle(Form("%s %1.1f<p_{T}<%1.1f, 20-60 Centrality",histlabel[i].c_str(),vmsa::pt_low[energy][2],vmsa::pt_up[energy][5]));
    h_mRcMass_pt[KEY_pt]->GetXaxis()->SetTitle("p_{T} GeV/c");
    h_mRcMass_pt[KEY_pt]->GetYaxis()->SetTitle("Counts");
    h_mRcMass_pt[KEY_pt]->SetMarkerStyle(20);
    h_mRcMass_pt[KEY_pt]->SetMarkerColor(kOrange+7);
    h_mRcMass_pt[KEY_pt]->SetLineColor(kOrange+7);

    int min = h_mRcMass_pt[KEY_pt]->GetMinimum();
    int max = h_mRcMass_pt[KEY_pt]->GetMaximum();
    if(h_mRcMass_pt[KEY_RC_pt]->GetMinimum() < min) min = h_mRcMass_pt[KEY_RC_pt]->GetMinimum();
    if(h_mRcMass_pt[KEY_RC_pt]->GetMaximum() > max) max = h_mRcMass_pt[KEY_RC_pt]->GetMaximum();
    h_mRcMass_pt[KEY_pt]->GetYaxis()->SetRangeUser(1/*min*0.9*/,max*2.0);

    h_mRcMass_pt[KEY_pt]->Draw("pE");

    h_mRcMass_pt[KEY_RC_pt]->SetMarkerStyle(24);
    h_mRcMass_pt[KEY_RC_pt]->SetMarkerColor(kBlack);
    h_mRcMass_pt[KEY_RC_pt]->SetLineColor(kBlack);
    h_mRcMass_pt[KEY_RC_pt]->Draw("pE same");

    if(i == 0)
    {
      legpt->AddEntry(h_mRcMass_pt[KEY_pt],"Pythia","p");
      legpt->AddEntry(h_mRcMass_pt[KEY_RC_pt],"Embedding","p");
    }
    legpt->Draw("same");
  }
  c1->SaveAs(Form("figures/%s/%s/pTstudy/%s/Cos2PhiStarPhi_KaonDistribution_ptonly_PythiaEmbedding_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),folder.c_str(),sim_py.c_str(),sim_em.c_str()));
  //File_OutPut->Close();
}

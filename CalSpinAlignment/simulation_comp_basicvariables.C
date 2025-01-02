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

void simulation_comp_basicvariables(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 40, string sim_py = "Rc5", string sim_em = "Rc8", bool dataon = true)//defaultF = 0 is BESII, defaultF = 1 is BESI
{

  const int nopts = 2; 
  const int offset = 8;
  string spectra[nopts] = { "../../../../data",
                            "PhiEmbedding_noweights_nomccuts"
                            };
 
  string filename[nopts] = {"Yields_Phi_SE_19GeV_20240828_cent2060_kaonspectra_basicinfo_fpoly_fixedcharge.root",
                            "PhiEmbedding_noweights_nomccuts.root"};
 
  string filenameME = "Yields_Phi_ME_19GeV_20240828_cent2060_kaonspectra_basicinfo_fpoly_fixedcharge.root";

  string label[nopts] = {"Data","Embed"};  
  string filelabel[nopts] = {"Data","Embedding"};
  
  //const int nopts = 11; 
  //string spectra[nopts] = { "DataVariables",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1"};
 
  //string filename[nopts] = {"DataVariables_19GeV.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root"};
 
  //string label[nopts] = {"Data","MC #phi-meson","MC (K^{+}+K^{-})","MC with RC","RC","+rapidity cut","+|#eta| cut","+TPC Eff","+ToF Match","+m^{2} PID Eff","+n#sigma_{K} PID Eff"};  
  //string filelabel[nopts] = {"DATA","RC0","RC1","RC2","RC3","RC4","RC5","RC6","RC7","RC8","RC9"};
  
  int color[14] = {kBlack,kOrange+7,kBlue,kGray+2,kRed,kGreen,kViolet,kBlack,kOrange+7,kBlue,kGray+2,kRed,kGreen,kViolet};
  int marker[14] = {20,20,20,20,24,24,24,21,21,21,21,25,25,25};
 
  float integralevent[2] = {0.0};
  float integral[2][nopts] = {0.0};

  // Vertex
  TH3F *h_mVertex[2];  
  TH1F *h_mVx[2];
  TH1F *h_mVy[2];
  TH1F *h_mVz[2];
  TH1F *h_mVx_Ratio;
  TH1F *h_mVy_Ratio;
  TH1F *h_mVz_Ratio;
  
  // ToF Match
  TH1F *h_mNToFMatch[2];
  TH1F *h_mNToFMatch_Ratio;
  
  // RefMult
  TH1F *h_mRefMult[2];
  TH1F *h_mRefMult_Ratio;

  // Dca
  TH1F *h_mDca[2][nopts];
  TH1F *h_mDca_Ratio[2][nopts];

  // NHits
  TH1F *h_mNHits[2][nopts];
  TH1F *h_mNHits_Ratio[2][nopts];

  // NHitsRatio
  TH1F *h_mNHitsRatio[2][nopts];
  TH1F *h_mNHitsRatio_Ratio[2][nopts];

  // NHitsRatio
  TH2F *h_mDEdx[2][nopts];
  TH2F *h_mDEdx_Ratio[2][nopts];
  TH1F *h_mDEdx_1D[2][nopts];
  TH1F *h_mDEdx_1D_Ratio[2][nopts];

  string charge[2] = {"plus","minus"};

  for(int ifile = 0; ifile < nopts; ifile++)
  {
    string inputfile = Form("effaccfiles/%s/%s/%s/%s",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[ifile].c_str(),filename[ifile].c_str());
    TFile *File_Input = TFile::Open(inputfile.c_str());
    
    string inputfileME = Form("effaccfiles/%s/%s/%s/%s",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[ifile].c_str(),filenameME.c_str());
    TFile *File_InputME;// = TFile::Open(inputfileME.c_str());
    if(ifile == 0 && dataon) File_InputME = TFile::Open(inputfileME.c_str());


    if(ifile < 2) 
    {
      string HistName = Form("h_mVertex_Cent_%d",i_cent);
      string HistNameRename = Form("h_mVertex_Cent_%d_%d",i_cent,ifile+offset);
      cout << HistName << endl;
      h_mVertex[ifile] = (TH3F*) ( (TH1F*) File_Input->Get(HistName.c_str()) )->Clone(HistNameRename.c_str());
      integralevent[ifile] = h_mVertex[ifile]->Integral(0,-1,0,-1,0,-1);
      if(ifile == 1) h_mVertex[ifile]->Scale(integralevent[0]/integralevent[1]);
      string NameVx = Form("vx_%d",ifile);
      string NameVy = Form("vy_%d",ifile);
      string NameVz = Form("vz_%d",ifile);
      h_mVx[ifile] = (TH1F*) h_mVertex[ifile]->ProjectionX(NameVx.c_str(),0,-1,0,-1,"e");
      h_mVy[ifile] = (TH1F*) h_mVertex[ifile]->ProjectionY(NameVy.c_str(),0,-1,0,-1,"e");
      h_mVz[ifile] = (TH1F*) h_mVertex[ifile]->ProjectionZ(NameVz.c_str(),0,-1,0,-1,"e");
    
   
      for(int i = 1; i <= h_mVx[ifile]->GetNbinsX(); i++)
      {
        cout << "ifile = " << ifile << ", bin = " << i << ", [" << h_mVx[ifile]->GetBinLowEdge(i) << "," << h_mVx[ifile]->GetBinLowEdge(i)+h_mVx[ifile]->GetBinWidth(i) << "}, nentries = " << h_mVx[ifile]->GetBinContent(i) << endl;
      } 
      for(int i = 1; i <= h_mVy[ifile]->GetNbinsX(); i++)
      {
        cout << "ifile = " << ifile << ", bin = " << i << ", [" << h_mVy[ifile]->GetBinLowEdge(i) << "," << h_mVy[ifile]->GetBinLowEdge(i)+h_mVy[ifile]->GetBinWidth(i) << "}, nentries = " << h_mVy[ifile]->GetBinContent(i) << endl;
      } 
      for(int i = 1; i <= h_mVz[ifile]->GetNbinsX(); i++)
      {
        cout << "ifile = " << ifile << ", bin = " << i << ", [" << h_mVz[ifile]->GetBinLowEdge(i) << "," << h_mVz[ifile]->GetBinLowEdge(i)+h_mVz[ifile]->GetBinWidth(i) << "}, nentries = " << h_mVz[ifile]->GetBinContent(i) << endl;
      } 


      h_mVz[ifile]->Rebin(4);

      if(ifile == 1) 
      { 
        h_mVx_Ratio = (TH1F*) h_mVx[0]->Clone();
        h_mVx_Ratio->Divide(h_mVx[1]);
        h_mVy_Ratio = (TH1F*) h_mVy[0]->Clone();
        h_mVy_Ratio->Divide(h_mVy[1]);
        h_mVz_Ratio = (TH1F*) h_mVz[0]->Clone();
        h_mVz_Ratio->Divide(h_mVz[1]);
      }     
 
      HistName = Form("h_mNToFMatch_Cent_%d",i_cent);
      cout << HistName << endl;
      HistNameRename = Form("h_mNToFMatch_Cent_%d_%d",i_cent,ifile+offset);
      h_mNToFMatch[ifile] = (TH1F*) ( (TH1F*) File_Input->Get(HistName.c_str()) )->Clone(HistNameRename.c_str());
      if(ifile == 1)
      {
        h_mNToFMatch[1]->Scale(integralevent[0]/integralevent[1]); 
        h_mNToFMatch_Ratio = (TH1F*) h_mNToFMatch[0]->Clone();
        h_mNToFMatch_Ratio->Divide(h_mNToFMatch[1]);
      }

      HistName = Form("h_mRefMult_Cent_%d",i_cent);
      cout << HistName << endl;
      HistNameRename = Form("h_mRefMult_Cent_%d_%d",i_cent,ifile+offset);
      h_mRefMult[ifile] = (TH1F*) ( (TH1F*) File_Input->Get(HistName.c_str()) )->Clone(HistNameRename.c_str());
      if(ifile == 1)
      {
        h_mRefMult[1]->Scale(integralevent[0]/integralevent[1]); 
        h_mRefMult_Ratio = (TH1F*) h_mRefMult[0]->Clone();
        h_mRefMult_Ratio->Divide(h_mRefMult[1]);
      }
    }
    for(int icharge = 0; icharge < 2; ++icharge)
    {
      string HistName = Form("h_mDca_K%s_Cent_9_%d",charge[icharge].c_str(),ifile+offset);
      string HistNameME;
      cout << HistName << endl;
      if(ifile == 0 && dataon) 
      {
        HistName = Form("h_mDca_K%s_Cent_9_SE",charge[icharge].c_str());
        HistNameME = Form("h_mDca_K%s_Cent_9_ME",charge[icharge].c_str());
      }
      h_mDca[icharge][ifile] = (TH1F*) ( (TH1F*) File_Input->Get(HistName.c_str()) )->Clone(); 
      if(ifile == 0 && dataon) h_mDca[icharge][ifile]->Add((TH1F*) File_InputME->Get(HistNameME.c_str()),-1.);
      integral[icharge][ifile] = h_mDca[icharge][ifile]->Integral(0,-1);
      if(ifile > 0)
      {
        h_mDca[icharge][ifile]->Scale(integral[icharge][0]/integral[icharge][ifile]); 
        h_mDca_Ratio[icharge][ifile] = (TH1F*) h_mDca[icharge][0]->Clone();
        h_mDca_Ratio[icharge][ifile]->Divide(h_mDca[icharge][ifile]);
        h_mDca_Ratio[icharge][ifile]->Print();
      }
    
      h_mDca[icharge][ifile]->Print();

      HistName = Form("h_mNHits_K%s_Cent_9_%d",charge[icharge].c_str(),ifile+offset);
      cout << HistName << endl;
      if(ifile == 0 && dataon) 
      {
        HistName = Form("h_mNHits_K%s_Cent_9_SE",charge[icharge].c_str());
        HistNameME = Form("h_mNHits_K%s_Cent_9_ME",charge[icharge].c_str());
      }
      h_mNHits[icharge][ifile] = (TH1F*) ( (TH1F*) File_Input->Get(HistName.c_str()) )->Clone();
      if(ifile == 0 && dataon) h_mNHits[icharge][ifile]->Add((TH1F*) File_InputME->Get(HistNameME.c_str()),-1.);
      if(ifile > 0)
      {
        h_mNHits[icharge][ifile]->Scale(integral[icharge][0]/integral[icharge][ifile]); 
        h_mNHits_Ratio[icharge][ifile] = (TH1F*) h_mNHits[icharge][0]->Clone();
        h_mNHits_Ratio[icharge][ifile]->Divide(h_mNHits[icharge][ifile]);
      }

      HistName = Form("h_mNHitsRatio_K%s_Cent_9_%d",charge[icharge].c_str(),ifile+offset);
      cout << HistName << endl;
      if(ifile == 0 && dataon) 
      {
        HistName = Form("h_mNHitsRatio_K%s_Cent_9_SE",charge[icharge].c_str());
        HistNameME = Form("h_mNHitsRatio_K%s_Cent_9_ME",charge[icharge].c_str());
      }
      h_mNHitsRatio[icharge][ifile] = (TH1F*) ( (TH1F*) File_Input->Get(HistName.c_str()) )->Clone();
      if(ifile == 0 && dataon) h_mNHitsRatio[icharge][ifile]->Add((TH1F*) File_InputME->Get(HistNameME.c_str()),-1.);
      if(ifile > 0)
      {
        h_mNHitsRatio[icharge][ifile]->Scale(integral[icharge][0]/integral[icharge][ifile]); 
        h_mNHitsRatio_Ratio[icharge][ifile] = (TH1F*) h_mNHitsRatio[icharge][0]->Clone();
        h_mNHitsRatio_Ratio[icharge][ifile]->Divide(h_mNHitsRatio[icharge][ifile]);
      }

      HistName = Form("h_mDEdx_K%s_Cent_9_%d",charge[icharge].c_str(),ifile+offset);
      cout << HistName << endl;
      if(ifile == 0 && dataon) 
      {
        HistName = Form("h_mDEdx_K%s_Cent_9_SE",charge[icharge].c_str());
        HistNameME = Form("h_mDEdx_K%s_Cent_9_ME",charge[icharge].c_str());
      }
      h_mDEdx[icharge][ifile] = (TH2F*) ( (TH1F*) File_Input->Get(HistName.c_str()) )->Clone();
      if(ifile == 0 && dataon) h_mDEdx[icharge][ifile]->Add((TH1F*) File_InputME->Get(HistNameME.c_str()),-1.);
      string Name1D = Form("h_mDEdx_1D_K%s_Cent_9_%d",charge[icharge].c_str(),ifile);
      h_mDEdx_1D[icharge][ifile] = (TH1F*) h_mDEdx[icharge][ifile]->ProjectionY(Name1D.c_str(),0,-1,"e");
      if(ifile > 0)
      {
        h_mDEdx_1D[icharge][ifile]->Scale(integral[icharge][0]/integral[icharge][ifile]); 
        h_mDEdx_1D_Ratio[icharge][ifile] = (TH1F*) h_mDEdx_1D[icharge][0]->Clone();
        h_mDEdx_1D_Ratio[icharge][ifile]->Divide(h_mDEdx_1D[icharge][ifile]);
      }
    }
  }

  string histpart[2] = {"K^{-}","K^{+}"};

  TCanvas *cev = new TCanvas("cev","cev",10,10,2000,800);
  cev->Divide(5,2);
  for(int i = 0; i < 10; i++)
  {
    cev->cd(i+1)->SetLeftMargin(0.15);
    cev->cd(i+1)->SetBottomMargin(0.15);
    cev->cd(i+1)->SetTicks(1,1);
    cev->cd(i+1)->SetGrid(0,0);
  } 
 
  cev->cd(1);
  cev->cd(1)->SetLogy();
  TLegend *legVx = new TLegend(0.6,0.6,0.8,0.8);
  double max = -9999999999;
  double min = 9999999999;
  for(int ifile = 0; ifile < 2; ifile++)
  {
    double tmax = h_mVx[ifile]->GetMaximum();
    double tmin = h_mVx[ifile]->GetMinimum();
    if(tmax > max) tmax = max;
    if(tmin < min) tmin = min;
    cout << "tmax = " << tmax << endl;
    cout << "tmin = " << tmin << endl;
  }  
  h_mVx[0]->GetXaxis()->SetTitle("v_{x} (cm)");
  h_mVx[0]->GetYaxis()->SetTitle("Counts");
  //h_mVx[0]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);
  h_mVx[0]->SetMarkerStyle(20);
  h_mVx[0]->SetMarkerColor(kOrange+7);
  h_mVx[0]->SetLineColor(kOrange+7);
  h_mVx[1]->SetMarkerStyle(20);
  h_mVx[1]->SetMarkerColor(kBlue);
  h_mVx[1]->SetLineColor(kBlue);
  legVx->AddEntry(h_mVx[0],"Data","p");
  legVx->AddEntry(h_mVx[1],"Embed","p");
  h_mVx[0]->Draw("pE");
  h_mVx[1]->Draw("pE same");
  legVx->Draw("same");  

  cev->cd(2);
  cev->cd(2)->SetLogy();
  TLegend *legVy = new TLegend(0.6,0.6,0.8,0.8);
  max = -9999999999;
  min = 9999999999;
  for(int ifile = 0; ifile < 2; ifile++)
  {
    double tmax = h_mVy[ifile]->GetMaximum();
    double tmin = h_mVy[ifile]->GetMinimum();
    if(tmax > max) tmax = max;
    if(tmin < min) tmin = min;
  }  
  h_mVy[0]->GetXaxis()->SetTitle("v_{y} (cm)");
  h_mVy[0]->GetYaxis()->SetTitle("Counts");
  //h_mVy[0]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);
  h_mVy[0]->SetMarkerStyle(20);
  h_mVy[0]->SetMarkerColor(kOrange+7);
  h_mVy[0]->SetLineColor(kOrange+7);
  h_mVy[1]->SetMarkerStyle(20);
  h_mVy[1]->SetMarkerColor(kBlue);
  h_mVy[1]->SetLineColor(kBlue);
  legVy->AddEntry(h_mVy[0],"Data","p");
  legVy->AddEntry(h_mVy[1],"Embed","p");
  h_mVy[0]->Draw("pE");
  h_mVy[1]->Draw("pE same");
  legVy->Draw("same");  

  cev->cd(3);
  TLegend *legVz = new TLegend(0.6,0.6,0.8,0.8);
  max = -9999999999;
  min = 9999999999;
  for(int ifile = 0; ifile < 2; ifile++)
  {
    double tmax = h_mVz[ifile]->GetMaximum();
    double tmin = h_mVz[ifile]->GetMinimum();
    if(tmax > max) tmax = max;
    if(tmin < min) tmin = min;
  }  
  h_mVz[0]->GetXaxis()->SetTitle("v_{z} (cm)");
  h_mVz[0]->GetYaxis()->SetTitle("Counts");
  //h_mVz[0]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);
  h_mVz[0]->SetMarkerStyle(20);
  h_mVz[0]->SetMarkerColor(kOrange+7);
  h_mVz[0]->SetLineColor(kOrange+7);
  h_mVz[1]->SetMarkerStyle(20);
  h_mVz[1]->SetMarkerColor(kBlue);
  h_mVz[1]->SetLineColor(kBlue);
  legVz->AddEntry(h_mVz[0],"Data","p");
  legVz->AddEntry(h_mVz[1],"Embed","p");
  h_mVz[0]->Draw("pE");
  h_mVz[1]->Draw("pE same");
  legVz->Draw("same");  

  cev->cd(4);
  TLegend *legNToFMatch = new TLegend(0.6,0.6,0.8,0.8);
  max = -9999999999;
  min = 9999999999;
  for(int ifile = 0; ifile < 2; ifile++)
  {
    double tmax = h_mNToFMatch[ifile]->GetMaximum();
    double tmin = h_mNToFMatch[ifile]->GetMinimum();
    if(tmax > max) tmax = max;
    if(tmin < min) tmin = min;
  }  
  h_mNToFMatch[0]->GetXaxis()->SetTitle("NToFMatch");
  h_mNToFMatch[0]->GetYaxis()->SetTitle("Counts");
  //h_mNToFMatch[0]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);
  h_mNToFMatch[0]->SetMarkerStyle(20);
  h_mNToFMatch[0]->SetMarkerColor(kOrange+7);
  h_mNToFMatch[0]->SetLineColor(kOrange+7);
  h_mNToFMatch[1]->SetMarkerStyle(20);
  h_mNToFMatch[1]->SetMarkerColor(kBlue);
  h_mNToFMatch[1]->SetLineColor(kBlue);
  legNToFMatch->AddEntry(h_mNToFMatch[0],"Data","p");
  legNToFMatch->AddEntry(h_mNToFMatch[1],"Embed","p");
  h_mNToFMatch[0]->Draw("pE");
  h_mNToFMatch[1]->Draw("pE same");
  legNToFMatch->Draw("same");  

  cev->cd(5);
  TLegend *legRefMult = new TLegend(0.6,0.6,0.8,0.8);
  max = -9999999999;
  min = 9999999999;
  for(int ifile = 0; ifile < 2; ifile++)
  {
    double tmax = h_mRefMult[ifile]->GetMaximum();
    double tmin = h_mRefMult[ifile]->GetMinimum();
    if(tmax > max) tmax = max;
    if(tmin < min) tmin = min;
  }  
  h_mRefMult[0]->GetXaxis()->SetTitle("RefMult");
  h_mRefMult[0]->GetYaxis()->SetTitle("Counts");
  //h_mRefMult[0]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);
  h_mRefMult[0]->SetMarkerStyle(20);
  h_mRefMult[0]->SetMarkerColor(kOrange+7);
  h_mRefMult[0]->SetLineColor(kOrange+7);
  h_mRefMult[1]->SetMarkerStyle(20);
  h_mRefMult[1]->SetMarkerColor(kBlue);
  h_mRefMult[1]->SetLineColor(kBlue);
  legRefMult->AddEntry(h_mRefMult[0],"Data","p");
  legRefMult->AddEntry(h_mRefMult[1],"Embed","p");
  h_mRefMult[0]->Draw("pE");
  h_mRefMult[1]->Draw("pE same");
  legRefMult->Draw("same");  

  cev->cd(6);
  h_mVx_Ratio->GetXaxis()->SetTitle("v_{x} (cm)");
  h_mVx_Ratio->GetYaxis()->SetTitle("Data/Embed");
  h_mVx_Ratio->Draw("pE same");

  cev->cd(7);
  h_mVy_Ratio->GetXaxis()->SetTitle("v_{y} (cm)");
  h_mVy_Ratio->GetYaxis()->SetTitle("Data/Embed");
  h_mVy_Ratio->Draw("pE same");

  cev->cd(8);
  h_mVz_Ratio->GetXaxis()->SetTitle("v_{z} (cm)");
  h_mVz_Ratio->GetYaxis()->SetTitle("Data/Embed");
  h_mVz_Ratio->Draw("pE same");

  cev->cd(9);
  h_mNToFMatch_Ratio->GetXaxis()->SetTitle("NToFMatch");
  h_mNToFMatch_Ratio->GetYaxis()->SetTitle("Data/Embed");
  h_mNToFMatch_Ratio->Draw("pE same");

  cev->cd(10);
  h_mRefMult_Ratio->GetXaxis()->SetTitle("RefMult");
  h_mRefMult_Ratio->GetYaxis()->SetTitle("Data/Embed");
  h_mRefMult_Ratio->Draw("pE same");
  
  cev->SaveAs(Form("figures/%s/%s/pTstudy/BasicParameterComparisons/EventParameterComparison_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),i_cent));


  TCanvas *c = new TCanvas("c","c",10,10,800,800);
  c->Divide(2,2);
  for(int ipart = 0; ipart < 4; ipart++) 
  {
    c->cd(ipart+1)->SetLeftMargin(0.15);
    c->cd(ipart+1)->SetBottomMargin(0.15);
    c->cd(ipart+1)->SetTicks(1,1);
    c->cd(ipart+1)->SetGrid(0,0);
  }
  
  cout << "DCA" << endl;
  TLegend *legDca[2];
  TLegend *legDcaRatio[2];
  for(int ipart = 0; ipart < 2; ipart++)
  {
    legDca[ipart] = new TLegend(0.6,0.6,0.8,0.8);
    legDcaRatio[ipart] = new TLegend(0.6,0.6,0.8,0.8);

    max = -9999999999;
    min = 9999999999;
    for(int ifile = 0; ifile < nopts; ifile++)
    {
      double tmax = h_mDca[ipart][ifile]->GetMaximum();
      double tmin = h_mDca[ipart][ifile]->GetMinimum();
      if(tmax > max) tmax = max;
      if(tmin < min) tmin = min;
      cout << "tmax = " << tmax << endl;
      cout << "tmin = " << tmin << endl;
    }  

    c->cd(ipart+1);
    c->cd(ipart+1)->SetLogy();
    for(int ifile = 0; ifile < nopts; ifile++)
    {
      h_mDca[ipart][ifile]->SetTitle(Form("%s",histpart[ipart].c_str()));
      h_mDca[ipart][ifile]->GetXaxis()->SetTitle("|DCA| (cm)");
      h_mDca[ipart][ifile]->GetYaxis()->SetTitle("Counts");
      h_mDca[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mDca[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mDca[ipart][ifile]->SetLineColor(color[ifile]);
      //h_mDca[ipart][ifile]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);
      legDca[ipart]->AddEntry(h_mDca[ipart][ifile],label[ifile].c_str(),"p");
      if(ifile == 0) h_mDca[ipart][ifile]->Draw("pE"); 
      else           h_mDca[ipart][ifile]->Draw("pE same"); 
    }
    legDca[ipart]->Draw("same");
    cout << "Proccessed DCA" << endl;
 
    c->cd(ipart+3);
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mDca_Ratio[ipart][ifile]->SetTitle(Form("%s %s/%s",histpart[ipart].c_str(),label[0].c_str(),label[ifile].c_str()));
      h_mDca_Ratio[ipart][ifile]->GetXaxis()->SetTitle("|DCA| (cm)");
      h_mDca_Ratio[ipart][ifile]->GetYaxis()->SetTitle("Ratio");
      h_mDca_Ratio[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mDca_Ratio[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mDca_Ratio[ipart][ifile]->SetLineColor(color[ifile]);
      legDcaRatio[ipart]->AddEntry(h_mDca_Ratio[ipart][ifile],Form("%s/%s",label[0].c_str(),label[ifile].c_str()),"p");
      if(ifile == 1) h_mDca_Ratio[ipart][ifile]->Draw("pE"); 
      else           h_mDca_Ratio[ipart][ifile]->Draw("pE same"); 
    }
    legDcaRatio[ipart]->Draw("same");
    cout << "Proccessed DCA Ratio" << endl;
  }  
  c->SaveAs(Form("figures/%s/%s/pTstudy/BasicParameterComparisons/Dca_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  

  cout << "NHits" << endl;
  TLegend *legNHits[2];
  TLegend *legNHits_Ratio[2];
  for(int ipart = 0; ipart < 2; ipart++)
  {
    legNHits[ipart] = new TLegend(0.4,0.6,0.6,0.8);
    legNHits_Ratio[ipart] = new TLegend(0.4,0.3,0.6,0.5);

    max = -9999999999;
    min = 9999999999;
    for(int ifile = 0; ifile < nopts; ifile++)
    {
      double tmax = h_mNHits[ipart][ifile]->GetMaximum();
      double tmin = h_mNHits[ipart][ifile]->GetMinimum();
      if(tmax > max) tmax = max;
      if(tmin < min) tmin = min;
    }  

    c->cd(ipart+1);
    c->cd(ipart+1)->SetLogy(0);
    for(int ifile = 0; ifile < nopts; ifile++)
    {
      h_mNHits[ipart][ifile]->SetTitle(Form("%s",histpart[ipart].c_str()));
      h_mNHits[ipart][ifile]->GetXaxis()->SetTitle("NHitsFit");
      h_mNHits[ipart][ifile]->GetYaxis()->SetTitle("Counts");
      h_mNHits[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mNHits[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mNHits[ipart][ifile]->SetLineColor(color[ifile]);
      //h_mNHits[ipart][ifile]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);
      legNHits[ipart]->AddEntry(h_mNHits[ipart][ifile],label[ifile].c_str(),"p");
      if(ifile == 0) h_mNHits[ipart][ifile]->Draw("pE"); 
      else           h_mNHits[ipart][ifile]->Draw("pE same"); 
    }
    legNHits[ipart]->Draw("same");

    c->cd(ipart+3);
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mNHits_Ratio[ipart][ifile]->SetTitle(Form("%s %s/%s",histpart[ipart].c_str(),label[0].c_str(),label[ifile].c_str()));
      h_mNHits_Ratio[ipart][ifile]->GetXaxis()->SetTitle("NHitsFit");
      h_mNHits_Ratio[ipart][ifile]->GetYaxis()->SetTitle("Ratio");
      h_mNHits_Ratio[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mNHits_Ratio[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mNHits_Ratio[ipart][ifile]->SetLineColor(color[ifile]);
      legNHits_Ratio[ipart]->AddEntry(h_mNHits_Ratio[ipart][ifile],Form("%s/%s",label[0].c_str(),label[ifile].c_str()),"p");
      if(ifile == 1) h_mNHits_Ratio[ipart][ifile]->Draw("pE"); 
      else           h_mNHits_Ratio[ipart][ifile]->Draw("pE same"); 
    }
    legNHits_Ratio[ipart]->Draw("same");
  }  
  c->SaveAs(Form("figures/%s/%s/pTstudy/BasicParameterComparisons/NHits_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  

  cout << "NHitsRatio" << endl;
  TLegend *legNHitsRatio[2];
  TLegend *legNHitsRatioRatio[2];
  for(int ipart = 0; ipart < 2; ipart++)
  {
    legNHitsRatio[ipart] = new TLegend(0.2,0.6,0.4,0.8);
    legNHitsRatioRatio[ipart] = new TLegend(0.2,0.6,0.4,0.8);

    max = -9999999999;
    min = 9999999999;
    for(int ifile = 0; ifile < nopts; ifile++)
    {
      double tmax = h_mNHitsRatio[ipart][ifile]->GetMaximum();
      double tmin = h_mNHitsRatio[ipart][ifile]->GetMinimum();
      if(tmax > max) tmax = max;
      if(tmin < min) tmin = min;
    }  

    c->cd(ipart+1);
    c->cd(ipart+1)->SetLogy(0);
    for(int ifile = 0; ifile < nopts; ifile++)
    {
      h_mNHitsRatio[ipart][ifile]->SetTitle(Form("%s",histpart[ipart].c_str()));
      h_mNHitsRatio[ipart][ifile]->GetXaxis()->SetTitle("NHitsFit/NHitsMax");
      h_mNHitsRatio[ipart][ifile]->GetYaxis()->SetTitle("Counts");
      h_mNHitsRatio[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mNHitsRatio[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mNHitsRatio[ipart][ifile]->SetLineColor(color[ifile]);
      //h_mNHitsRatio[ipart][ifile]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);
      legNHitsRatio[ipart]->AddEntry(h_mNHitsRatio[ipart][ifile],label[ifile].c_str(),"p");
      if(ifile == 0) h_mNHitsRatio[ipart][ifile]->Draw("pE"); 
      else           h_mNHitsRatio[ipart][ifile]->Draw("pE same"); 
    }
    legNHitsRatio[ipart]->Draw("same");

    c->cd(ipart+3);
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mNHitsRatio_Ratio[ipart][ifile]->SetTitle(Form("%s %s/%s",histpart[ipart].c_str(),label[0].c_str(),label[ifile].c_str()));
      h_mNHitsRatio_Ratio[ipart][ifile]->GetXaxis()->SetTitle("NHitsFit/NHitsMax");
      h_mNHitsRatio_Ratio[ipart][ifile]->GetYaxis()->SetTitle("Ratio");
      h_mNHitsRatio_Ratio[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mNHitsRatio_Ratio[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mNHitsRatio_Ratio[ipart][ifile]->SetLineColor(color[ifile]);
      legNHitsRatioRatio[ipart]->AddEntry(h_mNHitsRatio_Ratio[ipart][ifile],Form("%s/%s",label[0].c_str(),label[ifile].c_str()),"p");
      if(ifile == 1) h_mNHitsRatio_Ratio[ipart][ifile]->Draw("pE"); 
      else           h_mNHitsRatio_Ratio[ipart][ifile]->Draw("pE same"); 
    }
    legNHitsRatioRatio[ipart]->Draw("same");
  }  
  c->SaveAs(Form("figures/%s/%s/pTstudy/BasicParameterComparisons/NHitsRatio_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  

  cout << "DEDx" << endl;
  TLegend *legDEdx_1D[2];
  TLegend *legDEdx_1DRatio[2];
  for(int ipart = 0; ipart < 2; ipart++)
  {
    legDEdx_1D[ipart] = new TLegend(0.4,0.6,0.6,0.8);
    legDEdx_1DRatio[ipart] = new TLegend(0.4,0.6,0.6,0.8);

    max = -9999999999;
    min = 9999999999;
    for(int ifile = 0; ifile < nopts; ifile++)
    {
      double tmax = h_mDEdx_1D[ipart][ifile]->GetMaximum();
      double tmin = h_mDEdx_1D[ipart][ifile]->GetMinimum();
      if(tmax > max) tmax = max;
      if(tmin < min) tmin = min;
    }  

    c->cd(ipart+1);
    c->cd(ipart+1)->SetLogy();
    for(int ifile = 0; ifile < nopts; ifile++)
    {
      h_mDEdx_1D[ipart][ifile]->SetTitle(Form("%s",histpart[ipart].c_str()));
      h_mDEdx_1D[ipart][ifile]->GetXaxis()->SetTitle("dE/dx");
      h_mDEdx_1D[ipart][ifile]->GetYaxis()->SetTitle("Counts");
      h_mDEdx_1D[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mDEdx_1D[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mDEdx_1D[ipart][ifile]->SetLineColor(color[ifile]);
      //h_mDEdx_1D[ipart][ifile]->GetYaxis()->SetRangeUser(0.9*min,1.1*max);
      legDEdx_1D[ipart]->AddEntry(h_mDEdx_1D[ipart][ifile],label[ifile].c_str(),"p");
      if(ifile == 0) h_mDEdx_1D[ipart][ifile]->Draw("pE"); 
      else           h_mDEdx_1D[ipart][ifile]->Draw("pE same"); 
    }
    legDEdx_1D[ipart]->Draw("same");

    c->cd(ipart+3);
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mDEdx_1D_Ratio[ipart][ifile]->SetTitle(Form("%s %s/%s",histpart[ipart].c_str(),label[0].c_str(),label[ifile].c_str()));
      h_mDEdx_1D_Ratio[ipart][ifile]->GetXaxis()->SetTitle("dE/dx");
      h_mDEdx_1D_Ratio[ipart][ifile]->GetYaxis()->SetTitle("Ratio");
      h_mDEdx_1D_Ratio[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mDEdx_1D_Ratio[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mDEdx_1D_Ratio[ipart][ifile]->SetLineColor(color[ifile]);
      legDEdx_1DRatio[ipart]->AddEntry(h_mDEdx_1D_Ratio[ipart][ifile],Form("%s/%s",label[0].c_str(),label[ifile].c_str()),"p");
      if(ifile == 1) h_mDEdx_1D_Ratio[ipart][ifile]->Draw("pE"); 
      else           h_mDEdx_1D_Ratio[ipart][ifile]->Draw("pE same"); 
    }
    legDEdx_1DRatio[ipart]->Draw("same");
  }  
  c->SaveAs(Form("figures/%s/%s/pTstudy/BasicParameterComparisons/DEdx_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  

}


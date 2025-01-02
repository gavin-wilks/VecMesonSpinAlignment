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

void simulation_comp_kaonspectra(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 40, string sim_py = "Rc5", string sim_em = "Rc8")//defaultF = 0 is BESII, defaultF = 1 is BESI
{
  const int nopts = 9; 
  string spectra[nopts] = { "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1"};
 
  string filename[nopts] = {"PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
                            "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root"};
 
//  string spectra[nopts] = { "PhiEmbedding_NoHelicity_20240822EP_phiweight1",
//                            "PhiEmbedding_NoHelicity_20240822EP_phiweight1",
//                            "PhiEmbedding_NoHelicity_20240822EP_phiweight1",
//                            "PhiEmbedding_NoHelicity_20240822EP_phiweight1",
//                            "PhiEmbedding_NoHelicity_20240822EP_phiweight1",
//                            "PhiEmbedding_NoHelicity_20240822EP_phiweight1",
//                            "PhiEmbedding_NoHelicity_20240822EP_phiweight1"};
// 
//  string filename[nopts] = {"Embedding_19GeV_20240820.root",
//                            "Embedding_19GeV_20240820.root",
//                            "Embedding_19GeV_20240820.root",
//                            "Embedding_19GeV_20240820.root",
//                            "Embedding_19GeV_20240820.root",
//                            "Embedding_19GeV_20240820.root",
//                            "Embedding_19GeV_20240820.root"};
//
  //string label[nopts] = {"MC","RapidityCut","TPC Eff","|#eta| cut","ToF Match","m^{2} PID Eff","n#sigma_{K} PID Eff"};  
  string label[nopts] = {"MC","MC with RC","RC","+rapidity cut","+|#eta| cut","+TPC Eff","+ToF Match","+m^{2} PID Eff","+n#sigma_{K} PID Eff"};  
  string filelabel[nopts] = {"RC0","RC1","RC2","RC3","RC4","RC5","RC6","RC7","RC8"};
 
  string rcnum[nopts] = {"rc0","rc1","rc2","rc3","rc4","rc5","rc6","rc7","rc8"}; 
  float integral[2][nopts] = {0.0};


  //const int nopts = 7; 
 
  //string spectra[nopts] = {"Pythia_20240808_EtaCutFirst",
  //                         "Pythia_20240808_EtaCutFirst",
  //                         "Pythia_20240808_EtaCutFirst",
  //                         "Pythia_20240808_EtaCutFirst",
  //                         "Pythia_20240808_EtaCutFirst",
  //                         "Pythia_20240808_EtaCutFirst",
  //                         "Pythia_20240808_EtaCutFirst"};
 
  //string filename[nopts] = {"Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",
  //                          "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",
  //                          "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",
  //                          "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",
  //                          "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",
  //                          "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",
  //                          "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"};

  ////string label[nopts] = {"MC","RapidityCut","TPC Eff","|#eta| cut","ToF Match","m^{2} PID Eff","n#sigma_{K} PID Eff"};  
  ////string label[nopts] = {"MC","RapidityCut","|#eta| cut","TPC Eff","ToF Match","m^{2} PID Eff","n#sigma_{K} PID Eff"};  
  ////string label[nopts] = {"S0)MC","S1)S0+RapidityCut","S2)S1+|#eta| cut","S3)S2+TPC Eff","S4)S3+ToF Match","S5)S4+m^{2} PID Eff","S6)S5+n#sigma_{K} PID Eff"};  
  //string label[nopts] = {"MC","RC+RapidityCut","+|#eta| cut","+TPC Eff","+ToF Match","+m^{2} PID Eff","+n#sigma_{K} PID Eff"};  
  //string filelabel[nopts] = {"MC","RC0","RC1","RC2","RC3","RC4","RC5"};
 
  //string rcnum[nopts] = {"mc","rc0","rc1","rc2","rc3","rc4","rc5"}; 
  //float integral[2][nopts] = {0.0};

  //const int nopts = 2; 
 
  //string spectra[nopts] = {"Pythia_NoWeights_20240731",
  //                         "Pythia_PIDEff_Helicity2D_20240729"};
 
  //string filename[nopts] = {"Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",
  //                          "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"};

  //string label[nopts] = {"Default","Helicity Frame SA"};  
  //string filelabel[nopts] = {"Default","HFSA"};
 
  //int rcnum[nopts] = {5,5}; 
  //float integral[2][nopts] = {0.0};

  TH2D *h_mCounts[2][nopts];  
  TH2D *h_mCounts_Ratio[2][nopts];  
  TH1D *h_mCounts_pt[2][nopts];  
  TH1D *h_mCounts_Ratio_pt[2][nopts];  
  TH1D *h_mCounts_ptbin[2][nopts][50];  
  TH1D *h_mCounts_Ratio_ptbin[2][nopts][50];  
  TH1D *h_mCounts_y[2][nopts];  
  TH1D *h_mCounts_Ratio_y[2][nopts];  
  TH1D *h_mCounts_ybin[2][nopts][50];  
  TH1D *h_mCounts_Ratio_ybin[2][nopts][50];  

  int nbinsx, nbinsy;

  string daughter[2] = {"kminus","kplus"};

  for(int ipart = 0; ipart < 2; ipart++)
  {
    for(int ifile = 0; ifile < nopts; ifile++)  
    {
      string inputfile = Form("effaccfiles/%s/%s/%s/%s",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[ifile].c_str(),filename[ifile].c_str());
      TFile *File_InPut = TFile::Open(inputfile.c_str());
      string KEY = Form("%s_%s_cent%d",rcnum[ifile].c_str(),daughter[ipart].c_str(),i_cent);
      cout << KEY << endl;
      string KEY_count = Form("%s_%s_cent%d_%d",rcnum[ifile].c_str(),daughter[ipart].c_str(),i_cent,ifile);
      cout << KEY_count << endl;
      h_mCounts[ipart][ifile] = (TH2D*)((TH2D*) File_InPut->Get(KEY.c_str()))->Clone(KEY_count.c_str());
      h_mCounts[ipart][ifile]->RebinY(2);
      h_mCounts[ipart][ifile]->RebinX(3);
      h_mCounts[ipart][ifile]->Print();
      nbinsx = h_mCounts[ipart][ifile]->GetNbinsX();
      nbinsy = h_mCounts[ipart][ifile]->GetNbinsY();
      integral[ipart][ifile] = h_mCounts[ipart][ifile]->Integral(1,nbinsx,1,nbinsy);
    }
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      //h_mCounts[ipart][ifile]->Scale(integral[ipart][0]/integral[ipart][ifile]);
      h_mCounts_Ratio[ipart][ifile] = (TH2D*) h_mCounts[ipart][ifile]->Clone();
      h_mCounts_Ratio[ipart][ifile]->Divide(h_mCounts[ipart][0]);
    }
    for(int ifile = 0; ifile < nopts; ifile++)  
    {
      ///////////// 1D pT ////////////
      h_mCounts_pt[ipart][ifile] = (TH1D*) h_mCounts[ipart][ifile]->ProjectionX(Form("pt%d",ifile),1,nbinsy,"e"); 
      for(int iy = 1; iy <= nbinsy; iy++)
      {
        h_mCounts_ptbin[ipart][ifile][iy-1] = (TH1D*) h_mCounts[ipart][ifile]->ProjectionX(Form("pt%d_%d",ifile,iy),iy,iy,"e"); 
      }
      ///////////// 1D y /////////////
      h_mCounts_y[ipart][ifile] = (TH1D*) h_mCounts[ipart][ifile]->ProjectionY(Form("y%d",ifile),1,nbinsx,"e"); 
      for(int ipt = 1; ipt <= nbinsx; ipt++)
      {
        h_mCounts_ybin[ipart][ifile][ipt-1] = (TH1D*) h_mCounts[ipart][ifile]->ProjectionY(Form("y%d_%d",ifile,ipt),ipt,ipt,"e"); 
      }
    }
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      ///////////// 1D pT ////////////
      h_mCounts_Ratio_pt[ipart][ifile] = (TH1D*) h_mCounts_pt[ipart][ifile]->Clone();
      h_mCounts_Ratio_pt[ipart][ifile]->Divide(h_mCounts_pt[ipart][0]);
      for(int iy = 0; iy < nbinsy; iy++)
      {
        h_mCounts_Ratio_ptbin[ipart][ifile][iy] = (TH1D*) h_mCounts_ptbin[ipart][ifile][iy]->Clone();
        h_mCounts_Ratio_ptbin[ipart][ifile][iy]->Divide(h_mCounts_ptbin[ipart][0][iy]);
      }

      ///////////// 1D y /////////////
      h_mCounts_Ratio_y[ipart][ifile] = (TH1D*) h_mCounts_y[ipart][ifile]->Clone();
      h_mCounts_Ratio_y[ipart][ifile]->Divide(h_mCounts_y[ipart][0]);
      for(int ipt = 0; ipt < nbinsx; ipt++)
      {
        h_mCounts_Ratio_ybin[ipart][ifile][ipt] = (TH1D*) h_mCounts_ybin[ipart][ifile][ipt]->Clone();
        h_mCounts_Ratio_ybin[ipart][ifile][ipt]->Divide(h_mCounts_ybin[ipart][0][ipt]);
      }
    }
  }   

  string histpart[2] = {"K^{-}","K^{+}"};

  TCanvas *c = new TCanvas("c","c",10,10,800,400);
  c->Divide(2,1);
  for(int ipart = 0; ipart < 2; ipart++) 
  {
    c->cd(ipart+1)->SetLeftMargin(0.15);
    c->cd(ipart+1)->SetBottomMargin(0.15);
    c->cd(ipart+1)->SetTicks(1,1);
    c->cd(ipart+1)->SetGrid(0,0);
  }
  
  //int color[4] = {kBlack,kOrange+7,kBlue,kGray+2};
  int color[9] = {kBlack,kOrange+7,kBlue,kGray+2,kRed,kGreen,kViolet,kBlack,kOrange+7};
  int marker[9] = {20,20,20,20,24,24,24,24,24};

  for(int ifile = 1; ifile < nopts; ifile++)
  {
    for(int ipart = 0; ipart < 2; ipart++)
    {
      c->cd(ipart+1);
      h_mCounts_Ratio[ipart][ifile]->SetTitle(Form("%s %s/%s",histpart[ipart].c_str(),label[0].c_str(),label[ifile].c_str()));
      h_mCounts_Ratio[ipart][ifile]->GetXaxis()->SetTitle("p_{T}");
      h_mCounts_Ratio[ipart][ifile]->GetYaxis()->SetTitle("y");
      h_mCounts_Ratio[ipart][ifile]->Draw("colz");
      c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_%s%sRatio_ptyspectra2D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),daughter[ipart].c_str(),filelabel[ifile].c_str(),filelabel[0].c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
    }
  }  
  for(int ifile = 0; ifile < nopts; ifile++)
  {
    for(int ipart = 0; ipart < 2; ipart++)
    {
      c->cd(ipart+1);
      h_mCounts[ipart][ifile]->SetTitle(Form("%s %s",histpart[ipart].c_str(),label[ifile].c_str()));
      h_mCounts[ipart][ifile]->GetXaxis()->SetTitle("p_{T}");
      h_mCounts[ipart][ifile]->GetYaxis()->SetTitle("y");
      h_mCounts[ipart][ifile]->Draw("colz");
      c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_%s_ptyspectra2D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),daughter[ipart].c_str(),filelabel[ifile].c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
    }
  }  

  ///////////// pT 1D //////////////////////
  //c->cd();
  TLegend *leg = new TLegend(0.5,0.2,0.7,0.4);

  string filefullstring = filelabel[0];

  for(int ipart = 0; ipart < 2; ipart++)
  {
    c->cd(ipart+1);
    c->cd(ipart+1)->SetLogy();
    h_mCounts_pt[ipart][0]->SetTitle(Form("%s Centrality %d",histpart[ipart].c_str(),i_cent));
    h_mCounts_pt[ipart][0]->GetXaxis()->SetTitle("p_{T}");
    h_mCounts_pt[ipart][0]->GetYaxis()->SetTitle("Counts");
    h_mCounts_pt[ipart][0]->SetMarkerStyle(marker[0]);
    h_mCounts_pt[ipart][0]->SetMarkerColor(color[0]);
    h_mCounts_pt[ipart][0]->SetLineColor(color[0]);
    double maxy = h_mCounts_pt[ipart][0]->GetMaximum();
    double miny = h_mCounts_pt[ipart][0]->GetMinimum();
    if(ipart == 0) leg->AddEntry(h_mCounts_pt[ipart][0],label[0].c_str(),"p");

    //string filefullstring = filelabel[0];

    for(int ifile = 1; ifile < nopts; ifile++)
    {
      if(h_mCounts_pt[ipart][ifile]->GetMaximum() > maxy) maxy = h_mCounts_pt[ipart][ifile]->GetMaximum();
      if(h_mCounts_pt[ipart][ifile]->GetMinimum() < miny) miny = h_mCounts_pt[ipart][ifile]->GetMinimum();
      if(ipart == 0) filefullstring += filelabel[ifile];
    }
    h_mCounts_pt[ipart][0]->GetYaxis()->SetRangeUser(0.1/*miny*0.9*/,maxy*2.0);
    h_mCounts_pt[ipart][0]->Draw("pE");
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mCounts_pt[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mCounts_pt[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mCounts_pt[ipart][ifile]->SetLineColor(color[ifile]);
      h_mCounts_pt[ipart][ifile]->Draw("pE same");
      if(ipart == 0) leg->AddEntry(h_mCounts_pt[ipart][ifile],label[ifile].c_str(),"p");
    }
    leg->Draw("same"); 
  }

  c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/kaons_%s_ptspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// pT 1D //////////////////////

  ///////////// pT 1D individual y bins //////////////////////
  TCanvas *c2 = new TCanvas("c2","c2",10,10,2500,2500);
  c2->Divide(5,5);
  for(int ipart = 0; ipart < 2; ipart++)
  {
    for(int iy = 0; iy < nbinsy; iy++)
    {
      c2->cd(iy+1);
      c2->cd(iy+1)->SetLogy();
      c2->cd(iy+1)->SetLeftMargin(0.15);
      c2->cd(iy+1)->SetBottomMargin(0.15);
      c2->cd(iy+1)->SetTicks(1,1);
      c2->cd(iy+1)->SetGrid(0,0);
      h_mCounts_ptbin[ipart][0][iy]->SetTitle(Form("%s Centrality %d, %d/%d<y<%d/%d",histpart[ipart].c_str(),i_cent,(iy-nbinsy/2),nbinsy/2,(iy-nbinsy/2+1),nbinsy/2));
      h_mCounts_ptbin[ipart][0][iy]->GetXaxis()->SetTitle("p_{T}");
      h_mCounts_ptbin[ipart][0][iy]->GetYaxis()->SetTitle("Counts");
      h_mCounts_ptbin[ipart][0][iy]->SetMarkerStyle(marker[0]);
      h_mCounts_ptbin[ipart][0][iy]->SetMarkerColor(color[0]);
      h_mCounts_ptbin[ipart][0][iy]->SetLineColor(color[0]);
      double maxbin = h_mCounts_ptbin[ipart][0][iy]->GetMaximum();
      double minbin = h_mCounts_ptbin[ipart][0][iy]->GetMinimum();

      for(int ifile = 1; ifile < nopts; ifile++)
      {
        if(h_mCounts_ptbin[ipart][ifile][iy]->GetMaximum() > maxbin) maxbin = h_mCounts_ptbin[ipart][ifile][iy]->GetMaximum();
        if(h_mCounts_ptbin[ipart][ifile][iy]->GetMinimum() < minbin) minbin = h_mCounts_ptbin[ipart][ifile][iy]->GetMinimum();
      }
      h_mCounts_ptbin[ipart][0][iy]->GetYaxis()->SetRangeUser(0.1/*minbin*0.9*/,maxbin*2.0);
      h_mCounts_ptbin[ipart][0][iy]->Draw("pE");
      for(int ifile = 1; ifile < nopts; ifile++)
      {
        h_mCounts_ptbin[ipart][ifile][iy]->SetMarkerStyle(marker[ifile]);
        h_mCounts_ptbin[ipart][ifile][iy]->SetMarkerColor(color[ifile]);
        h_mCounts_ptbin[ipart][ifile][iy]->SetLineColor(color[ifile]);
        h_mCounts_ptbin[ipart][ifile][iy]->Draw("pE same");
      }
      leg->Draw("same"); 
    }
    c2->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_%s_individualybins_ptspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),daughter[ipart].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  }
  ///////////// pT 1D individual y bins //////////////////////

  ///////////// pT 1D ratio //////////////////////
  //c->cd();
  TLegend *legratio = new TLegend(0.5,0.2,0.7,0.4);

  for(int ipart = 0; ipart < 2; ipart++)
  {
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      c->cd(ipart+1);
      c->cd(ipart+1)->SetLogy(0);
      h_mCounts_Ratio_pt[ipart][ifile]->SetTitle(Form("%s Centrality %d",histpart[ipart].c_str(),i_cent));
      h_mCounts_Ratio_pt[ipart][ifile]->GetXaxis()->SetTitle("p_{T}");
      h_mCounts_Ratio_pt[ipart][ifile]->GetYaxis()->SetTitle("Counts");
      h_mCounts_Ratio_pt[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mCounts_Ratio_pt[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mCounts_Ratio_pt[ipart][ifile]->SetLineColor(color[ifile]);
      h_mCounts_Ratio_pt[ipart][ifile]->GetYaxis()->SetRangeUser(0.0,1.2);
      if(ipart == 0) legratio->AddEntry(h_mCounts_Ratio_pt[ipart][ifile],Form("%s/%s",label[ifile].c_str(),label[0].c_str()),"p");
      if (ifile == 1) h_mCounts_Ratio_pt[ipart][ifile]->Draw("pE");
      if (ifile >  1) h_mCounts_Ratio_pt[ipart][ifile]->Draw("pE same");
    }
 
    legratio->Draw("same"); 
  }
  c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/kaons_%s_ratioptspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// pT 1D ratio //////////////////////

  ///////////// pT 1D ratio //////////////////////
  for(int ipart = 0; ipart < 2; ipart++)
  {
    for(int iy = 0; iy < nbinsy; iy++)
    {
      c2->cd(iy+1);
      c2->cd(iy+1)->SetLogy(0);
      for(int ifile = 1; ifile < nopts; ifile++)
      {
        h_mCounts_Ratio_ptbin[ipart][ifile][iy]->SetTitle(Form("%s Centrality %d, %d/%d<y<%d/%d",histpart[ipart].c_str(),i_cent,(iy-nbinsy/2),nbinsy/2,(iy-nbinsy/2+1),nbinsy/2));
        h_mCounts_Ratio_ptbin[ipart][ifile][iy]->GetXaxis()->SetTitle("p_{T}");
        h_mCounts_Ratio_ptbin[ipart][ifile][iy]->GetYaxis()->SetTitle("Counts");
        h_mCounts_Ratio_ptbin[ipart][ifile][iy]->SetMarkerStyle(marker[ifile]);
        h_mCounts_Ratio_ptbin[ipart][ifile][iy]->SetMarkerColor(color[ifile]);
        h_mCounts_Ratio_ptbin[ipart][ifile][iy]->SetLineColor(color[ifile]);
        h_mCounts_Ratio_ptbin[ipart][ifile][iy]->GetYaxis()->SetRangeUser(0.0,1.5);
        //legratio->AddEntry(h_mCounts_Ratio_ptbin[ifile][iy],Form("%s/%s",label[ifile].c_str(),label[0].c_str()),"p");
        if (ifile == 1) h_mCounts_Ratio_ptbin[ipart][ifile][iy]->Draw("pE");
        if (ifile >  1) h_mCounts_Ratio_ptbin[ipart][ifile][iy]->Draw("pE same");
      }
 
      legratio->Draw("same"); 
    }
    c2->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_%s_individualybins_ratioptspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),daughter[ipart].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  }
  ///////////// pT 1D ratio //////////////////////

  ///////////// y 1D ///////////////////////
  TLegend *legy = new TLegend(0.4,0.1,0.6,0.3);

  for(int ipart = 0; ipart < 2; ipart++)
  {
    c->cd(ipart+1)->SetLogy(0);
    h_mCounts_y[ipart][0]->SetTitle(Form("%s Centrality %d",histpart[ipart].c_str(),i_cent));
    h_mCounts_y[ipart][0]->GetXaxis()->SetTitle("y");
    h_mCounts_y[ipart][0]->GetYaxis()->SetTitle("Counts");
    h_mCounts_y[ipart][0]->SetMarkerStyle(marker[0]);
    h_mCounts_y[ipart][0]->SetMarkerColor(color[0]);
    h_mCounts_y[ipart][0]->SetLineColor(color[0]);
    double max = h_mCounts_y[ipart][0]->GetMaximum();
    double min = h_mCounts_y[ipart][0]->GetMinimum();
    if(ipart == 0)legy->AddEntry(h_mCounts_y[ipart][0],label[0].c_str(),"p");

    for(int ifile = 1; ifile < nopts; ifile++)
    {
      if(h_mCounts_y[ipart][ifile]->GetMaximum() > max) max = h_mCounts_y[ipart][ifile]->GetMaximum();
      if(h_mCounts_y[ipart][ifile]->GetMinimum() < min) min = h_mCounts_y[ipart][ifile]->GetMinimum();
    }
    h_mCounts_y[ipart][0]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);
    h_mCounts_y[ipart][0]->Draw("pE");
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mCounts_y[ipart][ifile]->SetMarkerStyle(marker[ifile]);
      h_mCounts_y[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mCounts_y[ipart][ifile]->SetLineColor(color[ifile]);
      if(ipart == 0) legy->AddEntry(h_mCounts_y[ipart][ifile],label[ifile].c_str(),"p");
      h_mCounts_y[ipart][ifile]->Draw("pE same");
    }
    legy->Draw("same"); 
  }
  c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/kaons_%s_yspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// y 1D ///////////////////////

  ///////////// y 1D individual y bins //////////////////////
  for(int ipart = 0; ipart < 2; ipart++) 
  {
    for(int ipt = 0; ipt < nbinsx; ipt++)
    {
      c2->cd(ipt+1);
      c2->cd(ipt+1)->SetLogy(0);
      c2->cd(ipt+1)->SetLeftMargin(0.15);
      c2->cd(ipt+1)->SetBottomMargin(0.15);
      c2->cd(ipt+1)->SetTicks(1,1);
      c2->cd(ipt+1)->SetGrid(0,0);
      h_mCounts_ybin[ipart][0][ipt]->SetTitle(Form("%s Centrality %d, %1.1f<p_{T}<%1.1f",histpart[ipart].c_str(),i_cent,vmsa::ptRawStart[ipt],vmsa::ptRawStop[ipt]));
      h_mCounts_ybin[ipart][0][ipt]->GetXaxis()->SetTitle("p_{T}");
      h_mCounts_ybin[ipart][0][ipt]->GetYaxis()->SetTitle("Counts");
      h_mCounts_ybin[ipart][0][ipt]->SetMarkerStyle(marker[0]);
      h_mCounts_ybin[ipart][0][ipt]->SetMarkerColor(color[0]);
      h_mCounts_ybin[ipart][0][ipt]->SetLineColor(color[0]);
      double maxbin = h_mCounts_ybin[ipart][0][ipt]->GetMaximum();
      double minbin = h_mCounts_ybin[ipart][0][ipt]->GetMinimum();

      for(int ifile = 1; ifile < nopts; ifile++)
      {
        if(h_mCounts_ybin[ipart][ifile][ipt]->GetMaximum() > maxbin) maxbin = h_mCounts_ybin[ipart][ifile][ipt]->GetMaximum();
        if(h_mCounts_ybin[ipart][ifile][ipt]->GetMinimum() < minbin) minbin = h_mCounts_ybin[ipart][ifile][ipt]->GetMinimum();
      }
      h_mCounts_ybin[ipart][0][ipt]->GetYaxis()->SetRangeUser(minbin*0.9,maxbin*1.1);
      h_mCounts_ybin[ipart][0][ipt]->Draw("pE");
      for(int ifile = 1; ifile < nopts; ifile++)
      {
        h_mCounts_ybin[ipart][ifile][ipt]->SetMarkerStyle(marker[ifile]);
        h_mCounts_ybin[ipart][ifile][ipt]->SetMarkerColor(color[ifile]);
        h_mCounts_ybin[ipart][ifile][ipt]->SetLineColor(color[ifile]);
        h_mCounts_ybin[ipart][ifile][ipt]->Draw("pE same");
      }
      legy->Draw("same"); 
    }
    c2->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_%s_individualptbins_yspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),daughter[ipart].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  }
  ///////////// y 1D individual y bins //////////////////////

  ///////////// y 1D ratio //////////////////////

  TLegend *legyratio = new TLegend(0.5,0.2,0.7,0.4);

  for(int ipart = 0; ipart < 2; ipart++)
  {
    c->cd(ipart+1)->SetLogy(0);
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mCounts_Ratio_y[ipart][ifile]->SetTitle(Form("%s Centrality %d",histpart[ipart].c_str(),i_cent));
      h_mCounts_Ratio_y[ipart][ifile]->GetXaxis()->SetTitle("y");
      h_mCounts_Ratio_y[ipart][ifile]->GetYaxis()->SetTitle("Counts");
      h_mCounts_Ratio_y[ipart][ifile]->SetMarkerStyle(marker[0]);
      h_mCounts_Ratio_y[ipart][ifile]->SetMarkerColor(color[ifile]);
      h_mCounts_Ratio_y[ipart][ifile]->SetLineColor(color[ifile]);
      h_mCounts_Ratio_y[ipart][ifile]->GetYaxis()->SetRangeUser(0.0,1.2);
      if(ipart == 0) legyratio->AddEntry(h_mCounts_Ratio_y[ipart][ifile],Form("%s/%s",label[ifile].c_str(),label[0].c_str()),"p");
      if (ifile == 1) h_mCounts_Ratio_y[ipart][ifile]->Draw("pE");
      if (ifile >  1) h_mCounts_Ratio_y[ipart][ifile]->Draw("pE same");
    }
    legyratio->Draw("same"); 
  }
  c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/kaons_%s_ratioyspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// y 1D ratio //////////////////////

  ///////////// pT 1D ratio //////////////////////
  for(int ipart = 0; ipart < 2; ipart++)
  {
    for(int ipt = 0; ipt < nbinsx; ipt++)
    {
      c2->cd(ipt+1);
      for(int ifile = 1; ifile < nopts; ifile++)
      {
        h_mCounts_Ratio_ybin[ipart][ifile][ipt]->SetTitle(Form("%s Centrality %d, %1.1f<p_{T}<%1.1f",histpart[ipart].c_str(),i_cent,vmsa::ptRawStart[ipt],vmsa::ptRawStop[ipt]));
        h_mCounts_Ratio_ybin[ipart][ifile][ipt]->GetXaxis()->SetTitle("y");
        h_mCounts_Ratio_ybin[ipart][ifile][ipt]->GetYaxis()->SetTitle("Counts");
        h_mCounts_Ratio_ybin[ipart][ifile][ipt]->SetMarkerStyle(marker[ifile]);
        h_mCounts_Ratio_ybin[ipart][ifile][ipt]->SetMarkerColor(color[ifile]);
        h_mCounts_Ratio_ybin[ipart][ifile][ipt]->SetLineColor(color[ifile]);
        h_mCounts_Ratio_ybin[ipart][ifile][ipt]->GetYaxis()->SetRangeUser(0.0,1.5);
        //legratio->AddEntry(h_mCounts_Ratio_ybin[ifile][ipt],Form("%s/%s",label[ifile].c_str(),label[0].c_str()),"p");
        if (ifile == 1) h_mCounts_Ratio_ybin[ipart][ifile][ipt]->Draw("pE");
        if (ifile >  1) h_mCounts_Ratio_ybin[ipart][ifile][ipt]->Draw("pE same");
      }
 
      legratio->Draw("same"); 
    }
    c2->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_%s_individualptbins_ratioyspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),daughter[ipart].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  }
  ///////////// pT 1D ratio //////////////////////


  //h_mCounts_Ratio_y->SetTitle(Form("Centrality %d",i_cent));
  //h_mCounts_Ratio_y->GetXaxis()->SetTitle("y");
  //h_mCounts_Ratio_y->GetYaxis()->SetTitle("Pythia/Embedding");
  //h_mCounts_Ratio_y->Draw("pE");
  
  //c->SaveAs(Form("figures/%s/%s/pTstudy/%s/PythiaEmbeddingRatio_yspectra1D_%s_Order%d_%s_Cent%d_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),  spectra.c_str(),correction.c_str(),order,etamode.c_str(),i_cent,sim_py.c_str(),sim_em.c_str()));  

}


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

void simulation_comp_phimesonspectra(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", int order = 2, string etamode = "eta1_eta1", int yspectra = 40, string sim_py = "Rc5", string sim_em = "Rc8")//defaultF = 0 is BESII, defaultF = 1 is BESI
{

  //const int nopts = 2; 
 
  //string spectra[nopts] = {"Pythia_NoWeights_20240731",
  //                         "Pythia_PIDEff_Helicity2D_20240729"};
 
  //string filename[nopts] = {"Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root",
  //                          "Eff_19GeV_SingleParticle_noToF_Mode0_EtaMode0_PtBins2_5.root"};

  //string label[nopts] = {"Default","Helicity Frame SA"};  
  //string filelabel[nopts] = {"Default","HFSA"};
 
  //string rcnum[nopts] = {"Rc5","Rc5"}; 
  //float integral[nopts] = {0.0};
  const int nopts = 10; 
  string spectra[nopts] = { 
                            "PhiEmbedding_noweights_nomccuts_rapidity",
                            "PhiEmbedding_noweights_nomccuts_rapidity",
                            "PhiEmbedding_noweights_nomccuts_rapidity",
                            "PhiEmbedding_noweights_nomccuts_rapidity",
                            "PhiEmbedding_noweights_nomccuts_rapidity",
                            "PhiEmbedding_noweights_nomccuts_rapidity",
                            "PhiEmbedding_noweights_nomccuts_rapidity",
                            "PhiEmbedding_noweights_nomccuts_rapidity",
                            "PhiEmbedding_noweights_nomccuts_rapidity",
                            "PhiEmbedding_noweights_nomccuts_rapidity"};
 
  string filename[nopts] = {
                            "PhiEmbedding_noweights_nomccuts_rapidity.root",
                            "PhiEmbedding_noweights_nomccuts_rapidity.root",
                            "PhiEmbedding_noweights_nomccuts_rapidity.root",
                            "PhiEmbedding_noweights_nomccuts_rapidity.root",
                            "PhiEmbedding_noweights_nomccuts_rapidity.root",
                            "PhiEmbedding_noweights_nomccuts_rapidity.root",
                            "PhiEmbedding_noweights_nomccuts_rapidity.root",
                            "PhiEmbedding_noweights_nomccuts_rapidity.root",
                            "PhiEmbedding_noweights_nomccuts_rapidity.root",
                            "PhiEmbedding_noweights_nomccuts_rapidity.root"};
 
  string label[nopts] = {"MC #phi-meson","MC (K^{+}+K^{-})","MC with RC","RC","+rapidity cut","+|#eta| cut","+TPC Eff","+ToF Match","+m^{2} PID Eff","+n#sigma_{K} PID Eff"};  
  string filelabel[nopts] = {"RC0","RC1","RC2","RC3","RC4","RC5","RC6","RC7","RC8","RC9"};
  
  int color[14] = {kBlack,kOrange+7,kBlue,kGray+2,kRed,kGreen,kViolet,kBlack,kOrange+7,kBlue,kGray+2,kRed,kGreen,kViolet};
  int marker[14] = {20,20,20,20,24,24,24,21,21,21,21,25,25,25};
 
  //const int nopts = 9; 
  //string spectra[nopts] = { "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1"};
 
  //string filename[nopts] = {"PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root",
  //                          "PhiEmbedding_FuncEmbedWeight_20240827_FinerToFEta_phiswitch1.root"};
 
  //int color[9] = {kBlack,kOrange+7,kBlue,kGray+2,kRed,kGreen,kViolet,kBlack,kOrange+7};
  //int marker[9] = {20,20,20,20,24,24,24,24,24};
  ////string label[nopts] = {"MC","RapidityCut","TPC Eff","|#eta| cut","ToF Match","m^{2} PID Eff","n#sigma_{K} PID Eff"};  
  //string label[nopts] = {"MC","MC with RC","RC","+rapidity cut","+|#eta| cut","+TPC Eff","+ToF Match","+m^{2} PID Eff","+n#sigma_{K} PID Eff"};  
  //string filelabel[nopts] = {"RC0","RC1","RC2","RC3","RC4","RC5","RC6","RC7","RC8"};
  //string label[nopts] = {"MC","RC+RapidityCut","+|#eta| cut","+TPC Eff","+ToF Match","+m^{2} PID Eff","+n#sigma_{K} PID Eff"};  
  //string filelabel[nopts] = {"RC2","RC3","RC4","RC5","RC6","RC7","RC8"};
 
  string rcnum[nopts] = {"Rc0","Rc1","Rc2","Rc3","Rc4","Rc5","Rc6","Rc7","Rc8","Rc9"}; 
  float integral[nopts] = {0.0};


  TH2D *h_mCounts[nopts];  
  TH2D *h_mCounts_Ratio[nopts];  
  TH1D *h_mCounts_pt[nopts];  
  TH1D *h_mCounts_Ratio_pt[nopts];  
  TH1D *h_mCounts_ptbin[nopts][50];  
  TH1D *h_mCounts_Ratio_ptbin[nopts][50];  
  TH1D *h_mCounts_y[nopts];  
  TH1D *h_mCounts_Ratio_y[nopts];  
  TH1D *h_mCounts_ybin[nopts][50];  
  TH1D *h_mCounts_Ratio_ybin[nopts][50];  

  int nbinsx, nbinsy;

  for(int ifile = 0; ifile < nopts; ifile++)  
  {
    string inputfile = Form("effaccfiles/%s/%s/%s/%s",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),spectra[ifile].c_str(),filename[ifile].c_str());
    TFile *File_InPut = TFile::Open(inputfile.c_str());
    string KEY = Form("h_m%sEffPtY_Cent_%d",rcnum[ifile].c_str(),i_cent);
    cout << KEY << endl;
    string KEY_count = Form("h_mRcEffPtY_Cent_%d_%d",i_cent,ifile);
    cout << KEY_count << endl;
    h_mCounts[ifile] = (TH2D*)((TH2D*) File_InPut->Get(KEY.c_str()))->Clone(KEY_count.c_str());
    h_mCounts[ifile]->Print();
    nbinsx = h_mCounts[ifile]->GetNbinsX();
    nbinsy = h_mCounts[ifile]->GetNbinsY();
    integral[ifile] = h_mCounts[ifile]->Integral(1,nbinsx,1,nbinsy);
  }
  for(int ifile = 1; ifile < nopts; ifile++)
  {
    //h_mCounts[ifile]->Scale(integral[0]/integral[ifile]);
    h_mCounts_Ratio[ifile] = (TH2D*) h_mCounts[ifile]->Clone();
    h_mCounts_Ratio[ifile]->Divide(h_mCounts[0]);
  }
  for(int ifile = 0; ifile < nopts; ifile++)  
  {
    ///////////// 1D pT ////////////
    h_mCounts_pt[ifile] = (TH1D*) h_mCounts[ifile]->ProjectionX(Form("pt%d",ifile),1,nbinsy,"e"); 
    for(int iy = 1; iy <= nbinsy; iy++)
    {
      h_mCounts_ptbin[ifile][iy-1] = (TH1D*) h_mCounts[ifile]->ProjectionX(Form("pt%d_%d",ifile,iy),iy,iy,"e"); 
    }
    ///////////// 1D y /////////////
    h_mCounts_y[ifile] = (TH1D*) h_mCounts[ifile]->ProjectionY(Form("y%d",ifile),1,nbinsx,"e"); 
    for(int ipt = 1; ipt <= nbinsx; ipt++)
    {
      h_mCounts_ybin[ifile][ipt-1] = (TH1D*) h_mCounts[ifile]->ProjectionY(Form("y%d_%d",ifile,ipt),ipt,ipt,"e"); 
    }
  }
  for(int ifile = 1; ifile < nopts; ifile++)
  {
    ///////////// 1D pT ////////////
    h_mCounts_Ratio_pt[ifile] = (TH1D*) h_mCounts_pt[ifile]->Clone();
    h_mCounts_Ratio_pt[ifile]->Divide(h_mCounts_pt[0]);
    for(int iy = 0; iy < nbinsy; iy++)
    {
      h_mCounts_Ratio_ptbin[ifile][iy] = (TH1D*) h_mCounts_ptbin[ifile][iy]->Clone();
      h_mCounts_Ratio_ptbin[ifile][iy]->Divide(h_mCounts_ptbin[0][iy]);
    }

    ///////////// 1D y /////////////
    h_mCounts_Ratio_y[ifile] = (TH1D*) h_mCounts_y[ifile]->Clone();
    h_mCounts_Ratio_y[ifile]->Divide(h_mCounts_y[0]);
    for(int ipt = 0; ipt < nbinsx; ipt++)
    {
      h_mCounts_Ratio_ybin[ifile][ipt] = (TH1D*) h_mCounts_ybin[ifile][ipt]->Clone();
      h_mCounts_Ratio_ybin[ifile][ipt]->Divide(h_mCounts_ybin[0][ipt]);
    }
  }

  TCanvas *c = new TCanvas("c","c",10,10,400,400);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.15);
  c->SetTicks(1,1);
  c->SetGrid(0,0);
  

  for(int ifile = 1; ifile < nopts; ifile++)
  {
    h_mCounts_Ratio[ifile]->SetTitle(Form("%s/%s",label[0].c_str(),label[ifile].c_str()));
    h_mCounts_Ratio[ifile]->GetXaxis()->SetTitle("p_{T}");
    h_mCounts_Ratio[ifile]->GetYaxis()->SetTitle("y");
    h_mCounts_Ratio[ifile]->Draw("colz");
    c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s%sRatio_ptyspectra2D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filelabel[ifile].c_str(),filelabel[0].c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  }  
  for(int ifile = 0; ifile < nopts; ifile++)
  {
    h_mCounts[ifile]->SetTitle(Form("%s",label[ifile].c_str()));
    h_mCounts[ifile]->GetXaxis()->SetTitle("p_{T}");
    h_mCounts[ifile]->GetYaxis()->SetTitle("y");
    h_mCounts[ifile]->Draw("colz");
    c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_ptyspectra2D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filelabel[ifile].c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  }  

  ///////////// pT 1D //////////////////////
  c->cd();
  c->cd()->SetLogy();
  TLegend *leg = new TLegend(0.5,0.2,0.7,0.4);

  h_mCounts_pt[0]->SetTitle(Form("Centrality %d",i_cent));
  h_mCounts_pt[0]->GetXaxis()->SetTitle("p_{T}");
  h_mCounts_pt[0]->GetYaxis()->SetTitle("Counts");
  h_mCounts_pt[0]->SetMarkerStyle(20);
  h_mCounts_pt[0]->SetMarkerColor(color[0]);
  h_mCounts_pt[0]->SetLineColor(color[0]);
  double maxy = h_mCounts_pt[0]->GetMaximum();
  double miny = h_mCounts_pt[0]->GetMinimum();
  leg->AddEntry(h_mCounts_pt[0],label[0].c_str(),"p");

  string filefullstring = filelabel[0];

  for(int ifile = 1; ifile < nopts; ifile++)
  {
    if(h_mCounts_pt[ifile]->GetMaximum() > maxy) maxy = h_mCounts_pt[ifile]->GetMaximum();
    if(h_mCounts_pt[ifile]->GetMinimum() < miny) miny = h_mCounts_pt[ifile]->GetMinimum();
    filefullstring += filelabel[ifile];
  }
  h_mCounts_pt[0]->GetYaxis()->SetRangeUser(0.1/*miny*0.9*/,maxy*2.0);
  h_mCounts_pt[0]->Draw("pE");
  for(int ifile = 1; ifile < nopts; ifile++)
  {
    h_mCounts_pt[ifile]->SetMarkerStyle(marker[ifile]);
    h_mCounts_pt[ifile]->SetMarkerColor(color[ifile]);
    h_mCounts_pt[ifile]->SetLineColor(color[ifile]);
    h_mCounts_pt[ifile]->Draw("pE same");
    leg->AddEntry(h_mCounts_pt[ifile],label[ifile].c_str(),"p");
  }
  leg->Draw("same"); 

  c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_ptspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// pT 1D //////////////////////

  ///////////// pT 1D individual y bins //////////////////////
  TCanvas *c2 = new TCanvas("c2","c2",10,10,2500,2500);
  c2->Divide(5,5);
  for(int iy = 0; iy < nbinsy; iy++)
  {
    c2->cd(iy+1);
    c2->cd(iy+1)->SetLogy();
    c2->cd(iy+1)->SetLeftMargin(0.15);
    c2->cd(iy+1)->SetBottomMargin(0.15);
    c2->cd(iy+1)->SetTicks(1,1);
    c2->cd(iy+1)->SetGrid(0,0);
    h_mCounts_ptbin[0][iy]->SetTitle(Form("Centrality %d, %d/%d<y<%d/%d",i_cent,(iy-nbinsy/2),nbinsy/2,(iy-nbinsy/2+1),nbinsy/2));
    h_mCounts_ptbin[0][iy]->GetXaxis()->SetTitle("p_{T}");
    h_mCounts_ptbin[0][iy]->GetYaxis()->SetTitle("Counts");
    h_mCounts_ptbin[0][iy]->SetMarkerStyle(20);
    h_mCounts_ptbin[0][iy]->SetMarkerColor(color[0]);
    h_mCounts_ptbin[0][iy]->SetLineColor(color[0]);
    double maxbin = h_mCounts_ptbin[0][iy]->GetMaximum();
    double minbin = h_mCounts_ptbin[0][iy]->GetMinimum();

    for(int ifile = 1; ifile < nopts; ifile++)
    {
      if(h_mCounts_ptbin[ifile][iy]->GetMaximum() > maxbin) maxbin = h_mCounts_ptbin[ifile][iy]->GetMaximum();
      if(h_mCounts_ptbin[ifile][iy]->GetMinimum() < minbin) minbin = h_mCounts_ptbin[ifile][iy]->GetMinimum();
    }
    h_mCounts_ptbin[0][iy]->GetYaxis()->SetRangeUser(0.1/*minbin*0.9*/,maxbin*2.0);
    h_mCounts_ptbin[0][iy]->Draw("pE");
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mCounts_ptbin[ifile][iy]->SetMarkerStyle(marker[ifile]);
      h_mCounts_ptbin[ifile][iy]->SetMarkerColor(color[ifile]);
      h_mCounts_ptbin[ifile][iy]->SetLineColor(color[ifile]);
      h_mCounts_ptbin[ifile][iy]->Draw("pE same");
    }
    leg->Draw("same"); 
  }
  c2->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_individualybins_ptspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// pT 1D individual y bins //////////////////////

  ///////////// pT 1D ratio //////////////////////
  c->cd();
  c->cd()->SetLogy(0);
  TLegend *legratio = new TLegend(0.5,0.2,0.7,0.4);

  for(int ifile = 1; ifile < nopts; ifile++)
  {
    h_mCounts_Ratio_pt[ifile]->SetTitle(Form("Centrality %d",i_cent));
    h_mCounts_Ratio_pt[ifile]->GetXaxis()->SetTitle("p_{T}");
    h_mCounts_Ratio_pt[ifile]->GetYaxis()->SetTitle("Counts");
    h_mCounts_Ratio_pt[ifile]->SetMarkerStyle(marker[ifile]);
    h_mCounts_Ratio_pt[ifile]->SetMarkerColor(color[ifile]);
    h_mCounts_Ratio_pt[ifile]->SetLineColor(color[ifile]);
    h_mCounts_Ratio_pt[ifile]->GetYaxis()->SetRangeUser(0.0,1.2);
    legratio->AddEntry(h_mCounts_Ratio_pt[ifile],Form("%s/%s",label[ifile].c_str(),label[0].c_str()),"p");
    if (ifile == 1) h_mCounts_Ratio_pt[ifile]->Draw("pE");
    if (ifile >  1) h_mCounts_Ratio_pt[ifile]->Draw("pE same");
  }
 
  legratio->Draw("same"); 

  c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_ratioptspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// pT 1D ratio //////////////////////

  ///////////// pT 1D ratio //////////////////////
  for(int iy = 0; iy < nbinsy; iy++)
  {
    c2->cd(iy+1);
    c2->cd(iy+1)->SetLogy(0);
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mCounts_Ratio_ptbin[ifile][iy]->SetTitle(Form("Centrality %d, %d/%d<y<%d/%d",i_cent,(iy-nbinsy/2),nbinsy/2,(iy-nbinsy/2+1),nbinsy/2));
      h_mCounts_Ratio_ptbin[ifile][iy]->GetXaxis()->SetTitle("p_{T}");
      h_mCounts_Ratio_ptbin[ifile][iy]->GetYaxis()->SetTitle("Counts");
      h_mCounts_Ratio_ptbin[ifile][iy]->SetMarkerStyle(marker[ifile]);
      h_mCounts_Ratio_ptbin[ifile][iy]->SetMarkerColor(color[ifile]);
      h_mCounts_Ratio_ptbin[ifile][iy]->SetLineColor(color[ifile]);
      h_mCounts_Ratio_ptbin[ifile][iy]->GetYaxis()->SetRangeUser(0.0,1.2);
      //legratio->AddEntry(h_mCounts_Ratio_ptbin[ifile][iy],Form("%s/%s",label[ifile].c_str(),label[0].c_str()),"p");
      if (ifile == 1) h_mCounts_Ratio_ptbin[ifile][iy]->Draw("pE");
      if (ifile >  1) h_mCounts_Ratio_ptbin[ifile][iy]->Draw("pE same");
    }
 
    legratio->Draw("same"); 
  }
  c2->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_individualybins_ratioptspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// pT 1D ratio //////////////////////

  ///////////// y 1D ///////////////////////
  c->cd();
  c->cd()->SetLogy(0);
  TLegend *legy = new TLegend(0.4,0.1,0.6,0.3);

  h_mCounts_y[0]->SetTitle(Form("Centrality %d",i_cent));
  h_mCounts_y[0]->GetXaxis()->SetTitle("y");
  h_mCounts_y[0]->GetYaxis()->SetTitle("Counts");
  h_mCounts_y[0]->SetMarkerStyle(20);
  h_mCounts_y[0]->SetMarkerColor(color[0]);
  h_mCounts_y[0]->SetLineColor(color[0]);
  double max = h_mCounts_y[0]->GetMaximum();
  double min = h_mCounts_y[0]->GetMinimum();
  legy->AddEntry(h_mCounts_y[0],label[0].c_str(),"p");

  for(int ifile = 1; ifile < nopts; ifile++)
  {
    if(h_mCounts_y[ifile]->GetMaximum() > max) max = h_mCounts_y[ifile]->GetMaximum();
    if(h_mCounts_y[ifile]->GetMinimum() < min) min = h_mCounts_y[ifile]->GetMinimum();
  }
  h_mCounts_y[0]->GetYaxis()->SetRangeUser(min*0.9,max*1.1);
  h_mCounts_y[0]->Draw("pE");
  for(int ifile = 1; ifile < nopts; ifile++)
  {
    h_mCounts_y[ifile]->SetMarkerStyle(marker[ifile]);
    h_mCounts_y[ifile]->SetMarkerColor(color[ifile]);
    h_mCounts_y[ifile]->SetLineColor(color[ifile]);
    legy->AddEntry(h_mCounts_y[ifile],label[ifile].c_str(),"p");
    h_mCounts_y[ifile]->Draw("pE same");
  }
  legy->Draw("same"); 

  c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_yspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// y 1D ///////////////////////

  ///////////// y 1D individual y bins //////////////////////
  for(int ipt = 0; ipt < nbinsx; ipt++)
  {
    c2->cd(ipt+1);
    c2->cd(ipt+1)->SetLogy(0);
    c2->cd(ipt+1)->SetLeftMargin(0.15);
    c2->cd(ipt+1)->SetBottomMargin(0.15);
    c2->cd(ipt+1)->SetTicks(1,1);
    c2->cd(ipt+1)->SetGrid(0,0);
    h_mCounts_ybin[0][ipt]->SetTitle(Form("Centrality %d, %1.1f<p_{T}<%1.1f",i_cent,vmsa::ptRawStart[ipt],vmsa::ptRawStop[ipt]));
    h_mCounts_ybin[0][ipt]->GetXaxis()->SetTitle("p_{T}");
    h_mCounts_ybin[0][ipt]->GetYaxis()->SetTitle("Counts");
    h_mCounts_ybin[0][ipt]->SetMarkerStyle(20);
    h_mCounts_ybin[0][ipt]->SetMarkerColor(color[0]);
    h_mCounts_ybin[0][ipt]->SetLineColor(color[0]);
    double maxbin = h_mCounts_ybin[0][ipt]->GetMaximum();
    double minbin = h_mCounts_ybin[0][ipt]->GetMinimum();

    for(int ifile = 1; ifile < nopts; ifile++)
    {
      if(h_mCounts_ybin[ifile][ipt]->GetMaximum() > maxbin) maxbin = h_mCounts_ybin[ifile][ipt]->GetMaximum();
      if(h_mCounts_ybin[ifile][ipt]->GetMinimum() < minbin) minbin = h_mCounts_ybin[ifile][ipt]->GetMinimum();
    }
    h_mCounts_ybin[0][ipt]->GetYaxis()->SetRangeUser(minbin*0.9,maxbin*1.1);
    h_mCounts_ybin[0][ipt]->Draw("pE");
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mCounts_ybin[ifile][ipt]->SetMarkerStyle(marker[ifile]);
      h_mCounts_ybin[ifile][ipt]->SetMarkerColor(color[ifile]);
      h_mCounts_ybin[ifile][ipt]->SetLineColor(color[ifile]);
      h_mCounts_ybin[ifile][ipt]->Draw("pE same");
    }
    legy->Draw("same"); 
  }
  c2->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_individualptbins_yspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// y 1D individual y bins //////////////////////

  ///////////// y 1D ratio //////////////////////
  c->cd();
  c->cd()->SetLogy(0);

  TLegend *legyratio = new TLegend(0.5,0.2,0.7,0.4);

  for(int ifile = 1; ifile < nopts; ifile++)
  {
    h_mCounts_Ratio_y[ifile]->SetTitle(Form("Centrality %d",i_cent));
    h_mCounts_Ratio_y[ifile]->GetXaxis()->SetTitle("y");
    h_mCounts_Ratio_y[ifile]->GetYaxis()->SetTitle("Counts");
    h_mCounts_Ratio_y[ifile]->SetMarkerStyle(marker[ifile]);
    h_mCounts_Ratio_y[ifile]->SetMarkerColor(color[ifile]);
    h_mCounts_Ratio_y[ifile]->SetLineColor(color[ifile]);
    h_mCounts_Ratio_y[ifile]->GetYaxis()->SetRangeUser(0.0,1.2);
    legyratio->AddEntry(h_mCounts_Ratio_y[ifile],Form("%s/%s",label[ifile].c_str(),label[0].c_str()),"p");
    if (ifile == 1) h_mCounts_Ratio_y[ifile]->Draw("pE");
    if (ifile >  1) h_mCounts_Ratio_y[ifile]->Draw("pE same");
  }
  legyratio->Draw("same"); 

  c->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_ratioyspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// y 1D ratio //////////////////////

  ///////////// pT 1D ratio //////////////////////
  for(int ipt = 0; ipt < nbinsx; ipt++)
  {
    c2->cd(ipt+1);
    for(int ifile = 1; ifile < nopts; ifile++)
    {
      h_mCounts_Ratio_ybin[ifile][ipt]->SetTitle(Form("Centrality %d, %1.1f<p_{T}<%1.1f",i_cent,vmsa::ptRawStart[ipt],vmsa::ptRawStop[ipt]));
      h_mCounts_Ratio_ybin[ifile][ipt]->GetXaxis()->SetTitle("y");
      h_mCounts_Ratio_ybin[ifile][ipt]->GetYaxis()->SetTitle("Counts");
      h_mCounts_Ratio_ybin[ifile][ipt]->SetMarkerStyle(marker[ifile]);
      h_mCounts_Ratio_ybin[ifile][ipt]->SetMarkerColor(color[ifile]);
      h_mCounts_Ratio_ybin[ifile][ipt]->SetLineColor(color[ifile]);
      h_mCounts_Ratio_ybin[ifile][ipt]->GetYaxis()->SetRangeUser(0.0,1.2);
      //legratio->AddEntry(h_mCounts_Ratio_ybin[ifile][ipt],Form("%s/%s",label[ifile].c_str(),label[0].c_str()),"p");
      if (ifile == 1) h_mCounts_Ratio_ybin[ifile][ipt]->Draw("pE");
      if (ifile >  1) h_mCounts_Ratio_ybin[ifile][ipt]->Draw("pE same");
    }
 
    legratio->Draw("same"); 
  }
  c2->SaveAs(Form("figures/%s/%s/pTstudy/SimulationComparisons/%s_individualptbins_ratioyspectra1D_%s_Order%d_%s_Cent%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),filefullstring.c_str(),correction.c_str(),order,etamode.c_str(),i_cent));  
  ///////////// pT 1D ratio //////////////////////


  //h_mCounts_Ratio_y->SetTitle(Form("Centrality %d",i_cent));
  //h_mCounts_Ratio_y->GetXaxis()->SetTitle("y");
  //h_mCounts_Ratio_y->GetYaxis()->SetTitle("Pythia/Embedding");
  //h_mCounts_Ratio_y->Draw("pE");
  
  //c->SaveAs(Form("figures/%s/%s/pTstudy/%s/PythiaEmbeddingRatio_yspectra1D_%s_Order%d_%s_Cent%d_Py%s_Em%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),  spectra.c_str(),correction.c_str(),order,etamode.c_str(),i_cent,sim_py.c_str(),sim_em.c_str()));  

}


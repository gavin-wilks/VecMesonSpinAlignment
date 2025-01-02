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

#ifndef _PlotQA_
#define _PlotQA_  1
#endif


void plotGlobal2DOffDiag(string correction = "raw", int energy = 4, int pid = 0, int year = 0, bool random3D = false, int order = 2, string etamode = "eta1_eta1")
{

  TGraMap g_mRho;
  TGraMap g_mReal;
  TGraMap g_mImag;
  TGraMap g_mReRho1n1;
  TGraMap g_mImRho1n1;
  
  std::string EP[2] = {"1st","2nd"};
  //string inputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/rho00/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //string inputfile = Form("../output/AuAu%s/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile = Form("../output/AuAu%s/%s/AccPhiPtSys_%s_PolySys_Global_2D_OffDiag.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  if( correction == "raw") inputfile = Form("../output/AuAu%s/%s/PhiPrime_RawPhiPtSys_%s_PolySys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  
  string inputfile1D = "../output/AuAu19GeV/Phi/Rho_RawSysErrors_F_0_eta1_eta1_PolySys.root";
  TFile *File_InPut1D = TFile::Open(inputfile1D.c_str());
  string StatErrorRho = Form("g_rho00_order%d_%s_%s_StatError",order,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());

  TCanvas *c_spin = new TCanvas("c_spin","c_spin",10,10,900,600);
  c_spin->Divide(3,2);
  for(int i = 0; i < 6; i++)
  {
    c_spin->cd(i+1);
    c_spin->cd(i+1)->SetLeftMargin(0.15);
    c_spin->cd(i+1)->SetBottomMargin(0.15);
    c_spin->cd(i+1)->SetTicks(1,1);
    c_spin->cd(i+1)->SetGrid(0,0);
  }
   
  string KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),          0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_real = Form("realRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),        0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_imag = Form("imagRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),        0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_rerho1n1 = Form("rerho1n1Raw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_imrho1n1 = Form("imrho1n1Raw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);

  c_spin->cd(4);
 
  TLegend *leg = new TLegend(0.2,0.2,0.5,0.4);

  g_mRho[KEY_rho] = (TGraphAsymmErrors*) File_InPut->Get(KEY_rho.c_str());
  //g_mRho[KEY_rho]->GetYaxis()->SetTitle("#Delta#rho_{00}");
  g_mRho[KEY_rho]->GetYaxis()->SetTitle("Raw #rho_{00}");
  g_mRho[KEY_rho]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_mRho[KEY_rho]->SetMarkerStyle(20); 
  g_mRho[KEY_rho]->SetMarkerSize(0.75); 
  leg->AddEntry(g_mRho[KEY_rho],"2D Method","p");
  g_mRho[KEY_rho]->Draw("APE");
  g_mRho[StatErrorRho] = (TGraphAsymmErrors*) File_InPut1D->Get(StatErrorRho.c_str());
  for(int i = 0; i < g_mRho[StatErrorRho]->GetN(); i++)
  { 
    double pt, rho; 
    g_mRho[StatErrorRho]->GetPoint(i, pt, rho);
    g_mRho[StatErrorRho]->SetPoint(i, pt+0.05, rho);
  }
  g_mRho[StatErrorRho]->SetMarkerSize(0.75); 
  g_mRho[StatErrorRho]->SetMarkerStyle(20); 
  g_mRho[StatErrorRho]->SetMarkerColor(kBlue); 
  g_mRho[StatErrorRho]->SetLineColor(kBlue); 
  leg->AddEntry(g_mRho[StatErrorRho],"1D Method","p");
  g_mRho[StatErrorRho]->Draw("PE");
  leg->Draw("same");

  c_spin->cd(2);
  g_mReal[KEY_real] = (TGraphAsymmErrors*) File_InPut->Get(KEY_real.c_str());
  g_mReal[KEY_real]->GetYaxis()->SetTitle("Re(#rho_{10})-Re(#rho_{0-1})");
  g_mReal[KEY_real]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_mReal[KEY_real]->SetMarkerStyle(20); 
  g_mReal[KEY_real]->SetMarkerSize(0.75); 
  g_mReal[KEY_real]->Draw("APE");

  c_spin->cd(3);
  g_mImag[KEY_imag] = (TGraphAsymmErrors*) File_InPut->Get(KEY_imag.c_str());
  g_mImag[KEY_imag]->GetYaxis()->SetTitle("Im(#rho_{10})-Im(#rho_{0-1})");
  g_mImag[KEY_imag]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_mImag[KEY_imag]->SetMarkerStyle(20); 
  g_mImag[KEY_imag]->SetMarkerSize(0.75); 
  g_mImag[KEY_imag]->Draw("APE");

  c_spin->cd(5);
  g_mReRho1n1[KEY_rerho1n1] = (TGraphAsymmErrors*) File_InPut->Get(KEY_rerho1n1.c_str());
  g_mReRho1n1[KEY_rerho1n1]->GetYaxis()->SetTitle("Re(#rho_{1-1})");
  g_mReRho1n1[KEY_rerho1n1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_mReRho1n1[KEY_rerho1n1]->SetMarkerStyle(20); 
  g_mReRho1n1[KEY_rerho1n1]->SetMarkerSize(0.75); 
  g_mReRho1n1[KEY_rerho1n1]->Draw("APE");

  c_spin->cd(6);
  g_mImRho1n1[KEY_imrho1n1] = (TGraphAsymmErrors*) File_InPut->Get(KEY_imrho1n1.c_str());
  g_mImRho1n1[KEY_imrho1n1]->GetYaxis()->SetTitle("Im(#rho_{1-1})");
  g_mImRho1n1[KEY_imrho1n1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_mImRho1n1[KEY_imrho1n1]->SetMarkerStyle(20); 
  g_mImRho1n1[KEY_imrho1n1]->SetMarkerSize(0.75); 
  g_mImRho1n1[KEY_imrho1n1]->Draw("APE");

  string outputname = Form("./figures/%s/%s/pTstudy/Global2D%s_OffDiag_ExtractedSpinDensityElements_%s_Order%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),vmsa::mBeamEnergy[energy].c_str(),order);
  c_spin->Print(outputname.c_str());

  TCanvas *c_spin2 = new TCanvas("c_spin2","c_spin2",10,10,600,300);
  c_spin2->Divide(2,1);
  for(int i = 0; i < 2; i++)
  {
    c_spin2->cd(i+1);
    c_spin2->cd(i+1)->SetLeftMargin(0.15);
    c_spin2->cd(i+1)->SetBottomMargin(0.15);
    c_spin2->cd(i+1)->SetTicks(1,1);
    c_spin2->cd(i+1)->SetGrid(0,0);
  }

  c_spin2->cd(1);
  g_mRho[KEY_rho]->Draw("APE");
  g_mRho[StatErrorRho]->Draw("PE");
  leg->Draw("same");

 
  TLegend *legoff = new TLegend(0.2,0.2,0.55,0.45);

  c_spin2->cd(2);
  g_mReal[KEY_real]->GetYaxis()->SetRangeUser(-0.065,0.055);
  g_mReal[KEY_real]->GetYaxis()->SetTitle("Raw Value");
  g_mReal[KEY_real]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g_mReal[KEY_real]->SetMarkerColor(kBlack); 
  g_mReal[KEY_real]->SetLineColor(  kBlack); 
  legoff->AddEntry(g_mReal[KEY_real],"Re(#rho_{10})-Re(#rho_{0-1})","p");
  g_mReal[KEY_real]->Draw("APE");

  for(int i = 0; i < g_mImag[KEY_imag]->GetN(); i++)
  { 
    double pt, rho; 
    g_mImag[KEY_imag]->GetPoint(i, pt, rho);
    g_mImag[KEY_imag]->SetPoint(i, pt+0.05, rho);
  }
  g_mImag[KEY_imag]->SetMarkerColor(kBlue);
  g_mImag[KEY_imag]->SetLineColor(  kBlue);
  legoff->AddEntry(g_mImag[KEY_imag],"Im(#rho_{10})-Im(#rho_{0-1})","p");
  g_mImag[KEY_imag]->Draw("PE");

  for(int i = 0; i < g_mReRho1n1[KEY_rerho1n1]->GetN(); i++)
  { 
    double pt, rho; 
    g_mReRho1n1[KEY_rerho1n1]->GetPoint(i, pt, rho);
    g_mReRho1n1[KEY_rerho1n1]->SetPoint(i, pt+0.1, rho);
  }
  g_mReRho1n1[KEY_rerho1n1]->SetMarkerColor(kOrange+7);
  g_mReRho1n1[KEY_rerho1n1]->SetLineColor(  kOrange+7);
  legoff->AddEntry(g_mReRho1n1[KEY_rerho1n1],"Re(#rho_{1-1})","p");
  g_mReRho1n1[KEY_rerho1n1]->Draw("PE");

  for(int i = 0; i < g_mImRho1n1[KEY_imrho1n1]->GetN(); i++)
  { 
    double pt, rho; 
    g_mImRho1n1[KEY_imrho1n1]->GetPoint(i, pt, rho);
    g_mImRho1n1[KEY_imrho1n1]->SetPoint(i, pt+0.15, rho);
  }
  g_mImRho1n1[KEY_imrho1n1]->SetMarkerColor(kGray+2);
  g_mImRho1n1[KEY_imrho1n1]->SetLineColor(  kGray+2);
  legoff->AddEntry(g_mImRho1n1[KEY_imrho1n1],"Im(#rho_{1-1})","p");
  g_mImRho1n1[KEY_imrho1n1]->Draw("PE");

  legoff->Draw("same");

  outputname = Form("./figures/%s/%s/pTstudy/Global2D%s_AllOnOnePlot_OffDiag_ExtractedSpinDensityElements_%s_Order%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),vmsa::mBeamEnergy[energy].c_str(),order);
  c_spin2->Print(outputname.c_str());
 
  
  
}

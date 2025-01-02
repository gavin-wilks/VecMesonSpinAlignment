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


void plotHelicity2DOffDiag_CompSim(int energy = 4, int pid = 0, int year = 0, bool random3D = false, int order = 2, string etamode = "eta1_eta1", string option = "_Embed")
{

  TGraMap g_mRho;
  TGraMap g_mReal;
  TGraMap g_mImag;
  TGraMap g_mReRho1n1;
  TGraMap g_mImRho1n1;

  TGraMap g2_mRho;
  TGraMap g2_mReal;
  TGraMap g2_mImag;
  TGraMap g2_mReRho1n1;
  TGraMap g2_mImRho1n1;
  
  std::string EP[2] = {"1st","2nd"};
  //string inputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/rho00/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //string inputfile = Form("../output/AuAu%s/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile = Form("../output/AuAu%s/%s/AccPhiPtSys_%s_PolySys_Helicity_2D_OffDiag.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  string inputfile2 = Form("../output/AuAu%s/%s/AccPhiPtSys_%s_PolySys_Helicity_2D_OffDiag%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str(),option.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TFile *File_InPut2 = TFile::Open(inputfile2.c_str());

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

  c_spin->cd(1);
  g_mRho[KEY_rho] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut->Get(KEY_rho.c_str()))->Clone();
  g2_mRho[KEY_rho] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut2->Get(KEY_rho.c_str()))->Clone();
  //g_mRho[KEY_rho]->GetYaxis()->SetTitle("#Delta#rho_{00}");
  g_mRho[KEY_rho]->GetYaxis()->SetTitle("#rho_{00}");
  g_mRho[KEY_rho]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g2_mRho[KEY_rho]->SetMarkerStyle(20); 
  g2_mRho[KEY_rho]->SetMarkerColor(kRed); 
  g2_mRho[KEY_rho]->SetLineColor(kRed); 
  g_mRho[KEY_rho]->SetMarkerStyle(20); 
  g_mRho[KEY_rho]->Draw("APE");
  g2_mRho[KEY_rho]->Draw("PE same");

  c_spin->cd(2);
  g_mReal[KEY_real] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut->Get(KEY_real.c_str()))->Clone();
  g2_mReal[KEY_real] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut2->Get(KEY_real.c_str()))->Clone();
  g_mReal[KEY_real]->GetYaxis()->SetTitle("Re(#rho_{10})-Re(#rho_{0-1})");
  g_mReal[KEY_real]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g2_mReal[KEY_real]->SetMarkerStyle(20); 
  g2_mReal[KEY_real]->SetMarkerColor(kRed); 
  g2_mReal[KEY_real]->SetLineColor(kRed); 
  g_mReal[KEY_real]->SetMarkerStyle(20); 
  g_mReal[KEY_real]->Draw("APE");
  g2_mReal[KEY_real]->Draw("PE same");

  c_spin->cd(3);
  g_mImag[KEY_imag] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut->Get(KEY_imag.c_str()))->Clone();
  g2_mImag[KEY_imag] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut2->Get(KEY_imag.c_str()))->Clone();
  g_mImag[KEY_imag]->GetYaxis()->SetTitle("Im(#rho_{10})-Im(#rho_{0-1})");
  g_mImag[KEY_imag]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g2_mImag[KEY_imag]->SetMarkerStyle(20); 
  g2_mImag[KEY_imag]->SetMarkerColor(kRed); 
  g2_mImag[KEY_imag]->SetLineColor(kRed); 
  g_mImag[KEY_imag]->SetMarkerStyle(20); 
  g_mImag[KEY_imag]->Draw("APE");
  g2_mImag[KEY_imag]->Draw("PE same");

  c_spin->cd(4);
  g_mReRho1n1[KEY_rerho1n1] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut->Get(KEY_rerho1n1.c_str()))->Clone();
  g2_mReRho1n1[KEY_rerho1n1] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut2->Get(KEY_rerho1n1.c_str()))->Clone();
  g_mReRho1n1[KEY_rerho1n1]->GetYaxis()->SetTitle("Re(#rho_{1-1})");
  g_mReRho1n1[KEY_rerho1n1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g2_mReRho1n1[KEY_rerho1n1]->SetMarkerStyle(20); 
  g2_mReRho1n1[KEY_rerho1n1]->SetMarkerColor(kRed); 
  g2_mReRho1n1[KEY_rerho1n1]->SetLineColor(kRed); 
  g_mReRho1n1[KEY_rerho1n1]->SetMarkerStyle(20); 
  g_mReRho1n1[KEY_rerho1n1]->Draw("APE");
  g2_mReRho1n1[KEY_rerho1n1]->Draw("PE same");

  c_spin->cd(5);
  g_mImRho1n1[KEY_imrho1n1] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut->Get(KEY_imrho1n1.c_str()))->Clone();
  g2_mImRho1n1[KEY_imrho1n1] = (TGraphAsymmErrors*)((TGraphAsymmErrors*) File_InPut2->Get(KEY_imrho1n1.c_str()))->Clone();
  g_mImRho1n1[KEY_imrho1n1]->GetYaxis()->SetTitle("Im(#rho_{1-1})");
  g_mImRho1n1[KEY_imrho1n1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  g2_mImRho1n1[KEY_imrho1n1]->SetMarkerStyle(20); 
  g2_mImRho1n1[KEY_imrho1n1]->SetMarkerColor(kRed); 
  g2_mImRho1n1[KEY_imrho1n1]->SetLineColor(kRed); 
  g_mImRho1n1[KEY_imrho1n1]->SetMarkerStyle(20); 
  g_mImRho1n1[KEY_imrho1n1]->Draw("APE");
  g2_mImRho1n1[KEY_imrho1n1]->Draw("PE same");

  string outputname = Form("./figures/%s/%s/pTstudy/SimulationComparison_Helicity2DCorrected_OffDiag_ExtractedSpinDensityElements_%s_Order%d%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,option.c_str());
  c_spin->Print(outputname.c_str());
  
}

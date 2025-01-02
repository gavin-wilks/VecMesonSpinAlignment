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


void plot2DOffDiagCent_Comparison(string correction = "Eff", int energy = 4, int pid = 0, int year = 0, bool random3D = false, int order = 2, string etamode = "eta1_eta1", int frameopt = 0)
{

  std::string frame = "Global";
  if(frameopt == 1) frame = "Helicity";

  TGraMap g_mRho;
  
  std::string EP[2] = {"1st","2nd"};
  //string inputfile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/%s/rho00/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //string inputfile = Form("../output/AuAu%s/%s/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile    = Form("../output/AuAu%s/%s/%s2DCent_RawPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),frame.c_str(),etamode.c_str());
  string inputfileEff = Form("../output/AuAu%s/%s/%s2DCent_EffPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),frame.c_str(),etamode.c_str());

  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TFile *File_InPutEff = TFile::Open(inputfileEff.c_str());
  
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
 
  string KEYS[5] = {"rho","real","imag","rerho1n1","imrho1n1"}; 

  for(int i = 0; i < 5; i++)
  {
    string KEY     = Form("%sRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",KEYS[i].c_str(),9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    string KEY_Raw = Form("Raw_%sRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",KEYS[i].c_str(),9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    string KEY_Eff = Form("Eff_%sRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",KEYS[i].c_str(),9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);

    g_mRho[KEY_Raw] = (TGraphAsymmErrors*) File_InPut->Get(KEY.c_str());
    g_mRho[KEY_Eff] = (TGraphAsymmErrors*) File_InPutEff->Get(KEY.c_str());

    for(int i = 0; i < g_mRho[KEY_Eff]->GetN(); i++)
    { 
      double pt, rho; 
      g_mRho[KEY_Eff]->GetPoint(i, pt, rho);
      g_mRho[KEY_Eff]->SetPoint(i, pt+1.5, rho);
    }

  }  

  c_spin->cd(4);
 
  TLegend *leg = new TLegend(0.2,0.2,0.5,0.4);

  string xtitle[5] = {"#Delta#rho_{00}","Re(#rho_{10})-Re(#rho_{0-1})","Im(#rho_{10})-Im(#rho_{0-1})","Re(#rho_{1-1})","Im(#rho_{1-1})"};

  for(int i = 0; i < 5; i++)
  {
    string KEY_Raw = Form("Raw_%sRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",KEYS[i].c_str(),9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    string KEY_Eff = Form("Eff_%sRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",KEYS[i].c_str(),9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);

    c_spin->cd(i+1);

    g_mRho[KEY_Raw]->GetYaxis()->SetTitle(xtitle[i].c_str());
    g_mRho[KEY_Raw]->GetXaxis()->SetTitle("Centrality (%)");
    g_mRho[KEY_Raw]->SetMarkerStyle(20); 
    g_mRho[KEY_Raw]->SetMarkerSize(0.75); 
    g_mRho[KEY_Eff]->SetMarkerSize(0.75); 
    g_mRho[KEY_Eff]->SetMarkerStyle(20); 
    g_mRho[KEY_Eff]->SetMarkerColor(kBlue); 
    g_mRho[KEY_Eff]->SetLineColor(kBlue); 
    if(i == 0) leg->AddEntry(g_mRho[KEY_Raw],"Raw","p");
    if(i == 0) leg->AddEntry(g_mRho[KEY_Eff],"Eff Corrected","p");
    g_mRho[KEY_Raw]->Draw("APE");
    g_mRho[KEY_Eff]->Draw("PE");
    leg->Draw("same");

  }

  string outputname = Form("./figures/%s/%s/centralitystudy/%s2DCent_%s_OffDiag_ExtractedSpinDensityElements_%s_Order%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),correction.c_str(),vmsa::mBeamEnergy[energy].c_str(),order);
  c_spin->Print(outputname.c_str());

  //TCanvas *c_spin2 = new TCanvas("c_spin2","c_spin2",10,10,600,300);
  //c_spin2->Divide(2,1);
  //for(int i = 0; i < 2; i++)
  //{
  //  c_spin2->cd(i+1);
  //  c_spin2->cd(i+1)->SetLeftMargin(0.15);
  //  c_spin2->cd(i+1)->SetBottomMargin(0.15);
  //  c_spin2->cd(i+1)->SetTicks(1,1);
  //  c_spin2->cd(i+1)->SetGrid(0,0);
  //}

  //c_spin2->cd(1);
  //g_mRho[KEY_rho]->Draw("APE");
  //g_mRho[StatErrorRho]->Draw("PE");
  //leg->Draw("same");

 
  //TLegend *legoff = new TLegend(0.2,0.2,0.55,0.45);

  //c_spin2->cd(2);
  //g_mReal[KEY_real]->GetYaxis()->SetRangeUser(-0.065,0.055);
  //g_mReal[KEY_real]->GetYaxis()->SetTitle("Raw Value");
  //g_mReal[KEY_real]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //g_mReal[KEY_real]->SetMarkerColor(kBlack); 
  //g_mReal[KEY_real]->SetLineColor(  kBlack); 
  //legoff->AddEntry(g_mReal[KEY_real],"Re(#rho_{10})-Re(#rho_{0-1})","p");
  //g_mReal[KEY_real]->Draw("APE");

  //for(int i = 0; i < g_mImag[KEY_imag]->GetN(); i++)
  //{ 
  //  double pt, rho; 
  //  g_mImag[KEY_imag]->GetPoint(i, pt, rho);
  //  g_mImag[KEY_imag]->SetPoint(i, pt+0.05, rho);
  //}
  //g_mImag[KEY_imag]->SetMarkerColor(kBlue);
  //g_mImag[KEY_imag]->SetLineColor(  kBlue);
  //legoff->AddEntry(g_mImag[KEY_imag],"Im(#rho_{10})-Im(#rho_{0-1})","p");
  //g_mImag[KEY_imag]->Draw("PE");

  //for(int i = 0; i < g_mReRho1n1[KEY_rerho1n1]->GetN(); i++)
  //{ 
  //  double pt, rho; 
  //  g_mReRho1n1[KEY_rerho1n1]->GetPoint(i, pt, rho);
  //  g_mReRho1n1[KEY_rerho1n1]->SetPoint(i, pt+0.1, rho);
  //}
  //g_mReRho1n1[KEY_rerho1n1]->SetMarkerColor(kOrange+7);
  //g_mReRho1n1[KEY_rerho1n1]->SetLineColor(  kOrange+7);
  //legoff->AddEntry(g_mReRho1n1[KEY_rerho1n1],"Re(#rho_{1-1})","p");
  //g_mReRho1n1[KEY_rerho1n1]->Draw("PE");

  //for(int i = 0; i < g_mImRho1n1[KEY_imrho1n1]->GetN(); i++)
  //{ 
  //  double pt, rho; 
  //  g_mImRho1n1[KEY_imrho1n1]->GetPoint(i, pt, rho);
  //  g_mImRho1n1[KEY_imrho1n1]->SetPoint(i, pt+0.15, rho);
  //}
  //g_mImRho1n1[KEY_imrho1n1]->SetMarkerColor(kGray+2);
  //g_mImRho1n1[KEY_imrho1n1]->SetLineColor(  kGray+2);
  //legoff->AddEntry(g_mImRho1n1[KEY_imrho1n1],"Im(#rho_{1-1})","p");
  //g_mImRho1n1[KEY_imrho1n1]->Draw("PE");

  //legoff->Draw("same");

  //outputname = Form("./figures/%s/%s/pTstudy/Global2D%s_AllOnOnePlot_OffDiag_ExtractedSpinDensityElements_%s_Order%d.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),vmsa::mBeamEnergy[energy].c_str(),order);
  //c_spin2->Print(outputname.c_str());
 
  
  
}

#include <string>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/type.h"
#include "../Utility/draw.h"
#include "../Utility/functions.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

double FuncAD(double *x_val, double *par);

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE);

void calSysError(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Eff", bool random3D = true)
{
  string inputfileHframe = Form("../output/AuAu%s/%s/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  if(random3D) inputfile = Form("../output/AuAu%s/%s/3DRandom/%sPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  TFile *File_InPutHframe = TFile::Open(inputfileHframe.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TGraMap g_mRho;
  TH1DMap h_mRho;  
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
	    string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	    g_mRho[KEY_rho] = (TGraphAsymmErrors*)File_InPut->Get(KEY_rho.c_str())->Clone();
	    string KEY_rho_ptbin = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	    h_mRho[KEY_rho_ptbin] = (TH1D*)File_InPut->Get(KEY_rho_ptbin.c_str());
	  }
	}
      }
    }
  }
  TH1F *h_frame = (TH1F*)File_InPutHframe->Get("h_frame");

#if _PlotQA_
  /*if(correction == "AccRes")
  {
    TCanvas *c_rhocorr = new TCanvas("c_rhocorr","c_rhocorr",10,10,800,800);
    c_rhocorr->cd();
    c_rhocorr->cd()->SetLeftMargin(0.15);
    c_rhocorr->cd()->SetBottomMargin(0.15);
    c_rhocorr->cd()->SetTicks(1,1);
    c_rhocorr->cd()->SetGrid(0,0);

    TF1 *Func_rho = new TF1("Func_rho",FuncAD,0,1,4);
    string key = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",2,9,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
    TH1F *PtCos = (TH1F*)File_InPut->Get(key.c_str())->Clone();
        PtCos->GetXaxis()->SetTitleOffset(1.2);
    PtCos->GetXaxis()->SetTitle("cos#theta*");
    PtCos->GetYaxis()->SetTitle("yield");
    PtCos->GetYaxis()->SetTitleOffset(1.0);
    PtCos->SetMarkerColor(2);
    PtCos->SetMarkerSize(1.8);
    PtCos->SetMarkerStyle(21);
    PtCos->Draw("pE");
    Func_rho->SetParameter(0, PtCos->GetBinContent(5));
    Func_rho->SetParameter(1, 0.3);
    Func_rho->FixParameter(2, 0.00663655); //hardcoded
    Func_rho->FixParameter(3, 0.331673); //hardcoded
    PtCos->Fit(Func_rho, "NMI"); 
    Func_rho->Draw("same");
  }*/
  
  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  h_frame->DrawCopy("pE");
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

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
	    string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,1.1);
	  }
	}
      }
    }
  }
#endif
 
  string KEY_Default = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
  cout << "DEFAULT: " << KEY_Default << endl;
  TGraphAsymmErrors *g_SysErrors = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_StatErrors = new TGraphAsymmErrors();
    
  double sysErr[6][5]; // N points with 5 sources of systematics for each

  for(Int_t i_point = 0; i_point < 6; ++i_point)
  {
    double sysDca[vmsa::Dca_stop];
    double sysNSig[vmsa::nSigKaon_stop];
    double sysNorm[vmsa::Norm_stop];
    double sysSig[vmsa::Sig_stop];
    double sysMeth[vmsa::Method_stop];

    double pt_def, rho_def;
    g_mRho[KEY_Default]->GetPoint(i_point,pt_def,rho_def); 

    cout << "pt: " << pt_def << "   rho00: " << rho_def;

    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    { 
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysDca[i_dca] = rho_sys;
    }

    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
    {
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,i_sig,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysNSig[i_sig] = rho_sys;
    }

    for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
    {
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),i_norm,0,vmsa::mInteMethod[1].c_str());
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysNorm[i_norm] = rho_sys;
    }	 

    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[1].c_str());
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysSig[i_sigma] = rho_sys;
    }

    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
    {
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str());
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysMeth[i_method] = rho_sys;
    }

    Double_t rho_min[5] = { TMath::MinElement(vmsa::Dca_stop,sysDca),
                            TMath::MinElement(vmsa::nSigKaon_stop,sysNSig),
                            TMath::MinElement(vmsa::Norm_stop,sysNorm),
                            TMath::MinElement(vmsa::Sig_start,sysSig),
                            TMath::MinElement(vmsa::Method_stop,sysMeth)    };

    Double_t rho_max[5] = { TMath::MaxElement(vmsa::Dca_stop,sysDca),
                            TMath::MaxElement(vmsa::nSigKaon_stop,sysNSig),
                            TMath::MaxElement(vmsa::Norm_stop,sysNorm),
                            TMath::MaxElement(vmsa::Sig_start,sysSig),
                            TMath::MaxElement(vmsa::Method_stop,sysMeth)    };
  
    double SysError_rho = 0.0;
    for(int i = 0; i < 5; i++)
    {
      double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
      cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
      SysError_rho += sourcei;
    }
 
    SysError_rho = TMath::Sqrt(SysError_rho);

    Double_t pt, rho;
    g_mRho[KEY_Default]->GetPoint(i_point,pt,rho);

    //float mean_rho = total_rho/(float)counter;
    g_SysErrors->SetPoint(i_point,pt,rho);
    g_SysErrors->SetPointError(i_point,0.0,0.0,SysError_rho,SysError_rho);

    double StatError_rho = g_mRho[KEY_Default]->GetErrorYhigh(i_point);
    g_StatErrors->SetPoint(i_point,pt,rho);
    g_StatErrors->SetPointError(i_point,0.0,0.0,StatError_rho,StatError_rho);
  }
 
  cout << "Finished pT bin calculations" << endl;
  if(correction == "AccRes")
  { 
  string KEY_DefaultW = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());

  double sysDca[vmsa::Dca_stop];
  double sysNSig[vmsa::nSigKaon_stop];
  double sysNorm[vmsa::Norm_stop];
  double sysSig[vmsa::Sig_stop];
  double sysMeth[vmsa::Method_stop];

  double rho_def = h_mRho[KEY_DefaultW]->GetBinContent(7); 
  double rho_defErr = h_mRho[KEY_DefaultW]->GetBinError(7); 

  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  { 
    string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
    sysDca[i_dca] = h_mRho[KEY]->GetBinContent(7);
  }

  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
  {
    string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,i_sig,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
    sysNSig[i_sig] = h_mRho[KEY]->GetBinContent(7);
  }

  for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
  {
    string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),i_norm,0,vmsa::mInteMethod[1].c_str());
    sysNorm[i_norm] = h_mRho[KEY]->GetBinContent(7);
  }	 

  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
  {
    string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[1].c_str());
    sysSig[i_sigma] = h_mRho[KEY]->GetBinContent(7);
  }

  for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
  {
    string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str());
    sysMeth[i_method] = h_mRho[KEY]->GetBinContent(7);
  }

  Double_t rho_min[5] = { TMath::MinElement(vmsa::Dca_stop,sysDca),
                          TMath::MinElement(vmsa::nSigKaon_stop,sysNSig),
                          TMath::MinElement(vmsa::Norm_stop,sysNorm),
                          TMath::MinElement(vmsa::Sig_start,sysSig),
                          TMath::MinElement(vmsa::Method_stop,sysMeth)    };

  Double_t rho_max[5] = { TMath::MaxElement(vmsa::Dca_stop,sysDca),
                          TMath::MaxElement(vmsa::nSigKaon_stop,sysNSig),
                          TMath::MaxElement(vmsa::Norm_stop,sysNorm),
                          TMath::MaxElement(vmsa::Sig_start,sysSig),
                          TMath::MaxElement(vmsa::Method_stop,sysMeth)    };
  
  double SysError_rho = 0.0;
  for(int i = 0; i < 5; i++)
  {
    double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
    cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
    SysError_rho += sourcei;
  }
 
  SysError_rho = TMath::Sqrt(SysError_rho);

  cout << "Final corrected rho00 = " << rho_def << " +/- " << rho_defErr << " (stat) +/- " << SysError_rho << " (sys)" << endl;
  
  double val[6];
  double sysE[6];
  double statE[6];

  for(int ipt = 0; ipt < 6; ipt++)
  {
    double pt,rho;
    g_StatErrors->GetPoint(ipt,pt,rho);
    val[ipt] = rho;
    statE[ipt] = g_StatErrors->GetErrorYhigh(ipt);  
    sysE[ipt] = g_SysErrors->GetErrorYhigh(ipt);  
  }

  cout << "double ptbincenter={0.6,0.8,1.5,2.1,2.7,3.6};" << endl;
  cout << "rho00={"; for (int i = 0; i < 6; i++) {if (i<5) cout << val[i] << ","; else cout << val[i] << "};" << endl;}
  cout << "stat={"; for (int i = 0; i < 6; i++) {if (i<5) cout << statE[i] << ","; else cout << statE[i] << "};" << endl;}
  cout << "sys={"; for (int i = 0; i < 6; i++) {if (i<5) cout << sysE[i] << ","; else cout << sysE[i] << "};" << endl;}
  }
  TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
  TGraphAsymmErrors *besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");

  //cout << "All good" << endl;

  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,800,800);
  c_rho_SysError->cd();
  c_rho_SysError->cd()->SetLeftMargin(0.15);
  c_rho_SysError->cd()->SetBottomMargin(0.15);
  c_rho_SysError->cd()->SetTicks(1,1);
  c_rho_SysError->cd()->SetGrid(0,0);
  //h_frame->GetXaxis()->SetNdivisions(505,"I");
  h_frame->GetYaxis()->SetRangeUser(0.23,0.51);
  //h_frame->GetXaxis()->SetLimits(0,9.95);
  //h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  //h_frame->GetXaxis()->SetLimits(0,10);
  h_frame->Draw("pE");
  //g_StatErrors->SetMarkerStyle(20);
  //g_StatErrors->SetMarkerColor(kGray+2);
  //g_StatErrors->SetLineColor(2);
  //g_StatErrors->Draw("pE same");
  //g_SysErrors->SetMarkerStyle(20);
  //g_SysErrors->SetMarkerColor(kGray+2);
  //g_SysErrors->SetLineColor(kGray+2);
  //g_SysErrors->Draw("pE same");
//  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  for(int i = 0; i < 4; i++)
  {
    besi19_sys->SetPointEXhigh(i,0.0);
    besi19_sys->SetPointEXlow(i,0.0);
    besi19->SetPointEXhigh(i,0.0);
    besi19->SetPointEXlow(i,0.0);
  }
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);
  gStyle->SetEndErrorSize(8);
  //gStyle->SetEndErrorWidth(10);
  besi19->SetLineColor(kBlack);
  besi19->SetLineWidth(1.0);
  besi19->SetMarkerStyle(20);
  besi19->SetMarkerSize(1.3);
  besi19->SetMarkerColor(kBlack);
  besi19->SetLineColor(kBlack);
  besi19->Draw("pE Z same");


  besi19_sys->SetLineColor(kBlack);
  besi19_sys->SetMarkerStyle(1);
  besi19_sys->SetMarkerSize(10);
  besi19_sys->SetMarkerColor(kBlack);
  besi19_sys->SetLineColor(kBlack);
  besi19_sys->Draw("pE [] same");

  Draw_TGAE_Point_new_Symbol(0.5,0.46,0.0,0.0,0.0,0.0,20,kRed,1.3); 
  string leg_count;
  if(correction == "Eff") leg_count = "Efficiency Corrected #phi BES-II";
  if(correction == "Raw") leg_count = "Raw #phi BES-II";
  if(correction == "AccRes") leg_count = "#phi BES-II";
  plotTopLegend((char*)leg_count.c_str(),0.65,0.457,0.03,1,0.0,42,0);

  string leg_besi = "#phi BES-I";
  Draw_TGAE_Point_new_Symbol(0.5,0.44,0.0,0.0,0.0,0.0,20,kBlack,1.3);
  plotTopLegend((char*)leg_besi.c_str(),0.65,0.437,0.03,1,0.0,42,0);
 
  //string leg_sp = "STAR Preliminary";
  //plotTopLegend((char*)leg_sp.c_str(),0.3,0.48,0.03,1,0.0,42,0);

  PlotLine(0.25,0.6,0.418,0.418,1,2,2);
  string leg_line = "#rho_{00} = 1/3";
  plotTopLegend((char*)leg_line.c_str(),0.65,0.417,0.03,1,0.0,42,0);

  string leg_energy = Form("Au+Au %s %d", vmsa::mBeamEnergy[energy].c_str(), vmsa::mBeamYear[energy]);
  plotTopLegend((char*)leg_energy.c_str(),1.5,0.265,0.04,1,0.0,42,0);
  plotTopLegend((char*)"20%-60%",2.0,0.245,0.04,1,0.0,42,0);

  g_StatErrors->SetMarkerStyle(20);
  g_StatErrors->SetMarkerColor(kRed);
  g_StatErrors->SetLineColor(kRed);
  g_StatErrors->SetMarkerSize(1.3);

  g_SysErrors->SetMarkerStyle(20);
  g_SysErrors->SetMarkerSize(1.3);
  g_SysErrors->SetMarkerColor(kRed);
  g_SysErrors->SetLineColor(kRed);
  //g_SysErrors->Draw("pE || same");
  g_StatErrors->SetPoint(0,-1.0,-1.0);
  g_SysErrors->SetPoint(0,-1.0,-1.0);
  g_StatErrors->SetPoint(1,-1.0,-1.0);
  g_SysErrors->SetPoint(1,-1.0,-1.0);
  g_StatErrors->Draw("pE Z same");
  //g_SysErrors->Draw("pE [] same");

  plotSysErrorsBox(g_SysErrors,kRed,energy);

  for (int ipt = 0; ipt < 4; ipt++)
  {
    double pt;
    double b1val, b1stat, b1sys;
    double b2val, b2stat, b2sys;
    besi19->GetPoint(ipt,pt,b1val);
    b1stat = besi19->GetErrorYhigh(ipt);
    b1sys = besi19_sys->GetErrorYhigh(ipt);
    g_StatErrors->GetPoint(ipt+2,pt,b2val);
    b2stat = g_StatErrors->GetErrorYhigh(ipt+2);
    b2sys = g_SysErrors->GetErrorYhigh(ipt+2);

    double sigma = TMath::Abs((b2val-b1val)/TMath::Sqrt(b1stat*b1stat + b2stat*b2stat + b1sys*b1sys + b2sys*b2sys));
    std::cout << std::fixed << std::setprecision(4);
    cout << "pt = "<<  pt << ":    BESI = " << b1val << " +/- " << b1stat << " (stat) +/- " << b1sys << "(sys)    BESII = " << b2val << " +/- " << b2stat << " (stat) +/- " << b2sys << "(sys)    sigma = " << std::setprecision(2) << sigma << endl; 
  }

  string OutPutFile = Form("../output/AuAu%s/%s/Rho_%sSysErrors.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  if(random3D) OutPutFile = Form("../output/AuAu%s/%s/3DRandom/Rho_%sSysErrors.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  cout << "OutPutFile set to: " << OutPutFile.c_str() << endl;
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  TString StatErrorRho = Form("g_rho00_%s_%s_StatError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  g_StatErrors->SetName(StatErrorRho.Data());
  g_StatErrors->Write();
  TString SysErrorRho = Form("g_rho00_%s_%s_SysError",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  g_SysErrors->SetName(SysErrorRho.Data());
  g_SysErrors->Write();
  File_OutPut->Close();
}


void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE)
{
  const int nEnergy = 6;
  TBox *bSys[nEnergy];
  for(int i_energy = 0; i_energy < 6; ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);
    
    //bSys[i_energy] = new TBox(energy/1.08,rho-err,energy*1.08,rho+err);
    bSys[i_energy] = new TBox(vmsa::pt_low[beamE][i_energy],rho-err,vmsa::pt_up[beamE][i_energy],rho+err);
    bSys[i_energy]->SetFillColor(0);
    bSys[i_energy]->SetFillStyle(0);
    bSys[i_energy]->SetLineStyle(1);
    bSys[i_energy]->SetLineWidth(1);
    bSys[i_energy]->SetLineColor(plot_color);
    bSys[i_energy]->Draw("l Same");
  }
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

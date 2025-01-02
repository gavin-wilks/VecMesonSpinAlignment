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

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

double FuncAD(double *x_val, double *par);

//void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void calSysError_FixedSys(Int_t energy = 4, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", bool random3D = false, int order = 1, int defaultF = 0, string etamode = "eta1_eta1")//defaultF = 0 is BESII, defaultF = 1 is BESI
{

  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;



  string inputfileHframe = Form("../output/AuAu%s/%s/RawPhiPtSys_eta1_eta1.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_FINALFORPRELIM.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/PtDependence_v2off_pToff/%sPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/PtDependenceEffThenAcc_EPRes/%sPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  if(energy == 4 && correction == "AccRes" && order == 2) inputfile = Form("../output/AuAu%s/%s/FixedResolutionFinal/%sPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  if(energy == 3 && correction == "AccRes" && order == 2) inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_%s_noSmear_EffAcc.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_sub_noDelta_2ndSubEPRes.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  //if(order == 1) inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_Psi1.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  //inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_4th.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  //inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_2ndUpdated.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  if(random3D) inputfile = Form("../output/AuAu%s/%s/3DRandom/%sPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  TFile *File_InPutHframe = TFile::Open(inputfileHframe.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TGraMap g_mRho;
  TH1DMap h_mRho;  
  TH1DMap h_mCounts;  
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
            for(int i_F = 0; i_F < 1; i_F++)
            {
              for(int i_eff = 0; i_eff < 1; i_eff++)
              {
                if(correction != "AccRes")
                {
                  if(i_F > 0 || i_eff > 0) continue;
	          string KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	          g_mRho[KEY_rho] = (TGraphAsymmErrors*)File_InPut->Get(KEY_rho.c_str())->Clone();
	          string KEY_rho_ptbin = Form("rhoFinalWeighted_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	          h_mRho[KEY_rho_ptbin] = (TH1D*)File_InPut->Get(KEY_rho_ptbin.c_str());
	          for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt) // pt loop
	          {
	            string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	            h_mCounts[KEY_counts] = (TH1D*) File_InPut->Get(KEY_counts.c_str());
	          }
                }
                else
                {
	          string KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F,i_eff);
	          g_mRho[KEY_rho] = (TGraphAsymmErrors*)File_InPut->Get(KEY_rho.c_str())->Clone();
	          string KEY_rho_ptbin = Form("rhoFinalWeighted_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F,i_eff);
	          h_mRho[KEY_rho_ptbin] = (TH1D*)File_InPut->Get(KEY_rho_ptbin.c_str());
                }
              }
            }
	  }
	}
      }
    }
  }
  TH1F *h_frame = (TH1F*)File_InPutHframe->Get("h_frame");

#if _PlotQA_
  if(correction == "AccRes")
  {
    TCanvas *c_rhocorr = new TCanvas("c_rhocorr","c_rhocorr",10,10,800,600);
    c_rhocorr->cd();
    c_rhocorr->cd()->SetLeftMargin(0.15);
    c_rhocorr->cd()->SetBottomMargin(0.15);
    c_rhocorr->cd()->SetTicks(1,1);
    c_rhocorr->cd()->SetGrid(0,0);

    TH1F *h_frame_cos = new TH1F("h_frame_cos","h_frame_cos",100,-0.05,9.95);
    for(int i_bin = 0; i_bin < 100; ++i_bin)
    {
      h_frame_cos->SetBinContent(i_bin+1,-10.0);
      h_frame_cos->SetBinError(i_bin+1,1.0);
    }
    h_frame_cos->SetTitle("");
    h_frame_cos->SetStats(0);
    h_frame_cos->GetXaxis()->SetRangeUser(0.0,1.0);
    h_frame_cos->GetXaxis()->SetNdivisions(505,'N');
    h_frame_cos->GetXaxis()->SetLabelSize(0.03);
    h_frame_cos->GetXaxis()->SetTitle("|cos#theta*|");
    h_frame_cos->GetXaxis()->SetTitleSize(0.05);
    h_frame_cos->GetXaxis()->SetTitleOffset(1.2);
    h_frame_cos->GetXaxis()->CenterTitle();

    h_frame_cos->GetYaxis()->SetRangeUser(0.00335,0.00360);
    h_frame_cos->GetYaxis()->SetNdivisions(505,'N');
    h_frame_cos->GetYaxis()->SetTitle("#frac{dN^{corr}_{#phi}}{d(cos#theta*)} (arb. units)");
    h_frame_cos->GetYaxis()->SetTitleSize(0.05);
    h_frame_cos->GetYaxis()->SetMaxDigits(2);
    h_frame_cos->GetYaxis()->SetLabelSize(0.03);
    h_frame_cos->GetYaxis()->CenterTitle();
    h_frame_cos->DrawCopy("pE");

    Draw_TGAE_Point_new_Symbol(0.07,0.0035325,0.0,0.0,0.0,0.0,20,kBlack,1.3); 
    string leg_count = "#phi Data";
    plotTopLegend((char*)leg_count.c_str(),0.1,0.00353,0.04,1,0.0,42,0);
    string leg_pt = "|y| < 1.0, 1.2 < p_{T} < 1.8 GeV/c";
    plotTopLegend((char*)leg_pt.c_str(),0.1,0.00351,0.04,1,0.0,42,0);

    string leg_sp = "STAR Preliminary";
    plotTopLegend((char*)leg_sp.c_str(),0.1,0.00358,0.04,2,0.0,42,0);

    PlotLine(0.04,0.085,0.0035525,0.0035525,kRed,2,1);
    string leg_line = "N_{0}[(1-#rho_{00}^{obs})+(3#rho_{00}^{obs}-1)cos^{2}#theta*]";
    plotTopLegend((char*)leg_line.c_str(),0.1,0.00355,0.04,1,0.0,42,0);

    string leg_energy = Form("Au+Au %s %d", vmsa::mBeamEnergyText[energy].c_str(), vmsa::mBeamYear[energy]);
    plotTopLegend((char*)leg_energy.c_str(),0.55,0.003375,0.04,1,0.0,42,0);
    plotTopLegend((char*)"20%-60%",0.65,0.00336,0.04,1,0.0,42,0);

    TF1 *Func_rho = new TF1("Func_rho","[0]*((1.-[1])+(3.*[1]-1.)*x*x)",0,1);
    string key = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",2,9,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),defaultF);
    TH1F *PtCos = (TH1F*)File_InPut->Get(key.c_str())->Clone();
    for(int i = 1; i <= 7; i++)
    {
      double y = PtCos->GetBinContent(i);
      double yerr = PtCos->GetBinError(i);
      double yscaled = y/459885291.;
      double yscalederr = yerr/459885291.;
      PtCos->SetBinContent(i,yscaled);
      PtCos->SetBinError(i,yscalederr);
    }
    //PtCos->GetXaxis()->SetTitleOffset(1.2);
    //PtCos->GetXaxis()->SetTitle("cos#theta*");
    //PtCos->GetYaxis()->SetTitle("yield");
    //PtCos->GetYaxis()->SetTitleOffset(1.0);
    PtCos->SetMarkerColor(kBlack);
    PtCos->SetMarkerSize(1.3);
    PtCos->SetMarkerStyle(20);
    cout << "Minimum = " << PtCos->GetMinimum() << endl;
    cout << "Maximum = " << PtCos->GetMaximum() << endl;
    PtCos->Draw("pEX0 same");
    Func_rho->SetParameter(0, PtCos->GetBinContent(5));
    Func_rho->SetParameter(1, 0.333);
    //Func_rho->FixParameter(2, 0.0160525); //hardcoded
    //Func_rho->FixParameter(3, 0.331673); //hardcoded
    PtCos->Fit(Func_rho, "NMI"); 
    Func_rho->SetLineColor(kRed);
    Func_rho->Draw("same");
    c_rhocorr->SaveAs("figures/BESII_19p6GeV_2060_CorrectedYieldsCos.eps");
    c_rhocorr->SaveAs("figures/BESII_19p6GeV_2060_CorrectedYieldsCos.pdf");
  }
  
  /*TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
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
            for(int i_F = 0; i_F < 2; i_F++)
            {
              if(correction != "AccRes")
              {
                if(i_F > 0) continue;
	        string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
	        Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,1.1);
              }
              else
              {
	        string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F);
	        Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,1.1);
              }
            }
	  }
	}
      }
    }
  }*/
#endif
 
  string KEY_Default = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
  if(correction == "AccRes") KEY_Default = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),defaultF);
  cout << "DEFAULT: " << KEY_Default << endl;
  TGraphAsymmErrors *g_SysErrors = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_StatErrors = new TGraphAsymmErrors();
    
  double sysErr[6][6]; // N points with 5 sources of systematics for each

  double weight_rho00 = 0.0;
  double weight_error_stat_rho00 = 0.0;
  double weight_all = 0.0;

  for(Int_t i_point = 2; i_point < 6; ++i_point)
  {
    double sysDca[9];
    double sysNSig[9];
    double sysNorm[9];
    //double sysSig[vmsa::Sig_stop];
    //double sysMeth[vmsa::Method_stop];
    double sysF[5];
    double sysEff[5];

    double pt_def, rho_def;
    g_mRho[KEY_Default]->GetPoint(i_point,pt_def,rho_def); 
    double rho_defErr = g_mRho[KEY_Default]->GetErrorYhigh(i_point); 
    string KEY_Default_yield = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_point,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
    double weight;
   // TH1F *PtCos = (TH1F*)h_mCounts[KEY_Default_yield]->Clone();
   // weight = PtCos->Integral(1,7);
   
    weight_rho00 += rho_def/rho_defErr/rho_defErr;
    weight_all += 1./rho_defErr/rho_defErr;

    cout << "pt: " << pt_def << "   rho00: " << rho_def << endl;
    
    //sysDca[0] = rho_def;   
    //sysNSig[0] = rho_def;   
    //sysNorm[0] = rho_def;   
    //sysF[0] = rho_def;   
  
    int idx = 0;
    for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
    { 
      for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
      {
        for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
        {
          if(i_dca == 0 && (i_sigma != 0 || i_method == 0)) continue;
          if(i_dca != 0 && i_sigma != 0 && i_method == 1) continue;
          string KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
          if(correction == "AccRes") KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_cent,EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str(),defaultF);
          double pt_sys, rho_sys;
          g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
          sysDca[idx] = rho_sys;
          if(random3D) sysDca[idx] = 4.*(sysDca[idx]-1./3.)+1./3.;
          idx++;
        }
      }
    }
    idx = 0;

    for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
    {
      for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
      {
        for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
        {
          if(i_sig == 0 && (i_sigma != 0 || i_method == 0)) continue;
          if(i_sig != 0 && i_sigma != 0 && i_method == 1) continue;
          string KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
          if(correction == "AccRes") KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_cent,EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
          double pt_sys, rho_sys;
          g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
          sysNSig[idx] = rho_sys;
          if(random3D) sysNSig[idx] = 4.*(sysNSig[idx]-1./3.)+1./3.;
          idx++;
        }
      }
    }
    idx = 0;

    for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
    {
      for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
      {
        for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
        {
          if(i_norm == 0 && (i_sigma != 0 || i_method == 0)) continue;
          if(i_norm != 0 && i_sigma != 0 && i_method == 1) continue;
          string KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
          if(correction == "AccRes") KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
          double pt_sys, rho_sys;
          g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
          sysNorm[idx] = rho_sys;
          if(random3D) sysNorm[idx] = 4.*(sysNorm[idx]-1./3.)+1./3.;
          idx++;
        }
      }
    }	 
    idx = 0;

    /*for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[1].c_str());
      if(correction == "AccRes") KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[1].c_str(),defaultF);
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysSig[i_sigma] = rho_sys;
      if(random3D) sysSig[i_sigma] = 4.*(sysSig[i_sigma]-1./3.)+1./3.;
    }

    for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
    {
      string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str());
      if(correction == "AccRes") KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str(),defaultF);
      double pt_sys, rho_sys;
      g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
      sysMeth[i_method] = rho_sys;
      if(random3D) sysMeth[i_method] = 4.*(sysMeth[i_method]-1./3.)+1./3.;
    }*/

    if(correction == "AccRes")
    {
      /*for(int i_F = 0; i_F < 2; ++i_F)
      {
        for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
        {
          for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
          {
            if(i_F == 0 && (i_sigma != 0 || i_method == 0)) continue;
            if(i_F != 0 && i_sigma != 0 && i_method == 1) continue;
            string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F);
            double pt_sys, rho_sys;
            g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
            sysF[idx] = rho_sys;
            if(random3D) sysF[idx] = 4.*(sysF[idx]-1./3.)+1./3.;
            idx++;
          }
        }
      }*/
      idx = 0;
 
      /*for(int i_eff = 0; i_eff < 2; ++i_eff)
      {
        for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
        {
          for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
          {
            if(i_eff == 0 && (i_sigma != 0 || i_method == 0)) continue;
            if(i_eff != 0 && i_sigma != 0 && i_method == 1) continue;
            string KEY = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_%d",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF,i_eff);
            double pt_sys, rho_sys;
            g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
            sysEff[idx] = rho_sys;
            if(random3D) sysEff[idx] = 4.*(sysEff[idx]-1./3.)+1./3.;
            idx++;
          }
        }
      }*/
    }
    idx = 0;

    Double_t rho_min[5] = { TMath::MinElement(9,sysDca),
                            TMath::MinElement(9,sysNSig),
                            TMath::MinElement(9,sysNorm)
                            //(correction == "AccRes")? TMath::MinElement(5,sysF) : 0.0,
                            //(correction == "AccRes")? TMath::MinElement(5,sysEff) : 0.0
                          };

    Double_t rho_max[5] = { TMath::MaxElement(9,sysDca),
                            TMath::MaxElement(9,sysNSig),
                            TMath::MaxElement(9,sysNorm)
                            //(correction == "AccRes")? TMath::MaxElement(5,sysF) : 0.0,
                            //(correction == "AccRes")? TMath::MaxElement(5,sysEff) : 0.0
                          };
  
    double SysError_rho = 0.0;
    for(int i = 0; i < 3; i++)
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

  weight_rho00 /= weight_all;
  weight_error_stat_rho00 = TMath::Sqrt(1./weight_all);
  cout << std::setprecision(4);
  cout << "Final rho00 = " << weight_rho00 << " +/- " << weight_error_stat_rho00 << " (stat.) " << endl;
 
  cout << "Finished pT bin calculations" << endl;
 
  vecFMap vRho;
  vecFMap vRhoErr;
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
            for(int i_F = 0; i_F < 1; i_F++)
            {
              for(int i_eff = 0; i_eff < 1; i_eff++)
              {
                string KEY_rho;
	        if(correction == "AccRes") KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F,i_eff);
	        if(correction == "Raw") KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
                //cout << KEY_rho << endl;
                double weight_rho00f = 0.0;
                double weight_error_stat_rho00f = 0.0;
                double weight_allf = 0.0;

                for(Int_t i_point = 2; i_point < 6; ++i_point)
                {
                  //cout << "point = " << i_point << endl;
                  Double_t pt, rho;
                  g_mRho[KEY_rho]->GetPoint(i_point,pt,rho);
                  //cout << "pt = " << pt << ", rho = " << rho << endl;
                  double StatError_rho = g_mRho[KEY_rho]->GetErrorYhigh(i_point);
                  weight_rho00f += rho/StatError_rho/StatError_rho;
                  weight_allf += 1./StatError_rho/StatError_rho; 
                }
                weight_rho00f /= weight_allf;
                weight_error_stat_rho00f = TMath::Sqrt(1./weight_allf);
                vRho[KEY_rho].push_back(weight_rho00f);
                vRhoErr[KEY_rho].push_back(weight_error_stat_rho00f); 
              }
            }
          }
        }
      }
    }
  }

  double rhofinal, rhofinalstat, rhofinalsys;

  std::ofstream fw_sources(Form("ptDependenceSources_Integrated2060_Order%d_%s.txt",order,vmsa::mBeamEnergy[energy].c_str()), std::ofstream::out);
  //fw.open();

  fw_sources << "sqrt{s_{nn}}      EP_Order      drho_{dca}        drho_{nσ_{kaon}}       drho_{NormalizationRange}" << "\n";
  std::ofstream fw(Form("ptDependence_Integrated2060_Order%d_%s.txt",order,vmsa::mBeamEnergy[energy].c_str()), std::ofstream::out);
  //fw.open();

  fw << "sqrt{s_{nn}}      EP_Order      dca(cm)        nσ_{kaon}        NormalizationRange         IntegrationMethod          IntegrationRange             ρ_{00}" << "\n";
  if(correction == "AccRes")
  { 
  string KEY_DefaultW = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),defaultF);

  double sysDca[9];
  double sysNSig[9];
  double sysNorm[9];
  double sysF[5];
  double sysEff[5];

  double rho_def = vRho[KEY_DefaultW][0]; 
  double rho_defErr = vRhoErr[KEY_DefaultW][0];//h_mRho[KEY_DefaultW]->GetBinError(7); 

  rhofinal = rho_def;
  rhofinalstat = rho_defErr;  

  int idx = 0;
  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  { 
    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
      {
        if(i_dca == 0 && (i_sigma != 0 || i_method == 0)) continue;
        if(i_dca != 0 && i_sigma != 0 && i_method == 1) continue;
        string KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_cent,EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
        sysDca[idx] = vRho[KEY][0];
        fw << vmsa::mBeamEnergy[energy] << "     " << order << "     " << std::fixed << std::setprecision(1) << vmsa::mDcaSys[i_dca] << "     " << vmsa::mNSigmaKaonSys[0] << "      " << "[" << std::setprecision(2) << vmsa::Norm_Start[0][0] << "," << vmsa::Norm_Stop[0][0] << "]" << "     " << vmsa::mInteMethod[i_method] << "     " << "[" << std::setprecision(1) << -vmsa::nSigVecSys[i_sigma] << "Γ," << vmsa::nSigVecSys[i_sigma] << "Γ]" << "       " << std::setprecision(4) << sysDca[idx] << "\n"; 
        idx++;
      }
    }
  }
  idx = 0;

  for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < vmsa::nSigKaon_stop; i_sig++)
  {
    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
      {
        if(i_sig == 0 && (i_sigma != 0 || i_method == 0)) continue;
        if(i_sig != 0 && i_sigma != 0 && i_method == 1) continue;
        string KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_cent,EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
        sysNSig[idx] = vRho[KEY][0];//h_mRho[KEY]->GetBinContent(7);
        fw << vmsa::mBeamEnergy[energy] << "     " << order << "     " << std::fixed << std::setprecision(1) << vmsa::mDcaSys[0] << "     " << vmsa::mNSigmaKaonSys[i_sig] << "      " << "[" << std::setprecision(2) << vmsa::Norm_Start[0][0] << "," << vmsa::Norm_Stop[0][0] << "]" << "     " << vmsa::mInteMethod[i_method] << "     " << "[" << std::setprecision(1) << -vmsa::nSigVecSys[i_sigma] << "Γ," << vmsa::nSigVecSys[i_sigma] << "Γ]" << "       " << std::setprecision(4) << sysNSig[idx] << "\n"; 
        idx++;
      }
    }
  }
  idx = 0;

  for(int i_norm = vmsa::Norm_start; i_norm < vmsa::Norm_stop; ++i_norm)
  {
    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
      {
        if(i_norm == 0 && (i_sigma != 0 || i_method == 0)) continue;
        if(i_norm != 0 && i_sigma != 0 && i_method == 1) continue;
        cout << idx << "   i_norm = " << i_norm << "   i_sigma = " << i_sigma << "   i_method = " << i_method << endl;
        string KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
        sysNorm[idx] = vRho[KEY][0];//h_mRho[KEY]->GetBinContent(7);
        if(i_norm < 2) fw << vmsa::mBeamEnergy[energy] << "     " << order << "     " << std::fixed << std::setprecision(1) << vmsa::mDcaSys[0] << "     " << vmsa::mNSigmaKaonSys[0] << "      " << "[" << std::setprecision(2) << vmsa::Norm_Start[0][i_sigma] << "," << vmsa::Norm_Stop[0][i_sigma] << "]" << "     " << vmsa::mInteMethod[i_method] << "     " << "[" << std::setprecision(1) << -vmsa::nSigVecSys[i_sigma] << "Γ," << vmsa::nSigVecSys[i_sigma] << "Γ]" << "       " << std::setprecision(4) << sysNorm[idx] << "\n"; 
        if(i_norm == 2) fw << vmsa::mBeamEnergy[energy] << "     " << order << "     " << std::fixed << std::setprecision(1) << vmsa::mDcaSys[0] << "     " << vmsa::mNSigmaKaonSys[0] << "      " << "[" << std::setprecision(2) << vmsa::Norm_Start[0][0] << "," << vmsa::Norm_Stop[0][0] << "]" << "," << "[" << std::setprecision(2) << vmsa::Norm_Start[0][1] << "," << vmsa::Norm_Stop[0][1] << "]" << "     " << vmsa::mInteMethod[i_method] << "     " << "[" << std::setprecision(1) << -vmsa::nSigVecSys[i_sigma] << "Γ," << vmsa::nSigVecSys[i_sigma] << "Γ]" << "       " << std::setprecision(4) << sysNorm[idx] << "\n"; 
        idx++;
      }
    }
  }	 
  idx = 0;

  /*for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
  {
    string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[1].c_str(),defaultF);
    sysSig[i_sigma] = h_mRho[KEY]->GetBinContent(7);
  }

  for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
  {
    string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str(),defaultF);
    sysMeth[i_method] = h_mRho[KEY]->GetBinContent(7);
  }
*/

  /*for(int i_F = 0; i_F < 2; ++i_F)
  {
    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
      {
        if(i_F == 0 && (i_sigma != 0 || i_method == 0)) continue;
        if(i_F != 0 && i_sigma != 0 && i_method == 1) continue;
        string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F);
        sysF[idx] = h_mRho[KEY]->GetBinContent(7);
        idx++;
      }
    }
  }
  idx = 0;*/
/*
  for(int i_eff = 0; i_eff < 2; ++i_eff)
  {
    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
      {
        if(i_eff == 0 && (i_sigma != 0 || i_method == 0)) continue;
        if(i_eff != 0 && i_sigma != 0 && i_method == 1) continue;
        string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_%d",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF,i_eff);
        sysEff[idx] = h_mRho[KEY]->GetBinContent(7);
        idx++;
      }
    }
  }
  idx = 0;
*/
  Double_t rho_min[5] = { TMath::MinElement(9,sysDca),
                          TMath::MinElement(9,sysNSig),
                          TMath::MinElement(9,sysNorm)
                          //TMath::MinElement(5,sysF),    
                          //TMath::MinElement(5,sysEff)    
                        };

  Double_t rho_max[5] = { TMath::MaxElement(9,sysDca),
                          TMath::MaxElement(9,sysNSig),
                          TMath::MaxElement(9,sysNorm)
                          //TMath::MaxElement(5,sysF),    
                          //TMath::MaxElement(5,sysEff)    
                        };
  fw.close();
  
  fw_sources << vmsa::mBeamEnergy[energy] << "     " << order << "     " << std::fixed << std::setprecision(5);
  double SysError_rho = 0.0;
  for(int i = 0; i < 3; i++)
  {
    double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
    cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
    SysError_rho += sourcei;
    fw_sources << TMath::Sqrt(sourcei) << "     ";
  }
  fw_sources << "\n"; 
  fw_sources.close(); 

  SysError_rho = TMath::Sqrt(SysError_rho);
  
  rhofinalsys = SysError_rho;

  cout << "Final corrected rho00 = " << std::fixed << std::setprecision(4) << rho_def << " +/- " << rho_defErr << " (stat) +/- " << SysError_rho << " (sys)" << endl;
  cout << "sigma from 1/3 = " << std::fixed << std::setprecision(2) << TMath::Abs((rho_def-1./3.)/TMath::Sqrt(rho_defErr*rho_defErr + SysError_rho*SysError_rho)) << endl;
  //cout << "sigma from BESI = " << TMath::Abs((rho_def-0.37)/TMath::Sqrt(0.008*0.008 + 0.007*0.007 + rho_defErr*rho_defErr + SysError_rho*SysError_rho)) << endl;
  cout << "sigma from BESII 2nd order = " << TMath::Abs((rho_def-0.3505)/TMath::Sqrt(0.0024*0.0024 + 0.0025*0.0025 + rho_defErr*rho_defErr + SysError_rho*SysError_rho)) << endl;
  cout << "sigma from BESII 2nd order (stat only) = " << TMath::Abs((rho_def-0.3505)/TMath::Sqrt(0.0024*0.0024 + rho_defErr*rho_defErr)) << endl;
  }
  double val[6];
  double sysE[6];
  double statE[6];

  for(int ipt = 2; ipt < 6; ipt++)
  {
    double pt,rho;
    g_StatErrors->GetPoint(ipt,pt,rho);
    val[ipt] = rho;
    statE[ipt] = g_StatErrors->GetErrorYhigh(ipt);  
    sysE[ipt] = g_SysErrors->GetErrorYhigh(ipt); 
    
    if(random3D)
    {
      val[ipt] = 4.*(rho-1./3.)+1./3.;
      statE[ipt] = 4.*g_StatErrors->GetErrorYhigh(ipt);  
      sysE[ipt] = g_SysErrors->GetErrorYhigh(ipt); 
      g_StatErrors->SetPoint(ipt,pt,val[ipt]);
      g_SysErrors->SetPoint(ipt,pt,val[ipt]);
      g_StatErrors->SetPointError(ipt,0.0,0.0,statE[ipt],statE[ipt]); 
    }
  }

  cout << "double ptbincenter={0.6,0.8,1.5,2.1,2.7,3.6};" << endl;
  cout << "rho00={"; for (int i = 2; i < 6; i++) {if (i<5) cout << val[i] << ","; else cout << val[i] << "};" << endl;}
  cout << "stat={"; for (int i = 2; i < 6; i++) {if (i<5) cout << statE[i] << ","; else cout << statE[i] << "};" << endl;}
  cout << "sys={"; for (int i = 2; i < 6; i++) {if (i<5) cout << sysE[i] << ","; else cout << sysE[i] << "};" << endl;}
  
  TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *besi19;  
  TGraphAsymmErrors *besi19_sys;
  if(order == 2) 
  {
    besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
    besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");
  }
  if(order == 1) 
  {
    besi19 = (TGraphAsymmErrors*)besi->Get("rho00_1stEP_pt_stat_19;1");
    besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_1stEP_pt_sys_19;1");
  }


  TGraph *besi_stat;
  TGraph *besi_sys; 
  if(order == 1)
  {
    besi_stat = (TGraph*)besi->Get("rho00_1stEP_energy_stat;1");
    besi_sys = (TGraph*)besi->Get("rho00_1stEP_energy_sys;1");
  }
  if(order == 2)
  {
    besi_stat = (TGraph*)besi->Get("rho00_2ndEP_energy_stat;1");
    besi_sys = (TGraph*)besi->Get("rho00_2ndEP_energy_sys;1");
  }


  for(int i = 0; i < besi_stat->GetN(); i++)
  {
    double en1,val1,err1;
    double en2,val2,err2;
    besi_stat->GetPoint(i,en1,val1);
    err1 = besi_stat->GetErrorY(i);
    besi_sys->GetPoint(i,en2,val2);
    err2 = besi_sys->GetErrorY(i);
    cout << std::setprecision(6);
    cout << "energy = " << en1 << ", value = " << val1 << ", stat = " << err1 << endl;
    cout << "energy = " << en2 << ", value = " << val2 << ",  sys = " << err2 << endl;
    if(en1 == 19.6 || i == 1) 
    {
      double sigma = TMath::Abs((rhofinal-val1)/TMath::Sqrt(err1*err1 + err2*err2 + rhofinalstat*rhofinalstat + rhofinalsys*rhofinalsys));
      double sigma13 = TMath::Abs((1./3.-val1)/TMath::Sqrt(err1*err1 + err2*err2));
      cout << "energy = " << en1 << ", value = " << rhofinal << ", stat = " << rhofinalstat << ", sys = " << rhofinalsys << "    BESII" << endl;
      cout << "sigma from BES-I = " << sigma << endl;
      cout << "BES-I sigma from 1/3 = " << sigma13 << endl;
    }
  }

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

  if(energy != 3){
  for(int i = 0; i < 4; i++)
  {
    besi19_sys->SetPointEXhigh(i,0.0);
    besi19_sys->SetPointEXlow(i,0.0);
    besi19->SetPointEXhigh(i,0.0);
    besi19->SetPointEXlow(i,0.0);
  }
  }
  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);
  gStyle->SetEndErrorSize(3);
  //gStyle->SetEndErrorWidth(10);
  if(energy != 3){
  besi19->SetLineColor(kBlack);
  besi19->SetLineWidth(1.0);
  besi19->SetMarkerStyle(20);
  besi19->SetMarkerSize(1.3);
  besi19->SetMarkerColor(kBlack);
  besi19->SetLineColor(kBlack);
  //besi19->Draw("pE Z same");


  besi19_sys->SetLineColor(kBlack);
  besi19_sys->SetMarkerStyle(1);
  besi19_sys->SetMarkerSize(10);
  besi19_sys->SetMarkerColor(kBlack);
  besi19_sys->SetLineColor(kBlack);
  //besi19_sys->Draw("pE [] same");
  }
  Draw_TGAE_Point_new_Symbol(0.5,0.46,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2); 
  string leg_count;
  if(correction == "Eff") leg_count = "Efficiency Corrected #phi BES-II";
  if(correction == "Raw") leg_count = "Raw #phi BES-II (|y| < 1.0)";
  if(correction == "AccRes") leg_count = "#phi BES-II (|y| < 1.0)";
  plotTopLegend((char*)leg_count.c_str(),0.65,0.457,0.03,1,0.0,42,0);

  if(energy != 3){
  string leg_besi = "#phi BES-I  (|y| < 1.0)";
  Draw_TGAE_Point_new_Symbol(0.5,0.44,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
  plotTopLegend((char*)leg_besi.c_str(),0.65,0.437,0.03,1,0.0,42,0);
  }
  string leg_sp = "STAR Preliminary";
  plotTopLegend((char*)leg_sp.c_str(),0.3,0.48,0.03,2,0.0,42,0);

  PlotLine(0.25,0.6,0.418,0.418,1,2,2);
  string leg_line = "#rho_{00} = 1/3";
  plotTopLegend((char*)leg_line.c_str(),0.65,0.417,0.03,1,0.0,42,0);

  string leg_energy = Form("Au+Au %s", vmsa::mBeamEnergyText[energy].c_str());
  plotTopLegend((char*)leg_energy.c_str(),1.65,0.265,0.04,1,0.0,42,0);
  plotTopLegend((char*)"20%-60%",2.0,0.245,0.04,1,0.0,42,0);

  g_StatErrors->SetMarkerStyle(20);
  g_StatErrors->SetMarkerColor(kRed);
  g_StatErrors->SetLineColor(kRed);
  g_StatErrors->SetMarkerSize(1.3);

  g_SysErrors->SetMarkerStyle(20);
  g_SysErrors->SetMarkerSize(1.3);
  g_SysErrors->SetMarkerColor(kRed);
  g_SysErrors->SetLineColor(kRed);
  //g_SysErrors->Draw("pE [] same");
  //if(correction == "AccRes") 
  //{
  g_StatErrors->SetPoint(0,-1.0,-1.0);
  g_SysErrors->SetPoint(0,-1.0,-1.0);
  g_StatErrors->SetPoint(1,-1.0,-1.0);
  g_SysErrors->SetPoint(1,-1.0,-1.0);
  //g_SysErrors->Draw("pE [] same");
  //}

  //g_StatErrors->SetLineWidth(2);
  //besi19->SetLineWidth(2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_StatErrors,style_phi_2nd,color_phi_2nd,colorDiff_phi,size_marker+0.2);
  if(energy != 3){
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)besi19,style_phi_1st,color_phi_1st,colorDiff_phi,size_marker-0.2);
  }
  // plotSysErrors(g_rho_2nd_sys[total_pad-1],color_phi_2nd);
  plotSysErrorsBox(g_SysErrors,color_phi_2nd);
  if(energy != 3){
  plotSysErrorsBox(besi19_sys,color_phi_1st);
  }
  //g_StatErrors->Draw("pE same");
  //plotSysErrorsBox(g_SysErrors,kRed,energy);
  c_rho_SysError->SaveAs(Form("figures/%s/%s/pTstudy/finaloutputrho_%s_Order%d_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str()));  
  c_rho_SysError->SaveAs(Form("figures/%s/%s/pTstudy/finaloutputrho_%s_Order%d_%s.eps",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str()));  
  c_rho_SysError->SaveAs(Form("figures/%s/%s/pTstudy/finaloutputrho_%s_Order%d_%s.png",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),order,etamode.c_str()));  

  /*for (int ipt = 0; ipt < 4; ipt++)
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
  }*/
  if(random3D)
  {
    for (int ipt = 0; ipt < 6; ipt++)
    {
      double pt;
      double b2val, b2stat, b2sys;
      g_StatErrors->GetPoint(ipt,pt,b2val);
      b2stat = g_StatErrors->GetErrorYhigh(ipt);
      b2sys = g_SysErrors->GetErrorYhigh(ipt);

      double sigma = TMath::Abs((b2val-1./3.)/TMath::Sqrt(b2stat*b2stat + b2sys*b2sys));
      std::cout << std::fixed << std::setprecision(4);
      cout << "pt = "<<  pt << ":   BESII = " << b2val << " +/- " << b2stat << " (stat) +/- " << b2sys << "(sys)    sigma = " << std::setprecision(2) << sigma << endl; 
    }
  }

  string OutPutFile = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),defaultF,etamode.c_str());
  if(order == 1) OutPutFile = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_%d_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),defaultF,etamode.c_str());
  //string OutPutFile = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),defaultF,etamode.c_str());
  //OutPutFile = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_%d_2ndUpdated.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),defaultF);
  //string OutPutFile = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_%d_4th.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),defaultF);
  if(random3D) OutPutFile = Form("../output/AuAu%s/%s/3DRandom/Rho_%sSysErrors.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  cout << "OutPutFile set to: " << OutPutFile.c_str() << endl;
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  TString StatErrorRho = Form("g_rho00_order%d_%s_%s_StatError",order,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  g_StatErrors->SetName(StatErrorRho.Data());
  g_StatErrors->Write();
  TString SysErrorRho = Form("g_rho00_order%d_%s_%s_SysError",order,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  g_SysErrors->SetName(SysErrorRho.Data());
  g_SysErrors->Write();
  File_OutPut->Close();
}

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nPt = g_rho->GetN();
  TBox *bSys[nPt];
  for(int i_pt = 0; i_pt < g_rho->GetN(); ++i_pt) // plot sys errors
  {
    double pt, rho;
    g_rho->GetPoint(i_pt,pt,rho);
    double err = g_rho->GetErrorYhigh(i_pt);

    bSys[i_pt] = new TBox(pt-0.08,rho-err,pt+0.08,rho+err);
    bSys[i_pt]->SetFillColor(0);
    bSys[i_pt]->SetFillStyle(0);
    bSys[i_pt]->SetLineStyle(1);
    bSys[i_pt]->SetLineWidth(1);
    bSys[i_pt]->SetLineColor(plot_color);
    bSys[i_pt]->Draw("l Same");
  }
}

/*void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE)
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
}*/

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

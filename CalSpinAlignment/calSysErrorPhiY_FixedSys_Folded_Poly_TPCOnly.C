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
//#include "../Utility/draw.h"
#include "../../../draw.h"
#include "../Utility/functions.h"
#include "resolution_y.h"

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

double FuncAD(double *x_val, double *par);

//void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void calSysErrorPhiY_FixedSys_Folded_Poly_TPCOnly(Int_t energy = 4, Int_t pid = 0, string correction = "Raw", bool random3D = false, int ptQA = 0, int centQA = 2, int order = 2, std::string etamode = "eta1_eta1", std::string fileoption = "", int yspectra = -1)//defaultF = 0 is BESII, defaultF = 1 is BESI
{

  double Res_12[3] = {0.0};
  double Res_12_err[3] = {0.0};

  if(energy == 3 && order == 1)
  {
    for(int i = 0; i < 3; i++)
    {
      Res_12[i] = r14GeVep1[i][0];
      Res_12_err[i] = r14GeVep1[i][1];
    }
  }
  if(energy == 3 && order == 2)
  {
    for(int i = 0; i < 3; i++)
    {
      Res_12[i] = r14GeVep2[i][0];
      Res_12_err[i] = r14GeVep2[i][1];
    }
  }
  if(energy == 4 && order == 1)
  {
    for(int i = 0; i < 3; i++)
    {
      Res_12[i] = r19GeVep1[i][0];
      Res_12_err[i] = r19GeVep1[i][1];
    }
  }
  if(energy == 4 && order == 2)
  {
    for(int i = 0; i < 3; i++)
    {
      Res_12[i] = r19GeVep2[i][0];
      Res_12_err[i] = r19GeVep2[i][1];
    }
  }
  
  std::string spectra = "";
  if(yspectra == 0) spectra = "_NoRapiditySpectra";
  if(yspectra == 1) spectra = "_NoRapiditySpectra_EP";
  if(yspectra == 2) spectra = "_NoRapiditySpectra_EP_noV2";
  if(yspectra == 3) spectra = "_WithRapiditySpectra";
  if(yspectra == 4) spectra = "_WithRapiditySpectra_HalfSigma";
  if(yspectra == 5) spectra = "_NoRapiditySpectra_FixedFirstEP";
  if(yspectra == 6) spectra = "_NoRapiditySpectra_FixedFirstEP_PhiPsi";
  if(yspectra == 7) spectra = "_NoRapiditySpectra_PhiPsi_EP";

  std::string EP[2] = {"1st","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  int ypadding = 0;
  if(etamode == "eta1_eta1" && (correction == "AccRes" || correction == "Eff" || correction == "Acc")) ypadding = 2;
  if(etamode == "eta0p4" && (correction == "AccRes" || correction == "Eff" || correction == "Acc")) ypadding = 5;
  if(etamode == "eta0p6" && (correction == "AccRes" || correction == "Eff" || correction == "Acc")) ypadding = 4;
  if(etamode == "eta0p8" && (correction == "AccRes" || correction == "Eff" || correction == "Acc")) ypadding = 3;


  string inputfileHframe = Form("../output/AuAu%s/%s/TPCOnly_RawRhoEtaSys_%s_Poly.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/%sRhoEtaSys_%s%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),fileoption.c_str());
  string inputfile = Form("../output/AuAu%s/%s/TPCOnly_%sRhoEtaSys_%s_Poly.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/%s/%sRhoEtaSys_%s_Poly_FirstOrder%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),fileoption.c_str(),correction.c_str(),etamode.c_str(),spectra.c_str());
  if(correction == "RawRes") inputfile = Form("../output/AuAu%s/%s/%s/RawRhoEtaSys_%s_Poly%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),fileoption.c_str(),etamode.c_str(),spectra.c_str());
  if(correction == "RawRes" && order == 1) inputfile = Form("../output/AuAu%s/%s/%s/RawRhoEtaSys_%s_Poly_FirstOrder%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),fileoption.c_str(),etamode.c_str(),spectra.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/AccEffY/%sRhoEtaSys_%s%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),fileoption.c_str());
  if(order == 1 && random3D) inputfile = Form("../output/AuAu%s/%s/3DRandom/%sRhoEtaSys_%s_Poly_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  TFile *File_InPutHframe = TFile::Open(inputfileHframe.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TGraMap g_mRho;
  TH1DMap h_mRho; 
  TH1DMap h_mCounts;
 
  for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
  {           
    //for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; i_cent++) // Centrality loop
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2060; i_cent++) // Centrality loop
    {
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
                for(int i_poly = 0; i_poly < 3; i_poly++)
                {
                  string KEY_rho = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  g_mRho[KEY_rho] = (TGraphAsymmErrors*)File_InPut->Get(KEY_rho.c_str())->Clone();
                  
                  if(correction == "RawRes")
                  {
                    for(int i = 0; i < g_mRho[KEY_rho]->GetN(); i++)
                    {      
                      double pt, rho_obs; 
                      g_mRho[KEY_rho]->GetPoint(i,pt,rho_obs);
                      double rho_obs_err = g_mRho[KEY_rho]->GetErrorYhigh(i);          
       
                      double drhodobs = 4./(1.+3.*Res_12[i_cent]);  // calculate d(rho)/d(rho_obs)
                      double drhodR = -12.*(rho_obs - 1./3.)/(1.+3.*Res_12[i_cent])/(1.+3.*Res_12[i_cent]); // calculate d(rho)/d(R)

                      double real_rho = 1./3. + 4./(1.+3.*Res_12[i_cent])*(rho_obs - 1./3.);
                      double real_rho_error = TMath::Sqrt((rho_obs_err*rho_obs_err)*(drhodobs*drhodobs) + (Res_12_err[i_cent]*Res_12_err[i_cent])*(drhodR*drhodR));                 
                      
                      g_mRho[KEY_rho]->SetPoint(i,pt,real_rho);
                      g_mRho[KEY_rho]->SetPointError(i,0.0,0.0,real_rho_error,real_rho_error);
                    }
                  }
                  for(int i_eta = vmsa::eta_start+ypadding; i_eta < vmsa::eta_stop-ypadding; i_eta++) // Rapidity bin loop
                  {
                    string KEY_counts = Form("eta_%d_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_eta,i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                    h_mCounts[KEY_counts] = (TH1D*)File_InPut->Get(KEY_counts.c_str())->Clone();
                  }
                }
              } 
            }
          }
        }
      }
    }
  }
  //TH1F *h_frame = (TH1F*)File_InPutHframe->Get("h_frame");

  
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-1.51,1.51);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(-1.5,1.5);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  if(etamode == "eta1_eta1") h_frame->GetYaxis()->SetRangeUser(0.24,0.55);
  if(etamode == "eta0p4") h_frame->GetYaxis()->SetRangeUser(0.24,0.55);
  if(etamode == "eta0p6") h_frame->GetYaxis()->SetRangeUser(0.24,0.55);
  if(etamode == "eta0p8") h_frame->GetYaxis()->SetRangeUser(0.24,0.55);
  if(etamode == "eta1_eta1p5") h_frame->GetYaxis()->SetRangeUser(0.0,0.6);
  if(etamode == "eta1p5_eta1p5") h_frame->GetYaxis()->SetRangeUser(0.0,0.6);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();

#if _PlotQA_
  /*if(correction == "AccRes")
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

    h_frame_cos->GetYaxis()->SetRangeUser(0.00256,0.00274);
    h_frame_cos->GetYaxis()->SetNdivisions(505,'N');
    h_frame_cos->GetYaxis()->SetTitle("#frac{dN^{corr}_{#phi}}{d(cos#theta*)} (arb. units)");
    h_frame_cos->GetYaxis()->SetTitleSize(0.05);
    h_frame_cos->GetYaxis()->SetMaxDigits(2);
    h_frame_cos->GetYaxis()->SetLabelSize(0.03);
    h_frame_cos->GetYaxis()->CenterTitle();
    h_frame_cos->DrawCopy("pE");

    Draw_TGAE_Point_new_Symbol(0.07,0.0027025,0.0,0.0,0.0,0.0,20,kBlack,1.3); 
    string leg_count = "#phi Data";
    plotTopLegend((char*)leg_count.c_str(),0.1,0.0027,0.04,1,0.0,42,0);
    string leg_pt = "|y| < 1.0, 1.2 < p_{T} < 1.8 GeV/c";
    plotTopLegend((char*)leg_pt.c_str(),0.1,0.002685,0.04,1,0.0,42,0);

    string leg_sp = "STAR Preliminary";
    plotTopLegend((char*)leg_sp.c_str(),0.1,0.0027425,0.04,2,0.0,42,0);

    PlotLine(0.04,0.085,0.0027225,0.0027225,kRed,2,1);
    string leg_line = "N_{0}[(1+ #frac{B'F}{2})+(A'+F)cos^{2}#theta*+(A'F- #frac{B'F}{2})cos^{4}#theta*]";
    plotTopLegend((char*)leg_line.c_str(),0.1,0.00272,0.04,1,0.0,42,0);

    string leg_energy = Form("Au+Au %s %d", vmsa::mBeamEnergyText[energy].c_str(), vmsa::mBeamYear[energy]);
    plotTopLegend((char*)leg_energy.c_str(),0.55,0.00259,0.04,1,0.0,42,0);
    plotTopLegend((char*)"20%-60%",0.65,0.00258,0.04,1,0.0,42,0);

    TF1 *Func_rho = new TF1("Func_rho",FuncAD,0,1,4);
    string key = Form("pt_%d_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",2,9,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
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
    PtCos->Draw("pEX0 same");
    Func_rho->SetParameter(0, PtCos->GetBinContent(5));
    Func_rho->SetParameter(1, 0.3);
    Func_rho->FixParameter(2, 0.0160525); //hardcoded
    Func_rho->FixParameter(3, 0.331673); //hardcoded
    PtCos->Fit(Func_rho, "NMI"); 
    Func_rho->SetLineColor(kRed);
    Func_rho->Draw("same");
    c_rhocorr->SaveAs("figures/BESII_19p6GeV_2060_CorrectedYieldsCos.eps");
    c_rhocorr->SaveAs("figures/BESII_19p6GeV_2060_CorrectedYieldsCos.pdf");
  }*/
 

 
  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
  c_rho->cd();
  c_rho->cd()->SetLeftMargin(0.15);
  c_rho->cd()->SetBottomMargin(0.15);
  c_rho->cd()->SetTicks(1,1);
  c_rho->cd()->SetGrid(0,0);
  //h_frame->GetYaxis()->SetRangeUser(0.26,0.4);
  h_frame->DrawCopy("pE");
  PlotLine(0.0,80.0,1.0/3.0,1.0/3.0,1,2,2);

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
              for(int i_poly = 0; i_poly < 3; i_poly++)
              {
             // for(int i_F = 0; i_F < 2; i_F++)
             // {
             //   if(correction != "AccRes")
             //   {
             //     if(i_F > 0) continue;
                  string KEY_rho = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ptQA,centQA,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,0,1.1);
             //   }
             //   else
             //   {
             //     string KEY_rho = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",i_cent,i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F);
             //     Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRho[KEY_rho],24,i_sigma+10*i_method+1,1.1);
             //   }
             // }
              }
            }
          }
        }
      }
    }
#endif
  
  TGraMap g_mSysErrors; 
  TGraMap g_mStatErrors; 
  //TGraphAsymmErrors *g_SysErrors = new TGraphAsymmErrors();
  //TGraphAsymmErrors *g_StatErrors = new TGraphAsymmErrors();

  for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
  {           
    //for(int i_cent = vmsa::centStart; i_cent < vmsa::centStop; i_cent++) // Centrality loop
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2060; i_cent++) // Centrality loop
    {
      string KEY_Default = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
      g_mSysErrors[KEY_Default] = new TGraphAsymmErrors();
      g_mStatErrors[KEY_Default] = new TGraphAsymmErrors();
      //string KEY_Default = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
      //cout << "DEFAULT: " << KEY_Default << endl;
      //double sysErr[9][20]; // N points with 5 sources of systematics for each

      for(Int_t i_point = vmsa::eta_start+ypadding; i_point < vmsa::eta_stop-ypadding; ++i_point)
      {
        double sysDca[9];
        double sysNSig[9];
        double sysNorm[9];
        double sysPoly[5];
        //double sysSig[vmsa::Sig_stop];
        //double sysMeth[vmsa::Method_stop];
        double sysF[5];
        double sysEff[5];

        double pt_def, rho_def;
        g_mRho[KEY_Default]->GetPoint(i_point-ypadding,pt_def,rho_def); 

        //cout << "pt: " << pt_def << "   rho00: " << rho_def << endl;
        
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
              string KEY = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",i_pt,i_cent,EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
      //        if(correction == "AccRes") KEY = Form("rhoRaw_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_dca,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str(),defaultF);
              double pt_sys, rho_sys;
              g_mRho[KEY]->GetPoint(i_point-ypadding,pt_sys,rho_sys);
              sysDca[idx] = rho_sys;
              //if(random3D) sysDca[idx] = 4.*(sysDca[idx]-1./3.)+1./3.;
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
              string KEY = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",i_pt,i_cent,EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
        //      if(correction == "AccRes") KEY = Form("rhoRaw_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
              double pt_sys, rho_sys;
              g_mRho[KEY]->GetPoint(i_point-ypadding,pt_sys,rho_sys);
              sysNSig[idx] = rho_sys;
              //if(random3D) sysNSig[idx] = 4.*(sysNSig[idx]-1./3.)+1./3.;
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
              string KEY = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
              //if(correction == "AccRes") KEY = Form("rhoRaw_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
              double pt_sys, rho_sys;
              g_mRho[KEY]->GetPoint(i_point-ypadding,pt_sys,rho_sys);
              sysNorm[idx] = rho_sys;
              //if(random3D) sysNorm[idx] = 4.*(sysNorm[idx]-1./3.)+1./3.;
              idx++;
            }
          }
        }	 
        idx = 0;

        for(int i_poly = 0; i_poly < 2; ++i_poly)
        {
          for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
          {
            for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
            {
              if(i_poly == 0 && (i_sigma != 0 || i_method == 0)) continue;
              if(i_poly != 0 && i_sigma != 0 && i_method == 1) continue;
              string KEY = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
              //if(correction == "AccRes") KEY = Form("rhoRaw_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
              double pt_sys, rho_sys;
              g_mRho[KEY]->GetPoint(i_point-ypadding,pt_sys,rho_sys);
              sysPoly[idx] = rho_sys;
              //if(random3D) sysNorm[idx] = 4.*(sysNorm[idx]-1./3.)+1./3.;
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

       /* if(correction == "AccRes")
        {
          for(int i_F = 0; i_F < 2; ++i_F)
          {
            for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
            {
              for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
              {
                if(i_F == 0 && (i_sigma != 0 || i_method == 0)) continue;
                if(i_F != 0 && i_sigma != 0 && i_method == 1) continue;
                string KEY = Form("rhoRaw_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F);
                double pt_sys, rho_sys;
                g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
                sysF[idx] = rho_sys;
                if(random3D) sysF[idx] = 4.*(sysF[idx]-1./3.)+1./3.;
                idx++;
              }
            }
          }
          idx = 0;

          for(int i_eff = 0; i_eff < 2; ++i_eff)
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
          }
          idx = 0;
        }*/

        Double_t rho_min[4] = { TMath::MinElement(9,sysDca),
                                TMath::MinElement(9,sysNSig),
                                TMath::MinElement(9,sysNorm),
                                TMath::MinElement(5,sysPoly)
                                //(correction == "AccRes")? TMath::MinElement(5,sysF) : 0.0,
                                //(correction == "AccRes")? TMath::MinElement(5,sysEff) : 0.0
                              };

        Double_t rho_max[4] = { TMath::MaxElement(9,sysDca),
                                TMath::MaxElement(9,sysNSig),
                                TMath::MaxElement(9,sysNorm),
                                TMath::MaxElement(5,sysPoly)
                                //(correction == "AccRes")? TMath::MaxElement(5,sysF) : 0.0,
                                //(correction == "AccRes")? TMath::MaxElement(5,sysEff) : 0.0
                              };
      
        double SysError_rho = 0.0;
        for(int i = 0; i < 4; i++)
        {
          double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
          //cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
          SysError_rho += sourcei;
        }
 
        SysError_rho = TMath::Sqrt(SysError_rho);

        Double_t pt, rho;
        g_mRho[KEY_Default]->GetPoint(i_point-ypadding,pt,rho);

        //float mean_rho = total_rho/(float)counter;
        g_mSysErrors[KEY_Default]->SetPoint(i_point-ypadding,pt,rho);
        g_mSysErrors[KEY_Default]->SetPointError(i_point-ypadding,0.0,0.0,SysError_rho,SysError_rho);

        double StatError_rho = g_mRho[KEY_Default]->GetErrorYhigh(i_point-ypadding);
        g_mStatErrors[KEY_Default]->SetPoint(i_point-ypadding,pt,rho);
        g_mStatErrors[KEY_Default]->SetPointError(i_point-ypadding,0.0,0.0,StatError_rho,StatError_rho);
      }
    }
  } 
  //cout << "Finished pT bin calculations" << endl;


  TGraMap g_mStatErrorsInt;
  vecDMap vRho;
  vecDMap vRhoErr;
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
            for(int i_poly = 0; i_poly < 3; i_poly++)
            {
              string KEY_rho = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
              g_mStatErrorsInt[KEY_rho] = new TGraphAsymmErrors();

              double weight_y = 0.0;
              double weight_error_stat_y = 0.0;
              double weight_all_y = 0.0;
              double weight_y_statweight = 0;
              double weight_error_stat_y_statweight = 0;
              double weight_all_y_statweight = 0;

              for(int i_eta = vmsa::eta_start+ypadding; i_eta < vmsa::eta_stop-ypadding; i_eta++) // Rapidity bin loop
              {
                double weight_rho00 = 0;
                double weight_error_stat_rho00 = 0;
                double weight_all = 0;
                double weight_rho00_sw = 0;
                double weight_error_stat_rho00_sw = 0;
                double weight_all_sw = 0;

                for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
                {           
                  for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2060; i_cent++) // Centrality loop
                  { 
                    string KEY_Default = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                    string KEY_counts = Form("eta_%d_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_eta,i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                    h_mCounts[KEY_counts];
                    TH1F *PtCos = (TH1F*)h_mCounts[KEY_counts]->Clone();
                    double weight = PtCos->Integral(1,7); 
                     
                    double etaStat, rhoStat, rhoErrStat;
                    g_mRho[KEY_Default]->GetPoint(i_eta-ypadding,etaStat,rhoStat);   
                    if(rhoStat == -999.0) { /*cout << "skip this rapidity" << endl*/; continue;}
                    rhoErrStat = g_mRho[KEY_Default]->GetErrorYhigh(i_eta-ypadding);   

                    weight_rho00 += rhoStat*weight;
                    //cout << "eta = " << i_eta << ", Cent = " << i_cent << ", pT = " << i_pt << ",    rho00 = " << rhoStat << ", weight = " << weight << endl;
                    weight_error_stat_rho00 += rhoErrStat*rhoErrStat*weight*weight;
                    //cout << "eta = " << i_eta << ", Cent = " << i_cent << ", pT = " << i_pt << ",    rho00 stat error = " << rhoErrStat << ", weight = " << weight << endl;
                    weight_all += weight;

                    weight_rho00_sw += rhoStat/rhoErrStat/rhoErrStat;
                    weight_all_sw += 1./rhoErrStat/rhoErrStat;
 
                    weight_y += rhoStat*weight;
                    weight_error_stat_y += rhoErrStat*rhoErrStat*weight*weight;
                    weight_all_y += weight;       
                    weight_y_statweight += rhoStat/rhoErrStat/rhoErrStat;
                    weight_all_y_statweight += 1./rhoErrStat/rhoErrStat; 

                  } 
                }
                weight_rho00 /= weight_all;
                //cout << "eta = " << i_eta << ",    weight_rho00 = " << weight_rho00 << endl;
                weight_error_stat_rho00 = TMath::Sqrt(weight_error_stat_rho00)/weight_all;
  
                weight_rho00_sw /= weight_all_sw;
                //cout << "eta = " << i_eta << ",    weight_rho00 = " << weight_rho00_sw << endl;
                weight_error_stat_rho00_sw = TMath::Sqrt(1./weight_all_sw);
 
                //double eta_mean = (vmsa::eta_bin[i_eta] + vmsa::eta_bin[i_eta+1])/2.0;
                //g_mStatErrorsInt->SetPoint(i_eta-ypadding, eta_mean, weight_rho00);
                //g_mStatErrorsInt->SetPointError(i_eta-ypadding, 0.0, 0.0, weight_error_stat_rho00, weight_error_stat_rho00);
                //g_mSysErrorsInt->SetPoint(i_eta-ypadding, eta_mean, weight_rho00);
                //g_mSysErrorsInt->SetPointError(i_eta-ypadding, 0.0, 0.0, weight_error_sys_rho00, weight_error_sys_rho00);  
                double eta_mean = (vmsa::eta_bin[i_eta] + vmsa::eta_bin[i_eta+1])/2.0;
                g_mStatErrorsInt[KEY_rho]->SetPoint(i_eta-ypadding, eta_mean, weight_rho00_sw);
                g_mStatErrorsInt[KEY_rho]->SetPointError(i_eta-ypadding, 0.0, 0.0, weight_error_stat_rho00_sw, weight_error_stat_rho00_sw);
              }

             // cout << "KEY_rho: " << KEY_rho << endl;
              weight_y /= weight_all_y;
              weight_error_stat_y = TMath::Sqrt(weight_error_stat_y)/weight_all_y;

              weight_y_statweight /= weight_all_y_statweight;
              weight_error_stat_y_statweight = TMath::Sqrt(1./weight_all_y_statweight);
              //cout << std::setprecision(4);
              //cout << "Total Integrated rho00 yield weight = " << weight_y << " +/- " << weight_error_stat_y << " stat. " << endl;//" +/- " << weight_error_sys_y << "sys." << endl; 
              //cout << "Total Integrated rho00 stat  weight = " << weight_y_statweight << " +/- " << weight_error_stat_y_statweight << " stat. " << endl;//" +/- " << weight_error_sys_y_statweight << "sys." << endl; 
              vRho[KEY_rho].push_back(weight_y_statweight);
              vRhoErr[KEY_rho].push_back(weight_error_stat_y_statweight);
            }
          }
        }
      }
    }
  }

  vecDMap vFoldedRho;
  vecDMap vFoldedRhoErr;
  vecDMap vFoldedDiff;
  TGraMap g_mStatErrorsIntFolded;
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
            for(int i_poly = 0; i_poly < 3; i_poly++)
            {
              string KEY_rho = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
              g_mStatErrorsIntFolded[KEY_rho] = new TGraphAsymmErrors();

              double weight_y_statweight = 0;
              double weight_error_stat_y_statweight = 0;
              double weight_all_y_statweight = 0;
              double negweight_y_statweight = 0;
              double negweight_error_stat_y_statweight = 0;
              double negweight_all_y_statweight = 0;
              double posweight_y_statweight = 0;
              double posweight_error_stat_y_statweight = 0;
              double posweight_all_y_statweight = 0;
              int offset = 5;
              if(random3D && correction == "Raw") offset = 7;
              for(int i_eta = vmsa::eta_start+ypadding; i_eta < vmsa::eta_stop-ypadding-offset; i_eta++) // Rapidity bin loop, the -5 is because of folding
              {
                int yneg = i_eta;
                int ypos = vmsa::eta_stop-1-i_eta;
                //cout << "rapidity indeces: yneg = " << yneg << ", ypos = " << ypos << endl;

                double weight_rho00_sw = 0;
                double weight_error_stat_rho00_sw = 0;
                double weight_all_sw = 0;

                double rapidity, rhoyneg, rhoypos;
                g_mStatErrorsInt[KEY_rho]->GetPoint(yneg-ypadding, rapidity, rhoyneg);
                g_mStatErrorsInt[KEY_rho]->GetPoint(ypos-ypadding, rapidity, rhoypos);
                double rhoynegErr = g_mStatErrorsInt[KEY_rho]->GetErrorYhigh(yneg-ypadding);
                double rhoyposErr = g_mStatErrorsInt[KEY_rho]->GetErrorYhigh(ypos-ypadding);

                //if(rhoStat == -999.0) { cout << "skip this rapidity" << endl; continue;}

                weight_rho00_sw += rhoyneg/rhoynegErr/rhoynegErr;
                weight_all_sw += 1./rhoynegErr/rhoynegErr;
                weight_rho00_sw += rhoypos/rhoyposErr/rhoyposErr;
                weight_all_sw += 1./rhoyposErr/rhoyposErr;
 
                weight_y_statweight += rhoyneg/rhoynegErr/rhoynegErr;
                weight_all_y_statweight += 1./rhoynegErr/rhoynegErr; 
                weight_y_statweight += rhoypos/rhoyposErr/rhoyposErr;
                weight_all_y_statweight += 1./rhoyposErr/rhoyposErr; 

                negweight_y_statweight += rhoyneg/rhoynegErr/rhoynegErr;
                negweight_all_y_statweight += 1./rhoynegErr/rhoynegErr; 

                posweight_y_statweight += rhoypos/rhoyposErr/rhoyposErr;
                posweight_all_y_statweight += 1./rhoyposErr/rhoyposErr; 

                weight_rho00_sw /= weight_all_sw;
                //cout << "eta = " << i_eta << ",    weight_rho00 = " << weight_rho00_sw << endl;
                weight_error_stat_rho00_sw = TMath::Sqrt(1./weight_all_sw);
 
                double eta_mean = fabs(vmsa::eta_bin[i_eta] + vmsa::eta_bin[i_eta+1])/2.0;
                g_mStatErrorsIntFolded[KEY_rho]->SetPoint(yneg-ypadding, eta_mean, weight_rho00_sw);
                g_mStatErrorsIntFolded[KEY_rho]->SetPointError(yneg-ypadding, 0.0, 0.0, weight_error_stat_rho00_sw, weight_error_stat_rho00_sw);
                vFoldedDiff[KEY_rho].push_back(fabs(rhoyneg-rhoypos));
              }

              //cout << "KEY_rho: " << KEY_rho << endl;

              weight_y_statweight /= weight_all_y_statweight;
              weight_error_stat_y_statweight = TMath::Sqrt(1./weight_all_y_statweight);
              negweight_y_statweight /= negweight_all_y_statweight;
              negweight_error_stat_y_statweight = TMath::Sqrt(1./negweight_all_y_statweight);
              posweight_y_statweight /= posweight_all_y_statweight;
              posweight_error_stat_y_statweight = TMath::Sqrt(1./posweight_all_y_statweight);
              //cout << std::setprecision(4);
              //cout << "Total Integrated rho00 stat  weight = " << weight_y_statweight << " +/- " << weight_error_stat_y_statweight << " stat. " << endl;//" +/- " << weight_error_sys_y_statweight << "sys." << endl; 
              vFoldedRho[KEY_rho].push_back(negweight_y_statweight);
              vFoldedRhoErr[KEY_rho].push_back(negweight_error_stat_y_statweight);
              vFoldedRho[KEY_rho].push_back(posweight_y_statweight);
              vFoldedRhoErr[KEY_rho].push_back(posweight_error_stat_y_statweight);
              vFoldedRho[KEY_rho].push_back(weight_y_statweight);
              vFoldedRhoErr[KEY_rho].push_back(weight_error_stat_y_statweight);
              vFoldedDiff[KEY_rho].push_back(fabs(negweight_y_statweight-posweight_y_statweight));
            }
          }
        }
      }
    }
  }

  //cout << "Total Integrated Value" << endl;
  {
    string KEY_Default = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
    double sysDca[9];
    double sysNSig[9];
    double sysNorm[9];
    double sysPoly[5];

    double pt_def, rho_def;
    rho_def = vFoldedRho[KEY_Default][2];
    //cout << "rho00: " << rho_def << endl;
    
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
          string KEY = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
          sysDca[idx] = vFoldedRho[KEY][2];
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
          string KEY = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
          sysNSig[idx] = vFoldedRho[KEY][2];
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
          string KEY = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
          sysNorm[idx] = vFoldedRho[KEY][2];
          idx++;
        }
      }
    }	 
    idx = 0;

    for(int i_poly = 0; i_poly < 2; ++i_poly)
    {
      for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
      {
        for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
        {
          if(i_poly == 0 && (i_sigma != 0 || i_method == 0)) continue;
          if(i_poly != 0 && i_sigma != 0 && i_method == 1) continue;
          string KEY = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          sysPoly[idx] = vFoldedRho[KEY][2];
          idx++;
        }
      }
    }	 
    idx = 0;

    Double_t rho_min[4] = { TMath::MinElement(9,sysDca),
                            TMath::MinElement(9,sysNSig),
                            TMath::MinElement(9,sysNorm),
                            TMath::MinElement(5,sysPoly)
                          };

    Double_t rho_max[4] = { TMath::MaxElement(9,sysDca),
                            TMath::MaxElement(9,sysNSig),
                            TMath::MaxElement(9,sysPoly),
                            TMath::MaxElement(5,sysNorm)
                          };
  
    double SysError_rho = 0.0;
    for(int i = 0; i < 4; i++)
    {
      double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
      //cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
      SysError_rho += sourcei;
    }
    SysError_rho += TMath::Power(vFoldedDiff[KEY_Default][5]/TMath::Sqrt(12.0),2);
 
    SysError_rho = TMath::Sqrt(SysError_rho);
    
    //cout << "FINAL ACTUAL RHO = " << vFoldedRho[KEY_Default][2] << " +/- " << vFoldedRhoErr[KEY_Default][2] << " stat. +/- " << SysError_rho << endl;

  }
  


  TGraphAsymmErrors *g_mSysErrorsFinal = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_mStatErrorsFinal = new TGraphAsymmErrors();

  string KEY_Default = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
  //string KEY_Default = Form("rhoRaw_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
  //cout << "DEFAULT: " << KEY_Default << endl;
  //double sysErr[9][20]; // N points with 5 sources of systematics for each
  int offset = 5; 
  if(random3D && correction == "Raw") offset = 7;
  for(Int_t i_point = vmsa::eta_start+ypadding; i_point < vmsa::eta_stop-ypadding-offset; ++i_point)
  {
    double sysDca[9];
    double sysNSig[9];
    double sysNorm[9];
    double sysPoly[5];
    //double sysSig[vmsa::Sig_stop];
    //double sysMeth[vmsa::Method_stop];
    double sysF[5];
    double sysEff[5];

    double pt_def, rho_def;
    g_mStatErrorsIntFolded[KEY_Default]->GetPoint(i_point-ypadding,pt_def,rho_def); 

    //cout << "pt: " << pt_def << "   rho00: " << rho_def << endl;
    
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
          string KEY = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
  //        if(correction == "AccRes") KEY = Form("rhoRaw_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_dca,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str(),defaultF);
          double pt_sys, rho_sys;
          g_mStatErrorsIntFolded[KEY]->GetPoint(i_point-ypadding,pt_sys,rho_sys);
          sysDca[idx] = rho_sys;
          //if(random3D) sysDca[idx] = 4.*(sysDca[idx]-1./3.)+1./3.;
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
          string KEY = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
    //      if(correction == "AccRes") KEY = Form("rhoRaw_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
          double pt_sys, rho_sys;
          g_mStatErrorsIntFolded[KEY]->GetPoint(i_point-ypadding,pt_sys,rho_sys);
          sysNSig[idx] = rho_sys;
          //if(random3D) sysNSig[idx] = 4.*(sysNSig[idx]-1./3.)+1./3.;
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
          string KEY = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
          //if(correction == "AccRes") KEY = Form("rhoRaw_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
          double pt_sys, rho_sys;
          g_mStatErrorsIntFolded[KEY]->GetPoint(i_point-ypadding,pt_sys,rho_sys);
          sysNorm[idx] = rho_sys;
          //if(random3D) sysNorm[idx] = 4.*(sysNorm[idx]-1./3.)+1./3.;
          idx++;
        }
      }
    }	 
    idx = 0;

    for(int i_poly = 0; i_poly < 2; ++i_poly)
    {
      for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
      {
        for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
        {
          if(i_poly == 0 && (i_sigma != 0 || i_method == 0)) continue;
          if(i_poly != 0 && i_sigma != 0 && i_method == 1) continue;
          string KEY = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          //if(correction == "AccRes") KEY = Form("rhoRaw_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
          double pt_sys, rho_sys;
          g_mStatErrorsIntFolded[KEY]->GetPoint(i_point-ypadding,pt_sys,rho_sys);
          sysPoly[idx] = rho_sys;
          //if(random3D) sysNorm[idx] = 4.*(sysNorm[idx]-1./3.)+1./3.;
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

   /* if(correction == "AccRes")
    {
      for(int i_F = 0; i_F < 2; ++i_F)
      {
        for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
        {
          for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
          {
            if(i_F == 0 && (i_sigma != 0 || i_method == 0)) continue;
            if(i_F != 0 && i_sigma != 0 && i_method == 1) continue;
            string KEY = Form("rhoRaw_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F);
            double pt_sys, rho_sys;
            g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
            sysF[idx] = rho_sys;
            if(random3D) sysF[idx] = 4.*(sysF[idx]-1./3.)+1./3.;
            idx++;
          }
        }
      }
      idx = 0;

      for(int i_eff = 0; i_eff < 2; ++i_eff)
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
      }
      idx = 0;
    }*/

    Double_t rho_min[4] = { TMath::MinElement(9,sysDca),
                            TMath::MinElement(9,sysNSig),
                            TMath::MinElement(9,sysNorm),
                            TMath::MinElement(5,sysPoly)
                            //(correction == "AccRes")? TMath::MinElement(5,sysF) : 0.0,
                            //(correction == "AccRes")? TMath::MinElement(5,sysEff) : 0.0
                          };

    Double_t rho_max[4] = { TMath::MaxElement(9,sysDca),
                            TMath::MaxElement(9,sysNSig),
                            TMath::MaxElement(9,sysNorm),
                            TMath::MaxElement(5,sysPoly)
                            //(correction == "AccRes")? TMath::MaxElement(5,sysF) : 0.0,
                            //(correction == "AccRes")? TMath::MaxElement(5,sysEff) : 0.0
                          };
  
    double SysError_rho = 0.0;
    for(int i = 0; i < 4; i++)
    {
      double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
      //cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
      SysError_rho += sourcei;
    }
    SysError_rho += TMath::Power(vFoldedDiff[KEY_Default][i_point-ypadding]/TMath::Sqrt(12.0),2);
 
    SysError_rho = TMath::Sqrt(SysError_rho);

    Double_t pt, rho;
    g_mStatErrorsIntFolded[KEY_Default]->GetPoint(i_point-ypadding,pt,rho);

    //float mean_rho = total_rho/(float)counter;
    g_mSysErrorsFinal->SetPoint(i_point-ypadding,pt,rho);
    g_mSysErrorsFinal->SetPointError(i_point-ypadding,0.0,0.0,SysError_rho,SysError_rho);

    double StatError_rho = g_mStatErrorsIntFolded[KEY_Default]->GetErrorYhigh(i_point-ypadding);
    g_mStatErrorsFinal->SetPoint(i_point-ypadding,pt,rho);
    g_mStatErrorsFinal->SetPointError(i_point-ypadding,0.0,0.0,StatError_rho,StatError_rho);
  }
   
   
  /*if(correction == "AccRes")
  { 
  string KEY_DefaultW = Form("rhoFinalWeighted_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),defaultF);

  double sysDca[9];
  double sysNSig[9];
  double sysNorm[9];
  double sysF[5];
  double sysEff[5];

  double rho_def = h_mRho[KEY_DefaultW]->GetBinContent(7); 
  double rho_defErr = h_mRho[KEY_DefaultW]->GetBinError(7); 

  int idx = 0;
  for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
  { 
    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
      {
        if(i_dca == 0 && (i_sigma != 0 || i_method == 0)) continue;
        if(i_dca != 0 && i_sigma != 0 && i_method == 1) continue;
        string KEY = Form("rhoFinalWeighted_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
        sysDca[idx] = h_mRho[KEY]->GetBinContent(7);
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
        string KEY = Form("rhoFinalWeighted_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
        sysNSig[idx] = h_mRho[KEY]->GetBinContent(7);
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
        string KEY = Form("rhoFinalWeighted_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF);
        sysNorm[idx] = h_mRho[KEY]->GetBinContent(7);
        idx++;
      }
    }
  }	 
  idx = 0;

  for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
  {
    string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",i_cent,0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[1].c_str(),defaultF);
    sysSig[i_sigma] = h_mRho[KEY]->GetBinContent(7);
  }

  for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
  {
    string KEY = Form("rhoFinalWeighted_Centrality_%d_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d",i_cent,0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str(),defaultF);
    sysMeth[i_method] = h_mRho[KEY]->GetBinContent(7);
  }


  for(int i_F = 0; i_F < 2; ++i_F)
  {
    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
      {
        if(i_F == 0 && (i_sigma != 0 || i_method == 0)) continue;
        if(i_F != 0 && i_sigma != 0 && i_method == 1) continue;
        string KEY = Form("rhoFinalWeighted_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_0",0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_F);
        sysF[idx] = h_mRho[KEY]->GetBinContent(7);
        idx++;
      }
    }
  }
  idx = 0;

  for(int i_eff = 0; i_eff < 2; ++i_eff)
  {
    for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
    {
      for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
      {
        if(i_eff == 0 && (i_sigma != 0 || i_method == 0)) continue;
        if(i_eff != 0 && i_sigma != 0 && i_method == 1) continue;
        string KEY = Form("rhoFinalWeighted_2nd_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_F_%d_Eff_%d",0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),defaultF,i_eff);
        sysEff[idx] = h_mRho[KEY]->GetBinContent(7);
        idx++;
      }
    }
  }
  idx = 0;

  Double_t rho_min[5] = { TMath::MinElement(9,sysDca),
                          TMath::MinElement(9,sysNSig),
                          TMath::MinElement(9,sysNorm),
                          TMath::MinElement(5,sysF),    
                          TMath::MinElement(5,sysEff)    
                        };

  Double_t rho_max[5] = { TMath::MaxElement(9,sysDca),
                          TMath::MaxElement(9,sysNSig),
                          TMath::MaxElement(9,sysNorm),
                          TMath::MaxElement(5,sysF),    
                          TMath::MaxElement(5,sysEff)    
                        };
  
  double SysError_rho = 0.0;
  for(int i = 0; i < 5; i++)
  {
    double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
    cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
    SysError_rho += sourcei;
  }
 
  SysError_rho = TMath::Sqrt(SysError_rho);

  cout << "Final corrected rho00 = " << std::fixed << std::setprecision(4) << rho_def << " +/- " << rho_defErr << " (stat) +/- " << SysError_rho << " (sys)" << endl;
  cout << "sigma from 1/3 = " << std::fixed << std::setprecision(2) << TMath::Abs((rho_def-1./3.)/TMath::Sqrt(rho_defErr*rho_defErr + SysError_rho*SysError_rho)) << endl;
  cout << "sigma from BESI = " << TMath::Abs((rho_def-0.37)/TMath::Sqrt(0.008*0.008 + 0.007*0.007 + rho_defErr*rho_defErr + SysError_rho*SysError_rho)) << endl;
  }*/
  /*double val[6];
  double sysE[6];
  double statE[6];

  for(int ipt = 0; ipt < 6; ipt++)
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

  cout << "rho00={"; for (int i = 0; i < 9; i++) {if (i<8) cout << val[i] << ","; else cout << val[i] << "};" << endl;}
  cout << "stat={"; for (int i = 0; i < 9; i++) {if (i<8) cout << statE[i] << ","; else cout << statE[i] << "};" << endl;}
  cout << "sys={"; for (int i = 0; i < 9; i++) {if (i<8) cout << sysE[i] << ","; else cout << sysE[i] << "};" << endl;}
  */
  //TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  //TGraphAsymmErrors *besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_19;1");
  //TGraphAsymmErrors *besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_19;1");

  //cout << "All good" << endl;

  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,800,800);
  c_rho_SysError->cd();
  c_rho_SysError->cd()->SetLeftMargin(0.15);
  c_rho_SysError->cd()->SetBottomMargin(0.15);
  c_rho_SysError->cd()->SetTicks(1,1);
  c_rho_SysError->cd()->SetGrid(0,0);
  //h_frame->GetXaxis()->SetNdivisions(505,"I");
  //h_frame->GetYaxis()->SetRangeUser(0.28,0.5);
  //h_frame->GetXaxis()->SetLimits(0,9.95);
  h_frame->GetXaxis()->SetRangeUser(0.0,1.0);
  if(random3D) h_frame->GetYaxis()->SetRangeUser(0.3,0.37);
  //h_frame->GetXaxis()->SetLimits(0,10);
  h_frame->GetXaxis()->SetTitle("|y|");
  h_frame->DrawCopy("pE");
  //g_StatErrors->SetMarkerStyle(20);
  //g_StatErrors->SetMarkerColor(kGray+2);
  //g_StatErrors->SetLineColor(2);
  //g_StatErrors->Draw("pE same");
  //g_SysErrors->SetMarkerStyle(20);
  //g_SysErrors->SetMarkerColor(kGray+2);
  //g_SysErrors->SetLineColor(kGray+2);
  //g_SysErrors->Draw("pE same");
//  PlotLine(0.0,5.0,1.0/3.0,1.0/3.0,1,2,2);

  //for(int i = 0; i < 4; i++)
  //{
  //  besi19_sys->SetPointEXhigh(i,0.0);
  //  besi19_sys->SetPointEXlow(i,0.0);
  //  besi19->SetPointEXhigh(i,0.0);
  //  besi19->SetPointEXlow(i,0.0);
  //}
  PlotLine(0.0,1.0,1.0/3.0,1.0/3.0,1,2,2);
  //gStyle->SetEndErrorSize(3);
  //gStyle->SetEndErrorWidth(10);
  //besi19->SetLineColor(kBlack);
  //besi19->SetLineWidth(1.0);
  //besi19->SetMarkerStyle(20);
  //besi19->SetMarkerSize(1.3);
  //besi19->SetMarkerColor(kBlack);
  //besi19->SetLineColor(kBlack);
  //besi19->Draw("pE Z same");


  //besi19_sys->SetLineColor(kBlack);
  //besi19_sys->SetMarkerStyle(1);
  //besi19_sys->SetMarkerSize(10);
  //besi19_sys->SetMarkerColor(kBlack);
  //besi19_sys->SetLineColor(kBlack);
  //besi19_sys->Draw("pE [] same");

  Draw_TGAE_Point_new_Symbol(0.1,0.455,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2); 
  string leg_count;
  if(correction == "Eff") leg_count = "Efficiency Corrected #phi BES-II";
  if(correction == "Raw") leg_count = "Raw #phi BES-II";
  if(correction == "AccRes") leg_count = "#phi BES-II (|y| < 1.0)";
  plotTopLegend((char*)leg_count.c_str(),+0.15,0.454,0.03,1,0.0,42,0);

  //string leg_besi = "#phi BES-I  (|y| < 1.0)";
  //Draw_TGAE_Point_new_Symbol(0.5,0.44,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
  //plotTopLegend((char*)leg_besi.c_str(),0.65,0.437,0.03,1,0.0,42,0);
 
  //string leg_sp = "STAR Preliminary";
  //plotTopLegend((char*)leg_sp.c_str(),0.3,0.48,0.03,2,0.0,42,0);

  PlotLine(0.05,0.125,0.435,0.435,1,2,2);
  string leg_line = "#rho_{00} = 1/3";
  plotTopLegend((char*)leg_line.c_str(),0.15,0.434,0.03,1,0.0,42,0);

  string leg_energy = Form("Au+Au %s", vmsa::mBeamEnergyText[energy].c_str());
  plotTopLegend((char*)leg_energy.c_str(),0.35,0.29,0.04,1,0.0,42,0);
  plotTopLegend((char*)"0%-80%",0.45,27,0.04,1,0.0,42,0);

  g_mStatErrorsFinal->SetMarkerStyle(20);
  g_mStatErrorsFinal->SetMarkerColor(kRed);
  g_mStatErrorsFinal->SetLineColor(kRed);
  g_mStatErrorsFinal->SetMarkerSize(1.3);

  g_mSysErrorsFinal->SetMarkerStyle(20);
  g_mSysErrorsFinal->SetMarkerSize(1.3);
  g_mSysErrorsFinal->SetMarkerColor(kRed);
  g_mSysErrorsFinal->SetLineColor(kRed);
  //g_SysErrors->Draw("pE [] same");
  //if(correction == "AccRes") 
  //{
  //g_StatErrors->SetPoint(0,-1.0,-1.0);
  //g_SysErrors->SetPoint(0,-1.0,-1.0);
  //g_StatErrors->SetPoint(1,-1.0,-1.0);
  //g_SysErrors->SetPoint(1,-1.0,-1.0);
  //g_SysErrors->Draw("pE [] same");
  //}

  //g_StatErrors->SetLineWidth(2);
  //besi19->SetLineWidth(2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mStatErrorsFinal,style_phi_2nd,color_phi_2nd,colorDiff_phi,size_marker+0.2);
  //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)besi19,style_phi_1st,color_phi_1st,colorDiff_phi,size_marker-0.2);
  // plotSysErrors(g_rho_2nd_sys[total_pad-1],color_phi_2nd);
  plotSysErrorsBox(g_mSysErrorsFinal,color_phi_2nd);
  //plotSysErrorsBox(besi19_sys,color_phi_1st);
  //g_StatErrors->Draw("pE same");
  //plotSysErrorsBox(g_SysErrors,kRed,energy);
  c_rho_SysError->SaveAs(Form("./figures/%s/%s/rapiditystudy/TPCOnly_finaloutputrho_rapiditydependence_%s_%s_%s_order%d_%s_PolySys.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,etamode.c_str()));  

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
  }
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
  }*/

  
  string outputname = Form("./figures/%s/%s/rapiditystudy/TPCOnly_finaloutputrho_rapiditydependence_%s_ptandcentbins_%s_%s_order%d_%s_PolySys.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),correction.c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,etamode.c_str());
  string output_start = Form("%s[",outputname.c_str());

  TCanvas *c_pt = new TCanvas("c_pt","c_pt",10,10,1200,800);
  c_pt->Divide(3,2);

  c_pt->Print(output_start.c_str());

  std::string centStrings[3] = {"40-80","10-40","0-10"};

  h_frame->GetXaxis()->SetRangeUser(-1.5,1.5);

  for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
  {
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2060; i_cent++) // Centrality loop
    {
      c_pt->cd(i_pt*3+i_cent+1);
      c_pt->cd(i_pt*3+i_cent+1)->SetLeftMargin(0.15);
      c_pt->cd(i_pt*3+i_cent+1)->SetBottomMargin(0.15);
      c_pt->cd(i_pt*3+i_cent+1)->SetTicks(1,1);
      c_pt->cd(i_pt*3+i_cent+1)->SetGrid(0,0);
      h_frame->SetTitle(Form("%.2f<p_{T}<%.2f GeV/c, Cent (%s)",vmsa::pt_low_y[energy][i_pt],vmsa::pt_up_y[energy][i_pt],centStrings[i_cent].c_str()));
      h_frame->DrawCopy("pE");
      PlotLine(-1.5,1.5,1.0/3.0,1.0/3.0,1,2,2);
      Draw_TGAE_Point_new_Symbol(0.1,0.46,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2); 

      plotTopLegend((char*)leg_count.c_str(),0.15,0.457,0.04,1,0.0,42,0);

      //PlotLine(-0.25,-0.05,0.438,0.438,1,2,2);
      //plotTopLegend((char*)leg_line.c_str(),0.0,0.437,0.03,1,0.0,42,0);

      plotTopLegend((char*)leg_energy.c_str(),0.5,0.305,0.04,1,0.0,42,0);
      //plotTopLegend((char*)centStrings[i_cent].c_str(),-0.1,0.295,0.04,1,0.0,42,0);
      //plotTopLegend((char*)Form("%.2f<p_{T}<%.2f GeV/c",vmsa::pt_low_y[energy][i_pt],vmsa::pt_up_y[energy][i_pt]),-0.3,0.285,0.04,1,0.0,42,0);

      string KEY_Default = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mStatErrors[KEY_Default],style_phi_2nd,color_phi_2nd,colorDiff_phi,size_marker+0.2);

      plotSysErrorsBox(g_mSysErrors[KEY_Default],color_phi_2nd);     

    }
  }
  c_pt->Update();
  c_pt->Print(outputname.c_str());
  string output_stop = Form("%s]",outputname.c_str());
  c_pt->Print(output_stop.c_str()); // close pdf file




  string OutPutFile = Form("../output/AuAu%s/%s/TPCOnly_RhoEta_%sSysErrors_%s_PolySys%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),spectra.c_str());
  if(order == 1) OutPutFile = Form("../output/AuAu%s/%s/RhoEta_%sSysErrors_%s_PolySys_FirstOrder%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),spectra.c_str());
  if(order == 1 && random3D) OutPutFile = Form("../output/AuAu%s/%s/3DRandom/RhoEta_%sSysErrors_%s_PolySys_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  cout << "OutPutFile set to: " << OutPutFile.c_str() << endl;
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  TString StatErrorRho = Form("g_rho00_order%d_%s_%s_StatError",order,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  g_mStatErrorsFinal->SetName(StatErrorRho.Data());
  g_mStatErrorsFinal->Write();
  TString SysErrorRho = Form("g_rho00_order%d_%s_%s_SysError",order,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  g_mSysErrorsFinal->SetName(SysErrorRho.Data());
  g_mSysErrorsFinal->Write();

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
            for(int i_poly = 0; i_poly < 3; i_poly++)
            {
              string KEY_rho = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
              string KEY_rho_unfolded = Form("rhoRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_Unfolded",EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
              g_mStatErrorsIntFolded[KEY_rho]->SetName(KEY_rho.c_str());
              g_mStatErrorsIntFolded[KEY_rho]->Write();
              g_mStatErrorsInt[KEY_rho]->SetName(KEY_rho_unfolded.c_str()); 
              g_mStatErrorsInt[KEY_rho]->Write();
            }
          }
        }
      }
    }
  }
  for(int i_pt = vmsa::pt_rebin_first_y[energy]; i_pt <= vmsa::pt_rebin_last_y[energy]; ++i_pt) // pt loop
  {           
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2060; i_cent++) // Centrality loop
    {
      string KEY_rho = Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
      string KEY_rhoStat = Form("rhoRawStat_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
      string KEY_rhoSys  = Form("rhoRawSys_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
      g_mStatErrors[KEY_rho]->SetName(KEY_rhoStat.c_str());    
      g_mStatErrors[KEY_rho]->Write();    
      g_mSysErrors[KEY_rho]->SetName(KEY_rhoSys.c_str());    
      g_mSysErrors[KEY_rho]->Write();    
    }
  }

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

    bSys[i_pt] = new TBox(pt-0.02,rho-err,pt+0.02,rho+err);
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

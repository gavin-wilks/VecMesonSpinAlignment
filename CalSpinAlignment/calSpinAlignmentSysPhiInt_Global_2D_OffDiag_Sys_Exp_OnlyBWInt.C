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


void calSpinAlignmentSysPhiInt_Global_2D_OffDiag_Sys_Exp_OnlyBWInt(int energy = 4, int pid = 0, int year = 0, bool random3D = false, int order = 2, string etamode = "eta1_eta1", int frameopt = 0, int sysidx = 0, string flag = "")
{
  // energy is the index for which data set we are working with: 4 = 19.6 GeV
  // pid is the particle index that we set: 0 = phi-meson
  // year is obsolete 
  // random3D is a flag to switch which set of histograms is read in: false for main analysis
  // order = 1, 2 for which harmonic order of event plane we are interested in
  // etamode is an extra string that was used when exploring different eta cuts, default is "eta1_eta1"
  // frameopt = 0 corresponds to global frame, frameopt = 1 corresponds to helicity frame
  // sysidx = 0,1,2,...,9 and sets which file is loaded in for different systematics
  //       0: default cuts
  //       1: Dca < 1.0 cm
  //       2: Dca < 3.0 cm
  //       3: nsigma_kaon < 1.5
  //       4: nsigma_kaon < 2.0
  //       5: NhitsFit > 16
  //       6: NhitsFit > 18
  //       7: NHitsRatio >= 0.54
  //       8: NHitsRatio >= 0.56
  //       9: 0.18 < m^2 < 0.34 GeV/c^2  
  // flag sets an extra string that is appended to the end of the pdfs and files to distinguish small setting changes

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(2000000);

  // default case for sysidx = 1,...,9
  // says to only loop over the default initial value for normalization, signal extraction range
  // yield integration range and residual background function
  int stop_norm = 1;
  int sig_stop = 1;
  int start_method = 1;
  int stop_poly = 1;

  // here we set the index ranges for the 4 systematic variables that are used in the fitting and yield extraction process
  // these are NOT data level cuts like sysidx = 1,...,9 
  if(sysidx == 0) 
  {
    stop_norm = vmsa::Norm_stop;        // normalization 
    sig_stop = vmsa::Sig_stop;          // integration range
    start_method = vmsa::Method_start;  // yield extraction method, note that the default is 1, so here the starting index is changed to 0
    stop_poly = 2;                      // residual background function, this starts at 1
  }

  std::string frame = "Global";
  if(frameopt == 1) frame = "Helicity";
  
  std::string EP[2] = {"1st","2nd"};
  string inputfile = Form("../output/AuAu%s/%s/%s2DInt_InvMassSubBg_%s_SysIdx%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),frame.c_str(),etamode.c_str(),sysidx);
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/%s2DInt_InvMassSubBg_%s_FirstOrder_SysIdx%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),frame.c_str(),etamode.c_str(),sysidx);
  if(random3D) inputfile = Form("../output/AuAu%s/%s/3DRandom/InvMassSubBg.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  File_InPut->cd();
  TH1FMap h_mMass, h_mMass_InteTheta;
  vecFMap Par_InteTheta;
  double Centers;
  double Widths;
  double WidthsErr;
 
  TGraMap g_mChiNDF;
  TGraMap g_mPValue;
  // Here we read in histograms
  // integrated over cos(theta*) and do breit wiger fit to extract common fit parameter
  // the common fit parameters of mass and width will be fixed in future inidividual bin fits
  for(int i_norm = vmsa::Norm_start; i_norm < stop_norm; ++i_norm)
  {
    string KEY_InteTheta = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
    for(int i_theta = vmsa::CTS_start;  i_theta < vmsa::Cos_Theta_Star_Bins; i_theta++) // cos(theta*) loop
    {
      for(int i_phipsi = 0;  i_phipsi < vmsa::Beta_Bins; i_phipsi++) // cos(theta*) loop
      {
        string KEY = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
        h_mMass[KEY] = (TH1F*)File_InPut->Get(KEY.c_str());
    
        if(i_theta == vmsa::CTS_start && i_phipsi == 0) h_mMass_InteTheta[KEY_InteTheta] = (TH1F*)h_mMass[KEY]->Clone(KEY_InteTheta.c_str());
        else h_mMass_InteTheta[KEY_InteTheta]->Add(h_mMass[KEY],1.0);
      }
      for(int i_poly = 0; i_poly < stop_poly; i_poly++)
      {
        //if(i_poly != 0 && i_norm != 0) continue;
        string KEY_InteTheta_Poly = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
        ///string KEY_Poly = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
        //if(i_pt == vmsa::pt_rebin_first_2D[energy]) g_mChiNDF[KEY_Poly] = new TGraphAsymmErrors();             
        //if(i_pt == vmsa::pt_rebin_first_2D[energy]) g_mPValue[KEY_Poly] = new TGraphAsymmErrors();             
    
        TF1 *f_bw; 
        if(i_poly == 0) f_bw = new TF1("f_bw",Poly1MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
        if(i_poly == 1) f_bw = new TF1("f_bw",Poly2MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
        if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
        f_bw->SetParameter(0,vmsa::InvMass[pid]);
        //f_bw->SetParLimits(0,vmsa::InvMass[pid]-0.003,vmsa::InvMass[pid]+0.003);
        f_bw->SetParameter(1,vmsa::Width[pid]);
        f_bw->SetParameter(2,1.0);
        float norm = h_mMass_InteTheta[KEY_InteTheta]->GetMaximum()/f_bw->GetMaximum();
        f_bw->SetParameter(2,norm);
        f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
        h_mMass_InteTheta[KEY_InteTheta]->Fit(f_bw,"MNRI");
        Par_InteTheta[KEY_InteTheta_Poly].clear();
        Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(0)));
        Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(1)));
        Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(2)));
        Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(3)));
        Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(4)));
        if(i_poly >= 1) Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(5)));
        if(i_poly >= 2) Par_InteTheta[KEY_InteTheta_Poly].push_back(static_cast<float>(f_bw->GetParameter(6)));
          
        float pt_mean = 0.5;//(vmsa::pt_low_2D[energy][i_pt]+vmsa::pt_up_2D[energy][i_pt])/2.0;
        float chi2NDF = float(f_bw->GetChisquare())/float(f_bw->GetNDF()); 
        float pvalue = TMath::Prob(f_bw->GetChisquare(),f_bw->GetNDF());
        cout << "Poly " << i_poly + 1 << ", pT = " << pt_mean << ", chi2 = " << f_bw->GetChisquare() << ", NDF = " << f_bw->GetNDF() << ", pvalue = " << pvalue << endl;
        //if(i_pt >= 2) g_mChiNDF[KEY_Poly]->SetPoint(i_pt,pt_mean,chi2NDF);
        //if(i_pt >= 2) g_mChiNDF[KEY_Poly]->SetPoint(i_pt,pt_mean,pvalue);
        if(i_poly == 0 && sysidx == 0 && i_norm == 0)
        { 
          //saving the centers, widths, and widths errors for future use in the analysis pipeline
          Centers = f_bw->GetParameter(0);
          Widths = f_bw->GetParameter(1);
          WidthsErr = f_bw->GetParError(1);
        }
      }
    }
  }

  cout << "Finished first bit" << endl;

  // extract counts vs. pT with diffenretial integration ranges and methods
  vecFMap Par;
  TH2FMap h_mCounts;
  TH1FMap h_mCounts1D;
  vecFMap Par_rhoFit;
  TGraMap g_mRho;
  TGraMap g_mReal;
  TGraMap g_mImag;
  TGraMap g_mReRho1n1;
  TGraMap g_mImRho1n1;

  TH2FMap g_reducedchi2;
  TH2FMap h_fitstatus;
  TGraMap g_pull;


  // This is the main fitting loop
  for(int i_norm = vmsa::Norm_start; i_norm < stop_norm; ++i_norm)
  {
    for(int i_sigma = vmsa::Sig_start; i_sigma < sig_stop; ++i_sigma)
    {
      for(int i_method = start_method; i_method < vmsa::Method_stop; ++i_method)
      {
        for(int i_poly = 0; i_poly < stop_poly; i_poly++)
        {
          // initialize the TGraphs for all 5 SDMEs
          string KEY_rho = Form("rhoRaw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          string KEY_real = Form("realRaw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          string KEY_imag = Form("imagRaw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          string KEY_rerho1n1 = Form("rerho1n1Raw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          string KEY_imrho1n1 = Form("imrho1n1Raw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          g_mRho[KEY_rho] = new TGraphAsymmErrors();
          g_mReal[KEY_real] = new TGraphAsymmErrors();
          g_mImag[KEY_imag] = new TGraphAsymmErrors();
          g_mReRho1n1[KEY_rerho1n1] = new TGraphAsymmErrors();
          g_mImRho1n1[KEY_imrho1n1] = new TGraphAsymmErrors();

          // initialize different counting histograms
          string KEY_counts = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          string KEY_countsbg = Form("bg_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          string KEY_countssigbg = Form("sigbg_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          // This is the main signal histogram
          h_mCounts[KEY_counts] = new TH2F(KEY_counts.c_str(),KEY_counts.c_str(),vmsa::Cos_Theta_Star_Bins,-1.0,1.0,vmsa::Beta_Bins,0.0,2.0*TMath::Pi());
          // This counts just the background
          h_mCounts[KEY_countsbg] = new TH2F(KEY_countsbg.c_str(),KEY_countsbg.c_str(),vmsa::Cos_Theta_Star_Bins,-1.0,1.0,vmsa::Beta_Bins,0.0,2.0*TMath::Pi());
          // This counts the signal + background
          h_mCounts[KEY_countssigbg] = new TH2F(KEY_countssigbg.c_str(),KEY_countssigbg.c_str(),vmsa::Cos_Theta_Star_Bins,-1.0,1.0,vmsa::Beta_Bins,0.0,2.0*TMath::Pi());
  
          string KEY_InteTheta = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
          string KEY_InteTheta_Poly = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);

          // tracking other variables of the fit
          g_mChiNDF[KEY_InteTheta_Poly] = new TGraphAsymmErrors();             
          g_mPValue[KEY_InteTheta_Poly] = new TGraphAsymmErrors();             

          string KEY_reducedchi2 = Form("chi2_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          g_reducedchi2[KEY_reducedchi2] = new TH2F(KEY_reducedchi2.c_str(),KEY_reducedchi2.c_str(),vmsa::Cos_Theta_Star_Bins,-1.0,1.0,vmsa::Beta_Bins,0.0,2.0*TMath::Pi());                
 
          string KEY_fitstatus = Form("status_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          h_fitstatus[KEY_fitstatus] = new TH2F(KEY_fitstatus.c_str(),KEY_fitstatus.c_str(),vmsa::Cos_Theta_Star_Bins,-1.0,1.0,vmsa::Beta_Bins,0.0,2.0*TMath::Pi());                
 
          // loop over 2D cos(theta*) and beta bins
          for(int i_theta = vmsa::CTS_start;  i_theta < vmsa::Cos_Theta_Star_Bins; ++i_theta) // cos(theta*) loop
          {
            for(int i_phipsi = 0;  i_phipsi < vmsa::Beta_Bins; ++i_phipsi) // beta loop
            {
              string KEY = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
              string KEY_Poly = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
              string KEY_Poly_diff = Form("diff_CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
              string KEY_Poly_diffcum = Form("diffcum_CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);

              if(i_method == 1)
              {
                g_pull[KEY_Poly] = new TGraphAsymmErrors();         // this graph tracks the pull (data-fit)/(data error) vs M
                g_pull[KEY_Poly_diff] = new TGraphAsymmErrors();    // this graph tracks the data-fit vs M
                g_pull[KEY_Poly_diffcum] = new TGraphAsymmErrors(); // this graph tracks the cumulative (data-fit) across M vs M
              }

              // define full BW+background function
              TF1 *f_bw;
              if(i_poly == 0) f_bw = new TF1("f_bw",Poly1MBreitWigner, vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
              if(i_poly == 1) f_bw = new TF1("f_bw",Poly2MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
              if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
              //f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_Poly][0]);
              //f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_Poly][1]);
              f_bw->FixParameter(0,Par_InteTheta[KEY_InteTheta_Poly][0]); // fix mass 
              f_bw->FixParameter(1,Par_InteTheta[KEY_InteTheta_Poly][1]); // fix width
              f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_Poly][2]/100.0);
              f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_Poly][3]/100.0);
              f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_Poly][4]/100.0);
              if(i_poly >= 1) f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_Poly][5]);
              if(i_poly >= 2) f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta_Poly][6]);
              //f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta][5]);
              f_bw->SetRange(vmsa::BW_Start[pid],vmsa::BW_Stop[pid]);
              //if(h_mMass[KEY]->GetEntries() <= 1) continue;
              //cout << h_mMass[KEY]->GetEntries() << endl;
              TFitResultPtr result = h_mMass[KEY]->Fit(f_bw,"NMRIS");
            
              cout << "IsValid? = " << result->IsValid() << endl;
              cout << "Status = " << result->Status() << endl;
              cout << "Total calls = " << result->NCalls() << endl;
              cout << "EDM = " << result->Edm() << endl;
              cout << "CovMatrixStatus = " << result->CovMatrixStatus() << endl;
              h_fitstatus[KEY_fitstatus]->SetBinContent(i_theta+1,i_phipsi+1,result->IsValid()); 
 
              result->GetCovarianceMatrix().Print();

              //cout << Form("[0]+[1]*(x-%1.10f)",Par_InteTheta[KEY_InteTheta_Poly][0]) << endl;
              //cout << Form("[0]+[1]*(x-%1.10f)+[2]*pow(x-%1.10f,2)",Par_InteTheta[KEY_InteTheta_Poly][0],Par_InteTheta[KEY_InteTheta_Poly][0]) << endl;

              // BW without residual background
              TF1 *f_bw_only = new TF1("f_bw_only",BreitWigner, vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
              //f_bw_only->SetParameter(0,f_bw->GetParameter(0));
              //f_bw_only->SetParameter(1,f_bw->GetParameter(1));
              f_bw_only->FixParameter(0,f_bw->GetParameter(0));
              f_bw_only->FixParameter(1,f_bw->GetParameter(1));
              f_bw_only->SetParameter(2,f_bw->GetParameter(2));

            
              // propagate variables and covariance from full fit to the BW only function
              double paramsBW[3] = {0.0};//,result->GetParams()[5];
              TMatrixDSym covArrBW(3);

              covArrBW(0,0) = result->GetCovarianceMatrix()(0,0);
              covArrBW(0,1) = result->GetCovarianceMatrix()(0,1);
              covArrBW(0,2) = result->GetCovarianceMatrix()(0,2);
              covArrBW(1,0) = result->GetCovarianceMatrix()(1,0);
              covArrBW(1,1) = result->GetCovarianceMatrix()(1,1);
              covArrBW(1,2) = result->GetCovarianceMatrix()(1,2);
              covArrBW(2,0) = result->GetCovarianceMatrix()(2,0);
              covArrBW(2,1) = result->GetCovarianceMatrix()(2,1);
              covArrBW(2,2) = result->GetCovarianceMatrix()(2,2);

              cout << "covArrBW" << endl;
              covArrBW.Print();

              paramsBW[0] = result->GetParams()[0];
              paramsBW[1] = result->GetParams()[1];
              paramsBW[2] = result->GetParams()[2];

              // define residual background only function
              TF1 *f_bg;
              if(i_poly == 0) f_bg = new TF1("f_bg",Poly1M, vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
              if(i_poly == 1) f_bg = new TF1("f_bg",Poly2M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
              //if(i_poly == 2) f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
              //f_bg->SetParameter(0,f_bw->GetParameter(0));
              f_bg->FixParameter(0,f_bw->GetParameter(0));
              f_bg->SetParameter(1,f_bw->GetParameter(3));
              f_bg->SetParameter(2,f_bw->GetParameter(4));
              if(i_poly >= 1) f_bg->SetParameter(3,f_bw->GetParameter(5));
              //if(i_poly >= 2) f_bg->SetParameter(4,f_bw->GetParameter(6));
              //f_bg->SetParError(0,f_bw->GetParError(0));
              f_bg->SetParError(1,f_bw->GetParError(3));
              f_bg->SetParError(2,f_bw->GetParError(4));
              if(i_poly >= 1) f_bg->SetParError(3,f_bw->GetParError(5));
              //if(i_poly >= 2) f_bg->SetParError(4,f_bw->GetParError(6));
  
              // propagate variables andn covariance from full fit to the residual background only function 
              double params1[3] = {0.0};//,result->GetParams()[5];
              if(i_poly == 0) 
              {
                params1[0] = result->GetParams()[0];
                params1[1] = result->GetParams()[3];
                params1[2] = result->GetParams()[4];
              }
              double params2[4] = {0.0};
              if(i_poly == 1) 
              {
                params2[0] = result->GetParams()[0];
                params2[1] = result->GetParams()[3];
                params2[2] = result->GetParams()[4];
                params2[3] = result->GetParams()[5];
              }
              //double params3[4] = {0.0};
              //if(i_poly == 2) 
              //{
              //  params3[0] = result->GetParams()[3];
              //  params3[1] = result->GetParams()[4];
              //  params3[2] = result->GetParams()[5];
              //  params3[2] = result->GetParams()[6];
              //}
  
              TMatrixDSym covArr1(3);
              TMatrixDSym covArr2(4);
              //TMatrixDSym covArr3(4);
              if(i_poly == 0)
              {
                covArr1(0,0) = result->GetCovarianceMatrix()(0,0);
                covArr1(0,1) = result->GetCovarianceMatrix()(0,3);
                covArr1(0,2) = result->GetCovarianceMatrix()(0,4);
                covArr1(1,0) = result->GetCovarianceMatrix()(3,0);
                covArr1(1,1) = result->GetCovarianceMatrix()(3,3);
                covArr1(1,2) = result->GetCovarianceMatrix()(3,4);
                covArr1(2,0) = result->GetCovarianceMatrix()(4,0);
                covArr1(2,1) = result->GetCovarianceMatrix()(4,3);
                covArr1(2,2) = result->GetCovarianceMatrix()(4,4);
              }
              if(i_poly == 1)
              {
                covArr2(0,0) = result->GetCovarianceMatrix()(0,0);
                covArr2(0,1) = result->GetCovarianceMatrix()(0,3);
                covArr2(0,2) = result->GetCovarianceMatrix()(0,4);
                covArr2(0,3) = result->GetCovarianceMatrix()(0,5);
                covArr2(1,0) = result->GetCovarianceMatrix()(3,0);
                covArr2(1,1) = result->GetCovarianceMatrix()(3,3);
                covArr2(1,2) = result->GetCovarianceMatrix()(3,4);
                covArr2(1,3) = result->GetCovarianceMatrix()(3,5);
                covArr2(2,0) = result->GetCovarianceMatrix()(4,0);
                covArr2(2,1) = result->GetCovarianceMatrix()(4,3);
                covArr2(2,2) = result->GetCovarianceMatrix()(4,4);
                covArr2(2,3) = result->GetCovarianceMatrix()(4,5);
                covArr2(3,0) = result->GetCovarianceMatrix()(5,0);
                covArr2(3,1) = result->GetCovarianceMatrix()(5,3);
                covArr2(3,2) = result->GetCovarianceMatrix()(5,4);
                covArr2(3,3) = result->GetCovarianceMatrix()(5,5);
              }
              //if(i_poly == 2)
              //{
              //  covArr3(0,0) = result->GetCovarianceMatrix()(3,3);
              //  covArr3(0,1) = result->GetCovarianceMatrix()(3,4);
              //  covArr3(0,2) = result->GetCovarianceMatrix()(3,5);
              //  covArr3(0,3) = result->GetCovarianceMatrix()(3,6);
              //  covArr3(1,0) = result->GetCovarianceMatrix()(4,3);
              //  covArr3(1,1) = result->GetCovarianceMatrix()(4,4);
              //  covArr3(1,2) = result->GetCovarianceMatrix()(4,5);
              //  covArr3(1,3) = result->GetCovarianceMatrix()(4,6);
              //  covArr3(2,0) = result->GetCovarianceMatrix()(5,3);
              //  covArr3(2,1) = result->GetCovarianceMatrix()(5,4);
              //  covArr3(2,2) = result->GetCovarianceMatrix()(5,5);
              //  covArr3(2,3) = result->GetCovarianceMatrix()(5,6);
              //  covArr3(3,0) = result->GetCovarianceMatrix()(6,3);
              //  covArr3(3,1) = result->GetCovarianceMatrix()(6,4);
              //  covArr3(3,2) = result->GetCovarianceMatrix()(6,5);
              //  covArr3(3,3) = result->GetCovarianceMatrix()(6,6);
              //}
  
              float bin_width = h_mMass[KEY]->GetBinWidth(1);
              //float Inte_start = Par_InteTheta[KEY_InteTheta_Poly][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1]-0.5*bin_width;
              //float Inte_stop  = Par_InteTheta[KEY_InteTheta_Poly][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1]+0.5*bin_width;

              float Inte_start_true = Par_InteTheta[KEY_InteTheta_Poly][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1];
              float Inte_stop_true  = Par_InteTheta[KEY_InteTheta_Poly][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1];
              //float Inte_start_true = vmsa::BW_Start[pid];//Par_InteTheta[KEY_InteTheta_Poly][0]-vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1];
              //float Inte_stop_true  = vmsa::BW_Stop[pid];//Par_InteTheta[KEY_InteTheta_Poly][0]+vmsa::nSigVecSys[i_sigma]*Par_InteTheta[KEY_InteTheta_Poly][1];

              int bin_start = h_mMass[KEY]->FindBin(Inte_start_true);
              int bin_stop  = h_mMass[KEY]->FindBin(Inte_stop_true);

              // set integration length to the edge of the bins corresponding to m +/- N*Gamma             
              float Inte_start = h_mMass[KEY]->GetBinLowEdge(bin_start);
              float Inte_stop  = h_mMass[KEY]->GetBinLowEdge(bin_stop) + bin_width;

              float counts_bg = f_bg->Integral(Inte_start,Inte_stop)/bin_width;
              float errors_bg; 
              if(i_poly == 0) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params1,covArr1.GetMatrixArray())/bin_width;
              if(i_poly == 1) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params2,covArr2.GetMatrixArray())/bin_width;
              //if(i_poly == 2) errors_bg = f_bg->IntegralError(Inte_start,Inte_stop,params3,covArr3.GetMatrixArray())/bin_width;
              cout << KEY_Poly << endl;
              cout << "Background counts = " << counts_bg << " +/- " << errors_bg << endl;    


              float bin_center = float(i_theta-4.5)/5.;
              float bin_center_phipsi = float(i_phipsi)*TMath::Pi()/5.+TMath::Pi()/10.;

              h_mCounts[KEY_countsbg]->SetBinContent(h_mCounts[KEY_countsbg]->FindBin(bin_center,bin_center_phipsi),counts_bg);
              h_mCounts[KEY_countsbg]->SetBinError(h_mCounts[KEY_countsbg]->FindBin(bin_center,bin_center_phipsi),errors_bg);
              

              g_reducedchi2[KEY_reducedchi2]->SetBinContent(g_reducedchi2[KEY_reducedchi2]->FindBin(bin_center,bin_center_phipsi),f_bw->GetChisquare()/f_bw->GetNDF());

              if(i_method == 0) // bin counting 
              {
                float counts = 0.0;
                float errors = 0.0;
                for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
                {
                  counts += h_mMass[KEY]->GetBinContent(i_bin);
                  errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
                }
                h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center,bin_center_phipsi),counts-counts_bg);
                h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center,bin_center_phipsi),TMath::Sqrt(errors+errors_bg*errors_bg));
                h_mCounts[KEY_countssigbg]->SetBinContent(h_mCounts[KEY_countssigbg]->FindBin(bin_center,bin_center_phipsi),counts);
                h_mCounts[KEY_countssigbg]->SetBinError(h_mCounts[KEY_countssigbg]->FindBin(bin_center,bin_center_phipsi),TMath::Sqrt(errors));
              }
              if(i_method == 1) // default: BW integration
              {
                //float counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
                //float errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
                //h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center,bin_center_phipsi),counts_bw-counts_bg);
                //h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center,bin_center_phipsi),TMath::Sqrt(errors_bw*errors_bw+errors_bg*errors_bg));
       
                // main signal counts
                float counts_bw = f_bw_only->Integral(Inte_start,Inte_stop)/bin_width;
                float errors_bw = f_bw_only->IntegralError(Inte_start,Inte_stop,paramsBW,covArrBW.GetMatrixArray())/bin_width;
                h_mCounts[KEY_counts]->SetBinContent(h_mCounts[KEY_counts]->FindBin(bin_center,bin_center_phipsi),counts_bw);
                h_mCounts[KEY_counts]->SetBinError(h_mCounts[KEY_counts]->FindBin(bin_center,bin_center_phipsi),errors_bw);

                // signal + background counts
                float counts_sigbg = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
                float errors_sigbg = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
                h_mCounts[KEY_countssigbg]->SetBinContent(h_mCounts[KEY_countssigbg]->FindBin(bin_center,bin_center_phipsi),counts_sigbg);
                h_mCounts[KEY_countssigbg]->SetBinError(h_mCounts[KEY_countssigbg]->FindBin(bin_center,bin_center_phipsi),errors_sigbg);

                //for(int i_bin = 1; i_bin <= h_mMass[KEY]->GetNbinsX(); i_bin++)
                for(int i_bin = bin_start; i_bin <= bin_stop; i_bin++)
                {
                  double binlowedge = h_mMass[KEY]->GetBinLowEdge(i_bin);
                  double integral = f_bw->Integral(binlowedge,binlowedge+bin_width)/bin_width;
                  //g_pull[KEY_Poly]->SetPoint(i_bin-1,h_mMass[KEY]->GetBinCenter(i_bin),(h_mMass[KEY]->GetBinContent(i_bin)-integral)/h_mMass[KEY]->GetBinError(i_bin));
                  //g_pull[KEY_Poly_diff]->SetPoint(i_bin-1,h_mMass[KEY]->GetBinCenter(i_bin),(h_mMass[KEY]->GetBinContent(i_bin)-integral)/h_mMass[KEY]->GetBinError(i_bin));
                  //g_pull[KEY_Poly_diffcum]->SetPoint(i_bin-1,h_mMass[KEY]->GetBinCenter(i_bin),(h_mMass[KEY]->GetBinContent(i_bin)-integral)/h_mMass[KEY]->GetBinError(i_bin));
                  g_pull[KEY_Poly]->SetPoint(i_bin-bin_start,h_mMass[KEY]->GetBinCenter(i_bin),(h_mMass[KEY]->GetBinContent(i_bin)-integral)/h_mMass[KEY]->GetBinError(i_bin));
                  g_pull[KEY_Poly_diff]->SetPoint(i_bin-bin_start,h_mMass[KEY]->GetBinCenter(i_bin),(h_mMass[KEY]->GetBinContent(i_bin)-integral));

                  // this is what is meant by cumulative data-fit 
                  if(i_bin == bin_start) g_pull[KEY_Poly_diffcum]->SetPoint(i_bin-bin_start,h_mMass[KEY]->GetBinCenter(i_bin),(h_mMass[KEY]->GetBinContent(i_bin)-integral));
                  else   
                  {
                    double xval, yval; 
                    g_pull[KEY_Poly_diffcum]->GetPoint(i_bin-bin_start-1,xval,yval);
                    g_pull[KEY_Poly_diffcum]->SetPoint(i_bin-bin_start,h_mMass[KEY]->GetBinCenter(i_bin),(h_mMass[KEY]->GetBinContent(i_bin)-integral)+yval);
                  }
                }
                //push back all fit parameters for plotting 
                Par[KEY_Poly].clear();
                Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(0)));
                Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(1)));
                Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(2)));
                Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(3)));
                Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(4)));
                if(i_poly >= 1) Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(5)));
                if(i_poly >= 2) Par[KEY_Poly].push_back(static_cast<float>(f_bw->GetParameter(6)));
              }
            }
          }
          float pt_mean = 0.5;//(vmsa::pt_low_2D[energy][i_pt]+vmsa::pt_up_2D[energy][i_pt])/2.0;
  
          // define 2D fit function for the SDMEs
          TF2 *f_rho = new TF2("f_rho",SpinDensity2Dcos,0.0,1.0,0,2.0*TMath::Pi(),6);
          f_rho->SetParameter(0,0.33333);
          f_rho->SetParameter(1,0.0);
          f_rho->SetParameter(2,0.0);
          f_rho->SetParameter(3,0.0);
          f_rho->SetParameter(4,0.0);
          f_rho->SetParameter(5,h_mCounts[KEY_counts]->GetMaximum());
          h_mCounts[KEY_counts]->Fit(f_rho,"NMRI");
          Par_rhoFit[KEY_counts].clear();
          Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(0)));
          Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(1)));
          Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(2)));
          Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(3)));
          Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(4)));
          Par_rhoFit[KEY_counts].push_back(static_cast<float>(f_rho->GetParameter(5)));

          // add raw values to TGraphs
          g_mRho[KEY_rho]->SetPoint(0,pt_mean,f_rho->GetParameter(0));
          g_mRho[KEY_rho]->SetPointError(0,0.0,0.0,f_rho->GetParError(0),f_rho->GetParError(0));
          g_mReal[KEY_real]->SetPoint(0,pt_mean,f_rho->GetParameter(1));
          g_mReal[KEY_real]->SetPointError(0,0.0,0.0,f_rho->GetParError(1),f_rho->GetParError(1));
          g_mImag[KEY_imag]->SetPoint(0,pt_mean,f_rho->GetParameter(2));
          g_mImag[KEY_imag]->SetPointError(0,0.0,0.0,f_rho->GetParError(2),f_rho->GetParError(2));
          g_mReRho1n1[KEY_rerho1n1]->SetPoint(0,pt_mean,f_rho->GetParameter(3));
          g_mReRho1n1[KEY_rerho1n1]->SetPointError(0,0.0,0.0,f_rho->GetParError(3),f_rho->GetParError(3));
          g_mImRho1n1[KEY_imrho1n1]->SetPoint(0,pt_mean,f_rho->GetParameter(4));
          g_mImRho1n1[KEY_imrho1n1]->SetPointError(0,0.0,0.0,f_rho->GetParError(4),f_rho->GetParError(4));
    
          for(int ipar = 0; ipar < 5; ipar++)
          {
            cout << "ipar = " << ipar << " = " << f_rho->GetParameter(0) << " +/- " << f_rho->GetParError(0) << endl;
          }
        }
      }
    }
  }


  // This plots the reduced chi2 for the fits in a 2D grid of cos(theta*) and beta
  {
    string outputname = Form("./figures/%s/%s/pTstudy/%s2D_Chi2D_%s_%s_Order%d_%s_pTdependence_SysIdx%d%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,etamode.c_str(),sysidx,flag.c_str());
    string output_start = Form("%s[",outputname.c_str());
    string output_stop = Form("%s]",outputname.c_str());
    TCanvas *c_chi = new TCanvas("c_chi","c_chi",10,10,1600,1200);
    c_chi->Print(output_start.c_str());
    c_chi->Divide(4,3);

    for(int i_norm = vmsa::Norm_start; i_norm < stop_norm; ++i_norm)
    {
      for(int i_sigma = vmsa::Sig_start; i_sigma < sig_stop; ++i_sigma)
      {
        for(int i_method = start_method; i_method < vmsa::Method_stop; ++i_method)
        {
          for(int i_poly = 0; i_poly < stop_poly; i_poly++)
          {
            c_chi->cd(1);
            c_chi->cd(1)->SetLeftMargin(0.15);
            c_chi->cd(1)->SetBottomMargin(0.15);
            c_chi->cd(1)->SetTicks(1,1);
            c_chi->cd(1)->SetGrid(0,0);

    
            string KEY_reducedchi2 = Form("chi2_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
            g_reducedchi2[KEY_reducedchi2]->SetStats(0);
            g_reducedchi2[KEY_reducedchi2]->SetTitle(Form("inorm=%d,isig=%d,imethod=%d,poly=%d",i_norm,i_sigma,i_method,i_poly+1));
            g_reducedchi2[KEY_reducedchi2]->GetXaxis()->SetTitle("cos(#theta*)");
            g_reducedchi2[KEY_reducedchi2]->GetYaxis()->SetTitle("#beta");
            g_reducedchi2[KEY_reducedchi2]->Draw("colz");
            
            c_chi->Update();
            c_chi->Print(outputname.c_str());
          }
        }
      }
    }
       
    c_chi->Print(output_stop.c_str());
  }  


  // this plots the fit status in 2D for all plots (4000) is good 
  {
    string outputname = Form("./figures/%s/%s/pTstudy/%s2D_FitStatus_%s_%s_Order%d_%s_pTdependence_SysIdx%d%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,etamode.c_str(),sysidx,flag.c_str());
    string output_start = Form("%s[",outputname.c_str());
    string output_stop = Form("%s]",outputname.c_str());
    TCanvas *c_chi = new TCanvas("c_chi","c_chi",10,10,1600,1200);
    c_chi->Print(output_start.c_str());
    c_chi->Divide(4,3);

    for(int i_norm = vmsa::Norm_start; i_norm < stop_norm; ++i_norm)
    {
      for(int i_sigma = vmsa::Sig_start; i_sigma < sig_stop; ++i_sigma)
      {
        for(int i_method = start_method; i_method < vmsa::Method_stop; ++i_method)
        {
          for(int i_poly = 0; i_poly < stop_poly; i_poly++)
          {
            c_chi->cd(1);
            c_chi->cd(1)->SetLeftMargin(0.15);
            c_chi->cd(1)->SetBottomMargin(0.15);
            c_chi->cd(1)->SetTicks(1,1);
            c_chi->cd(1)->SetGrid(0,0);

    
            string KEY_fitstatus = Form("status_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
            h_fitstatus[KEY_fitstatus]->SetStats(0);
            h_fitstatus[KEY_fitstatus]->SetTitle(Form("inorm=%d,isig=%d,imethod=%d,poly=%d",i_norm,i_sigma,i_method,i_poly+1));
            h_fitstatus[KEY_fitstatus]->GetXaxis()->SetTitle("cos(#theta*)");
            h_fitstatus[KEY_fitstatus]->GetYaxis()->SetTitle("#beta");
            h_fitstatus[KEY_fitstatus]->Draw("colz");
            
            c_chi->Update();
            c_chi->Print(outputname.c_str());
          }
        }
      }
    }
    c_chi->Print(output_stop.c_str());
  }  
#if _PlotQA_

  // QA plots which show all of the cos(theta*) bin fits for a single beta bin in one pdf page
  for(int i_poly = 0; i_poly < stop_poly; i_poly++)
  {
    for(int i_norm = 0; i_norm < stop_norm; i_norm++)
    {
      string outputname = Form("./figures/%s/%s/pTstudy/%s2D_allThetaYields%s_%s_Order%d_%s_Integrated_Norm%d_Poly%d_SysIdx%d%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,etamode.c_str(),i_norm,i_poly+1,sysidx,flag.c_str());
      string output_start = Form("%s[",outputname.c_str());
      
      //TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1200,900);
      TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1600,1200);
      TCanvas *c_diff2d = new TCanvas("c_diff2d","c_diff2d",10,10,500,500);
      c_diff->Print(output_start.c_str());
      for(int i_phipsi = 0;  i_phipsi < vmsa::Beta_Bins; i_phipsi++)
      {
        c_diff->Clear();
        c_diff->Divide(4,3);
        for(int i_theta = 0;  i_theta < vmsa::Cos_Theta_Star_Bins+1; ++i_theta)
        {
          c_diff->cd(i_theta+1);
          c_diff->cd(i_theta+1)->SetLeftMargin(0.15);
          c_diff->cd(i_theta+1)->SetBottomMargin(0.15);
          c_diff->cd(i_theta+1)->SetTicks(1,1);
          c_diff->cd(i_theta+1)->SetGrid(0,0);
          if( i_theta < vmsa::Cos_Theta_Star_Bins)
          { // draw data for each bin
            string KEY_QA = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
            h_mMass[KEY_QA]->SetTitle("");
            h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
            h_mMass[KEY_QA]->SetTitle(Form("%d#pi/5<#beta<%d#pi/5, %d/5<cos(#theta*)<%d/5",i_phipsi,(i_phipsi+1),i_theta-5,i_theta-4));
            h_mMass[KEY_QA]->GetXaxis()->SetTitleSize(0.05);
            h_mMass[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
            h_mMass[KEY_QA]->GetXaxis()->CenterTitle();
  
            h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(h_mMass[KEY_QA]->GetMinimum(),1.1*h_mMass[KEY_QA]->GetMaximum());
            h_mMass[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
            h_mMass[KEY_QA]->GetYaxis()->SetTitle("Yields");
            h_mMass[KEY_QA]->GetYaxis()->SetTitleSize(0.05);
            h_mMass[KEY_QA]->GetYaxis()->SetLabelSize(0.03);
            h_mMass[KEY_QA]->GetYaxis()->CenterTitle();
  
            h_mMass[KEY_QA]->SetMarkerStyle(24);
            h_mMass[KEY_QA]->SetMarkerColor(kGray+2);
            h_mMass[KEY_QA]->SetMarkerSize(1.2);
            h_mMass[KEY_QA]->Draw("pE");
            PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
  
            //string KEY_InteTheta = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",EP[order-1].c_str(),vmsa::Dca_start,1/*vmsa::nSigKaon_start*/,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
            //float x1_0 = Par_InteTheta[KEY_InteTheta][0] - 2.0*Par_InteTheta[KEY_InteTheta][1];
            //float x2_0 = Par_InteTheta[KEY_InteTheta][0] + 2.0*Par_InteTheta[KEY_InteTheta][1];
            //float x1_1 = Par_InteTheta[KEY_InteTheta][0] - 2.5*Par_InteTheta[KEY_InteTheta][1];
            //float x2_1 = Par_InteTheta[KEY_InteTheta][0] + 2.5*Par_InteTheta[KEY_InteTheta][1];
            //float x1_2 = Par_InteTheta[KEY_InteTheta][0] - 3.0*Par_InteTheta[KEY_InteTheta][1];
            //float x2_2 = Par_InteTheta[KEY_InteTheta][0] + 3.0*Par_InteTheta[KEY_InteTheta][1];
            //float y = h_mMass[KEY_QA]->GetBinContent(h_mMass[KEY_QA]->FindBin(Par_InteTheta[KEY_InteTheta][0]));
            //h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(-0.3*y,1.2*y);
            //float ymin = h_mMass[KEY_QA]->GetMinimum();
            //PlotLine(x1_0,x1_0,-0.3*y,y,4,2,2);
            //PlotLine(x2_0,x2_0,-0.3*y,y,4,2,2);
            //PlotLine(x1_1,x1_1,-0.3*y,y,1,2,2);
            //PlotLine(x2_1,x2_1,-0.3*y,y,1,2,2);
            //PlotLine(x1_2,x1_2,-0.3*y,y,4,2,2);
            //PlotLine(x2_2,x2_2,-0.3*y,y,4,2,2);
          }
          if(i_theta == vmsa::Cos_Theta_Star_Bins) // this is the integrated fit
          {
            string KEY_InteTheta_QA = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
            string KEY_InteTheta_QA_Poly = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);

            cout << KEY_InteTheta_QA << endl;
            cout << KEY_InteTheta_QA_Poly << endl;
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle("");
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle(Form("-1<cos(#theta*)<1"));
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-})");
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleSize(0.05);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleOffset(1.2);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->CenterTitle();
  
            if(i_phipsi == 0) h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetRangeUser(h_mMass_InteTheta[KEY_InteTheta_QA]->GetMinimum(),1.1*h_mMass_InteTheta[KEY_InteTheta_QA]->GetMaximum());
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetNdivisions(505,'N');
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitle("Yields");
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitleSize(0.05);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetLabelSize(0.03);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->CenterTitle();
  
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerStyle(24);
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerColor(kGray+2);
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerSize(1.2);
            h_mMass_InteTheta[KEY_InteTheta_QA]->Draw("pE");
            PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
            cout << "Draw Plot" << endl;

            // draw the functions
            TF1 *f_bw;
            if(i_poly == 0) f_bw = new TF1("f_bw",Poly1MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
            if(i_poly == 1) f_bw = new TF1("f_bw",Poly2MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
            if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
            f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA_Poly][0]);
            f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA_Poly][1]);
            f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA_Poly][2]);
            f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA_Poly][3]);
            f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_QA_Poly][4]);
            if(i_poly >= 1) f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_QA_Poly][5]);
            if(i_poly >= 2) f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta_QA_Poly][6]);
            f_bw->SetLineColor(kOrange+7);
            f_bw->SetLineStyle(1);
            f_bw->SetLineWidth(2);
            f_bw->Draw("l same");

            cout << "Draw BW" << endl;

            TF1 *f_bg;
            if(i_poly == 0) f_bg = new TF1("f_bg",Poly1M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
            if(i_poly == 1) f_bg = new TF1("f_bg",Poly2M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
            if(i_poly == 2) f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
            f_bg->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA_Poly][0]);
            f_bg->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA_Poly][3]);
            f_bg->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA_Poly][4]);
            if(i_poly >= 1) f_bg->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA_Poly][5]);
            if(i_poly >= 2) f_bg->SetParameter(4,Par_InteTheta[KEY_InteTheta_QA_Poly][6]);
            f_bg->SetLineColor(kBlue);
            f_bg->SetLineStyle(2);
            f_bg->SetLineWidth(2);
            f_bg->Draw("l same");
            cout << "Draw Background" << endl;
          }
        }
 

        for(int i_theta = 0;  i_theta < vmsa::Cos_Theta_Star_Bins; ++i_theta)
        { // now  draw the functions over the data
          c_diff->cd(i_theta+1);
          string KEY_QA = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
          string KEY_QA_Poly = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
          //cout << KEY_QA << endl;
          TF1 *f_bw;
          if(i_poly == 0) f_bw = new TF1("f_bw",Poly1MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
          if(i_poly == 1) f_bw = new TF1("f_bw",Poly2MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
          if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
          f_bw->SetParameter(0,Par[KEY_QA_Poly][0]);
          f_bw->SetParameter(1,Par[KEY_QA_Poly][1]);
          f_bw->SetParameter(2,Par[KEY_QA_Poly][2]);
          f_bw->SetParameter(3,Par[KEY_QA_Poly][3]);
          f_bw->SetParameter(4,Par[KEY_QA_Poly][4]);
          if(i_poly >= 1) f_bw->SetParameter(5,Par[KEY_QA_Poly][5]);
          if(i_poly >= 2) f_bw->SetParameter(6,Par[KEY_QA_Poly][6]);
          f_bw->SetLineColor(kOrange+7);
          f_bw->SetLineStyle(1);
          f_bw->SetLineWidth(2);
          f_bw->Draw("l same");
  
          TF1 *f_bg;
          if(i_poly == 0) f_bg = new TF1("f_bg",Poly1M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
          if(i_poly == 1) f_bg = new TF1("f_bg",Poly2M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
          if(i_poly == 2) f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
          f_bg->SetParameter(0,Par[KEY_QA_Poly][0]);
          f_bg->SetParameter(1,Par[KEY_QA_Poly][3]);
          f_bg->SetParameter(2,Par[KEY_QA_Poly][4]);
          if(i_poly >= 1) f_bg->SetParameter(3,Par[KEY_QA_Poly][5]);
          if(i_poly >= 2) f_bg->SetParameter(4,Par[KEY_QA_Poly][6]);
          f_bg->SetLineColor(kBlue);
          f_bg->SetLineStyle(2);
          f_bg->SetLineWidth(2);
          f_bg->Draw("l same");
  
          TLegend *leg1 = new TLegend(0.2,0.6,0.4,0.8);
          leg1->AddEntry(h_mMass[KEY_QA],"data","p");
          leg1->AddEntry(f_bw,"sig+res","l");
          leg1->AddEntry(f_bg,"res","l");
          leg1->Draw("same");
        }
        c_diff->Update();
        c_diff->Print(outputname.c_str());
 

        TH1F* hist = new TH1F("hist","hist",100,0.98,1.08);
        for(int i = 1; i <= hist->GetNbinsX(); i++) hist->SetBinContent(i,-999);

        hist->GetXaxis()->SetNdivisions(505,'N');
        hist->GetXaxis()->SetLabelSize(0.03);
        hist->GetXaxis()->SetTitle("M(K^{+},K^{-})");
        hist->GetXaxis()->SetTitleSize(0.05);
        hist->GetXaxis()->SetTitleOffset(1.2);
        hist->GetXaxis()->CenterTitle();
  
        hist->GetYaxis()->SetNdivisions(505,'N');
        hist->GetYaxis()->SetTitle("Pull");
        hist->GetYaxis()->SetTitleSize(0.05);
        hist->GetYaxis()->SetLabelSize(0.03);
        hist->GetYaxis()->CenterTitle();
  


        //pulls for each fit as a function of invariant mass
        for(int i_theta = 0;  i_theta < vmsa::Cos_Theta_Star_Bins; ++i_theta)
        {
          c_diff->cd(i_theta+1);
          string KEY_QA = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
          string KEY_QA_Poly = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
          //cout << KEY_QA << endl;
  
          hist->SetTitle(Form("%d#pi/5<#beta<%d#pi/5, %d/5<cos(#theta*)<%d/5",i_phipsi,(i_phipsi+1),i_theta-5,i_theta-4));
          hist->GetYaxis()->SetRangeUser(-3.0,3.0);
          hist->DrawCopy();                

          g_pull[KEY_QA_Poly]->SetMarkerStyle(20);
          g_pull[KEY_QA_Poly]->Draw("P same");             

        }
        c_diff->Update();
        c_diff->Print(outputname.c_str());



 
      }
      c_diff2d->cd();
      c_diff2d->cd()->SetLeftMargin(0.15);
      c_diff2d->cd()->SetRightMargin(0.15);
      c_diff2d->cd()->SetBottomMargin(0.15);
      c_diff2d->cd()->SetTicks(1,1);
      c_diff2d->cd()->SetGrid(0,0);

      // this plots the raw yields for each of the below variations in yield extraction parameters
      for(int i_sigma = 0; i_sigma < sig_stop; ++i_sigma)
      {
        for(int i_method = start_method; i_method < vmsa::Method_stop; ++i_method)
        {
          string KEY_counts_QA = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          //cout << KEY_counts_QA << endl;
          //if(i_sigma == vmsa::Sig_start && i_method == start_method)
          //{
            h_mCounts[KEY_counts_QA]->SetStats(0);
            h_mCounts[KEY_counts_QA]->SetTitle(Form("Raw Yields, isig=%d, imethod=%d",i_sigma,i_method));
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitle("cos(#theta)");
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleSize(0.05);
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleOffset(1.2);
            h_mCounts[KEY_counts_QA]->GetXaxis()->CenterTitle();
  
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA]->GetMinimum(),1.2*h_mCounts[KEY_counts_QA]->GetMaximum());
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetNdivisions(505,'N');
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitle("#beta");
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitleSize(0.05);
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetLabelSize(0.03);
            h_mCounts[KEY_counts_QA]->GetYaxis()->CenterTitle();
  
            h_mCounts[KEY_counts_QA]->SetMarkerStyle(20);
            h_mCounts[KEY_counts_QA]->SetMarkerColor(1);
            h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
            h_mCounts[KEY_counts_QA]->Draw("Colz");
          //}
          //else
          //{
          //  h_mCounts1D[KEY_counts_QA]->SetMarkerStyle(24);
          //  h_mCounts1D[KEY_counts_QA]->SetMarkerColor(i_sigma+10*i_method+1);
          //  h_mCounts1D[KEY_counts_QA]->SetMarkerSize(1.2);
          //  h_mCounts1D[KEY_counts_QA]->Draw("pE same");
          //}
          //TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
          //f_rho->SetParameter(0,Par_rhoFit[KEY_counts_QA][0]);
          //f_rho->SetParameter(1,Par_rhoFit[KEY_counts_QA][1]);
          //f_rho->SetLineColor(i_sigma+10*i_method+1);
          //f_rho->SetLineWidth(2);
          //f_rho->SetLineStyle(2);
          //f_rho->Draw("l same");
          c_diff2d->Update();
          c_diff2d->Print(outputname.c_str());
        }
      }
      
      string output_stop = Form("%s]",outputname.c_str());
      c_diff->Print(output_stop.c_str()); // close pdf file
      
    }
     // if(!random3D) c_diff->SaveAs("../figures/c_diff_2_SysIdx%d.pdf");
     // if(random3D) c_diff->SaveAs("../figures/3DRandom/c_diff_2_SysIdx%d.pdf");
  }



  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // This QA plotting contains more information.
  // Here we plot the following:
  //      Invariant mass distributions and fits for each beta bin in a given cos(theta*) bin on one page 
  //      Pulls vs invariant mass
  //      data-fit vs invariant mass
  //      cumulative data-fit vs invariant mass
  //      2D signal distributions
  //      Signal ratios for counts/integration
  //          - in 1D and 2D 
  //      2D signal, background, and signal + background distributions along with ratios for (bin counting) / integration  
  //          - 1D projections are also included
  for(int i_poly = 0; i_poly < stop_poly; i_poly++)
  {
    for(int i_norm = 0; i_norm < stop_norm; i_norm++)
    {
      string outputname = Form("./figures/%s/%s/pTstudy/%s2D_phibins_allThetaYields%s_%s_Order%d_%s_Integrated_Norm%d_Poly%d_SysIdx%d%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,etamode.c_str(),i_norm,i_poly+1,sysidx,flag.c_str());
      string output_start = Form("%s[",outputname.c_str());
      
      //TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1200,900);
      TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1600,1200);
      //TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1200,1600);
      //TCanvas *c_diff = new TCanvas("c_diff","c_diff",10,10,1200,1600);
      TCanvas *c_diff2d = new TCanvas("c_diff2d","c_diff2d",10,10,500,500);
      c_diff->Print(output_start.c_str());
      for(int i_theta = 0;  i_theta < vmsa::Cos_Theta_Star_Bins; ++i_theta)
      {
        c_diff->Clear();
        c_diff->Divide(4,3);
        //c_diff->Divide(3,4);
        for(int i_phipsi = 0;  i_phipsi < vmsa::Beta_Bins+1; i_phipsi++)
        {
          c_diff->cd(i_phipsi+1);
          c_diff->cd(i_phipsi+1)->SetLeftMargin(0.15);
          c_diff->cd(i_phipsi+1)->SetBottomMargin(0.15);
          c_diff->cd(i_phipsi+1)->SetTicks(1,1);
          c_diff->cd(i_phipsi+1)->SetGrid(0,0);
          if( i_phipsi < vmsa::Cos_Theta_Star_Bins)
          { // plot data for individual bins
            string KEY_QA = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
            h_mMass[KEY_QA]->SetStats(0);
            h_mMass[KEY_QA]->SetTitle("");
            h_mMass[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mMass[KEY_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mMass[KEY_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
            h_mMass[KEY_QA]->SetTitle(Form("%d#pi/5<#beta<%d#pi/5, %d/5<cos(#theta*)<%d/5",i_phipsi,(i_phipsi+1),i_theta-5,i_theta-4));
            h_mMass[KEY_QA]->GetXaxis()->SetTitleSize(0.05);
            h_mMass[KEY_QA]->GetXaxis()->SetTitleOffset(1.2);
            h_mMass[KEY_QA]->GetYaxis()->SetTitleOffset(1.2);
            h_mMass[KEY_QA]->GetXaxis()->CenterTitle();
  
            h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(h_mMass[KEY_QA]->GetMinimum(),1.1*h_mMass[KEY_QA]->GetMaximum());
            h_mMass[KEY_QA]->GetYaxis()->SetNdivisions(505,'N');
            h_mMass[KEY_QA]->GetYaxis()->SetTitle("Counts");
            h_mMass[KEY_QA]->GetYaxis()->SetTitleSize(0.05);
            h_mMass[KEY_QA]->GetYaxis()->SetLabelSize(0.03);
            h_mMass[KEY_QA]->GetYaxis()->CenterTitle();
  
            h_mMass[KEY_QA]->SetMarkerStyle(24);
            h_mMass[KEY_QA]->SetMarkerColor(kGray+2);
            h_mMass[KEY_QA]->SetMarkerSize(1.2);
            h_mMass[KEY_QA]->Draw("pE");
            PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
  
            //string KEY_InteTheta = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",i_pt,i_cent,EP[order-1].c_str(),vmsa::Dca_start,1/*vmsa::nSigKaon_start*/,vmsa::mPID[pid].c_str(),vmsa::Norm_QA);
            //float x1_0 = Par_InteTheta[KEY_InteTheta][0] - 2.0*Par_InteTheta[KEY_InteTheta][1];
            //float x2_0 = Par_InteTheta[KEY_InteTheta][0] + 2.0*Par_InteTheta[KEY_InteTheta][1];
            //float x1_1 = Par_InteTheta[KEY_InteTheta][0] - 2.5*Par_InteTheta[KEY_InteTheta][1];
            //float x2_1 = Par_InteTheta[KEY_InteTheta][0] + 2.5*Par_InteTheta[KEY_InteTheta][1];
            //float x1_2 = Par_InteTheta[KEY_InteTheta][0] - 3.0*Par_InteTheta[KEY_InteTheta][1];
            //float x2_2 = Par_InteTheta[KEY_InteTheta][0] + 3.0*Par_InteTheta[KEY_InteTheta][1];
            //float y = h_mMass[KEY_QA]->GetBinContent(h_mMass[KEY_QA]->FindBin(Par_InteTheta[KEY_InteTheta][0]));
            //h_mMass[KEY_QA]->GetYaxis()->SetRangeUser(-0.3*y,1.2*y);
            //float ymin = h_mMass[KEY_QA]->GetMinimum();
            //PlotLine(x1_0,x1_0,-0.3*y,y,4,2,2);
            //PlotLine(x2_0,x2_0,-0.3*y,y,4,2,2);
            //PlotLine(x1_1,x1_1,-0.3*y,y,1,2,2);
            //PlotLine(x2_1,x2_1,-0.3*y,y,1,2,2);
            //PlotLine(x1_2,x1_2,-0.3*y,y,4,2,2);
            //PlotLine(x2_2,x2_2,-0.3*y,y,4,2,2);
          }
          if(i_phipsi == vmsa::Beta_Bins)
          {  // plot the data and functions for the integrated bin
            string KEY_InteTheta_QA = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
            string KEY_InteTheta_QA_Poly = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);

            cout << KEY_InteTheta_QA << endl;
            cout << KEY_InteTheta_QA_Poly << endl;
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetStats(0);
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle("");
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetTitle(Form("-1<cos(#theta*)<1"));
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleSize(0.05);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->SetTitleOffset(1.2);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitleOffset(1.2);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetXaxis()->CenterTitle();
  
            if(i_theta == 0) h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetRangeUser(h_mMass_InteTheta[KEY_InteTheta_QA]->GetMinimum(),1.1*h_mMass_InteTheta[KEY_InteTheta_QA]->GetMaximum());
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetNdivisions(505,'N');
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitle("Counts");
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetTitleSize(0.05);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->SetLabelSize(0.03);
            h_mMass_InteTheta[KEY_InteTheta_QA]->GetYaxis()->CenterTitle();
  
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerStyle(24);
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerColor(kGray+2);
            h_mMass_InteTheta[KEY_InteTheta_QA]->SetMarkerSize(1.2);
            h_mMass_InteTheta[KEY_InteTheta_QA]->Draw("pE");
            PlotLine(vmsa::InvMass_low[pid],vmsa::InvMass_high[pid],0.0,0.0,1,2,2);
            cout << "Draw Plot" << endl;


            TF1 *f_bw;
            if(i_poly == 0) f_bw = new TF1("f_bw",Poly1MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
            if(i_poly == 1) f_bw = new TF1("f_bw",Poly2MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
            if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
            f_bw->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA_Poly][0]);
            f_bw->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA_Poly][1]);
            f_bw->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA_Poly][2]);
            f_bw->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA_Poly][3]);
            f_bw->SetParameter(4,Par_InteTheta[KEY_InteTheta_QA_Poly][4]);
            if(i_poly >= 1) f_bw->SetParameter(5,Par_InteTheta[KEY_InteTheta_QA_Poly][5]);
            if(i_poly >= 2) f_bw->SetParameter(6,Par_InteTheta[KEY_InteTheta_QA_Poly][6]);
            f_bw->SetLineColor(kOrange+7);
            f_bw->SetLineStyle(1);
            f_bw->SetLineWidth(2);
            f_bw->Draw("l same");

            cout << "Draw BW" << endl;

            TF1 *f_bg;
            if(i_poly == 0) f_bg = new TF1("f_bg",Poly1M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
            if(i_poly == 1) f_bg = new TF1("f_bg",Poly2M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
            if(i_poly == 2) f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
            f_bg->SetParameter(0,Par_InteTheta[KEY_InteTheta_QA_Poly][0]);
            f_bg->SetParameter(1,Par_InteTheta[KEY_InteTheta_QA_Poly][3]);
            f_bg->SetParameter(2,Par_InteTheta[KEY_InteTheta_QA_Poly][4]);
            if(i_poly >= 1) f_bg->SetParameter(3,Par_InteTheta[KEY_InteTheta_QA_Poly][5]);
            if(i_poly >= 2) f_bg->SetParameter(4,Par_InteTheta[KEY_InteTheta_QA_Poly][6]);
            f_bg->SetLineColor(kBlue);
            f_bg->SetLineStyle(2);
            f_bg->SetLineWidth(2);
            f_bg->Draw("l same");
            cout << "Draw Background" << endl;
            TLegend *leg1 = new TLegend(0.2,0.6,0.4,0.8);
            leg1->AddEntry(h_mMass_InteTheta[KEY_InteTheta_QA],"data","p");
            leg1->AddEntry(f_bw,"sig+res","l");
            leg1->AddEntry(f_bg,"res","l");
            leg1->Draw("same");
          }
        }
 

        //for(int i_theta = 0;  i_theta < vmsa::Cos_Theta_Star_Bins; ++i_theta)
        for(int i_phipsi = 0;  i_phipsi < vmsa::Beta_Bins; i_phipsi++)
        { // plot the BW + residual background functions for individual bins 
          c_diff->cd(i_phipsi+1);
          string KEY_QA = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
          string KEY_QA_Poly = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
          //cout << KEY_QA << endl;
          TF1 *f_bw;
          if(i_poly == 0) f_bw = new TF1("f_bw",Poly1MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],5);
          if(i_poly == 1) f_bw = new TF1("f_bw",Poly2MBreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],6);
          if(i_poly == 2) f_bw = new TF1("f_bw",Poly3BreitWigner,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],7);
          f_bw->SetParameter(0,Par[KEY_QA_Poly][0]);
          f_bw->SetParameter(1,Par[KEY_QA_Poly][1]);
          f_bw->SetParameter(2,Par[KEY_QA_Poly][2]);
          f_bw->SetParameter(3,Par[KEY_QA_Poly][3]);
          f_bw->SetParameter(4,Par[KEY_QA_Poly][4]);
          if(i_poly >= 1) f_bw->SetParameter(5,Par[KEY_QA_Poly][5]);
          if(i_poly >= 2) f_bw->SetParameter(6,Par[KEY_QA_Poly][6]);
          f_bw->SetLineColor(kOrange+7);
          f_bw->SetLineStyle(1);
          f_bw->SetLineWidth(2);
          f_bw->Draw("l same");
  
          TF1 *f_bg;
          if(i_poly == 0) f_bg = new TF1("f_bg",Poly1M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],3);
          if(i_poly == 1) f_bg = new TF1("f_bg",Poly2M,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
          if(i_poly == 2) f_bg = new TF1("f_bg",Poly3,vmsa::BW_Start[pid],vmsa::BW_Stop[pid],4);
          f_bg->SetParameter(0,Par[KEY_QA_Poly][0]);
          f_bg->SetParameter(1,Par[KEY_QA_Poly][3]);
          f_bg->SetParameter(2,Par[KEY_QA_Poly][4]);
          if(i_poly >= 1) f_bg->SetParameter(3,Par[KEY_QA_Poly][5]);
          if(i_poly >= 2) f_bg->SetParameter(4,Par[KEY_QA_Poly][6]);
          f_bg->SetLineColor(kBlue);
          f_bg->SetLineStyle(2);
          f_bg->SetLineWidth(2);
          f_bg->Draw("l same");
  
          TLegend *leg1 = new TLegend(0.2,0.6,0.4,0.8);
          leg1->AddEntry(h_mMass[KEY_QA],"data","p");
          leg1->AddEntry(f_bw,"sig+res","l");
          leg1->AddEntry(f_bg,"res","l");
          leg1->Draw("same");
        }
        c_diff->Update();
        c_diff->Print(outputname.c_str());
 

        TH1F* hist = new TH1F("hist","hist",100,0.98,1.08);
        for(int i = 1; i <= hist->GetNbinsX(); i++) hist->SetBinContent(i,-999);

        hist->GetXaxis()->SetNdivisions(505,'N');
        hist->GetXaxis()->SetLabelSize(0.03);
        hist->GetXaxis()->SetTitle("M(K^{+},K^{-})");
        hist->GetXaxis()->SetTitleSize(0.05);
        hist->GetXaxis()->SetTitleOffset(1.2);
        hist->GetXaxis()->CenterTitle();
  
        hist->GetYaxis()->SetNdivisions(505,'N');
        hist->GetYaxis()->SetTitle("Pull");
        hist->GetYaxis()->SetTitleSize(0.05);
        hist->GetYaxis()->SetLabelSize(0.03);
        hist->GetYaxis()->CenterTitle();
  


        // Plotting the pulls
        for(int i_phipsi = 0;  i_phipsi < vmsa::Beta_Bins; i_phipsi++)
        {
          c_diff->cd(i_phipsi+1);
          string KEY_QA = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm);
          string KEY_QA_Poly = Form("CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
          //cout << KEY_QA << endl;
  
          hist->SetTitle(Form("%d#pi/5<#beta<%d#pi/5, %d/5<cos(#theta*)<%d/5",i_phipsi,(i_phipsi+1),i_theta-5,i_theta-4));
          hist->GetYaxis()->SetRangeUser(-3.0,3.0);
          hist->DrawCopy();                

          g_pull[KEY_QA_Poly]->SetMarkerStyle(20);
          g_pull[KEY_QA_Poly]->Draw("P same");             

        }
        c_diff->Update();
        c_diff->Print(outputname.c_str());

        // Plotting the data-fit
        for(int i_phipsi = 0;  i_phipsi < vmsa::Beta_Bins; i_phipsi++)
        {
          c_diff->cd(i_phipsi+1);
          string KEY_QA_Poly = Form("diff_CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
          //cout << KEY_QA << endl;
  

          g_pull[KEY_QA_Poly]->SetMarkerStyle(20);
          g_pull[KEY_QA_Poly]->Draw("APE");          
          g_pull[KEY_QA_Poly]->SetTitle(Form("%d#pi/5<#beta<%d#pi/5, %d/5<cos(#theta*)<%d/5",i_phipsi,(i_phipsi+1),i_theta-5,i_theta-4));
          g_pull[KEY_QA_Poly]->GetXaxis()->SetTitle("M");
          g_pull[KEY_QA_Poly]->GetYaxis()->SetTitle("data-fit");
          g_pull[KEY_QA_Poly]->GetXaxis()->SetRangeUser(0.98,1.08);  

        }
        c_diff->Update();
        c_diff->Print(outputname.c_str());

        // Plotting the cumulative data-fit
        for(int i_phipsi = 0;  i_phipsi < vmsa::Beta_Bins; i_phipsi++)
        {
          c_diff->cd(i_phipsi+1);
          string KEY_QA_Poly = Form("diffcum_CosThetaStar_%d_Beta_%d_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Poly%d",i_theta,i_phipsi,EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_poly+1);
          //cout << KEY_QA << endl;
  

          g_pull[KEY_QA_Poly]->SetMarkerStyle(20);
          g_pull[KEY_QA_Poly]->Draw("APE");          
          g_pull[KEY_QA_Poly]->SetTitle(Form("%d#pi/5<#beta<%d#pi/5, %d/5<cos(#theta*)<%d/5",i_phipsi,(i_phipsi+1),i_theta-5,i_theta-4));
          g_pull[KEY_QA_Poly]->GetXaxis()->SetTitle("M");
          g_pull[KEY_QA_Poly]->GetYaxis()->SetTitle("cumulative data-fit");
          g_pull[KEY_QA_Poly]->GetXaxis()->SetRangeUser(0.98,1.08);  

        }
        c_diff->Update();
        c_diff->Print(outputname.c_str());

 
      }
 

      c_diff2d->cd();
      c_diff2d->cd()->SetLeftMargin(0.15);
      c_diff2d->cd()->SetRightMargin(0.15);
      c_diff2d->cd()->SetBottomMargin(0.15);
      c_diff2d->cd()->SetTicks(1,1);
      c_diff2d->cd()->SetGrid(0,0);









      // raw yield plots

      for(int i_sigma = 0; i_sigma < sig_stop; ++i_sigma)
      {
        for(int i_method = start_method; i_method < vmsa::Method_stop; ++i_method)
        {
          string KEY_counts_QA = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          //cout << KEY_counts_QA << endl;
          //if(i_sigma == vmsa::Sig_start && i_method == start_method)
          //{
            h_mCounts[KEY_counts_QA]->SetStats(0);
            h_mCounts[KEY_counts_QA]->SetTitle(Form("Raw Yields isig=%d, imethod=%d",i_sigma,i_method));
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitle("cos(#theta)");
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleSize(0.05);
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleOffset(1.2);
            h_mCounts[KEY_counts_QA]->GetXaxis()->CenterTitle();
  
            //h_mCounts[KEY_counts_QA]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA]->GetMinimum(),1.2*h_mCounts[KEY_counts_QA]->GetMaximum());
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetNdivisions(505,'N');
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitle("#beta");
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitleSize(0.05);
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetLabelSize(0.03);
            h_mCounts[KEY_counts_QA]->GetYaxis()->CenterTitle();
  
            h_mCounts[KEY_counts_QA]->SetMarkerStyle(20);
            h_mCounts[KEY_counts_QA]->SetMarkerColor(1);
            h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
            h_mCounts[KEY_counts_QA]->Draw("Colz");
          //}
          //else
          //{
          //  h_mCounts1D[KEY_counts_QA]->SetMarkerStyle(24);
          //  h_mCounts1D[KEY_counts_QA]->SetMarkerColor(i_sigma+10*i_method+1);
          //  h_mCounts1D[KEY_counts_QA]->SetMarkerSize(1.2);
          //  h_mCounts1D[KEY_counts_QA]->Draw("pE same");
          //}
          //TF1 *f_rho = new TF1("f_rho",SpinDensity,0.0,1.0,2);
          //f_rho->SetParameter(0,Par_rhoFit[KEY_counts_QA][0]);
          //f_rho->SetParameter(1,Par_rhoFit[KEY_counts_QA][1]);
          //f_rho->SetLineColor(i_sigma+10*i_method+1);
          //f_rho->SetLineWidth(2);
          //f_rho->SetLineStyle(2);
          //f_rho->Draw("l same");
          c_diff2d->Update();
          c_diff2d->Print(outputname.c_str());
        }
      }


      // here are the ratio plots for (bin counting)/Integration
      // these are only performed for sysidx == 0 since this has all default data level cuts
      if(sysidx == 0)
      {
        for(int i_sigma = 0; i_sigma < sig_stop; ++i_sigma)
        {
          for(int i_method = start_method; i_method < vmsa::Method_stop; ++i_method)
          {
            string KEY_counts_QA = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
            string KEY_counts_QA_copy = Form("ratio_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[0].c_str(),i_poly+1);

            cout << "At the top" << endl;
            cout << KEY_counts_QA << endl;
            if(i_method == 0) h_mCounts[KEY_counts_QA_copy] = (TH2F*) h_mCounts[KEY_counts_QA]->Clone(KEY_counts_QA_copy.c_str());
            if(i_method == 1) 
            { 

              // 2D ratio plot of signals
              c_diff2d->cd();
              h_mCounts[KEY_counts_QA_copy]->Divide(h_mCounts[KEY_counts_QA]);
            
              h_mCounts[KEY_counts_QA_copy]->SetStats(0);
              h_mCounts[KEY_counts_QA_copy]->SetTitle(Form("Counts/Inte, isig=%d",i_sigma));
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetNdivisions(505,'N');
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetLabelSize(0.03);
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetTitle("cos(#theta)");
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetTitleSize(0.05);
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetTitleOffset(1.2);
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->CenterTitle();
  
              //h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA_copy]->GetMinimum(),1.2*h_mCounts[KEY_counts_QA_copy]->GetMaximum());
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetNdivisions(505,'N');
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetTitle("#beta");
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetTitleSize(0.05);
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetLabelSize(0.03);
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->CenterTitle();
  
              h_mCounts[KEY_counts_QA_copy]->SetMarkerStyle(20);
              h_mCounts[KEY_counts_QA_copy]->SetMarkerColor(1);
              h_mCounts[KEY_counts_QA_copy]->SetMarkerSize(1.2);
              h_mCounts[KEY_counts_QA_copy]->Draw("Colz");
              c_diff2d->Update();
              c_diff2d->Print(outputname.c_str());
            }
           
            cout << "about to plot individual plots" << endl; 
            // cos(theta*) dependent ratios for each beta bin
            for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)
            {
              string histname1D = Form("%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
              string histname1Dr = Form("r_%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
              h_mCounts1D[histname1D] = new TH1F(histname1D.c_str(),histname1D.c_str(),vmsa::Cos_Theta_Star_Bins,-1,1);
              h_mCounts1D[histname1Dr] = new TH1F(histname1Dr.c_str(),histname1Dr.c_str(),vmsa::Cos_Theta_Star_Bins,-1,1);
           
              for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
              {
                double value = h_mCounts[KEY_counts_QA]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double error = h_mCounts[KEY_counts_QA]->GetBinError(i_costheta+1,i_phipsi+1);
                double rvalue = h_mCounts[KEY_counts_QA_copy]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double rerror = h_mCounts[KEY_counts_QA_copy]->GetBinError(i_costheta+1,i_phipsi+1);
             
                h_mCounts1D[histname1D.c_str()]->SetBinContent(i_costheta+1,value);
                h_mCounts1D[histname1D.c_str()]->SetBinError(i_costheta+1,error);
                h_mCounts1D[histname1Dr.c_str()]->SetBinContent(i_costheta+1,rvalue);
                h_mCounts1D[histname1Dr.c_str()]->SetBinError(i_costheta+1,rerror);
              }
            }
            // beta dependent ratios for each cos(theta*) bin
            for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
            {
              string histname1D = Form("%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
              string histname1Dr = Form("r_%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
              h_mCounts1D[histname1D] = new TH1F(histname1D.c_str(),histname1D.c_str(),vmsa::Beta_Bins,0.0,2.0*TMath::Pi());
              h_mCounts1D[histname1Dr] = new TH1F(histname1Dr.c_str(),histname1Dr.c_str(),vmsa::Beta_Bins,0.0,2.0*TMath::Pi());

              for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)         
              {
                double value = h_mCounts[KEY_counts_QA]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double error = h_mCounts[KEY_counts_QA]->GetBinError(i_costheta+1,i_phipsi+1);
                double rvalue = h_mCounts[KEY_counts_QA_copy]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double rerror = h_mCounts[KEY_counts_QA_copy]->GetBinError(i_costheta+1,i_phipsi+1);
             
                h_mCounts1D[histname1D.c_str()]->SetBinContent(i_phipsi+1,value);
                h_mCounts1D[histname1D.c_str()]->SetBinError(i_phipsi+1,error);
                h_mCounts1D[histname1Dr.c_str()]->SetBinContent(i_phipsi+1,rvalue);
                h_mCounts1D[histname1Dr.c_str()]->SetBinError(i_phipsi+1,rerror);
              }
            }

            c_diff->cd(11)->Clear();
            for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)
            {
              string histname1D = Form("%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
              
              c_diff->cd(i_phipsi+1);
              h_mCounts1D[histname1D]->SetTitle(Form("%s, i_sig=%d, %d#pi/5<#beta<%d#pi/5",vmsa::mInteMethod[i_method].c_str(),i_sigma,i_phipsi,i_phipsi+1));
              h_mCounts1D[histname1D]->GetXaxis()->SetTitle("cos#theta*");
              h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts");
              h_mCounts1D[histname1D]->Draw("PE");
            }
            c_diff->Update();
            c_diff->Print(outputname.c_str());
            for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
            {
              string histname1D = Form("%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
              
              c_diff->cd(i_costheta+1);
              h_mCounts1D[histname1D]->SetTitle(Form("%s, i_sig=%d, %d/5<cos#theta*<%d/5",vmsa::mInteMethod[i_method].c_str(),i_sigma,i_costheta-5,i_costheta-4));
              h_mCounts1D[histname1D]->GetXaxis()->SetTitle("#beta");
              h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts");
              h_mCounts1D[histname1D]->Draw("PE");
            }
            c_diff->Update();
            c_diff->Print(outputname.c_str());
            if(i_method == 1)
            { // only when we are on the second yield integration method can we plot these ratios since they do not exist before
              for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)
              {
                string histname1D = Form("r_%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
                
                c_diff->cd(i_phipsi+1);
                h_mCounts1D[histname1D]->SetTitle(Form("Counts/Inte, i_sig=%d, %d#pi/5<#beta<%d#pi/5",i_sigma,i_phipsi,i_phipsi+1));
                h_mCounts1D[histname1D]->GetXaxis()->SetTitle("cos#theta*");
                h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts/Inte");
                h_mCounts1D[histname1D]->Draw("PE");
              }
              c_diff->Update();
              c_diff->Print(outputname.c_str());
              for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
              {
                string histname1D = Form("r_%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
                
                c_diff->cd(i_costheta+1);
                h_mCounts1D[histname1D]->SetTitle(Form("Counts/Inte, i_sig=%d, %d/5<cos#theta*<%d/5",i_sigma,i_costheta-5,i_costheta-4));
                h_mCounts1D[histname1D]->GetXaxis()->SetTitle("#beta");
                h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts/Inte");
                h_mCounts1D[histname1D]->Draw("PE");
              }
              c_diff->Update();
              c_diff->Print(outputname.c_str());
            }
          }
        }


        // Here we have the same set of the plots but for signal + background
        for(int i_sigma = 0; i_sigma < sig_stop; ++i_sigma)
        {
          for(int i_method = start_method; i_method < vmsa::Method_stop; ++i_method)
          {
            string KEY_counts_QA = Form("sigbg_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
            string KEY_counts_QA_copy = Form("ratio_sigbg_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[0].c_str(),i_poly+1);

            h_mCounts[KEY_counts_QA]->SetStats(0);
            h_mCounts[KEY_counts_QA]->SetTitle(Form("Sig+Bg isig=%d, imethod=%d",i_sigma,i_method));
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitle("cos(#theta)");
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleSize(0.05);
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleOffset(1.2);
            h_mCounts[KEY_counts_QA]->GetXaxis()->CenterTitle();
  
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetNdivisions(505,'N');
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitle("#beta");
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitleSize(0.05);
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetLabelSize(0.03);
            h_mCounts[KEY_counts_QA]->GetYaxis()->CenterTitle();
  
            h_mCounts[KEY_counts_QA]->SetMarkerStyle(20);
            h_mCounts[KEY_counts_QA]->SetMarkerColor(1);
            h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
            h_mCounts[KEY_counts_QA]->Draw("Colz");

            c_diff2d->Update();
            c_diff2d->Print(outputname.c_str());

            if(i_method == 0) h_mCounts[KEY_counts_QA_copy] = (TH2F*) h_mCounts[KEY_counts_QA]->Clone(KEY_counts_QA_copy.c_str());
            if(i_method == 1) 
            {
              c_diff2d->cd();
              h_mCounts[KEY_counts_QA_copy]->Divide(h_mCounts[KEY_counts_QA]);
              h_mCounts[KEY_counts_QA_copy]->SetStats(0);
              h_mCounts[KEY_counts_QA_copy]->SetTitle(Form("Sig+Bg Counts/Inte, isig=%d",i_sigma));
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetNdivisions(505,'N');
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetLabelSize(0.03);
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetTitle("cos(#theta)");
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetTitleSize(0.05);
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetTitleOffset(1.2);
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->CenterTitle();
  
              //h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA_copy]->GetMinimum(),1.2*h_mCounts[KEY_counts_QA_copy]->GetMaximum());
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetNdivisions(505,'N');
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetTitle("#beta");
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetTitleSize(0.05);
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetLabelSize(0.03);
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->CenterTitle();
  
              h_mCounts[KEY_counts_QA_copy]->SetMarkerStyle(20);
              h_mCounts[KEY_counts_QA_copy]->SetMarkerColor(1);
              h_mCounts[KEY_counts_QA_copy]->SetMarkerSize(1.2);
              h_mCounts[KEY_counts_QA_copy]->Draw("Colz");
              c_diff2d->Update();
              c_diff2d->Print(outputname.c_str());
            }
            

            for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)
            {
              string histname1D = Form("%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
              string histname1Dr = Form("r_%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
              h_mCounts1D[histname1D] = new TH1F(histname1D.c_str(),histname1D.c_str(),vmsa::Cos_Theta_Star_Bins,-1,1);
              h_mCounts1D[histname1Dr] = new TH1F(histname1Dr.c_str(),histname1Dr.c_str(),vmsa::Cos_Theta_Star_Bins,-1,1);
            
              for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
              {
                double value = h_mCounts[KEY_counts_QA]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double error = h_mCounts[KEY_counts_QA]->GetBinError(i_costheta+1,i_phipsi+1);
                double rvalue = h_mCounts[KEY_counts_QA_copy]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double rerror = h_mCounts[KEY_counts_QA_copy]->GetBinError(i_costheta+1,i_phipsi+1);
             
                h_mCounts1D[histname1D.c_str()]->SetBinContent(i_costheta+1,value);
                h_mCounts1D[histname1D.c_str()]->SetBinError(i_costheta+1,error);
                h_mCounts1D[histname1Dr.c_str()]->SetBinContent(i_costheta+1,rvalue);
                h_mCounts1D[histname1Dr.c_str()]->SetBinError(i_costheta+1,rerror);
              }
            }

            for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
            {
              string histname1D = Form("%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
              string histname1Dr = Form("r_%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
              h_mCounts1D[histname1D] = new TH1F(histname1D.c_str(),histname1D.c_str(),vmsa::Beta_Bins,0.0,2.0*TMath::Pi());
              h_mCounts1D[histname1Dr] = new TH1F(histname1Dr.c_str(),histname1Dr.c_str(),vmsa::Beta_Bins,0.0,2.0*TMath::Pi());

              for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)         
              {
                double value = h_mCounts[KEY_counts_QA]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double error = h_mCounts[KEY_counts_QA]->GetBinError(i_costheta+1,i_phipsi+1);
                double rvalue = h_mCounts[KEY_counts_QA_copy]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double rerror = h_mCounts[KEY_counts_QA_copy]->GetBinError(i_costheta+1,i_phipsi+1);
             
                h_mCounts1D[histname1D.c_str()]->SetBinContent(i_phipsi+1,value);
                h_mCounts1D[histname1D.c_str()]->SetBinError(i_phipsi+1,error);
                h_mCounts1D[histname1Dr.c_str()]->SetBinContent(i_phipsi+1,rvalue);
                h_mCounts1D[histname1Dr.c_str()]->SetBinError(i_phipsi+1,rerror);
              }
            }

            c_diff->cd(11)->Clear();
            for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)
            {
              string histname1D = Form("%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
              
              c_diff->cd(i_phipsi+1);
              h_mCounts1D[histname1D]->SetTitle(Form("Sig+Bg %s, i_sig=%d, %d#pi/5<#beta<%d#pi/5",vmsa::mInteMethod[i_method].c_str(),i_sigma,i_phipsi,i_phipsi+1));
              h_mCounts1D[histname1D]->GetXaxis()->SetTitle("cos#theta*");
              h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts");
              h_mCounts1D[histname1D]->Draw("PE");
            }
            c_diff->Update();
            c_diff->Print(outputname.c_str());
            for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
            {
              string histname1D = Form("%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
              
              c_diff->cd(i_costheta+1);
              h_mCounts1D[histname1D]->SetTitle(Form("Sig+Bg %s, i_sig=%d, %d/5<cos#theta*<%d/5",vmsa::mInteMethod[i_method].c_str(),i_sigma,i_costheta-5,i_costheta-4));
              h_mCounts1D[histname1D]->GetXaxis()->SetTitle("#beta");
              h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts");
              h_mCounts1D[histname1D]->Draw("PE");
            }
            c_diff->Update();
            c_diff->Print(outputname.c_str());
            if(i_method == 1)
            {
              for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)
              {
                string histname1D = Form("r_%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
                
                c_diff->cd(i_phipsi+1);
                h_mCounts1D[histname1D]->SetTitle(Form("Sig+Bg Counts/Inte, i_sig=%d, %d#pi/5<#beta<%d#pi/5",i_sigma,i_phipsi,i_phipsi+1));
                h_mCounts1D[histname1D]->GetXaxis()->SetTitle("cos#theta*");
                h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts/Inte");
                h_mCounts1D[histname1D]->Draw("PE");
              }
              c_diff->Update();
              c_diff->Print(outputname.c_str());
              for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
              {
                string histname1D = Form("r_%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
                
                c_diff->cd(i_costheta+1);
                h_mCounts1D[histname1D]->SetTitle(Form("Sig+Bg Counts/Inte, i_sig=%d, %d/5<cos#theta*<%d/5",i_sigma,i_costheta-5,i_costheta-4));
                h_mCounts1D[histname1D]->GetXaxis()->SetTitle("#beta");
                h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts/Inte");
                h_mCounts1D[histname1D]->Draw("PE");
              }
              c_diff->Update();
              c_diff->Print(outputname.c_str());
            }







            // here we have the same set of plots for background only (background is the same in both cases so this is redundant) 

            KEY_counts_QA = Form("bg_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
            KEY_counts_QA_copy = Form("ratio_bg_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[0].c_str(),i_poly+1);

            h_mCounts[KEY_counts_QA]->SetStats(0);
            h_mCounts[KEY_counts_QA]->SetTitle(Form("Bg isig=%d, imethod=%d",i_sigma,i_method));
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetNdivisions(505,'N');
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetLabelSize(0.03);
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitle("cos(#theta)");
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleSize(0.05);
            h_mCounts[KEY_counts_QA]->GetXaxis()->SetTitleOffset(1.2);
            h_mCounts[KEY_counts_QA]->GetXaxis()->CenterTitle();
  
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetNdivisions(505,'N');
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitle("#beta");
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetTitleSize(0.05);
            h_mCounts[KEY_counts_QA]->GetYaxis()->SetLabelSize(0.03);
            h_mCounts[KEY_counts_QA]->GetYaxis()->CenterTitle();
  
            h_mCounts[KEY_counts_QA]->SetMarkerStyle(20);
            h_mCounts[KEY_counts_QA]->SetMarkerColor(1);
            h_mCounts[KEY_counts_QA]->SetMarkerSize(1.2);
            h_mCounts[KEY_counts_QA]->Draw("Colz");

            c_diff2d->Update();
            c_diff2d->Print(outputname.c_str());


            if(i_method == 0) h_mCounts[KEY_counts_QA_copy] = (TH2F*) h_mCounts[KEY_counts_QA]->Clone(KEY_counts_QA_copy.c_str());
            if(i_method == 1) 
            {
              c_diff2d->cd();
              h_mCounts[KEY_counts_QA_copy]->Divide(h_mCounts[KEY_counts_QA]);
              h_mCounts[KEY_counts_QA_copy]->SetStats(0);
              h_mCounts[KEY_counts_QA_copy]->SetTitle(Form("Bg Counts/Inte, isig=%d",i_sigma));
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetNdivisions(505,'N');
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetLabelSize(0.03);
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetTitle("cos(#theta)");
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetTitleSize(0.05);
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->SetTitleOffset(1.2);
              h_mCounts[KEY_counts_QA_copy]->GetXaxis()->CenterTitle();
  
              //h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_counts_QA_copy]->GetMinimum(),1.2*h_mCounts[KEY_counts_QA_copy]->GetMaximum());
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetNdivisions(505,'N');
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetTitle("#beta");
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetTitleSize(0.05);
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->SetLabelSize(0.03);
              h_mCounts[KEY_counts_QA_copy]->GetYaxis()->CenterTitle();
  
              h_mCounts[KEY_counts_QA_copy]->SetMarkerStyle(20);
              h_mCounts[KEY_counts_QA_copy]->SetMarkerColor(1);
              h_mCounts[KEY_counts_QA_copy]->SetMarkerSize(1.2);
              h_mCounts[KEY_counts_QA_copy]->Draw("Colz");
              c_diff2d->Update();
              c_diff2d->Print(outputname.c_str());
            }

           
            for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)
            {
              string histname1D = Form("%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
              string histname1Dr = Form("r_%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
              h_mCounts1D[histname1D] = new TH1F(histname1D.c_str(),histname1D.c_str(),vmsa::Cos_Theta_Star_Bins,-1,1);
              h_mCounts1D[histname1Dr] = new TH1F(histname1Dr.c_str(),histname1Dr.c_str(),vmsa::Cos_Theta_Star_Bins,-1,1);
            
              for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
              {
                double value = h_mCounts[KEY_counts_QA]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double error = h_mCounts[KEY_counts_QA]->GetBinError(i_costheta+1,i_phipsi+1);
                double rvalue = h_mCounts[KEY_counts_QA_copy]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double rerror = h_mCounts[KEY_counts_QA_copy]->GetBinError(i_costheta+1,i_phipsi+1);
             
                h_mCounts1D[histname1D.c_str()]->SetBinContent(i_costheta+1,value);
                h_mCounts1D[histname1D.c_str()]->SetBinError(i_costheta+1,error);
                h_mCounts1D[histname1Dr.c_str()]->SetBinContent(i_costheta+1,rvalue);
                h_mCounts1D[histname1Dr.c_str()]->SetBinError(i_costheta+1,rerror);
              }
            }

            for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
            {
              string histname1D = Form("%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
              string histname1Dr = Form("r_%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
              h_mCounts1D[histname1D] = new TH1F(histname1D.c_str(),histname1D.c_str(),vmsa::Beta_Bins,0.0,2.0*TMath::Pi());
              h_mCounts1D[histname1Dr] = new TH1F(histname1Dr.c_str(),histname1Dr.c_str(),vmsa::Beta_Bins,0.0,2.0*TMath::Pi());

              for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)         
              {
                double value = h_mCounts[KEY_counts_QA]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double error = h_mCounts[KEY_counts_QA]->GetBinError(i_costheta+1,i_phipsi+1);
                double rvalue = h_mCounts[KEY_counts_QA_copy]->GetBinContent(i_costheta+1,i_phipsi+1);              
                double rerror = h_mCounts[KEY_counts_QA_copy]->GetBinError(i_costheta+1,i_phipsi+1);
             
                h_mCounts1D[histname1D.c_str()]->SetBinContent(i_phipsi+1,value);
                h_mCounts1D[histname1D.c_str()]->SetBinError(i_phipsi+1,error);
                h_mCounts1D[histname1Dr.c_str()]->SetBinContent(i_phipsi+1,rvalue);
                h_mCounts1D[histname1Dr.c_str()]->SetBinError(i_phipsi+1,rerror);
              }
            }

            c_diff->cd(11)->Clear();
            for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)
            {
              string histname1D = Form("%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
              
              c_diff->cd(i_phipsi+1);
              h_mCounts1D[histname1D]->SetTitle(Form("Bg %s, i_sig=%d, %d#pi/5<#beta<%d#pi/5",vmsa::mInteMethod[i_method].c_str(),i_sigma,i_phipsi,i_phipsi+1));
              h_mCounts1D[histname1D]->GetXaxis()->SetTitle("cos#theta*");
              h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts");
              h_mCounts1D[histname1D]->Draw("PE");
            }
            c_diff->Update();
            c_diff->Print(outputname.c_str());
            for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
            {
              string histname1D = Form("%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
              
              c_diff->cd(i_costheta+1);
              h_mCounts1D[histname1D]->SetTitle(Form("Bg %s, i_sig=%d, %d/5<cos#theta*<%d/5",vmsa::mInteMethod[i_method].c_str(),i_sigma,i_costheta-5,i_costheta-4));
              h_mCounts1D[histname1D]->GetXaxis()->SetTitle("#beta");
              h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts");
              h_mCounts1D[histname1D]->Draw("PE");
            }
            c_diff->Update();
            c_diff->Print(outputname.c_str());
            if(i_method == 1)
            {
              for(int i_phipsi = 0; i_phipsi < vmsa::Beta_Bins; i_phipsi++)
              {
                string histname1D = Form("r_%s_phi%d",KEY_counts_QA.c_str(),i_phipsi); 
                
                c_diff->cd(i_phipsi+1);
                h_mCounts1D[histname1D]->SetTitle(Form("Counts/Inte Bg, i_sig=%d, %d#pi/5<#beta<%d#pi/5",i_sigma,i_phipsi,i_phipsi+1));
                h_mCounts1D[histname1D]->GetXaxis()->SetTitle("cos#theta*");
                h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts/Inte");
                h_mCounts1D[histname1D]->Draw("PE");
              }
              c_diff->Update();
              c_diff->Print(outputname.c_str());
              for(int i_costheta = 0; i_costheta < vmsa::Cos_Theta_Star_Bins; i_costheta++)
              {
                string histname1D = Form("r_%s_cos%d",KEY_counts_QA.c_str(),i_costheta); 
                
                c_diff->cd(i_costheta+1);
                h_mCounts1D[histname1D]->SetTitle(Form("Counts/Inte Bg, i_sig=%d, %d/5<cos#theta*<%d/5",i_sigma,i_costheta-5,i_costheta-4));
                h_mCounts1D[histname1D]->GetXaxis()->SetTitle("#beta");
                h_mCounts1D[histname1D]->GetYaxis()->SetTitle("Counts/Inte");
                h_mCounts1D[histname1D]->Draw("PE");
              }
              c_diff->Update();
              c_diff->Print(outputname.c_str());
            }



          }
        }
      }
      string output_stop = Form("%s]",outputname.c_str());
      c_diff->Print(output_stop.c_str()); // close pdf file
    }
     // if(!random3D) c_diff->SaveAs("../figures/c_diff_2_SysIdx%d.pdf");
     // if(random3D) c_diff->SaveAs("../figures/3DRandom/c_diff_2_SysIdx%d.pdf");
  
  }
#endif
//
//  TCanvas *c_rho = new TCanvas("c_rho","c_rho",10,10,800,800);
//  c_rho->cd();
//  c_rho->cd()->SetLeftMargin(0.15);
//  c_rho->cd()->SetBottomMargin(0.15);
//  c_rho->cd()->SetTicks(1,1);
//  c_rho->cd()->SetGrid(0,0);
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.05,9.95);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,5.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.03);
  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.05);
  h_frame->GetXaxis()->SetTitleOffset(1.2);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.2,0.5);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.05);
  h_frame->GetYaxis()->SetLabelSize(0.03);
  h_frame->GetYaxis()->CenterTitle();

  // here we plot the raw SDMEs
  
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

  cout << "Created c_spin" << endl;   

  string KEY_rho = Form("rhoRaw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),          vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_real = Form("realRaw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),        vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_imag = Form("imagRaw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),        vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_rerho1n1 = Form("rerho1n1Raw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  string KEY_imrho1n1 = Form("imrho1n1Raw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);

  cout << KEY_rho << endl;
  cout << KEY_real << endl;
  cout << KEY_imag << endl;
  cout << KEY_rerho1n1 << endl;
  cout << KEY_imrho1n1 << endl;

  c_spin->cd(1);
  g_mRho[KEY_rho]->GetYaxis()->SetTitle("#rho_{00}");
  g_mRho[KEY_rho]->GetXaxis()->SetTitle("p_{T} GeV/c");
  g_mRho[KEY_rho]->SetMarkerStyle(20); 
  g_mRho[KEY_rho]->Draw("APE");

  c_spin->cd(2);
  g_mReal[KEY_real]->GetYaxis()->SetTitle("Re(#rho_{10}-#rho_{0-1})");
  g_mReal[KEY_real]->GetXaxis()->SetTitle("p_{T} GeV/c");
  g_mReal[KEY_real]->SetMarkerStyle(20); 
  g_mReal[KEY_real]->Draw("APE");

  c_spin->cd(3);
  g_mImag[KEY_imag]->GetYaxis()->SetTitle("Im(#rho_{10}-#rho_{0-1})");
  g_mImag[KEY_imag]->GetXaxis()->SetTitle("p_{T} GeV/c");
  g_mImag[KEY_imag]->SetMarkerStyle(20); 
  g_mImag[KEY_imag]->Draw("APE");

  c_spin->cd(4);
  g_mReRho1n1[KEY_rerho1n1]->GetYaxis()->SetTitle("Re(#rho_{1-1})");
  g_mReRho1n1[KEY_rerho1n1]->GetXaxis()->SetTitle("p_{T} GeV/c");
  g_mReRho1n1[KEY_rerho1n1]->SetMarkerStyle(20); 
  g_mReRho1n1[KEY_rerho1n1]->Draw("APE");

  c_spin->cd(5);
  g_mImRho1n1[KEY_imrho1n1]->GetYaxis()->SetTitle("Im(#rho_{1-1})");
  g_mImRho1n1[KEY_imrho1n1]->GetXaxis()->SetTitle("p_{T} GeV/c");
  g_mImRho1n1[KEY_imrho1n1]->SetMarkerStyle(20); 
  g_mImRho1n1[KEY_imrho1n1]->Draw("APE");


  string outputname = Form("./figures/%s/%s/pTstudy/INT_Global2D_ExtractedSpinDensityElements_%s_Order%d_SysIdx%d%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),order,sysidx,flag.c_str());
  c_spin->Print(outputname.c_str());


  // This writes the Centers and Widths file for the default case
  if(sysidx == 0) 
  {
    cout << "About to write file " << endl;
    // Create an output file stream object
    std::ofstream outFile(Form("CentersAndWidths_Int_Order%d_Method2.txt",order));

    // Check if the file opened successfully
    if (!outFile) {
        std::cerr << "Error opening file!" << std::endl;
        //return 1;
    }

    // Write data to the file
    outFile << "Pt    Centers    Widths" << endl;
    outFile << 0 << "    " << Centers << "    " << Widths <<endl;
    // Close the file
    outFile.close();

  // This writes the Widths and Widths error file for the default case
    std::ofstream outFileWE(Form("WidthsAndErrors_Int_Order%d_Method2.txt",order));

    // Check if the file opened successfully
    if (!outFileWE) {
        std::cerr << "Error opening file!" << std::endl;
        //return 1;
    }

    // Write data to the file
    outFileWE << "Pt    Widths    WidthsErr" << endl;
      outFileWE << 0 << "    " << Widths << "    " << WidthsErr <<endl;
    // Close the file
    outFileWE.close();
  }
 

  // Save the output TGraphs and the signal histograms needed in the correction code
  string outputfile = Form("../output/AuAu%s/%s/%s2DInt_RawPhiPtSys_%s_SysIdx%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),frame.c_str(),etamode.c_str(),sysidx);
  if(order == 1) outputfile = Form("../output/AuAu%s/%s/%s2DInt_RawPhiPtSys_%s_FirstOrder_SysIdx%d.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),frame.c_str(),etamode.c_str(),sysidx);
  if(random3D) outputfile = Form("../output/AuAu%s/%s/3DRandom/RawPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  for(int i_norm = vmsa::Norm_start; i_norm < stop_norm; ++i_norm)
  {
    for(int i_sigma = vmsa::Sig_start; i_sigma < sig_stop; ++i_sigma)
    {
      for(int i_method = start_method; i_method < vmsa::Method_stop; ++i_method)
      {
        for(int i_poly = 0; i_poly < stop_poly; i_poly++)
        {
          string KEY_rho = Form("rhoRaw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          string KEY_real = Form("realRaw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          string KEY_imag = Form("imagRaw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          string KEY_rerho1n1 = Form("rerho1n1Raw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          string KEY_imrho1n1 = Form("imrho1n1Raw_%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          g_mRho[KEY_rho]->SetName(KEY_rho.c_str());
          g_mRho[KEY_rho]->Write();
          g_mReal[KEY_real]->SetName(KEY_real.c_str());
          g_mReal[KEY_real]->Write();
          g_mImag[KEY_imag]->SetName(KEY_imag.c_str());
          g_mImag[KEY_imag]->Write();
          g_mReRho1n1[KEY_rerho1n1]->SetName(KEY_rerho1n1.c_str());
          g_mReRho1n1[KEY_rerho1n1]->Write();
          g_mImRho1n1[KEY_imrho1n1]->SetName(KEY_imrho1n1.c_str());
          g_mImRho1n1[KEY_imrho1n1]->Write();
          string KEY_counts = Form("%s_Dca_%d_Sig_%d_HF_%d_HR_%d_M2_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",EP[order-1].c_str(),vmsa::mSys[sysidx][0],vmsa::mSys[sysidx][1],vmsa::mSys[sysidx][2],vmsa::mSys[sysidx][3],vmsa::mSys[sysidx][4],vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          h_mCounts[KEY_counts]->Write();
          cout << KEY_counts << " has " << h_mCounts[KEY_counts]->Integral() << " phi-mesons" << endl;
          
        }
      }
    }
  }
  
  File_OutPut->Close();
}

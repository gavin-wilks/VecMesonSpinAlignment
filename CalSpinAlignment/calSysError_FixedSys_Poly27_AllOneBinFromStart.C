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

double FuncAD(double *x_val, double *par);

//void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color, int beamE);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

void calSysError_FixedSys_Poly27_AllOneBinFromStart(Int_t energy = 5, Int_t pid = 0, Int_t i_cent = 9, string correction = "Raw", bool random3D = false, int order = 2, int defaultF = 0, string etamode = "eta1_eta1", int yspectra = -1)//defaultF = 0 is BESII, defaultF = 1 is BESI
{
 
  std::string spectra = "";
  if(yspectra == 0) spectra = "_NoRapiditySpectra";
  if(yspectra == 1) spectra = "_NoRapiditySpectra_EP";
  if(yspectra == 2) spectra = "_NoRapiditySpectra_EP_noV2";
  if(yspectra == 3) spectra = "_WithRapiditySpectra";
  if(yspectra == 4) spectra = "_WithRapiditySpectra_HalfSigma";
  if(yspectra == 5) spectra = "_NoRapiditySpectra_FixedFirstEP";
  if(yspectra == 6) spectra = "_NoRapiditySpectra_FixedFirstEP_PhiPsi";
  if(yspectra == 7) spectra = "_NoRapiditySpectra_PhiPsi_EP";
  
  std::string EP[2] = {"","2nd"};
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int colorDiff_phi = 0;
  const float size_marker = 1.4;

  double Res_12 = 0.0;
  double Res_12_err = 0.0;

  if(energy == 3 && order == 1)
  {
    Res_12 = r14GeVep1[0];
    Res_12_err = r14GeVep1[1];
  }
  if(energy == 3 && order == 2)
  {
    Res_12 = r14GeVep2[0];
    Res_12_err = r14GeVep2[1];
  }
  if(energy == 4 && order == 1)
  {
    Res_12 = r19GeVep1[0];
    Res_12_err = r19GeVep1[1];
  }
  if(energy == 4 && order == 2)
  {
    Res_12 = r19GeVep2[0];
    Res_12_err = r19GeVep2[1];
  }


  string inputfileHframe = Form("../output/AuAu%s/%s/Poly/OneBigBinFromStart_RawPhiPtSys_eta1_eta1_Poly.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  //string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/%sPhiPtSys_FINALFORPRELIM.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/PtDependence_v2off_pToff/%sPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/PtDependenceEffThenAcc_EPRes/%sPhiPtSys_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str());
  //string inputfile = Form("../output/AuAu%s/%s/Poly/OneBigBinFromStart_%sPhiPtSys_%s_Poly%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),spectra.c_str());
  string inputfile = Form("../output/AuAu%s/%s/OneBigBinFromStart_%sPhiPtSys_%s_PolySys%s_NoPsi2Bin.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),spectra.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/Poly/%sPhiPtSys_%s_Poly_FirstOrder%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),spectra.c_str());
  if(correction == "RawRes") inputfile = Form("../output/AuAu%s/%s/Poly/RawPhiPtSys_%s_Poly%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str(),spectra.c_str());
  if(correction == "RawRes" && order == 1) inputfile = Form("../output/AuAu%s/%s/Poly/RawPhiPtSys_%s_Poly_FirstOrder%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str(),spectra.c_str());
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
            for(int i_poly = 0; i_poly < 3; i_poly++)
            {
              for(int i_F = 0; i_F < 1; i_F++)
              {
                for(int i_eff = 0; i_eff < 1; i_eff++)
                {
                  if(correction != "AccRes" && correction != "EffAcc")
                  {
                    if(i_F > 0 || i_eff > 0) continue;
	            string KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
	            g_mRho[KEY_rho] = (TGraphAsymmErrors*)File_InPut->Get(KEY_rho.c_str())->Clone();
                    if(correction == "RawRes")
                    {
                      for(int i = 0; i < g_mRho[KEY_rho]->GetN(); i++)
                      {      
                        double pt, rho_obs; 
                        g_mRho[KEY_rho]->GetPoint(i,pt,rho_obs);
                        double rho_obs_err = g_mRho[KEY_rho]->GetErrorYhigh(i);          
                        //cout << "rho_obs = " << rho_obs << " +/- " << rho_obs_err << endl;
                
                        double drhodobs = 4./(1.+3.*Res_12);  // calculate d(rho)/d(rho_obs)
                        double drhodR = -12.*(rho_obs - 1./3.)/(1.+3.*Res_12)/(1.+3.*Res_12); // calculate d(rho)/d(R)

                        double real_rho = 1./3. + 4./(1.+3.*Res_12)*(rho_obs - 1./3.);
                        double real_rho_error = TMath::Sqrt((rho_obs_err*rho_obs_err)*(drhodobs*drhodobs) + (Res_12_err*Res_12_err)*(drhodR*drhodR));                 
                        
                        //cout << "rho_real = " << real_rho << " +/- " << real_rho_error << endl;
                        g_mRho[KEY_rho]->SetPoint(i,pt,real_rho);
                        g_mRho[KEY_rho]->SetPointError(i,0.0,0.0,real_rho_error,real_rho_error);
                      }
                    }
	            //string KEY_rho_ptbin = Form("rhoFinalWeighted_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
	            //h_mRho[KEY_rho_ptbin] = (TH1D*)File_InPut->Get(KEY_rho_ptbin.c_str());
	            for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt) // pt loop
	            {
	              string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
	              h_mCounts[KEY_counts] = (TH1D*) File_InPut->Get(KEY_counts.c_str());
	            }
                  }
                  else
                  {
	            string KEY_rho = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_F_%d_Eff_%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1,i_F,i_eff);
	            g_mRho[KEY_rho] = (TGraphAsymmErrors*)File_InPut->Get(KEY_rho.c_str())->Clone();
	            string KEY_rho_ptbin = Form("rhoFinalWeighted_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_F_%d_Eff_%d",i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1,i_F,i_eff);
	            h_mRho[KEY_rho_ptbin] = (TH1D*)File_InPut->Get(KEY_rho_ptbin.c_str());
                  }
                }
              }
            }
	  }
	}
      }
    }
  }
  TH1F *h_frame = (TH1F*)File_InPutHframe->Get("h_frame");
 
  string KEY_Default = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
  if(correction == "AccRes" || correction == "EffAcc") KEY_Default = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_F_%d_Eff_0",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1,defaultF);
  cout << "DEFAULT: " << KEY_Default << endl;
  TGraphAsymmErrors *g_SysErrors = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_StatErrors = new TGraphAsymmErrors();
    
  double sysErr[vmsa::pt_rebin_last[energy]][6]; // N points with 5 sources of systematics for each

  double weight_rho00 = 0.0;
  double weight_error_stat_rho00 = 0.0;
  double weight_all = 0.0;

  for(Int_t i_point = 0; i_point < vmsa::pt_rebin_last[energy]; ++i_point)
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
    g_mRho[KEY_Default]->GetPoint(i_point,pt_def,rho_def); 
    double rho_defErr = g_mRho[KEY_Default]->GetErrorYhigh(i_point); 
    string KEY_Default_yield = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_point,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
    double weight;
   // TH1F *PtCos = (TH1F*)h_mCounts[KEY_Default_yield]->Clone();
   // weight = PtCos->Integral(1,7);
   
    weight_rho00 += rho_def/rho_defErr/rho_defErr;
    weight_all += 1./rho_defErr/rho_defErr;

    cout << "pt: " << pt_def << "   rho00: " << rho_def << " +/- " << rho_defErr << endl;
    
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
          string KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),1);
          if(correction == "AccRes" || correction == "EffAcc") KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_F_%d_Eff_0",i_cent,EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[i_method].c_str(),1,defaultF);
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
          string KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),1);
          if(correction == "AccRes" || correction == "EffAcc") KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_F_%d_Eff_0",i_cent,EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),1,defaultF);
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
          string KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),1);
          if(correction == "AccRes" || correction == "EffAcc") KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_F_%d_Eff_0",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),1,defaultF);
          double pt_sys, rho_sys;
          g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
          sysNorm[idx] = rho_sys;
          if(random3D) sysNorm[idx] = 4.*(sysNorm[idx]-1./3.)+1./3.;
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
          string KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          if(correction == "AccRes" || correction == "EffAcc") KEY = Form("rhoRaw_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d_F_%d_Eff_0",i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1,defaultF);
          double pt_sys, rho_sys;
          g_mRho[KEY]->GetPoint(i_point,pt_sys,rho_sys);
          sysPoly[idx] = rho_sys;
          if(random3D) sysPoly[idx] = 4.*(sysPoly[idx]-1./3.)+1./3.;
          idx++;
        }
      }
    }	 
    idx = 0;

    Double_t rho_min[5] = { TMath::MinElement(9,sysDca),
                            TMath::MinElement(9,sysNSig),
			    TMath::MinElement(9,sysNorm),
                            TMath::MinElement(5,sysPoly)
                            //(correction == "AccRes")? TMath::MinElement(5,sysF) : 0.0,
                            //(correction == "AccRes")? TMath::MinElement(5,sysEff) : 0.0
                          };

    Double_t rho_max[5] = { TMath::MaxElement(9,sysDca),
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

  const int nx = 4;
  char *label[nx] = {"CWR Fix Width","GAVIN Fix Width (Psi2 OFF)","GAVIN Fix Width (Psi2 ON)", "Nature"};
  float rho[nx] = {0.34430,0.344041,0.345318,0.346832};
  float rhoerr[nx] = {0.00080,0.00111084,0.00106476,0.00116909};

  TCanvas *c1 = new TCanvas("c1","",10,10,400,400);
  c1->SetGrid(0,0);
  c1->SetBottomMargin(0.25);
  c1->SetLeftMargin(0.20);
  TH1F *hint = new TH1F("hint","Au+Au 27 GeV p_{T} Integrated #rho_{00}, 1.2<p_{T}<5.4 GeV/c",nx,0,float(nx));
  hint->SetStats(0);
  //hint->SetFillColor(38);
  hint->GetYaxis()->SetRangeUser(0.340,0.349);
  hint->SetMarkerStyle(20);
  hint->GetYaxis()->SetTitle("#rho_{00}");
  for(int i = 1; i <= nx; i++)
  {
    hint->SetBinContent(i,rho[i-1]);
    hint->SetBinError(i,rhoerr[i-1]);
    hint->GetXaxis()->SetBinLabel(i,label[i-1]);
  }
  hint->Draw("pEX0");
  c1->SaveAs(Form("figures/%s/%s/pTstudy/OneBigBinFromStart_AUAU27GeV_Integrated_NoPsi2Bin.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));  

  //const int nx = 8;
  //char *label[nx] = {"(CWR) R18 Free Width","(CWR) R18 InvMass","(CWR) R18 Fixed Width - Nature Method","(Gavin) R18 - One Large p_{T} Bin","(Gavin) R18 - No VzVPD-VzTPC cut","(Gavin) R18 - Average Over p_{T}","(Nature) R18","(Nature) R11"};
  //float rho[nx] = {0.343486,0.344322,0.342628,0.345318,0.345082,0.3452,0.346832,0.344039};
  //float rhoerr[nx] = {0.0011847,0.00101282,0.000908083,0.00106476,0.00101743,0.001056,0.00116909,0.00220591};

  //TCanvas *c1 = new TCanvas("c1","",10,10,400,400);
  //c1->SetGrid(0,0);
  //c1->SetBottomMargin(0.25);
  //c1->SetLeftMargin(0.20);
  //TH1F *hint = new TH1F("hint","Au+Au 27 GeV p_{T} Integrated #rho_{00}, 1.2<p_{T}<5.4 GeV/c",nx,0,float(nx));
  //hint->SetStats(0);
  ////hint->SetFillColor(38);
  //hint->GetYaxis()->SetRangeUser(0.340,0.349);
  //hint->SetMarkerStyle(20);
  //hint->GetYaxis()->SetTitle("#rho_{00}");
  //for(int i = 1; i <= nx; i++)
  //{
  //  hint->SetBinContent(i,rho[i-1]);
  //  hint->SetBinError(i,rhoerr[i-1]);
  //  hint->GetXaxis()->SetBinLabel(i,label[i-1]);
  //}
  //hint->Draw("pEX0");
  //c1->SaveAs(Form("figures/%s/%s/pTstudy/OneBigBinFromStart_AUAU27GeV_Integrated.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));  

  const int nx2 = 4;
  char *label2[nx2] = {"CWR Data + CWR Code","CWR Data + GAVIN Code","GAVIN Data + GAVIN Code", "NATURE"};
  float rho2[nx2] = {0.34336,0.343631,0.345459,0.346832};
  float rhoerr2[nx2] = {0.00091,0.00123425,0.00108615,0.00116909};

  TCanvas *c12 = new TCanvas("c12","",10,10,400,400);
  c12->SetGrid(0,0);
  c12->SetBottomMargin(0.15);
  c12->SetLeftMargin(0.20);
  TH1F *hint2 = new TH1F("hint2","Au+Au 27GeV Run18 1.2<p_{T}<5.4 GeV/c, Fixed Width",nx2,0,float(nx2));
  hint2->SetStats(0);
  //hint->SetFillColor(38);
  hint2->GetYaxis()->SetRangeUser(0.340,0.349);
  hint2->SetMarkerStyle(20);
  hint2->GetYaxis()->SetTitle("#rho_{00}");
  for(int i = 1; i <= nx2; i++)
  {
    hint2->SetBinContent(i,rho2[i-1]);
    hint2->SetBinError(i,rhoerr2[i-1]);
    hint2->GetXaxis()->SetBinLabel(i,label2[i-1]);
  }
  hint2->Draw("pEX0");
  c12->SaveAs(Form("figures/%s/%s/pTstudy/OneBigBinFromStart_AUAU27GeV_Integrated_Run18Comp.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));  


/////////////////////////////////////////////////////////////////////////  
  const int nxr = 5;
  char *labelr[nxr] = {"(CWR) R18 Free Width","(CWR) R18 InvMass","(CWR) R18 Fixed Width - Nature Method","(Gavin) R18 - One Large p_{T} Bin","(Gavin) R18 - Average Over p_{T}"};
  float rhor[nxr+1] = {0.343486,0.344322,0.342628,0.3453,0.3452,0.346832};
  float rhoerrr[nxr+1] = {0.0011847,0.00101282,0.000908083,0.001065,0.001056,0.00116909};

  float rhod[nxr] = {0.0};  
  float rhoerrd[nxr] = {0.0};  

  for(int i = 0; i < nxr; i++)
  {
    rhod[i] = rhor[i]-rhor[5];
    rhoerrd[i] = rhod[i]*rhod[i]/(rhoerrr[i]*rhoerrr[i]+rhoerrr[5]*rhoerrr[5]);
    //rhoerrd[i] = rhod[i]*rhod[i]/TMath::Sqrt(rhoerrr[i]*rhoerrr[i]+rhoerrr[5]*rhoerrr[5]);
    //rhoerrd[i] = TMath::Sqrt(rhoerrr[i]*rhoerrr[i]+rhoerrr[5]*rhoerrr[5]-2*rhoerrr[i]*rhoerrr[5]);
    cout << "difference = " << rhod[i] << " +/- " << rhoerrd[i] << endl; 
  } 

  TH1F *hintd = new TH1F("hintd","Au+Au 27 Run 18 GeV p_{T} Integrated #rho_{00} Difference, 1.2<p_{T}<5.4 GeV/c",nxr,0,float(nxr));
  hintd->SetStats(0);
  //hint->SetFillColor(38);
  hintd->GetYaxis()->SetRangeUser(-0.02,0.02);
  hintd->GetYaxis()->SetTitle("#rho_{00}-#rho_{00,Nature}");
  hintd->SetMarkerStyle(20);
  for(int i = 1; i <= nxr; i++)
  {
    hintd->SetBinContent(i,rhod[i-1]);
    hintd->SetBinError(i,rhoerrd[i-1]);
    hintd->GetXaxis()->SetBinLabel(i,labelr[i-1]);
  }
  hintd->Draw("pEX0");
  //PlotLine(0.0,nx,0.0,0.0,1,2,2);
  c1->SaveAs(Form("figures/%s/%s/pTstudy/OneBigBinFromStart_AUAU27GeV_Integrated_Run18Difference.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str()));  
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

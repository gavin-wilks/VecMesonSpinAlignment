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

void calSysErrorRapidity_2D_OffDiag(Int_t energy = 4, Int_t pid = 0, string correction = "Raw", bool random3D = false, int order = 2, string etamode = "eta1_eta1", int yspectra = -1, int frameopt = 0, int y_padding = 2)//defaultF = 0 is BESII, defaultF = 1 is BESI
{

  std::string frame = "Global";
  if(frameopt == 1) frame = "Helicity";
 
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

  string parstring[5] = {"rho","real","imag","rerho1n1","imrho1n1"};

  string inputfileHframe = Form("../output/AuAu%s/%s/%s2DRapidity_RawPhiYSys_eta1_eta1.root",frame.c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
  string inputfile = Form("../output/AuAu%s/%s/%s2DRapidity_%sPhiYSys_%s%s.root",frame.c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),spectra.c_str());
  if(order == 1) inputfile = Form("../output/AuAu%s/%s/%s2DRapidity_%sPhiYSys_%s_FirstOrder%s.root",frame.c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),spectra.c_str());
  if(correction == "RawRes") inputfile = Form("../output/AuAu%s/%s/%s2DRapidity_RawPhiYSys_%s%s.root",frame.c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str(),spectra.c_str());
  if(correction == "RawRes" && order == 1) inputfile = Form("../output/AuAu%s/%s/%s2DRapidity_RawPhiYSys_%s_FirstOrder%s.root",frame.c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etamode.c_str(),spectra.c_str());
  if(random3D) inputfile = Form("../output/AuAu%s/%s/3DRandom/%sPhiPtSys.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  TFile *File_InPutHframe = TFile::Open(inputfileHframe.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TGraMap g_mPar;
  //TH2DMap h_mCounts;  
  for(int i_pt = vmsa::pt_rebin_first_y_2D[energy]; i_pt <= vmsa::pt_rebin_last_y_2D[energy]; ++i_pt) // pt loop
  {           
    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2D; i_cent++)
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
                  for(int i_par = 0; i_par < 5; i_par++)
                  {
                    string KEY = Form("%sRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                    g_mPar[KEY] = (TGraphAsymmErrors*)File_InPut->Get(KEY.c_str())->Clone();
                    //if(correction == "RawRes")
                    //{
                    //  for(int i = 0; i < g_mRho[KEY_rho]->GetN(); i++)
                    //  {      
                    //    double pt, rho_obs; 
                    //    g_mRho[KEY_rho]->GetPoint(i,pt,rho_obs);
                    //    double rho_obs_err = g_mRho[KEY_rho]->GetErrorYhigh(i);          
                    //    //cout << "rho_obs = " << rho_obs << " +/- " << rho_obs_err << endl;
                    //
                    //    double drhodobs = 4./(1.+3.*Res_12);  // calculate d(rho)/d(rho_obs)
                    //    double drhodR = -12.*(rho_obs - 1./3.)/(1.+3.*Res_12)/(1.+3.*Res_12); // calculate d(rho)/d(R)

                    //    double real_rho = 1./3. + 4./(1.+3.*Res_12)*(rho_obs - 1./3.);
                    //    double real_rho_error = TMath::Sqrt((rho_obs_err*rho_obs_err)*(drhodobs*drhodobs) + (Res_12_err*Res_12_err)*(drhodR*drhodR));                 
                    //    
                    //    //cout << "rho_real = " << real_rho << " +/- " << real_rho_error << endl;
                    //    g_mRho[KEY_rho]->SetPoint(i,pt,real_rho);
                    //    g_mRho[KEY_rho]->SetPointError(i,0.0,0.0,real_rho_error,real_rho_error);
                    //  }
                    //}
                  }
                  //for(int i_pt = vmsa::pt_rebin_first[energy]; i_pt < vmsa::pt_rebin_last[energy]; ++i_pt) // pt loop
                  //{
                  //  string KEY_counts = Form("pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                  //  h_mCounts[KEY_counts] = (TH2D*) File_InPut->Get(KEY_counts.c_str());
                  //}
                }
              }
            }
          }
        }
      }
    }
  }
  TH1F *h_frame = (TH1F*)File_InPutHframe->Get("h_frame");
  

  TGraMap g_SysErrors; 
  TGraMap g_StatErrors;

  double weight_par[5] = {0.0};
  double weight_error_stat_par[5] = {0.0};
  double weight_all[5] = {0.0};

  for(int i_par = 0; i_par < 5; i_par++)
  {
    for(int i_pt = vmsa::pt_rebin_first_y_2D[energy]; i_pt <= vmsa::pt_rebin_last_y_2D[energy]; ++i_pt) // pt loop
    {           
      for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2D; i_cent++)
      {
        string KEY_Default = Form("%sRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",parstring[i_par].c_str(),i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
        g_mSysErrors[KEY_Default] = new TGraphAsymmErrors();
        g_mStatErrors[KEY_Default] = new TGraphAsymmErrors();

        double sysErr[9][6]; // N points with 5 sources of systematics for each

        for(Int_t i_point = y_padding; i_point < vmsa::eta_total-y_padding; ++i_point)
        {
          double sysDca[9];
          double sysNSig[9];
          double sysNorm[9];
          double sysPoly[5];

          double y_def, rho_def;
          g_mPar[KEY_Default]->GetPoint(i_point,y_def,rho_def); 
          double rho_defErr = g_mPar[KEY_Default]->GetErrorYhigh(i_point); 
          double weight;
        
          weight_par[i_par] += rho_def/rho_defErr/rho_defErr;
          weight_all[i_par] += 1./rho_defErr/rho_defErr;

          cout << parstring[i_par] << ", y = " << y_def << ", value = " << rho_def << " +/- " << rho_defErr << endl;
          
          int idx = 0;
          for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
          { 
            for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
            {
              for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
              {
                if(i_dca == 0 && (i_sigma != 0 || i_method == 0)) continue;
                if(i_dca != 0 && i_sigma != 0 && i_method == 1) continue;
                string KEY = Form("%sRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),i_pt,i_cent,EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),1);
                double pt_sys, rho_sys;
                g_mPar[KEY]->GetPoint(i_point,pt_sys,rho_sys);
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
                string KEY = Form("%sRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),i_pt,i_cent,EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),1);
                double pt_sys, rho_sys;
                g_mPar[KEY]->GetPoint(i_point,pt_sys,rho_sys);
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
                string KEY = Form("%sRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),1);
                double pt_sys, rho_sys;
                g_mPar[KEY]->GetPoint(i_point,pt_sys,rho_sys);
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
                string KEY = Form("%sRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),i_pt,i_cent,EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                double pt_sys, rho_sys;
                g_mPar[KEY]->GetPoint(i_point,pt_sys,rho_sys);
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
                                };

          Double_t rho_max[5] = { TMath::MaxElement(9,sysDca),
                                  TMath::MaxElement(9,sysNSig),
                                  TMath::MaxElement(9,sysNorm),
                                  TMath::MaxElement(5,sysPoly)
                                };
        
          double SysError_rho = 0.0;
          for(int i = 0; i < 4; i++)
          {
            double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
            cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
            SysError_rho += sourcei;
          }
          SysError_rho = TMath::Sqrt(SysError_rho);

          g_SysErrors[KEY_Default]->SetPoint(i_point,y_def,rho_def);
          g_SysErrors[KEY_Default]->SetPointError(i_point,0.0,0.0,SysError_rho,SysError_rho);

          double StatError_rho = g_mPar[KEY_Default]->GetErrorYhigh(i_point);
          g_StatErrors[KEY_Default]->SetPoint(i_point,y_def,rho_def);
          g_StatErrors[KEY_Default]->SetPointError(i_point,0.0,0.0,StatError_rho,StatError_rho);
        }
      }
    }
    weight_par[i_par] /= weight_all[i_par];
    weight_error_stat_par[i_par][i_pt] = TMath::Sqrt(1./weight_all[i_par]);
    cout << std::setprecision(4);
    cout << "Final " << parstring[i_par] << " = " << weight_par[i_par] << " +/- " << weight_error_stat_par[i_par] << " (stat.) " << endl;
  }
      

  cout << "Finished pT bin calculations" << endl;
 
  vecFMap vRho;
  vecFMap vRhoErr;
  
  for(int i_par = 0; i_par < 5; i_par++)
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
                string KEY_Int = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                g_mStatErrors[KEY_Int] = new TGraphAsymmErrors();

                double weight_rho00f = 0.0;
                double weight_error_stat_rho00f = 0.0;
                double weight_allf = 0.0;

                for(Int_t i_point = y_padding; i_point < vmsa::eta_total-y_padding; ++i_point)
                {
                  double weight_rho00f_point = 0.0;
                  double weight_error_stat_rho00f_point = 0.0;
                  double weight_allf_point = 0.0;
                  
                  Double_t y, rho;
                  for(int i_pt = vmsa::pt_rebin_first_y_2D[energy]; i_pt <= vmsa::pt_rebin_last_y_2D[energy]; ++i_pt) // pt loop
                  {           
                    for(int i_cent = 0; i_cent < vmsa::cent_rebin_total_2D; i_cent++)
                    {
                      string KEY = Form("%sRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),i_pt,i_cent,EP[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      g_mPar[KEY]->GetPoint(i_point,y,rho);
                      double StatError_rho = g_mPar[KEY]->GetErrorYhigh(i_point);

                      weight_rho00f_point += rho/StatError_rho/StatError_rho;
                      weight_allf_point += 1./StatError_rho/StatError_rho; 

                      weight_rho00f += rho/StatError_rho/StatError_rho;
                      weight_allf += 1./StatError_rho/StatError_rho; 
                    }
                  }
                  weight_rho00f_point /= weight_allf_point;
                  weight_error_stat_rho00f_point = TMath::Sqrt(1./weight_allf_point);

                  g_mStatErrors[KEY_Int]->SetPoint(i_point,y,weight_rho00f);
                  g_mStatErrors[KEY_Int]->SetPointError(i_point,0.0,0.0,weight_error_stat_rho00f_point,weight_error_stat_rho00f_point);
                }
                weight_rho00f /= weight_allf;
                weight_error_stat_rho00f = TMath::Sqrt(1./weight_allf);
                
                vRho[KEY_Int].push_back(weight_rho00f);
                vRhoErr[KEY_Int].push_back(weight_error_stat_rho00f); 
              }
            }
          }
        }
      }
    }
  }

  //std::ofstream fw_sources(Form("ptDependenceSources_Integrated2060_Order%d_%s.txt",order,vmsa::mBeamEnergy[energy].c_str()), std::ofstream::out);
  ////fw.open();

  //fw_sources << "sqrt{s_{nn}}      EP_Order      drho_{dca}        drho_{nσ_{kaon}}       drho_{NormalizationRange}" << "\n";
  //std::ofstream fw(Form("ptDependence_Integrated2060_Order%d_%s.txt",order,vmsa::mBeamEnergy[energy].c_str()), std::ofstream::out);
  ////fw.open();

  //fw << "sqrt{s_{nn}}      EP_Order      dca(cm)        nσ_{kaon}        NormalizationRange         IntegrationMethod          IntegrationRange             ρ_{00}" << "\n";
  //if(correction == "AccRes" || correction == "EffAcc")
  //{ 
  for(int i_par = 0; i_par < 5; i_par++)
  {
  
    double rhofinal, rhofinalstat, rhofinalsys;

    string KEY_Default = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);

    double sysDca[9];
    double sysNSig[9];
    double sysNorm[9];
    double sysPoly[5];

    double rho_def = vRho[KEY_Default][0]; 
    double rho_defErr = vRhoErr[KEY_Default][0];//h_mRho[KEY_DefaultW]->GetBinError(7); 

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
          string KEY = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",parstring[i_par].c_str(),EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
          sysDca[idx] = vRho[KEY][0];
          //fw << vmsa::mBeamEnergy[energy] << "     " << order << "     " << std::fixed << std::setprecision(1) << vmsa::mDcaSys[i_dca] << "     " << vmsa::mNSigmaKaonSys[0] << "      " << "[" << std::setprecision(2) << vmsa::Norm_Start[0][0] << "," << vmsa::Norm_Stop[0][0] << "]" << "     " << vmsa::mInteMethod[i_method] << "     " << "[" << std::setprecision(1) << -vmsa::nSigVecSys[i_sigma] << "Γ," << vmsa::nSigVecSys[i_sigma] << "Γ]" << "       " << std::setprecision(4) << sysDca[idx] << "\n"; 
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
          string KEY = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",parstring[i_par].c_str(),EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str());
          sysNSig[idx] = vRho[KEY][0];//h_mRho[KEY]->GetBinContent(7);
          //fw << vmsa::mBeamEnergy[energy] << "     " << order << "     " << std::fixed << std::setprecision(1) << vmsa::mDcaSys[0] << "     " << vmsa::mNSigmaKaonSys[i_sig] << "      " << "[" << std::setprecision(2) << vmsa::Norm_Start[0][0] << "," << vmsa::Norm_Stop[0][0] << "]" << "     " << vmsa::mInteMethod[i_method] << "     " << "[" << std::setprecision(1) << -vmsa::nSigVecSys[i_sigma] << "Γ," << vmsa::nSigVecSys[i_sigma] << "Γ]" << "       " << std::setprecision(4) << sysNSig[idx] << "\n"; 
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
          //cout << idx << "   i_norm = " << i_norm << "   i_sigma = " << i_sigma << "   i_method = " << i_method << endl;
          string KEY = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",parstring[i_par].c_str(),EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str());
          sysNorm[idx] = vRho[KEY][0];//h_mRho[KEY]->GetBinContent(7);
          //if(i_norm < 2) fw << vmsa::mBeamEnergy[energy] << "     " << order << "     " << std::fixed << std::setprecision(1) << vmsa::mDcaSys[0] << "     " << vmsa::mNSigmaKaonSys[0] << "      " << "[" << std::setprecision(2) << vmsa::Norm_Start[0][i_sigma] << "," << vmsa::Norm_Stop[0][i_sigma] << "]" << "     " << vmsa::mInteMethod[i_method] << "     " << "[" << std::setprecision(1) << -vmsa::nSigVecSys[i_sigma] << "Γ," << vmsa::nSigVecSys[i_sigma] << "Γ]" << "       " << std::setprecision(4) << sysNorm[idx] << "\n"; 
          //if(i_norm == 2) fw << vmsa::mBeamEnergy[energy] << "     " << order << "     " << std::fixed << std::setprecision(1) << vmsa::mDcaSys[0] << "     " << vmsa::mNSigmaKaonSys[0] << "      " << "[" << std::setprecision(2) << vmsa::Norm_Start[0][0] << "," << vmsa::Norm_Stop[0][0] << "]" << "," << "[" << std::setprecision(2) << vmsa::Norm_Start[0][1] << "," << vmsa::Norm_Stop[0][1] << "]" << "     " << vmsa::mInteMethod[i_method] << "     " << "[" << std::setprecision(1) << -vmsa::nSigVecSys[i_sigma] << "Γ," << vmsa::nSigVecSys[i_sigma] << "Γ]" << "       " << std::setprecision(4) << sysNorm[idx] << "\n"; 
          idx++;
        }
      }
    }	 
    idx = 0;

    for(Int_t i_poly = 0; i_poly < 2; i_poly++)
    {
      for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
      {
        for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
        {
          if(i_poly == 0 && (i_sigma != 0 || i_method == 0)) continue;
          if(i_poly != 0 && i_sigma != 0 && i_method == 1) continue;
          string KEY = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
          sysPoly[idx] = vRho[KEY][0];//h_mRho[KEY]->GetBinContent(7);
          idx++;
        }
      }
    }
    idx = 0;

    Double_t rho_min[5] = { TMath::MinElement(9,sysDca),
                            TMath::MinElement(9,sysNSig),
                            TMath::MinElement(9,sysNorm),
                            TMath::MinElement(5,sysPoly)
                          };

    Double_t rho_max[5] = { TMath::MaxElement(9,sysDca),
                            TMath::MaxElement(9,sysNSig),
                            TMath::MaxElement(9,sysNorm),
                            TMath::MaxElement(5,sysPoly)
                          };
    //fw.close();
    
    //fw_sources << vmsa::mBeamEnergy[energy] << "     " << order << "     " << std::fixed << std::setprecision(5);
    double SysError_rho = 0.0;
    for(int i = 0; i < 4; i++)
    {
      double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
      //cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
      SysError_rho += sourcei;
      //fw_sources << TMath::Sqrt(sourcei) << "     ";
    }
    //fw_sources << "\n"; 
    //fw_sources.close(); 

    SysError_rho = TMath::Sqrt(SysError_rho);
    
    rhofinalsys = SysError_rho;

    cout << "Final corrected rho00 = " << std::fixed << std::setprecision(4) << rho_def << " +/- " << rho_defErr << " (stat) +/- " << SysError_rho << " (sys)" << endl;
    cout << "sigma from 0 = " << std::fixed << std::setprecision(2) << TMath::Abs((rho_def)/TMath::Sqrt(rho_defErr*rho_defErr + SysError_rho*SysError_rho)) << endl;
  }


  for(int i_par = 0; i_par < 5; i_par++)
  {
    string KEY_Default = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",parstring[i_par].c_str(),EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
    //g_mSysErrors[KEY_Default] = new TGraphAsymmErrors();
    g_mStatErrors[KEY_Default] = new TGraphAsymmErrors();

    double sysErr[9][6]; // N points with 5 sources of systematics for each

    for(Int_t i_point = y_padding; i_point < vmsa::eta_total_y_padding; ++i_point)
    {
      double sysDca[9];
      double sysNSig[9];
      double sysNorm[9];
      double sysPoly[5];

      double y_def, rho_def;
      g_mPar[KEY_Default]->GetPoint(i_point,y_def,rho_def); 
      double rho_defErr = g_mPar[KEY_Default]->GetErrorYhigh(i_point); 
      double weight;
    
      weight_par[i_par] += rho_def/rho_defErr/rho_defErr;
      weight_all[i_par] += 1./rho_defErr/rho_defErr;

      cout << parstring[i_par] << ", cent = " << y_def << ", value = " << rho_def << " +/- " << rho_defErr << endl;
      
      int idx = 0;
      for(Int_t i_dca = vmsa::Dca_start; i_dca < vmsa::Dca_stop; i_dca++)
      { 
        for(int i_sigma = vmsa::Sig_start; i_sigma < vmsa::Sig_stop; ++i_sigma)
        {
          for(int i_method = vmsa::Method_start; i_method < vmsa::Method_stop; ++i_method)
          {
            if(i_dca == 0 && (i_sigma != 0 || i_method == 0)) continue;
            if(i_dca != 0 && i_sigma != 0 && i_method == 1) continue;
            string KEY = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),EP[order-1].c_str(),i_dca,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),1);
            double pt_sys, rho_sys;
            g_mPar[KEY]->GetPoint(i_point,pt_sys,rho_sys);
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
            string KEY = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),EP[order-1].c_str(),0,i_sig,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),1);
            double pt_sys, rho_sys;
            g_mPar[KEY]->GetPoint(i_point,pt_sys,rho_sys);
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
            string KEY = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),1);
            double pt_sys, rho_sys;
            g_mPar[KEY]->GetPoint(i_point,pt_sys,rho_sys);
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
            string KEY = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",parstring[i_par].c_str(),EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
            double pt_sys, rho_sys;
            g_mPar[KEY]->GetPoint(i_point,pt_sys,rho_sys);
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
                            };

      Double_t rho_max[5] = { TMath::MaxElement(9,sysDca),
                              TMath::MaxElement(9,sysNSig),
                              TMath::MaxElement(9,sysNorm),
                              TMath::MaxElement(5,sysPoly)
                            };
    
      double SysError_rho = 0.0;
      for(int i = 0; i < 4; i++)
      {
        double sourcei = TMath::Power((rho_max[i] - rho_min[i])/TMath::Sqrt(12.0),2);
        cout << "rho_min = " << rho_min[i] << ", rho_max = " << rho_max[i] << endl;
        SysError_rho += sourcei;
      }
      SysError_rho = TMath::Sqrt(SysError_rho);

      g_SysErrors[KEY_Default]->SetPoint(i_point,y_def,rho_def);
      g_SysErrors[KEY_Default]->SetPointError(i_point,0.0,0.0,SysError_rho,SysError_rho);

      //double StatError_rho = g_mPar[KEY_Default]->GetErrorYhigh(i_point);
      //g_StatErrors[KEY_Default]->SetPoint(i_point,y_def,rho_def);
      //g_StatErrors[KEY_Default]->SetPointError(i_point,0.0,0.0,StatError_rho,StatError_rho);
    }
  }



  //cout << "sigma from BESI = " << TMath::Abs((rho_def-0.37)/TMath::Sqrt(0.008*0.008 + 0.007*0.007 + rho_defErr*rho_defErr + SysError_rho*SysError_rho)) << endl;
  //cout << "sigma from BESII 2nd order = " << TMath::Abs((rho_def-0.3505)/TMath::Sqrt(0.0024*0.0024 + 0.0025*0.0025 + rho_defErr*rho_defErr + SysError_rho*SysError_rho)) << endl;
  //cout << "sigma from BESII 2nd order (stat only) = " << TMath::Abs((rho_def-0.3505)/TMath::Sqrt(0.0024*0.0024 + rho_defErr*rho_defErr)) << endl;
  //}
  //double val[6];
  //double sysE[6];
  //double statE[6];

  //for(int ipt = 2; ipt < vmsa::pt_rebin_last[energy]; ipt++)
  //{
  //  double pt,rho;
  //  g_StatErrors->GetPoint(ipt,pt,rho);
  //  val[ipt] = rho;
  //  statE[ipt] = g_StatErrors->GetErrorYhigh(ipt);  
  //  sysE[ipt] = g_SysErrors->GetErrorYhigh(ipt); 
  //  
  //  if(random3D)
  //  {
  //    val[ipt] = 4.*(rho-1./3.)+1./3.;
  //    statE[ipt] = 4.*g_StatErrors->GetErrorYhigh(ipt);  
  //    sysE[ipt] = g_SysErrors->GetErrorYhigh(ipt); 
  //    g_StatErrors->SetPoint(ipt,pt,val[ipt]);
  //    g_SysErrors->SetPoint(ipt,pt,val[ipt]);
  //    g_StatErrors->SetPointError(ipt,0.0,0.0,statE[ipt],statE[ipt]); 
  //  }
  //}

  //cout << "double ptbincenter={0.6,0.8,1.5,2.1,2.7,3.6};" << endl;
  //cout << "rho00={"; for (int i = 2; i < 6; i++) {if (i<5) cout << val[i] << ","; else cout << val[i] << "};" << endl;}
  //cout << "stat={"; for (int i = 2; i < 6; i++) {if (i<5) cout << statE[i] << ","; else cout << statE[i] << "};" << endl;}
  //cout << "sys={"; for (int i = 2; i < 6; i++) {if (i<5) cout << sysE[i] << ","; else cout << sysE[i] << "};" << endl;}

  //TFile *bes27_11 = TFile::Open("../data/RawRhoPtSys_27GeVRun11.root");
  //TFile *bes27_18 = TFile::Open("../data/RawRhoPtSys_27GeVRun18.root");
 
  //TGraphAsymmErrors *g27_11 = (TGraphAsymmErrors*) ((TGraphAsymmErrors*)bes27_11->Get("rhoRaw_Centrality_9_2nd_Dca_0_Sig_1_Phi_Norm_0_Sigma_2_Inte"))->Clone();
  //TGraphAsymmErrors *g27_18 = (TGraphAsymmErrors*) ((TGraphAsymmErrors*)bes27_18->Get("rho00_2ndEP_pt_stat_27"))->Clone();

  cout << "Loaded Plots" << endl;
  
  //TFile *besi = TFile::Open("../data/rho00_stat_sys_Laxis.root");
  //TGraphAsymmErrors *besi19;  
  //TGraphAsymmErrors *besi19_sys;
  //if(order == 2 && energy == 4) 
  //{
  //  besi19 = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_stat_27;1");
  //  besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_2ndEP_pt_sys_27;1");
  //}
  //if(order == 1 && energy == 4) 
  //{
  //  besi19 = (TGraphAsymmErrors*)besi->Get("rho00_1stEP_pt_stat_27;1");
  //  besi19_sys = (TGraphAsymmErrors*)besi->Get("rho00_1stEP_pt_sys_27;1");
  //}


  //TGraph *besi_stat;
  //TGraph *besi_sys; 
  //if(order == 1 && energy == 4)
  //{
  //  besi_stat = (TGraph*)besi->Get("rho00_1stEP_energy_stat;1");
  //  besi_sys = (TGraph*)besi->Get("rho00_1stEP_energy_sys;1");
  //}
  //if(order == 2 && energy == 4)
  //{
  //  besi_stat = (TGraph*)besi->Get("rho00_2ndEP_energy_stat;1");
  //  besi_sys = (TGraph*)besi->Get("rho00_2ndEP_energy_sys;1");
  //}


  //for(int i = 0; i < besi_stat->GetN(); i++)
  //{
  //  double en1,val1,err1;
  //  double en2,val2,err2;
  //  besi_stat->GetPoint(i,en1,val1);
  //  err1 = besi_stat->GetErrorY(i);
  //  besi_sys->GetPoint(i,en2,val2);
  //  err2 = besi_sys->GetErrorY(i);
  //  cout << std::setprecision(6);
  //  cout << "energy = " << en1 << ", value = " << val1 << ", stat = " << err1 << endl;
  //  cout << "energy = " << en2 << ", value = " << val2 << ",  sys = " << err2 << endl;
  //  if(en1 == 19.6 || i == 1) 
  //  {
  //    double sigma = TMath::Abs((rhofinal-val1)/TMath::Sqrt(err1*err1 + err2*err2 + rhofinalstat*rhofinalstat + rhofinalsys*rhofinalsys));
  //    double sigma13 = TMath::Abs((1./3.-val1)/TMath::Sqrt(err1*err1 + err2*err2));
  //    cout << "energy = " << en1 << ", value = " << rhofinal << ", stat = " << rhofinalstat << ", sys = " << rhofinalsys << "    BESII" << endl;
  //    cout << "sigma from BES-I = " << sigma << endl;
  //    cout << "BES-I sigma from 1/3 = " << sigma13 << endl;
  //  }
  //}

  //cout << "All good" << endl;

  string partitle[5] = {"#Delta#rho_{00}","Re(#rho_{10}-#rho_{0-1})","Im(#rho_{10}-#rho_{0-1})","Re(#rho_{1-1}","Im(#rho_{1-1})"};

  TCanvas *c_rho_SysError = new TCanvas("c_rho_SysError","c_rho_SysError",600,10,1200,800);
  c_rho_SysError->Divide(3,2);
  for(int i_par = 0; i_par < 5; i_par++)
  {
    string KEY_Default = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",parstring[i_par].c_str(),EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());

    c_rho_SysError->cd(i_par+1);
    c_rho_SysError->cd(i_par+1)->SetLeftMargin(0.15);
    c_rho_SysError->cd(i_par+1)->SetBottomMargin(0.15);
    c_rho_SysError->cd(i_par+1)->SetTicks(1,1);
    c_rho_SysError->cd(i_par+1)->SetGrid(0,0);
    h_frame->GetYaxis()->SetRangeUser(-0.075,0.05);

    h_frame->GetYaxis()->SetTitle(partitle[i_par].c_str());
    h_frame->DrawCopy("pE");
    cout << "Draw a copy" << endl;

    //if(energy == 4){
    //for(int i = 0; i < 4; i++)
    //{
    //  besi19_sys->SetPointEXhigh(i,0.0);
    //  besi19_sys->SetPointEXlow(i,0.0);
    //  besi19->SetPointEXhigh(i,0.0);
    //  besi19->SetPointEXlow(i,0.0);
    //}
    //}
     PlotLine(0.0,80.0,0.0,0.0,1,2,2);
     gStyle->SetEndErrorSize(3);

     //if(energy == 4){
     //  besi19->SetLineColor(kBlack);
     //  besi19->SetLineWidth(1.0);
     //  besi19->SetMarkerStyle(20);
     //  besi19->SetMarkerSize(1.3);
     //  besi19->SetMarkerColor(kBlack);
     //  besi19->SetLineColor(kBlack);

     //  besi19_sys->SetLineColor(kBlack);
     //  besi19_sys->SetMarkerStyle(1);
     //  besi19_sys->SetMarkerSize(10);
     //  besi19_sys->SetMarkerColor(kBlack);
     //  besi19_sys->SetLineColor(kBlack);
     //}

     Draw_TGAE_Point_new_Symbol(0.5,0.46,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2); 
     string leg_count;
     if(correction == "Eff") leg_count = "Efficiency Corrected #phi BES-II";
     if(correction == "Raw") leg_count = "Raw #phi BES-II (|y| < 1.0)";
     if(correction == "AccRes") leg_count = "#phi BES-II (|y| < 1.0)";
     plotTopLegend((char*)leg_count.c_str(),0.65,0.457,0.03,1,0.0,42,0);

     //if(energy == 4){
     //  string leg_besi = "#phi BES-I  (|y| < 1.0)";
     //  Draw_TGAE_Point_new_Symbol(0.5,0.44,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
     //  plotTopLegend((char*)leg_besi.c_str(),0.65,0.437,0.03,1,0.0,42,0);
     //}
     //if(energy == 5){
     //string leg_bes11 = "#phi Run 11 (|y| < 1.0)";
     //Draw_TGAE_Point_new_Symbol(0.5,0.44,0.0,0.0,0.0,0.0,style_phi_1st,color_phi_1st,size_marker-0.2);
     //plotTopLegend((char*)leg_bes11.c_str(),0.65,0.437,0.03,1,0.0,42,0);
     //string leg_bes18 = "#phi Run 18 (|y| < 1.0)";
     //Draw_TGAE_Point_new_Symbol(0.5,0.42,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_1st,size_marker+0.2);
     //plotTopLegend((char*)leg_bes18.c_str(),0.65,0.417,0.03,1,0.0,42,0);
     //}
  //string leg_sp = "STAR Preliminary";
  //plotTopLegend((char*)leg_sp.c_str(),0.3,0.48,0.03,2,0.0,42,0);

    PlotLine(0.25,0.6,0.418,0.418,1,2,2);
    string leg_line = "#rho_{00} = 1/3";
    plotTopLegend((char*)leg_line.c_str(),0.65,0.417,0.03,1,0.0,42,0);

    string leg_energy = Form("Au+Au %s", vmsa::mBeamEnergyText[energy].c_str());
    plotTopLegend((char*)leg_energy.c_str(),1.65,0.265,0.04,1,0.0,42,0);
    plotTopLegend((char*)"20%-60%",2.0,0.245,0.04,1,0.0,42,0);

    g_StatErrors[KEY_Default]->SetMarkerStyle(20);
    g_StatErrors[KEY_Default]->SetMarkerColor(kRed);
    g_StatErrors[KEY_Default]->SetLineColor(kRed);
    g_StatErrors[KEY_Default]->SetMarkerSize(1.3);

    g_SysErrors[KEY_Default]->SetMarkerStyle(20);
    g_SysErrors[KEY_Default]->SetMarkerSize(1.3);
    g_SysErrors[KEY_Default]->SetMarkerColor(kRed);
    g_SysErrors[KEY_Default]->SetLineColor(kRed);
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
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_StatErrors[KEY_Default],style_phi_2nd,color_phi_2nd,colorDiff_phi,size_marker+0.2);
    //if(energy == 4){
    //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)besi19,style_phi_1st,color_phi_1st,colorDiff_phi,size_marker-0.2);
    //}
    //if(energy == 5){
    //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g27_18,style_phi_2nd,color_phi_1st,colorDiff_phi,size_marker+0.2);
    //Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g27_11,style_phi_1st,color_phi_1st,colorDiff_phi,size_marker-0.2);
    //}
    // plotSysErrors(g_rho_2nd_sys[total_pad-1],color_phi_2nd);
    plotSysErrorsBox(g_SysErrors[KEY_Default],color_phi_2nd);
    //if(energy == 4){
    //plotSysErrorsBox(besi19_sys,color_phi_1st);
   // }
  //g_StatErrors->Draw("pE same");
  //plotSysErrorsBox(g_SysErrors,kRed,energy);

  }
  c_rho_SysError->SaveAs(Form("figures/%s/%s/rapiditystudy/%s2D_finaloutputrho_%s_Order%d_%s.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),correction.c_str(),order,etamode.c_str()));  
  c_rho_SysError->SaveAs(Form("figures/%s/%s/rapiditystudy/%s2D_finaloutputrho_%s_Order%d_%s.eps",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),correction.c_str(),order,etamode.c_str()));  
  c_rho_SysError->SaveAs(Form("figures/%s/%s/rapiditystudy/%s2D_finaloutputrho_%s_Order%d_%s.png",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),frame.c_str(),correction.c_str(),order,etamode.c_str()));  

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
  //if(random3D)
  //{
  //  for (int ipt = 0; ipt < 6; ipt++)
  //  {
  //    double pt;
  //    double b2val, b2stat, b2sys;
  //    g_StatErrors->GetPoint(ipt,pt,b2val);
  //    b2stat = g_StatErrors->GetErrorYhigh(ipt);
  //    b2sys = g_SysErrors->GetErrorYhigh(ipt);

  //    double sigma = TMath::Abs((b2val-1./3.)/TMath::Sqrt(b2stat*b2stat + b2sys*b2sys));
  //    std::cout << std::fixed << std::setprecision(4);
  //    cout << "pt = "<<  pt << ":   BESII = " << b2val << " +/- " << b2stat << " (stat) +/- " << b2sys << "(sys)    sigma = " << std::setprecision(2) << sigma << endl; 
  //  }
  //}

  string OutPutFile = Form("../output/AuAu%s/%s/%s2DRapidity_Rho_%sSysErrors_%s_%s.root",frame.c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),spectra.c_str());
  if(order == 1) OutPutFile = Form("../output/AuAu%s/%s/%s2DRapidity_Rho_%sSysErrors_%s_FirstOrder%s.root",frame.c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),etamode.c_str(),spectra.c_str());
  //string OutPutFile = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_%d_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),defaultF,etamode.c_str());
  //OutPutFile = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_%d_2ndUpdated.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),defaultF);
  //string OutPutFile = Form("../output/AuAu%s/%s/Rho_%sSysErrors_F_%d_4th.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str(),defaultF);
  if(random3D) OutPutFile = Form("../output/AuAu%s/%s/3DRandom/Rho_%sSysErrors.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),correction.c_str());
  cout << "OutPutFile set to: " << OutPutFile.c_str() << endl;
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  for(int i_par = 0; i_par < 5; i_par++)
  {
    string KEY_Default = Form("%sRaw_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly1",parstring[i_par].c_str(),EP[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str());
    TString StatErrorRho = Form("g_%s_order%d_%s_%s_StatError",parstring[i_par].c_str(),order,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
    g_StatErrors[KEY_Default]->SetName(StatErrorRho.Data());
    g_StatErrors[KEY_Default]->Write();
    TString SysErrorRho = Form("g_%s_order%d_%s_%s_SysError",parstring[i_par].c_str(),order,vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
    g_SysErrors[KEY_Default]->SetName(SysErrorRho.Data());
    g_SysErrors[KEY_Default]->Write();
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

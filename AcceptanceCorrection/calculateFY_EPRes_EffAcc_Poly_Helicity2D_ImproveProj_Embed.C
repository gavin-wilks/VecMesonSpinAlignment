#include "../Utility/StSpinAlignmentCons.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "../Utility/type.h"
#include "../Utility/functions.h"
#include <string>
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH2.h"
#include "TH1.h"

void calculateFY_EPRes_EffAcc_Poly_Helicity2D_ImproveProj_Embed(const int energy = 4, const int pid = 0, bool doall = true, bool isBesI = false, bool random3D = false, int etamode = 0, int order = 2, int frameopt = 0) {

  const int ncos = 9;
  const int nphi = 12;

  const int y_padding = 2; 

  std::string frame = "Global";
  if(frameopt == 1) frame = "Helicity";

  std::string etastring;
  if(etamode == 0) etastring = "eta1_eta1";
  if(etamode == 3) etastring = "eta0p4";
  if(etamode == 4) etastring = "eta0p6";
  if(etamode == 5) etastring = "eta0p8";

  gStyle->SetOptDate(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(0.4);
  gStyle->SetLegendFont(22);
  gStyle->SetLabelFont(22);
  gStyle->SetTitleFont(22);

  Double_t cent_set[10] = {80,70,60,50,40,30,20,10,5,0};
  Double_t pt_set[8] = {0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 4.2, 5.4};
  Char_t Char[9][10] = {"70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","5-10%","0-5%"};
  Double_t centCent[9] = {75.0,65.0,55.0,45.0,35.0,25.0,15.0,7.5,2.5};

  if(doall)
  {
    double eff_plateau[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total][9][12];
    double eff_plateau_error[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total][9][12];
    
    TFile *eff_file_plateau;
    if(order == 1) eff_file_plateau = new TFile(Form("../Corrections/PhiEmbedding_Corrections_%s_FirstOrder.root",vmsa::mBeamEnergy[energy].c_str()),"READ");
    if(order == 2) eff_file_plateau = new TFile(Form("../Corrections/PhiEmbedding_Corrections_%s.root"           ,vmsa::mBeamEnergy[energy].c_str()),"READ");

    TH2D *eff_hist_plot[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total];
    TH2D        *histmc[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total];
    TH2D        *histrc[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total];
    
    TH1D *eff_hist_plot_cos[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total];
    TH1D        *histmc_cos[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total];
    TH1D        *histrc_cos[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total];

    TH1D *eff_hist_plot_phi[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total];
    TH1D        *histmc_phi[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total];
    TH1D        *histrc_phi[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total];

    TH1D *eff_hist_plot_cos_bin[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total][nphi];
    TH1D        *histmc_cos_bin[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total][nphi];
    TH1D        *histrc_cos_bin[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total][nphi];

    TH1D *eff_hist_plot_phi_bin[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total][ncos];
    TH1D        *histmc_phi_bin[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total][ncos];
    TH1D        *histrc_phi_bin[vmsa::cent_rebin_total_2D][vmsa::pt_rebin_y_2D][vmsa::eta_total][ncos];

    TH2DMap         yield_hist_plot;
    TH1DMap     yield_hist_plot_cos;
    TH1DMap     yield_hist_plot_phi;
    TH1DMap yield_hist_plot_cos_bin;
    TH1DMap yield_hist_plot_phi_bin;

    TH2DMap         rawyield_hist_plot;
    TH1DMap     rawyield_hist_plot_cos;
    TH1DMap     rawyield_hist_plot_phi;
    TH1DMap rawyield_hist_plot_cos_bin;
    TH1DMap rawyield_hist_plot_phi_bin;
    
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int iy = y_padding; iy < vmsa::eta_total - y_padding; iy++)
        {
          string histnamemc = Form("h_MC_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy);
          string histnamerc = Form("h_RC_Cent_%d_Pt_%d_Y_%d",icent,ipt,iy);
          histmc[icent][ipt][iy] = (TH2D*) ((TH2D*) eff_file_plateau->Get(histnamemc.c_str()))->Clone(Form("mc_%d_%d_%d",icent,ipt,iy));    
          histrc[icent][ipt][iy] = (TH2D*) ((TH2D*) eff_file_plateau->Get(histnamerc.c_str()))->Clone(Form("rc_%d_%d_%d",icent,ipt,iy));    

          histmc_cos[icent][ipt][iy] = (TH1D*) histmc[icent][ipt][iy]->ProjectionX(Form("mccos_%d_%d_%d",icent,ipt,iy),0,-1,"e");
          histrc_cos[icent][ipt][iy] = (TH1D*) histrc[icent][ipt][iy]->ProjectionX(Form("rccos_%d_%d_%d",icent,ipt,iy),0,-1,"e");

          histmc_phi[icent][ipt][iy] = (TH1D*) histmc[icent][ipt][iy]->ProjectionY(Form("mcphi_%d_%d_%d",icent,ipt,iy),0,-1,"e");
          histrc_phi[icent][ipt][iy] = (TH1D*) histrc[icent][ipt][iy]->ProjectionY(Form("rcphi_%d_%d_%d",icent,ipt,iy),0,-1,"e");

          eff_hist_plot[icent][ipt][iy] = (TH2D*) histrc[icent][ipt][iy]->Clone(Form("eff_%d_%d_%d",icent,ipt,iy));
          eff_hist_plot[icent][ipt][iy]->Divide(histrc[icent][ipt][iy],histmc[icent][ipt][iy],1,1,"B");
          eff_hist_plot[icent][ipt][iy]->Print();

          eff_hist_plot_cos[icent][ipt][iy] = (TH1D*) histrc_cos[icent][ipt][iy]->Clone(Form("effcos_%d_%d_%d",icent,ipt,iy));
          eff_hist_plot_cos[icent][ipt][iy]->Divide(histrc_cos[icent][ipt][iy],histmc_cos[icent][ipt][iy],1,1,"B");
          eff_hist_plot_cos[icent][ipt][iy]->Print();

          eff_hist_plot_phi[icent][ipt][iy] = (TH1D*) histrc_phi[icent][ipt][iy]->Clone(Form("effphi_%d_%d_%d",icent,ipt,iy));
          eff_hist_plot_phi[icent][ipt][iy]->Divide(histrc_phi[icent][ipt][iy],histmc_phi[icent][ipt][iy],1,1,"B");
          eff_hist_plot_phi[icent][ipt][iy]->Print();

          for(int k = 1; k <= ncos; k++)
          {
            histmc_phi_bin[icent][ipt][iy][k-1] = (TH1D*) histmc[icent][ipt][iy]->ProjectionY(Form("mcphi_%d_%d_%d_%d",icent,ipt,iy,k-1),k,k,"e");
            histrc_phi_bin[icent][ipt][iy][k-1] = (TH1D*) histrc[icent][ipt][iy]->ProjectionY(Form("rcphi_%d_%d_%d_%d",icent,ipt,iy,k-1),k,k,"e");

            eff_hist_plot_phi_bin[icent][ipt][iy][k-1] = (TH1D*) histrc_phi_bin[icent][ipt][iy][k-1]->Clone(Form("effphi_%d_%d_%d_%d",icent,ipt,iy,k-1));
            eff_hist_plot_phi_bin[icent][ipt][iy][k-1]->Divide(histrc_phi_bin[icent][ipt][iy][k-1],histmc_phi_bin[icent][ipt][iy][k-1],1,1,"B");
          }
          for(int l = 1; l <= nphi; l++)
          {
            histmc_cos_bin[icent][ipt][iy][l-1] = (TH1D*) histmc[icent][ipt][iy][j-1]->ProjectionX(Form("mccos_%d_%d_%d_%d",icent,ipt,iy,l-1),l,l,"e");
            histrc_cos_bin[icent][ipt][iy][l-1] = (TH1D*) histrc[icent][ipt][iy][j-1]->ProjectionX(Form("rccos_%d_%d_%d_%d",icent,ipt,iy,l-1),l,l,"e");

            eff_hist_plot_cos_bin[icent][ipt][iy][l-1] = (TH1D*) histrc_cos_bin[icent][ipt][iy][l-1]->Clone(Form("effcos_%d_%d_%d_%d",icent,ipt,iy,l-1));
            eff_hist_plot_cos_bin[icent][ipt][iy][l-1]->Divide(histrc_cos_bin[icent][ipt][iy][l-1],histmc_cos_bin[icent][ipt][iy][l-1],1,1,"B");
          }

          //for(int k = 1; k <= ncos; k++)
          //{
          //  for(int l = 1; l <= nphi; l++)
          //  {
          //    eff_plateau[icent][ipt][iy][k-1][l-1] = eff_hist_plot[icent][ipt][iy]->GetBinContent(k,l);
          //    eff_plateau_error[icent][ipt][iy][k-1][l-1] = eff_hist_plot[icent][ipt][iy]->GetBinError(k,l)/eff_plateau[icent][ipt][iy][k-1][l-1];
          //    cout << "efficiency = " << eff_plateau[icent][ipt][iy][k-1][l-1] << " +/- " << eff_plateau_error[icent][ipt][iy][k-1][l-1] << endl;
          //  }
          //}
        }
      }
    }

    TFile *input;
    if(order == 1) input = new TFile(Form("rho00/%s/%s/Poly/PhiPsi_Raw%sPtSys_%s_PolySys_FirstOrder.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str()),"READ");
    if(order == 2) input = new TFile(Form("rho00/%s/%s/Poly/RawPhiPtSys_%s_PolySys_Helicity_2D_OffDiag_Rapidity.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),etastring.c_str()),"READ");

    TString *Title;

    TGraMap g_rho00;
    TGraMap g_real;
    TGraMap g_imag;
    TGraMap g_rerho;
    TGraMap g_imrho;
    
    //cout << "Before the loop" << endl;

    vecFMap rho00;
    vecFMap real;
    vecFMap imag;
    vecFMap rerho;
    vecFMap imrho;

    TF2Map Func_obs;
    TH1DMap Func_obs_cos;
    TH1DMap Func_obs_phi;
    TH1DMap Func_obs_cos_bin;
    TH1DMap Func_obs_phi_bin;

    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; ++icent) // Centrality loop
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      { 
        for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
        {
          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
          {
            if( i_dca != 0 && i_sig != 0 ) continue;
            for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
            {
              for(int i_sigma = vmsa::Sig_start; i_sigma < 1/*vmsa::Sig_stop*/; ++i_sigma)
              {
                for(int i_method = 1/*vmsa::Method_start*/; i_method < vmsa::Method_stop; ++i_method)
                {
                  for(int i_poly = 0; i_poly < 1/*3*/; i_poly++)
                  {
                    for(Int_t iy = y_padding; iy < vmsa::eta_total - y_padding; iy++) 
                    {
                      string key = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      string key_clone = Form("clone_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);

                      rawyield_hist_plot[key] = (TH2F*)((TH2F*)input->Get(key.c_str()))->Clone(key_clone.c_str());
                      rawyield_hist_plot[key]->Print();
                       
                      string key_corr = Form("corrected_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      yield_hist_plot[key] = (TH2D*) rawyield_hist_plot[key]->Clone(key_corr.c_str());
                      //yield_hist_plot[key] = new TH2D(key_corr.c_str(), key_corr.c_str(), ncos, -1, 1, nphi, 0.0, TMath::Pi()*2.0);
                      yield_hist_plot[key]->Divide(eff_hist_plot[icent][ipt][iy]);
                     
                      //for(int i=1; i<=9; i++) 
                      //{
                      //  for(int j=1; j<=12; j++) 
                      //  {
                      //    float inte_mean = rawyield_hist_plot[PtBin-1]->GetBinContent(i,j);
                      //    float inte_mean_error = rawyield_hist_plot[PtBin-1]->GetBinError(i,j);
                      //    //cout << "Inte = " << inte_mean << " +/- " << inte_mean_error << endl;

                      //    double value = inte_mean/eff_plateau[i_cent][PtBin-1][i-1][j-1]; 
                      //    double value_err = inte_mean/eff_plateau[i_cent][PtBin-1][i-1][j-1]*TMath::Sqrt(inte_mean_error*inte_mean_error/inte_mean/inte_mean+eff_plateau_error[i_cent][PtBin-1][i-1][j-1]*eff_plateau_error[i_cent][PtBin-1][i-1][j-1]);
                      //    //cout << "corrected yield = " << value << " +/- " << value_err << endl;

                      //    yield_hist_plot[PtBin-1]->SetBinContent(i, j, value);
                      //    yield_hist_plot[PtBin-1]->SetBinError(i, j, value_err);
                     
                      //    //cout << "EfficiencyCalculation" << endl;
                      //  }
                      //}
                      //cout << "Right before 1D projections" << endl;

                      rawyield_hist_plot[key]->Print();
                      string key_cos = Form("cos_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      string key_phi = Form("phi_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      string key_cos_corr = Form("corrected_cos_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      string key_phi_corr = Form("corrected_phi_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      rawyield_hist_plot_cos[key] = (TH1F*) rawyield_hist_plot[key]->ProjectionX(key_cos,1,nphi,"e");
                      rawyield_hist_plot_phi[key] = (TH1F*) rawyield_hist_plot[key]->ProjectionY(key_phi,1,ncos,"e");
                      yield_hist_plot_cos[key] = (TH1D*) yield_hist_plot[key]->ProjectionX(key_cos_corr,1,nphi,"e");
                      yield_hist_plot_phi[key] = (TH1D*) yield_hist_plot[key]->ProjectionY(key_phi_corr,1,ncos,"e");

                      for(int iphi = 0; iphi < 12; iphi++)
                      {
                        string key_cos_phi = Form("cos_phi_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",iphi,ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                        string key_cos_phi_corr = Form("corrected_cos_phi_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",iphi,ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                        yield_hist_plot_cos_bin[key_cos_phi_corr] = (TH1D*) yield_hist_plot[key]->ProjectionX(key_cos_phi_corr.c_str(),iphi+1,iphi+1,"e");
                        rawyield_hist_plot_cos_bin[key_cos_phi] = (TH1F*) rawyield_hist_plot[key]->ProjectionX(key_cos_phi.c_str(),iphi+1,iphi+1,"e");
                      }
                      for(int icos = 0; icos < 9; icos++)
                      {
                        string key_phi_cos = Form("phi_cos_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",icos,ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                        string key_phi_cos_corr = Form("corrected_phi_cos_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",icos,ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                        yield_hist_plot_phi_bin[key_phi_cos_corr] = (TH1D*) yield_hist_plot[key]->ProjectionY(key_phi_cos_corr.c_str(),icos+1,icos+1,"e");
                        rawyield_hist_plot_phi_bin[key_phi_cos] = (TH1F*) rawyield_hist_plot[key]->ProjectionY(key_phi_cos.c_str(),icos+1,icos+1,"e");
                      }
                      //cout << "After 1D projections" << endl;
                      
                      string key_func = Form("func_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      Func_obs[key] = new TF2(key_func.c_str(),SpinDensity2Dcos,-1.0,1.0,0.0,2.0*TMath::Pi(),6);
                    
                      //cout << "Before fit" << endl;
                      Func_obs[key]->SetParameter(0,1./3.);
                      Func_obs[key]->SetParameter(1,0.0);
                      Func_obs[key]->SetParameter(2,0.0);
                      Func_obs[key]->SetParameter(3,0.0);
                      Func_obs[key]->SetParameter(4,0.0);
                      Func_obs[key]->SetParameter(5,yield_hist_plot[key]->GetMaximum());
                      //yield_hist_plot[key]->Print();
                      yield_hist_plot[key]->Fit(Func_obs[key],"NMRI"); // fit corrected distribution for observerved rho00
                      
                      // 1D projection in cos(theta*) for 2D fit
                      int nFineBinsCos = ncos;
                      int nFineBinsPhi = nphi;
                      //for(int iphi = 0; iphi < 12; iphi++)
                      //{
                      //  Func_obs_cos[iphi][PtBin-1] = new TF1(Form("Func_obs_cos%d_%d",PtBin-1,iphi),SpinDensity1Dphiset,-1.0,1.0,7);
                      //}
                      for(int icos = 0; icos < 9; icos++)
                      {
                        string key_func_phi = Form("func_phi_cos_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",icos,ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                        Func_obs_phi_bin[key_func_phi] = new TH1D(key_func_phi.c_str(),key_func_phi.c_str(),nFineBinsPhi,0.0,2.0*TMath::Pi());
                      }
                      for(int iphi = 0; iphi < 12; iphi++)
                      {
                        string key_func_cos = Form("func_cos_phi_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",iphi,ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                        Func_obs_cos_bin[key_func_cos] = new TH1D(key_func_cos.c_str(),key_func_cos.c_str(),nFineBinsCos,-1.0,1.0);
                      }
                      for(int icos = 1; icos <= nFineBinsCos; ++icos)
                      {
                        for(int iphi = 1; iphi <= nFineBinsPhi; ++iphi)
                        {
                          string key_func_cos = Form("func_cos_phi_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",iphi-1,ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                          string key_func_phi = Form("func_phi_cos_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",icos-1,ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);

                          double phi      = Func_obs_phi_bin[key_func_phi]->GetBinCenter(iphi-1);
                          double phiwidth = Func_obs_phi_bin[key_func_phi]->GetBinWidth(iphi-1);
                          double cos      = Func_obs_cos_bin[key_func_cos]->GetBinCenter(icos-1);
                          double coswidth = Func_obs_cos_bin[key_func_cos]->GetBinWidth(icos-1);

                          double integral = Func_obs[key]->Integral(cos - 0.5 * coswidth, 
                                                                    cos + 0.5 * coswidth,
                                                                    phi - 0.5 * phiwidth,
                                                                    phi + 0.5 * phiwidth);
                          integral /= coswidth;
                          integral /= phiwidth;

                          Func_obs_cos_bin[key_func_cos]->SetBinContent(icos, integral);
                          Func_obs_phi_bin[key_func_phi]->SetBinContent(iphi, integral);
                        }
                      }

                      string key_func_cos = Form("func_cos_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      string key_func_phi = Form("func_phi_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);

                      Func_obs_cos[key_func_cos] = new TH1D(key_func_cos.c_str(),key_func_cos.c_str(),nFineBinsCos,-1.0,1.0);
                      Func_obs_phi[key_func_phi] = new TH1D(key_func_phi.c_str(),key_func_phi.c_str(),nFineBinsPhi,0.0,2.0*TMath::Pi());
                      for(int icos = 1; icos <= nFineBinsCos; ++icos)
                      {
                        double cos = Func_obs_cos[key_func_cos]->GetBinCenter(icos);
                        double coswidth = Func_obs_cos[key_func_cos]->GetBinWidth(icos);
                        double totalvalue = 0.0;
                        for(int iphi = 1; iphi <= nFineBinsPhi; ++iphi)
                        {
                          double phi = Func_obs_phi[key_func_phi]->GetBinCenter(j);
                          double phiwidth = Func_obs_phi[key_func_phi]->GetBinWidth(j);
                          double cosintegral = Func_obs[key]->Integral(cos - 0.5 * coswidth, 
                                                                       cos + 0.5 * coswidth,
                                                                       phi - 0.5 * phiwidth,
                                                                       phi + 0.5 * phiwidth);
                          cosintegral /= coswidth;
                          cosintegral /= phiwidth;
                          totalvalue += cosintegral;
                        //cosintegral *= double(nFineBinsCos)/9;
                        }
                        Func_obs_cos[key_func_cos]->SetBinContent(icos, totalvalue);
                        //cout << "bin = " << icos << ", cos = " << cos << ", integral = " << totalvalue << endl;
                        //cout << "bin = " << icos << ", cos - 0.5 * coswidth = " << cos - 0.5 * coswidth << endl;
                        //cout << "bin = " << icos << ", cos + 0.5 * coswidth = " << cos + 0.5 * coswidth << endl;
                           

                        //cout << "Before scaling bin = " << i << ", cos = " << cos << ", integral = " << cosintegral << endl;
                      }
                      for(int iphi = 1; iphi <= nFineBinsPhi; ++iphi)
                      {
                        double phi = Func_obs_phi[key_func_phi]->GetBinCenter(iphi);
                        double phiwidth = Func_obs_phi[key_func_phi]->GetBinWidth(iphi);
                        double totalvalue = 0.0;
                        for(int icos = 1; icos <= nFineBinsCos; ++icos)
                        {
                          double cos = Func_obs_cos[key_func_cos]->GetBinCenter(icos);
                          double coswidth = Func_obs_cos[key_func_cos]->GetBinWidth(icos);
                          double phiintegral = Func_obs[key]->Integral(cos - 0.5 * coswidth,
                                                                       cos + 0.5 * coswidth, 
                                                                       phi - 0.5 * phiwidth,
                                                                       phi + 0.5 * phiwidth);
                          //cout << "Before scaling bin = " << i << ", phi = " << phi << ", integral = " << phiintegral << endl;
                          phiintegral /= phiwidth;
                          phiintegral /= coswidth;
                          totalvalue += phiintegral;
                        }
                        //phiintegral *= double(nFineBinsPhi)/12;
                        Func_obs_phi[key_func_phi]->SetBinContent(iphi, totalvalue);
                        //cout << "bin = " << iphi << ", phi = " << phi << ", integral = " << totalvalue << endl;
                        //cout << "bin = " << iphi << ", phi - 0.5 * phiwidth = " << phi - 0.5 * phiwidth << endl;
                        //cout << "bin = " << iphi << ", phi + 0.5 * phiwidth = " << phi + 0.5 * phiwidth << endl;
                      }


                      // 1D projection in cos(theta*) for 2D fit

                      //cout << "After fit" << endl;
                    
                      rho00[key].push_back(Func_obs[key]->GetParameter(0));
                       real[key].push_back(Func_obs[key]->GetParameter(1));
                       imag[key].push_back(Func_obs[key]->GetParameter(2));
                      rerho[key].push_back(Func_obs[key]->GetParameter(3));
                      imrho[key].push_back(Func_obs[key]->GetParameter(4));
      
                      rho00[key].push_back(Func_obs[key]->GetParError(0));
                       real[key].push_back(Func_obs[key]->GetParError(1));
                       imag[key].push_back(Func_obs[key]->GetParError(2));
                      rerho[key].push_back(Func_obs[key]->GetParError(3));
                      imrho[key].push_back(Func_obs[key]->GetParError(4));
                      
                      //cout << "rho00 = " << rho00[PtBin-1] << " +/- " << rho00err[PtBin-1] << endl;
                      //cout << "real = " << real[PtBin-1] << " +/- " << realerr[PtBin-1] << endl;
                      //cout << "imag = " << imag[PtBin-1] << " +/- " << imagerr[PtBin-1] << endl;
                      //cout << "rerho1n1 = " << rerho[PtBin-1] << " +/- " << rerhoerr[PtBin-1] << endl;
                      //cout << "imrho1n1 = " << imrho[PtBin-1] << " +/- " << imrhoerr[PtBin-1] << endl;

                      //Func_obs->SetLineColor(kRed);      
                      //Func_obs->DrawCopy("same");
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    TCanvas *c = new TCanvas("c","c",10,10,1600,1600);
    c->Divide(4,4);
    for(int i = 0; i < 16; i++)
    {   
      c->cd(i+1);
      c->cd(i+1)->SetLeftMargin(0.175);
      c->cd(i+1)->SetRightMargin(0.175);
      c->cd(i+1)->SetBottomMargin(0.175);
      c->cd(i+1)->SetTicks(1,1);
      c->cd(i+1)->SetGrid(0,0);
    }

    string outputpdf = Form("figures/%s/%s/RapidityEmbed2DHelicityEfficiency.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    string outputstartpdf = Form("%s[",outputpdf.c_str());
    string outputstoppdf = Form("%s]",outputpdf.c_str());

    c->Print(outputstartpdf.c_str());
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int i = 0; i < vmsa::eta_total; i++)
        {   
          c->cd(i+1);     

          eff_hist_plot[icent][ipt][i]->SetTitle(Form("Eff., cent %d-%d%%, %1.1f<p_{T}<%1.1f  %1.1f<y<%1.1f",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1]));
          eff_hist_plot[icent][ipt][i]->GetXaxis()->SetTitle("cos(#theta*)");
          eff_hist_plot[icent][ipt][i]->GetYaxis()->SetTitle("#phi'");
          eff_hist_plot[icent][ipt][i]->Draw("Colz");
        }
        c->Update();
        c->Print(outputpdf.c_str());
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   2D Efficiency ///////////////////////////
    string outputpdf = Form("figures/%s/%s/RapidityEmbed2DHelicityEfficiencyCorrectedYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    string outputstartpdf = Form("%s[",outputpdf.c_str());
    string outputstoppdf = Form("%s]",outputpdf.c_str());

    c->Print(outputstartpdf.c_str());
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int i = 0; i < vmsa::eta_total; i++)
        {   
          c->cd(i+1);
          
          string key = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
            
          yield_hist_plot[key]->SetTitle(Form("Corr. Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f  %1.1f<y<%1.1f",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1]));
          yield_hist_plot[key]->GetXaxis()->SetTitle("cos(#theta*)");
          yield_hist_plot[key]->GetYaxis()->SetTitle("#phi'");
          yield_hist_plot[key]->Draw("Colz");
          //Func_obs[i]->Draw("same");
        }
        c->Update();
        c->Print(outputpdf.c_str());
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   2D Efficiency ///////////////////////////
    //////////////////////   2D Corrected Yields ///////////////////////////
    string outputpdf = Form("figures/%s/%s/RapidityEmbed2DHelicity2DFunc.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    string outputstartpdf = Form("%s[",outputpdf.c_str());
    string outputstoppdf = Form("%s]",outputpdf.c_str());

    c->Print(outputstartpdf.c_str());
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int i = 0; i < vmsa::eta_total; i++)
        {   
          c->cd(i+1);

          string key = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
            
          yield_hist_plot[key]->SetTitle(Form("Corr. Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f  %1.1f<y<%1.1f",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1]));
          yield_hist_plot[key]->GetXaxis()->SetTitle("cos(#theta*)");
          yield_hist_plot[key]->GetYaxis()->SetTitle("#phi'");
          yield_hist_plot[key]->Draw("Colz");
          Func_obs[key]->Draw("Colz");
        }
        c->Update();
        c->Print(outputpdf.c_str());
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   2D Corrected Yields ///////////////////////////
    //////////////////////   2D Raw Yields ///////////////////////////
    string outputpdf = Form("figures/%s/%s/RapidityEmbed2DHelicityRawYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    string outputstartpdf = Form("%s[",outputpdf.c_str());
    string outputstoppdf = Form("%s]",outputpdf.c_str());

    c->Print(outputstartpdf.c_str());
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int i = 0; i < vmsa::eta_total; i++)
        {   
          c->cd(i+1);

          string key = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
        
          rawyield_hist_plot[key]->Print();
          rawyield_hist_plot[key]->SetTitle(Form("Raw Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f  %1.1f<y<%1.1f",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1]));
          rawyield_hist_plot[key]->GetXaxis()->SetTitle("cos(#theta*)");
          rawyield_hist_plot[key]->GetYaxis()->SetTitle("#phi'");
          rawyield_hist_plot[key]->Draw("Colz");
        }
        c->Update();
        c->Print(outputpdf.c_str());
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   2D Raw Yields ///////////////////////////
    //////////////////////   1D Cos Efficiency ///////////////////////////
    string outputpdf = Form("figures/%s/%s/RapidityEmbed1DCosHelicityEfficiency.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    string outputstartpdf = Form("%s[",outputpdf.c_str());
    string outputstoppdf = Form("%s]",outputpdf.c_str());

    c->Print(outputstartpdf.c_str());
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int i = 0; i < vmsa::eta_total; i++)
        {   
          c->cd(i+1);           

          //string key = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);

          eff_hist_plot_cos[icent][ipt][i]->SetTitle(Form("Eff., cent %d-%d%%, %1.1f<p_{T}<%1.1f  %1.1f<y<%1.1f",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1]));
          eff_hist_plot_cos[icent][ipt][i]->GetXaxis()->SetTitle("cos(#theta*)");
          eff_hist_plot_cos[icent][ipt][i]->GetYaxis()->SetTitle("Efficiency");
          eff_hist_plot_cos[icent][ipt][i]->Draw("pE");
        }
        c->Update();
        c->Print(outputpdf.c_str());
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Cos Efficiency ///////////////////////////
    //////////////////////   1D Cos Corrected Yields Small Bins ///////////////////////////
  
    outputpdf = Form("figures/%s/%s/RapidityEmbed1DCosHelicityEfficiencyCorrectedYields_SmallBins.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());
 
    c->Print(outputstartpdf.c_str());
    
    int color[3] = {kRed, kBlue, kGreen}; 
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int iphi = 0; iphi < 12; iphi++)
        {
          for(int i = 0; i < vmsa::eta_total; i++)
          {   
            c->cd(i+1);
       
            string key = Form("cos_phi_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",iphi,ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
      
            yield_hist_plot_cos_bin[key]->SetTitle(Form("Corr. Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f, %d#pi/6<#phi'<%d#pi/6",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1],iphi,iphi+1));
            yield_hist_plot_cos_bin[key]->GetXaxis()->SetTitle("cos(#theta*)");
            yield_hist_plot_cos_bin[key]->GetYaxis()->SetTitle("count");
            yield_hist_plot_cos_bin[key]->Draw("pE");

            string key_func = Form("func_cos_phi_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",iphi,ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
            
            Func_obs_cos_bin[key_func]->SetLineColor(color[0]);
            Func_obs_cos_bin[key_func]->Draw("C same");
          }
          c->Update();
          c->Print(outputpdf.c_str());
        }
      }    
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Cos Corrected Yields Small Bins ///////////////////////////
    //////////////////////   1D Cos Raw Yields Small Bins ///////////////////////////

    outputpdf = Form("figures/%s/%s/RapidityEmbed1DCosHelicityRawYields_SmallBins.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());
 
    c->Print(outputstartpdf.c_str());
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int iphi = 0; iphi < 12; iphi++)
        {
          for(int i = 0; i < vmsa::eta_total; i++)
          {   
            c->cd(i+1);
              
            string key = Form("cos_phi_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",iphi,ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
          
            rawyield_hist_plot_cos_bin[key]->SetTitle(Form("RawYields, cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f, %d#pi/6<#phi'<%d#pi/6",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1],iphi,iphi+1));
            rawyield_hist_plot_cos_bin[key]->GetXaxis()->SetTitle("cos(#theta*)");
            rawyield_hist_plot_cos_bin[key]->GetYaxis()->SetTitle("count");
            rawyield_hist_plot_cos_bin[key]->Draw("pE");
          }
          c->Update();
          c->Print(outputpdf.c_str());
        }
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Cos Raw Yields Small Bins ///////////////////////////
    //////////////////////   1D Cos Raw Yields Small Bins ///////////////////////////

    outputpdf = Form("figures/%s/%s/RapidityEmbed1DCosHelicityEfficiency_SmallBins.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());
 
    c->Print(outputstartpdf.c_str());
    
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int iphi = 0; iphi < 12; iphi++)
        {
          for(int i = 0; i < vmsa::eta_total; i++)
          {   
            c->cd(i+1);
              
            eff_hist_plot_cos_bin[icent][ipt][i][iphi]->SetTitle(Form("Eff., cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f, %d#pi/6<#phi'<%d#pi/6",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1],iphi,iphi+1));
            eff_hist_plot_cos_bin[icent][ipt][i][iphi]->GetXaxis()->SetTitle("cos(#theta*)");
            eff_hist_plot_cos_bin[icent][ipt][i][iphi]->GetYaxis()->SetTitle("count");
            eff_hist_plot_cos_bin[icent][ipt][i][iphi]->Draw("pE");
          }
          c->Update();
          c->Print(outputpdf.c_str());
        }
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Cos Efficiency Small Bins ///////////////////////////
    //////////////////////   1D Cos Corrected Yields ///////////////////////////

    outputpdf = Form("figures/%s/%s/RapidityEmbed1DCosHelicityEfficiencyCorrectedYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());
 
    c->Print(outputstartpdf.c_str());
    
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int i = 0; i < vmsa::eta_total; i++)
        {   
          c->cd(i+1);
            
          string key = Form("cos_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);

          yield_hist_plot_cos[key]->SetTitle(Form("Corr. Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1]));
          yield_hist_plot_cos[key]->GetXaxis()->SetTitle("cos(#theta*)");
          yield_hist_plot_cos[key]->GetYaxis()->SetTitle("count");
          yield_hist_plot_cos[key]->Draw("pE");

          string key_func = Form("func_cos_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
          Func_obs_cos[key_func]->SetLineColor(color[0]);
          Func_obs_cos[key_func]->Draw("C same");
        }
        c->Update();
        c->Print(outputpdf.c_str());
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Cos Corrected Yields ///////////////////////////
    //////////////////////   1D Cos Raw Yields ///////////////////////////
    outputpdf = Form("figures/%s/%s/RapidityEmbed1DCosHelicityRawYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());
 
    c->Print(outputstartpdf.c_str());
    
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int i = 0; i < vmsa::eta_total; i++)
        {   
          c->cd(i+1);

          string key = Form("cos_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
            
          rawyield_hist_plot_cos[key]->SetTitle(Form("Raw Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1]));
          rawyield_hist_plot_cos[key]->GetXaxis()->SetTitle("cos(#theta*)");
          rawyield_hist_plot_cos[key]->GetYaxis()->SetTitle("count");
          rawyield_hist_plot_cos[key]->Draw("pE");
        }
        c->Update();
        c->Print(outputpdf.c_str());
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Cos Raw Yields ///////////////////////////
    //////////////////////   1D Phi Efficiency  ///////////////////////////
    outputpdf = Form("figures/%s/%s/RapidityEmbed1DPhiHelicityEfficiency.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());
 
    c->Print(outputstartpdf.c_str());
    
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int i = 0; i < vmsa::eta_total; i++)
        {   
          c->cd(i+1);           
          eff_hist_plot_phi[icent][ipt][i]->SetTitle(Form("Eff., cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1]));
          eff_hist_plot_phi[icent][ipt][i]->GetXaxis()->SetTitle("phi'");
          eff_hist_plot_phi[icent][ipt][i]->GetYaxis()->SetTitle("Efficiency");
          eff_hist_plot_phi[icent][ipt][i]->Draw("pE");
        }
        c->Update();
        c->Print(outputpdf.c_str());
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Phi Efficiency  ///////////////////////////
    //////////////////////   1D Phi Raw Yields Small Bins  ///////////////////////////

    outputpdf = Form("figures/%s/%s/RapidityEmbed1DPhiHelicityRawYields_SmallBins.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());
 
    c->Print(outputstartpdf.c_str());
    
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int icos = 0; icos < 9; icos++)
        {
          for(int i = 0; i < vmsa::eta_total; i++)
          {   
            c->cd(i+1);

            string key = Form("phi_cos_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",icos,ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
              
            rawyield_hist_plot_phi_bin[key]->SetTitle(Form("Raw Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f, %1.1f/4.5<cos(#theta*)<%1.1f/4.5",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1],-4.5+icos,-3.5+icos));
            rawyield_hist_plot_phi_bin[key]->GetXaxis()->SetTitle("cos(#theta*)");
            rawyield_hist_plot_phi_bin[key]->GetYaxis()->SetTitle("count");
            rawyield_hist_plot_phi_bin[key]->Draw("pE");
          }
          c->Update();
          c->Print(outputpdf.c_str());
        }
      }
    }  
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Phi Raw Yields Small Bins  ///////////////////////////
    //////////////////////   1D Phi Efficiency Small Bins  ///////////////////////////

    outputpdf = Form("figures/%s/%s/RapidityEmbed1DPhiHelicityEfficiency_SmallBins.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());
 
    c->Print(outputstartpdf.c_str());
    
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int icos = 0; icos< 9; icos++)
        {
          for(int i = 0; i < vmsa::eta_total; i++)
          {   
            c->cd(i+1);
              
            eff_hist_plot_phi_bin[icent][ipt][i][icos]->SetTitle(Form("Eff., cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f, %1.1f/4.5<cos(#theta*)<%1.1f/4.5",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1],-4.5+icos,-3.5+icos));
            eff_hist_plot_phi_bin[icent][ipt][i][icos]->GetXaxis()->SetTitle("cos(#theta*)");
            eff_hist_plot_phi_bin[icent][ipt][i][icos]->GetYaxis()->SetTitle("count");
            eff_hist_plot_phi_bin[icent][ipt][i][icos]->Draw("pE");
          }
          c->Update();
          c->Print(outputpdf.c_str());
        }
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Phi Efficiency Small Bins  ///////////////////////////
    //////////////////////   1D Phi Corrected Yields Small Bins  ///////////////////////////


    outputpdf = Form("figures/%s/%s/RapidityEmbed1DPhiHelicityEfficiencyCorrectedYields_SmallBins.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());

    c->Print(outputstartpdf.c_str());
  
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int icos = 0; icos < 9; icos++)
        {
          for(int i = 0; i < vmsa::eta_total; i++)
          {   
            c->cd(i+1);
              
            string key = Form("phi_cos_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",icos,ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
            yield_hist_plot_phi_bin[key]->SetTitle(Form("Corr. Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f, %1.1f/4.5<cos(#theta*)<%1.1f/4.5",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1],-4.5+icos,-3.5+icos));
            yield_hist_plot_phi_bin[key]->GetXaxis()->SetTitle("phi'");
            yield_hist_plot_phi_bin[key]->GetYaxis()->SetTitle("count");
            yield_hist_plot_phi_bin[key]->Draw("pE");

            string key_func = Form("func_phi_cos_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",icos,ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
            Func_obs_phi_bin[key_func]->SetLineColor(color[0]);
            Func_obs_phi_bin[key_func]->Draw("C same");
          }
          c->Update();
          c->Print(outputpdf.c_str());
        }
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Phi Corrected Yields Small Bins  ///////////////////////////
    //////////////////////   1D Phi Corrected Yields  ///////////////////////////

    outputpdf = Form("figures/%s/%s/RapidityEmbed1DPhiHelicityEfficiencyCorrectedYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());
 
    c->Print(outputstartpdf.c_str());
    
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int i = 0; i < vmsa::eta_total; i++)
        {   
          c->cd(i+1);
            
          string key = Form("phi_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
          yield_hist_plot_phi[key]->SetTitle(Form("Corr. Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1]));
          yield_hist_plot_phi[key]->GetXaxis()->SetTitle("phi'");
          yield_hist_plot_phi[key]->GetYaxis()->SetTitle("count");
          yield_hist_plot_phi[key]->Draw("pE");

          string key = Form("func_phi_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
          Func_obs_phi[key_func]->SetLineColor(color[0]);
          Func_obs_phi[key_func]->Draw("C same");
        }  
        c->Update();
        c->Print(outputpdf.c_str());
      }    
    }
    c->Print(outputstoppdf.c_str());
    
    //////////////////////   1D Phi Corrected Yields  ///////////////////////////
    //////////////////////   1D Phi Raw Yields  ///////////////////////////
    
      
    outputpdf = Form("figures/%s/%s/RapidityEmbed1DPhiHelicityRawYields.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());
 
    c->Print(outputstartpdf.c_str());
    
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int i = 0; i < vmsa::eta_total; i++)
        {   
          c->cd(i+1);
            
          string key = Form("phi_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
          rawyield_hist_plot_phi[key]->SetTitle(Form("Raw Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1]));
          rawyield_hist_plot_phi[key]->GetXaxis()->SetTitle("phi'");
          rawyield_hist_plot_phi[key]->GetYaxis()->SetTitle("count");
          rawyield_hist_plot_phi[key]->Draw("pE");
        }
        c->Update();
        c->Print(outputpdf.c_str());
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Phi Raw Yields  ///////////////////////////
    //////////////////////   1D Phi Raw Yields Small Bins ///////////////////////////

    outputpdf = Form("figures/%s/%s/RapidityEmbed1DPhiHelicitRawYields_SmallBins.pdf",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str());
    outputstartpdf = Form("%s[",outputpdf.c_str());
    outputstoppdf = Form("%s]",outputpdf.c_str());

    c->Print(outputstartpdf.c_str());
  
    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; icent++)
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      {
        for(int icos = 0; icos < 9; icos++)
        {
          for(int i = 0; i < vmsa::eta_total; i++)
          {   
            c->cd(i+1);
              
            string key = Form("phi_cos_%d_pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",icos,ipt,i,icent,ordertext[order-1].c_str(),0,0,vmsa::mPID[pid].c_str(),0,0,vmsa::mInteMethod[1].c_str(),1);
            rawyield_hist_plot_phi_bin[key]->SetTitle(Form("Raw Yields, cent %d-%d%%, %1.1f<p_{T}<%1.1f %1.1f<y<%1.1f, %1.1f/4.5<cos(#theta*)<%1.1f/4.5",vmsa::centVal[vmsa::cent_rebin_2D[icent+1]],vmsa::centVal[vmsa::cent_rebin_2D[icent]],vmsa::pt_low_y_2D[energy][ipt],vmsa::pt_up_y_2D[energy][ipt],vmsa::eta_bin[i],vmsa::eta_bin[i+1],-4.5+icos,-3.5+icos));
            rawyield_hist_plot_phi_bin[key]->GetXaxis()->SetTitle("phi'");
            rawyield_hist_plot_phi_bin[key]->GetYaxis()->SetTitle("count");
            rawyield_hist_plot_phi_bin[key]->Draw("pE");
          }
          c->Update();
          c->Print(outputpdf.c_str());
        }
      }
    }
    c->Print(outputstoppdf.c_str());
    //////////////////////   1D Phi Raw Yields Small Bins ///////////////////////////

    string outputname;
    if(order == 1) outputname = Form("output/%s/%s/PhiPsi_AccRes%sPtSys_%s_Poly_FirstOrder%s%s%s.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str(),spectra.c_str(),eptext.c_str(),v2text.c_str());
    if(order == 2) outputname = Form("../output/AuAu%s/%s/AccPhiPtSys_%s_PolySys_Helicity_2D_OffDiag_Embed_Rapidity.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str(),etastring.c_str());
    if(isBesI) outputname = Form("output/%s/%s/BESI/AccRes%sPtSys.root",vmsa::mPID[pid].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::mPID[pid].c_str());
    TFile *output = new TFile(outputname.c_str(),"RECREATE");
    output->cd();

    for(int icent = 0; icent < vmsa::cent_rebin_total_2D; ++icent) // Centrality loop
    {
      for(int ipt = 0; ipt < vmsa::pt_rebin_y_2D; ipt++)
      { 
        for(Int_t i_dca = vmsa::Dca_start; i_dca < 1/*vmsa::Dca_stop*/; i_dca++)
        {
          for(Int_t i_sig = vmsa::nSigKaon_start; i_sig < 1/*vmsa::nSigKaon_stop*/; i_sig++)
          {
            if( i_dca != 0 && i_sig != 0 ) continue;
            for(int i_norm = vmsa::Norm_start; i_norm < 1/*vmsa::Norm_stop*/; ++i_norm)
            {
              for(int i_sigma = vmsa::Sig_start; i_sigma < 1/*vmsa::Sig_stop*/; ++i_sigma)
              {
                for(int i_method = 1/*vmsa::Method_start*/; i_method < vmsa::Method_stop; ++i_method)
                {
                  for(int i_poly = 0; i_poly < 1/*3*/; i_poly++)
                  {
                    Title = new TString(Form("rhoRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             
                    g_rho00 = new TGraphAsymmErrors();
                    g_rho00->SetName(Title->Data());
                    for(Int_t iy = y_padding; iy < vmsa::eta_total - y_padding; iy++) 
                    {
                      string key = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      double mean = (vmsa::eta_bin[iy]+vmsa::eta_bin[iy+1])/2.0; //(ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                      g_rho00->SetPoint(iy,mean,rho00[key][0]);
                      g_rho00->SetPointError(iy,0.0,0.0,rho00[key][1],rho00[key][1]);
                    }
                    g_rho00->Write();
                    delete Title;

                    Title = new TString(Form("realRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             
                    g_real = new TGraphAsymmErrors();
                    g_real->SetName(Title->Data());
                    for(Int_t iy = y_padding; iy < vmsa::eta_total - y_padding; iy++) 
                    {
                      string key = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      double mean = (vmsa::eta_bin[iy]+vmsa::eta_bin[iy+1])/2.0; //(ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                      g_real->SetPoint(iy,mean,real[key][0]);
                      g_real->SetPointError(iy,0.0,0.0,realerr[key][1],real[key][1]);
                    }
                    g_real->Write();
                    delete Title;

                    Title = new TString(Form("imagRaw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             
                    g_imag = new TGraphAsymmErrors();
                    g_imag->SetName(Title->Data());
                    for(Int_t iy = y_padding; iy < vmsa::eta_total - y_padding; iy++) 
                    {
                      string key = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      double mean = (vmsa::eta_bin[iy]+vmsa::eta_bin[iy+1])/2.0; //(ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                      g_imag->SetPoint(iy,mean,imag[key][0]);
                      g_imag->SetPointError(iy,0.0,0.0,imag[key][1],imag[key][1]);
                    }
                    g_imag->Write();
                    delete Title;

                    Title = new TString(Form("rerho1n1Raw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             
                    g_rerho = new TGraphAsymmErrors();
                    g_rerho->SetName(Title->Data());
                    for(Int_t iy = y_padding; iy < vmsa::eta_total - y_padding; iy++) 
                    {
                      string key = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      double mean = (vmsa::eta_bin[iy]+vmsa::eta_bin[iy+1])/2.0; //(ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                      g_rerho->SetPoint(iy,mean,rerho[key][0]);
                      g_rerho->SetPointError(iy,0.0,0.0,rerho[key][1],rerho[key][1]);
                    }
                    g_rerho->Write();
                    delete Title;

                    Title = new TString(Form("imrho1n1Raw_pt_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1));                             
                    g_imrho = new TGraphAsymmErrors();
                    g_imrho->SetName(Title->Data());
                    for(Int_t iy = y_padding; iy < vmsa::eta_total - y_padding; iy++) 
                    {
                      string key = Form("pt_%d_eta_%d_Centrality_%d_%s_Dca_%d_Sig_%d_%s_Norm_%d_Sigma_%d_%s_Poly%d",ipt,iy,icent,ordertext[order-1].c_str(),i_dca,i_sig,vmsa::mPID[pid].c_str(),i_norm,i_sigma,vmsa::mInteMethod[i_method].c_str(),i_poly+1);
                      double mean = (vmsa::eta_bin[iy]+vmsa::eta_bin[iy+1])/2.0; //(ptbin_rho00[PtBin]+ptbin_rho00[PtBin-1])/2.0;
                      g_imrho->SetPoint(iy,mean,imrho[key][0]);
                      g_imrho->SetPointError(iy,0.0,0.0,imrho[key][1],imrho[key][1]);
                    }
                    g_imrho->Write();
                    delete Title;

                    g_rho00->Print();
                    g_real->Print();
                    g_imag->Print();
                    g_rerho->Print();
                    g_imrho->Print();
                    
                    output->Close();
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


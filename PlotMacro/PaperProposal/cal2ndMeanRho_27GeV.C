#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"

void cal2ndMeanRho_27GeV()
{
  const int mEpTotal = 2;
  const int mEffTotal = 2;
  const int mDcaTotal = 3;
  const int mSigTotal = 3;
  const int mNormTotal = 3;
  const int mModeTotal = 4;
  const string mMode[4] = {"Sigma_2_Inte","Sigma_0_Count","Sigma_1_Count","Sigma_2_Count"};

  const int dca_default = 0;
  const int sig_default = 1;
  const int norm_default = 0;
  const int eff_default = 0;
  // const int mode_default = 0;


  string input_run11 = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/SysErrors/NewF_JHChen/rho00_27GeV_run11.root";
  TFile *File_Run11 = TFile::Open(input_run11.c_str());

  string input_run18 = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/SysErrors/NewF_JHChen/rho00_27GeV_run18.root";
  TFile *File_Run18 = TFile::Open(input_run18.c_str());


  //==============================================================
  //--------------------------------------------------------------
  // get default value
  string HistName_Default = Form("EP_2_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_Sigma_2_Inte",eff_default,dca_default,sig_default,norm_default);
  cout << "Default Histogram set to: " << HistName_Default.c_str() << endl;
  TH1F *h_rhoDef_Run11 = (TH1F*)File_Run11->Get(HistName_Default.c_str());
  TH1F *h_rhoDef_Run18 = (TH1F*)File_Run18->Get(HistName_Default.c_str());

  HistName_Default = Form("EP_2_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_Sigma_2_Inte",eff_default,dca_default,sig_default,norm_default);
  TH1F *h_rhoDef_Mean = (TH1F*)h_rhoDef_Run11->Clone(HistName_Default.c_str());
  h_rhoDef_Mean->Reset(); // get same histogram wo any bin content
  for(int i_pt = 0; i_pt < h_rhoDef_Run11->GetNbinsX(); ++i_pt) 
  {
    double rho_run11 = h_rhoDef_Run11->GetBinContent(i_pt+1);
    double err_run11 = h_rhoDef_Run11->GetBinError(i_pt+1);
    double rho_run18 = h_rhoDef_Run18->GetBinContent(i_pt+1);
    double err_run18 = h_rhoDef_Run18->GetBinError(i_pt+1);

    if(err_run11 > 0 && err_run18 > 0)
    {
      double rho_mean = (rho_run11/err_run11/err_run11+rho_run18/err_run18/err_run18)/(1./err_run11/err_run11+1./err_run18/err_run18);
      double err_mean = sqrt(1./(1./err_run11/err_run11+1./err_run18/err_run18));
      // cout << "rho_run11 = " << rho_run11 << ", rho_run18 = " << rho_run18 << ", rho_mean = " << rho_mean << endl;
      // cout << "err_run11 = " << err_run11 << ", err_run18 = " << err_run18 << ", err_mean = " << err_mean << endl;
      h_rhoDef_Mean->SetBinContent(i_pt+1,rho_mean);
      h_rhoDef_Mean->SetBinError(i_pt+1,err_mean);
    }
  }
  // get default value
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // eff sysmatic value
  TH1F *h_rhoSys_Eff_Run11[2][4]; // 0 for different Eff | 1 for different yields extraction
  TH1F *h_rhoSys_Eff_Run18[2][4]; // 0 for different Eff | 1 for different yields extraction
  TH1F *h_rhoSys_Eff_Mean[2][4]; // 0 for different Eff | 1 for different yields extraction
  for(int i_eff = 0; i_eff < 2; ++i_eff)
  {
    if(i_eff == eff_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string HistName = Form("EP_2_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",i_eff,dca_default,sig_default,norm_default,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from Eff: " << HistName.c_str() << endl;
      h_rhoSys_Eff_Run11[i_eff][i_mode] = (TH1F*)File_Run11->Get(HistName.c_str());
      h_rhoSys_Eff_Run18[i_eff][i_mode] = (TH1F*)File_Run18->Get(HistName.c_str());

      HistName = Form("EP_2_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",i_eff,dca_default,sig_default,norm_default,mMode[i_mode].c_str());
      h_rhoSys_Eff_Mean[i_eff][i_mode] = (TH1F*)h_rhoSys_Eff_Run11[i_eff][i_mode]->Clone(HistName.c_str());
      h_rhoSys_Eff_Mean[i_eff][i_mode]->Reset(); // get same histogram wo any bin content
      // h_rhoSys_Eff_Mean[i_eff][i_mode]->Sumw2(); 

      for(int i_pt = 0; i_pt < h_rhoSys_Eff_Mean[i_eff][i_mode]->GetNbinsX(); ++i_pt) 
      {
	double rho_run11 = h_rhoSys_Eff_Run11[i_eff][i_mode]->GetBinContent(i_pt+1);
	double err_run11 = h_rhoSys_Eff_Run11[i_eff][i_mode]->GetBinError(i_pt+1);
	double rho_run18 = h_rhoSys_Eff_Run18[i_eff][i_mode]->GetBinContent(i_pt+1);
	double err_run18 = h_rhoSys_Eff_Run18[i_eff][i_mode]->GetBinError(i_pt+1);

	if(err_run11 > 0 && err_run18 > 0)
	{
	  double rho_mean = (rho_run11/err_run11/err_run11+rho_run18/err_run18/err_run18)/(1./err_run11/err_run11+1./err_run18/err_run18);
	  double err_mean = sqrt(1./(1./err_run11/err_run11+1./err_run18/err_run18));
	  // cout << "rho_run11 = " << rho_run11 << ", rho_run18 = " << rho_run18 << ", rho_mean = " << rho_mean << endl;
	  // cout << "err_run11 = " << err_run11 << ", err_run18 = " << err_run18 << ", err_mean = " << err_mean << endl;
	  h_rhoSys_Eff_Mean[i_eff][i_mode]->SetBinContent(i_pt+1,rho_mean);
	  h_rhoSys_Eff_Mean[i_eff][i_mode]->SetBinError(i_pt+1,err_mean);
	}
      }
    }
  }
  // eff sysmatic value
  //--------------------------------------------------------------


  //--------------------------------------------------------------
  // dca sysmatic value
  TH1F *h_rhoSys_Dca_Run11[3][4]; // 0 for different Dca | 1 for different yields extraction
  TH1F *h_rhoSys_Dca_Run18[3][4]; // 0 for different Dca | 1 for different yields extraction
  TH1F *h_rhoSys_Dca_Mean[3][4]; // 0 for different Dca | 1 for different yields extraction
  for(int i_dca = 0; i_dca < 3; ++i_dca)
  {
    if(i_dca == dca_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string HistName = Form("EP_2_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,i_dca,sig_default,norm_default,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from Dca: " << HistName.c_str() << endl;
      h_rhoSys_Dca_Run11[i_dca][i_mode] = (TH1F*)File_Run11->Get(HistName.c_str());
      h_rhoSys_Dca_Run18[i_dca][i_mode] = (TH1F*)File_Run18->Get(HistName.c_str());

      HistName = Form("EP_2_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,i_dca,sig_default,norm_default,mMode[i_mode].c_str());
      h_rhoSys_Dca_Mean[i_dca][i_mode] = (TH1F*)h_rhoSys_Dca_Run11[i_dca][i_mode]->Clone(HistName.c_str());
      h_rhoSys_Dca_Mean[i_dca][i_mode]->Reset(); // get same histogram wo any bin content
      // h_rhoSys_Dca_Mean[i_dca][i_mode]->Sumw2(); 

      for(int i_pt = 0; i_pt < h_rhoSys_Dca_Mean[i_dca][i_mode]->GetNbinsX(); ++i_pt) 
      {
	double rho_run11 = h_rhoSys_Dca_Run11[i_dca][i_mode]->GetBinContent(i_pt+1);
	double err_run11 = h_rhoSys_Dca_Run11[i_dca][i_mode]->GetBinError(i_pt+1);
	double rho_run18 = h_rhoSys_Dca_Run18[i_dca][i_mode]->GetBinContent(i_pt+1);
	double err_run18 = h_rhoSys_Dca_Run18[i_dca][i_mode]->GetBinError(i_pt+1);

	if(err_run11 > 0 && err_run18 > 0)
	{
	  double rho_mean = (rho_run11/err_run11/err_run11+rho_run18/err_run18/err_run18)/(1./err_run11/err_run11+1./err_run18/err_run18);
	  double err_mean = sqrt(1./(1./err_run11/err_run11+1./err_run18/err_run18));
	  // cout << "rho_run11 = " << rho_run11 << ", rho_run18 = " << rho_run18 << ", rho_mean = " << rho_mean << endl;
	  // cout << "err_run11 = " << err_run11 << ", err_run18 = " << err_run18 << ", err_mean = " << err_mean << endl;
	  h_rhoSys_Dca_Mean[i_dca][i_mode]->SetBinContent(i_pt+1,rho_mean);
	  h_rhoSys_Dca_Mean[i_dca][i_mode]->SetBinError(i_pt+1,err_mean);
	}
      }
    }
  }
  // dca sysmatic value
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // sig sysmatic value
  TH1F *h_rhoSys_Sig_Run11[3][4]; // 0 for different Sig | 1 for different yields extraction
  TH1F *h_rhoSys_Sig_Run18[3][4]; // 0 for different Sig | 1 for different yields extraction
  TH1F *h_rhoSys_Sig_Mean[3][4]; // 0 for different Sig | 1 for different yields extraction
  for(int i_sig = 0; i_sig < 3; ++i_sig)
  {
    if(i_sig == sig_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string HistName = Form("EP_2_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,dca_default,i_sig,norm_default,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from Sig: " << HistName.c_str() << endl;
      h_rhoSys_Sig_Run11[i_sig][i_mode] = (TH1F*)File_Run11->Get(HistName.c_str());
      h_rhoSys_Sig_Run18[i_sig][i_mode] = (TH1F*)File_Run18->Get(HistName.c_str());

      HistName = Form("EP_2_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,dca_default,i_sig,norm_default,mMode[i_mode].c_str());
      h_rhoSys_Sig_Mean[i_sig][i_mode] = (TH1F*)h_rhoSys_Sig_Run11[i_sig][i_mode]->Clone(HistName.c_str());
      h_rhoSys_Sig_Mean[i_sig][i_mode]->Reset(); // get same histogram wo any bin content
      // h_rhoSys_Sig_Mean[i_sig][i_mode]->Sumw2(); 

      for(int i_pt = 0; i_pt < h_rhoSys_Sig_Mean[i_sig][i_mode]->GetNbinsX(); ++i_pt) 
      {
	double rho_run11 = h_rhoSys_Sig_Run11[i_sig][i_mode]->GetBinContent(i_pt+1);
	double err_run11 = h_rhoSys_Sig_Run11[i_sig][i_mode]->GetBinError(i_pt+1);
	double rho_run18 = h_rhoSys_Sig_Run18[i_sig][i_mode]->GetBinContent(i_pt+1);
	double err_run18 = h_rhoSys_Sig_Run18[i_sig][i_mode]->GetBinError(i_pt+1);

	if(err_run11 > 0 && err_run18 > 0)
	{
	  double rho_mean = (rho_run11/err_run11/err_run11+rho_run18/err_run18/err_run18)/(1./err_run11/err_run11+1./err_run18/err_run18);
	  double err_mean = sqrt(1./(1./err_run11/err_run11+1./err_run18/err_run18));
	  // cout << "rho_run11 = " << rho_run11 << ", rho_run18 = " << rho_run18 << ", rho_mean = " << rho_mean << endl;
	  // cout << "err_run11 = " << err_run11 << ", err_run18 = " << err_run18 << ", err_mean = " << err_mean << endl;
	  h_rhoSys_Sig_Mean[i_sig][i_mode]->SetBinContent(i_pt+1,rho_mean);
	  h_rhoSys_Sig_Mean[i_sig][i_mode]->SetBinError(i_pt+1,err_mean);
	}
      }
    }
  }
  // sig sysmatic value
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // norm sysmatic value
  TH1F *h_rhoSys_Norm_Run11[3][4]; // 0 for different Norm | 1 for different yields extraction
  TH1F *h_rhoSys_Norm_Run18[3][4]; // 0 for different Norm | 1 for different yields extraction
  TH1F *h_rhoSys_Norm_Mean[3][4]; // 0 for different Norm | 1 for different yields extraction
  for(int i_norm = 0; i_norm < 3; ++i_norm)
  {
    if(i_norm == norm_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string HistName = Form("EP_2_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,dca_default,sig_default,i_norm,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from Norm: " << HistName.c_str() << endl;
      h_rhoSys_Norm_Run11[i_norm][i_mode] = (TH1F*)File_Run11->Get(HistName.c_str());
      h_rhoSys_Norm_Run18[i_norm][i_mode] = (TH1F*)File_Run18->Get(HistName.c_str());

      HistName = Form("EP_2_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,dca_default,sig_default,i_norm,mMode[i_mode].c_str());
      h_rhoSys_Norm_Mean[i_norm][i_mode] = (TH1F*)h_rhoSys_Norm_Run11[i_norm][i_mode]->Clone(HistName.c_str());
      h_rhoSys_Norm_Mean[i_norm][i_mode]->Reset(); // get same histogram wo any bin content
      // h_rhoSys_Norm_Mean[i_norm][i_mode]->Sumw2(); 

      for(int i_pt = 0; i_pt < h_rhoSys_Norm_Mean[i_norm][i_mode]->GetNbinsX(); ++i_pt) 
      {
	double rho_run11 = h_rhoSys_Norm_Run11[i_norm][i_mode]->GetBinContent(i_pt+1);
	double err_run11 = h_rhoSys_Norm_Run11[i_norm][i_mode]->GetBinError(i_pt+1);
	double rho_run18 = h_rhoSys_Norm_Run18[i_norm][i_mode]->GetBinContent(i_pt+1);
	double err_run18 = h_rhoSys_Norm_Run18[i_norm][i_mode]->GetBinError(i_pt+1);

	if(err_run11 > 0 && err_run18 > 0)
	{
	  double rho_mean = (rho_run11/err_run11/err_run11+rho_run18/err_run18/err_run18)/(1./err_run11/err_run11+1./err_run18/err_run18);
	  double err_mean = sqrt(1./(1./err_run11/err_run11+1./err_run18/err_run18));
	  // cout << "rho_run11 = " << rho_run11 << ", rho_run18 = " << rho_run18 << ", rho_mean = " << rho_mean << endl;
	  // cout << "err_run11 = " << err_run11 << ", err_run18 = " << err_run18 << ", err_mean = " << err_mean << endl;
	  h_rhoSys_Norm_Mean[i_norm][i_mode]->SetBinContent(i_pt+1,rho_mean);
	  h_rhoSys_Norm_Mean[i_norm][i_mode]->SetBinError(i_pt+1,err_mean);
	}
      }
    }
  }
  // norm sysmatic value
  //--------------------------------------------------------------

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/SysErrors/NewF_JHChen/rho00_27GeV_2ndMean.root";
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();

  // save default rho00
  h_rhoDef_Mean->Write();

  // save dca cuts
  for(int i_dca = 0; i_dca < 3; ++i_dca)
  {
    if(i_dca == dca_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      h_rhoSys_Dca_Mean[i_dca][i_mode]->Write();
    }
  }

  // save sig cuts
  for(int i_sig = 0; i_sig < 3; ++i_sig)
  {
    if(i_sig == sig_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      h_rhoSys_Sig_Mean[i_sig][i_mode]->Write();
    }
  }

  // plot norm cuts
  for(int i_norm = 0; i_norm < 3; ++i_norm)
  {
    if(i_norm == norm_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      h_rhoSys_Norm_Mean[i_norm][i_mode]->Write();
    }
  }

  // plot eff cuts
  for(int i_eff = 0; i_eff < 2; ++i_eff)
  {
    if(i_eff == eff_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      h_rhoSys_Eff_Mean[i_eff][i_mode]->Write();
    }
  }
  File_OutPut->Close();
}

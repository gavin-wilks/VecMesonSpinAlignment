#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TBox.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"

void cal1stEnergyInteRho()
{
  const string mEnergy[5] = {"11GeV","19GeV","27GeV_1stMean","39GeV","62GeV"};
  const string mMode[4] = {"Sigma_2_Inte","Sigma_0_Count","Sigma_1_Count","Sigma_2_Count"};

  const int dca_default = 0;
  const int sig_default = 1;
  const int norm_default = 0;
  const int eff_default = 0;
  // const int mode_default = 0;

  cout << "Start to calculate Energy Integrated rho00 for all systematic combinations!" << endl;
  TFile *File_InPut[5];
  for(int i_energy = 0; i_energy < 5; ++i_energy)
  {
    string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/SysErrors/rho00_%s.root",mEnergy[i_energy].c_str());
    File_InPut[i_energy] = TFile::Open(inputfile.c_str());
    cout << "Open input files: " << inputfile.c_str() << endl;
  }

  //==============================================================
  //--------------------------------------------------------------
  // get default value
  string HistName_Default = Form("EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_Sigma_2_Inte",eff_default,dca_default,sig_default,norm_default);
  cout << "Default Histogram set to: " << HistName_Default.c_str() << endl;
  TH1F *h_rhoDef[5];

  string GraphName_Default = Form("g_EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_Sigma_2_Inte",eff_default,dca_default,sig_default,norm_default);
  TGraphAsymmErrors *g_rhoDef = new TGraphAsymmErrors(); // init TGraphAsymmErrors to store integrated rho00
  g_rhoDef->SetName(GraphName_Default.c_str());

  double rhoSumDef = 0.0;
  double errSumDef = 0.0;
  double wgtSumDef = 0.0;
  for(int i_energy = 0; i_energy < 5; ++i_energy)
  {
    h_rhoDef[i_energy] = (TH1F*)File_InPut[i_energy]->Get(HistName_Default.c_str());
    double rho = h_rhoDef[i_energy]->GetBinContent(1); // first bin stored pT averaged rho00
    double err = h_rhoDef[i_energy]->GetBinError(1);
    double weight = 1.0/(err*err);
    rhoSumDef += rho*weight;
    errSumDef += err*err*weight*weight;
    wgtSumDef += weight;
    cout << "Energy: " << mEnergy[i_energy].c_str() << ", rho = " << rho << " +/- " << err << endl;
  }
  double rhoDef = rhoSumDef/wgtSumDef;
  double errDef = sqrt(errSumDef)/wgtSumDef;
  g_rhoDef->SetPoint(0,0.0,rhoDef);
  g_rhoDef->SetPointError(0,0.0,0.0,errDef,errDef);
  cout << "Default Energy integrated rho00 = " << rhoDef << " +/- " << errDef << endl;
  cout << endl;
  // get default value
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // eff sysmatic value
  TH1F *h_rhoSys_Eff[5][2][4]; // 0 for energy | 1 for different Eff | 2 for different yields extraction
  TGraphAsymmErrors *g_rhoSys_Eff[2][4];
  for(int i_eff = 0; i_eff < 2; ++i_eff)
  {
    if(i_eff == eff_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string HistName = Form("EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",i_eff,dca_default,sig_default,norm_default,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from Eff: " << HistName.c_str() << endl;

      string GraphName = Form("g_EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",i_eff,dca_default,sig_default,norm_default,mMode[i_mode].c_str());
      g_rhoSys_Eff[i_eff][i_mode] = new TGraphAsymmErrors(); // init TGraphAsymmErrors to store integrated value
      g_rhoSys_Eff[i_eff][i_mode]->SetName(GraphName.c_str());
      double rhoSumEff = 0.0;
      double errSumEff = 0.0;
      double wgtSumEff = 0.0;
      for(int i_energy = 0; i_energy < 5; ++i_energy)
      {
	h_rhoSys_Eff[i_energy][i_eff][i_mode] = (TH1F*)File_InPut[i_energy]->Get(HistName.c_str());

	double rho = h_rhoSys_Eff[i_energy][i_eff][i_mode]->GetBinContent(1); // first bin stored pT averaged rho00
	double err = h_rhoSys_Eff[i_energy][i_eff][i_mode]->GetBinError(1);
	double weight = 1.0/(err*err);
	rhoSumEff += rho*weight;
	errSumEff += err*err*weight*weight;
	wgtSumEff += weight;
	cout << "Energy: " << mEnergy[i_energy].c_str() << ", rho = " << rho << " +/- " << err << endl;
      }
      double rhoEff = rhoSumEff/wgtSumEff;
      double errEff = sqrt(errSumEff)/wgtSumEff;
      g_rhoSys_Eff[i_eff][i_mode]->SetPoint(0,0.0,rhoEff);
      g_rhoSys_Eff[i_eff][i_mode]->SetPointError(0,0.0,0.0,errEff,errEff);
      cout << "Energy integrated rho00 with different Eff = " << rhoEff << " +/- " << errEff << endl;
      cout << endl;
    }
  }
  // eff sysmatic value
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // dca sysmatic value
  TH1F *h_rhoSys_Dca[5][3][4]; // 0 for energy | 1 for different Dca | 2 for different yields extraction
  TGraphAsymmErrors *g_rhoSys_Dca[3][4];
  for(int i_dca = 0; i_dca < 3; ++i_dca)
  {
    if(i_dca == dca_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string HistName = Form("EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,i_dca,sig_default,norm_default,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from Dca: " << HistName.c_str() << endl;

      string GraphName = Form("g_EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,i_dca,sig_default,norm_default,mMode[i_mode].c_str());
      g_rhoSys_Dca[i_dca][i_mode] = new TGraphAsymmErrors(); // init TGraphAsymmErrors to store integrated value
      g_rhoSys_Dca[i_dca][i_mode]->SetName(GraphName.c_str());
      double rhoSumDca = 0.0;
      double errSumDca = 0.0;
      double wgtSumDca = 0.0;
      for(int i_energy = 0; i_energy < 5; ++i_energy)
      {
	h_rhoSys_Dca[i_energy][i_dca][i_mode] = (TH1F*)File_InPut[i_energy]->Get(HistName.c_str());

	double rho = h_rhoSys_Dca[i_energy][i_dca][i_mode]->GetBinContent(1); // first bin stored pT averaged rho00
	double err = h_rhoSys_Dca[i_energy][i_dca][i_mode]->GetBinError(1);
	double weight = 1.0/(err*err);
	rhoSumDca += rho*weight;
	errSumDca += err*err*weight*weight;
	wgtSumDca += weight;
	cout << "Energy: " << mEnergy[i_energy].c_str() << ", rho = " << rho << " +/- " << err << endl;
      }
      double rhoDca = rhoSumDca/wgtSumDca;
      double errDca = sqrt(errSumDca)/wgtSumDca;
      g_rhoSys_Dca[i_dca][i_mode]->SetPoint(0,0.0,rhoDca);
      g_rhoSys_Dca[i_dca][i_mode]->SetPointError(0,0.0,0.0,errDca,errDca);
      cout << "Energy integrated rho00 with different Dca = " << rhoDca << " +/- " << errDca << endl;
      cout << endl;
    }
  }
  // dca sysmatic value
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // sig sysmatic value
  TH1F *h_rhoSys_Sig[5][3][4]; // 0 for energy | 1 for different Sig | 2 for different yields extraction
  TGraphAsymmErrors *g_rhoSys_Sig[3][4];
  for(int i_sig = 0; i_sig < 3; ++i_sig)
  {
    if(i_sig == sig_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string HistName = Form("EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,dca_default,i_sig,norm_default,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from Sig: " << HistName.c_str() << endl;

      string GraphName = Form("g_EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,dca_default,i_sig,norm_default,mMode[i_mode].c_str());
      g_rhoSys_Sig[i_sig][i_mode] = new TGraphAsymmErrors(); // init TGraphAsymmErrors to store integrated value
      g_rhoSys_Sig[i_sig][i_mode]->SetName(GraphName.c_str());
      double rhoSumSig = 0.0;
      double errSumSig = 0.0;
      double wgtSumSig = 0.0;
      for(int i_energy = 0; i_energy < 5; ++i_energy)
      {
	h_rhoSys_Sig[i_energy][i_sig][i_mode] = (TH1F*)File_InPut[i_energy]->Get(HistName.c_str());

	double rho = h_rhoSys_Sig[i_energy][i_sig][i_mode]->GetBinContent(1); // first bin stored pT averaged rho00
	double err = h_rhoSys_Sig[i_energy][i_sig][i_mode]->GetBinError(1);
	double weight = 1.0/(err*err);
	rhoSumSig += rho*weight;
	errSumSig += err*err*weight*weight;
	wgtSumSig += weight;
	cout << "Energy: " << mEnergy[i_energy].c_str() << ", rho = " << rho << " +/- " << err << endl;
      }
      double rhoSig = rhoSumSig/wgtSumSig;
      double errSig = sqrt(errSumSig)/wgtSumSig;
      g_rhoSys_Sig[i_sig][i_mode]->SetPoint(0,0.0,rhoSig);
      g_rhoSys_Sig[i_sig][i_mode]->SetPointError(0,0.0,0.0,errSig,errSig);
      cout << "Energy integrated rho00 with different Sig = " << rhoSig << " +/- " << errSig << endl;
      cout << endl;
    }
  }
  // sig sysmatic value
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // norm sysmatic value
  TH1F *h_rhoSys_Norm[5][3][4]; // 0 for Energy | 1 for different Norm | 2 for different yields extraction
  TGraphAsymmErrors *g_rhoSys_Norm[3][4];
  for(int i_norm = 0; i_norm < 3; ++i_norm)
  {
    if(i_norm == norm_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      string HistName = Form("EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,dca_default,sig_default,i_norm,mMode[i_mode].c_str());
      cout << "Read in Systematic Contribution from Norm: " << HistName.c_str() << endl;

      string GraphName = Form("g_EP_1_eff_%d_Dca_%d_Sig_%d_Phi_Norm_%d_%s",eff_default,dca_default,sig_default,i_norm,mMode[i_mode].c_str());
      g_rhoSys_Norm[i_norm][i_mode] = new TGraphAsymmErrors(); // init TGraphAsymmErrors to store integrated value
      g_rhoSys_Norm[i_norm][i_mode]->SetName(GraphName.c_str());
      double rhoSumNorm = 0.0;
      double errSumNorm = 0.0;
      double wgtSumNorm = 0.0;
      for(int i_energy = 0; i_energy < 5; ++i_energy)
      {
	h_rhoSys_Norm[i_energy][i_norm][i_mode] = (TH1F*)File_InPut[i_energy]->Get(HistName.c_str());

	double rho = h_rhoSys_Norm[i_energy][i_norm][i_mode]->GetBinContent(1); // first bin stored pT averaged rho00
	double err = h_rhoSys_Norm[i_energy][i_norm][i_mode]->GetBinError(1);
	double weight = 1.0/(err*err);
	rhoSumNorm += rho*weight;
	errSumNorm += err*err*weight*weight;
	wgtSumNorm += weight;
	cout << "Energy: " << mEnergy[i_energy].c_str() << ", rho = " << rho << " +/- " << err << endl;
      }
      double rhoNorm = rhoSumNorm/wgtSumNorm;
      double errNorm = sqrt(errSumNorm)/wgtSumNorm;
      g_rhoSys_Norm[i_norm][i_mode]->SetPoint(0,0.0,rhoNorm);
      g_rhoSys_Norm[i_norm][i_mode]->SetPointError(0,0.0,0.0,errNorm,errNorm);
      cout << "Energy integrated rho00 with different Norm = " << rhoNorm << " +/- " << errNorm << endl;
      cout << endl;
    }
  }
  // norm sysmatic value
  //--------------------------------------------------------------

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperProposal/SysErrors/rho00_1stEnergyInte.root";
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();

  // save default rho00
  g_rhoDef->Write();

  // save dca cuts
  for(int i_dca = 0; i_dca < 3; ++i_dca)
  {
    if(i_dca == dca_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      g_rhoSys_Dca[i_dca][i_mode]->Write();
    }
  }

  // save sig cuts
  for(int i_sig = 0; i_sig < 3; ++i_sig)
  {
    if(i_sig == sig_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      g_rhoSys_Sig[i_sig][i_mode]->Write();
    }
  }

  // plot norm cuts
  for(int i_norm = 0; i_norm < 3; ++i_norm)
  {
    if(i_norm == norm_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      g_rhoSys_Norm[i_norm][i_mode]->Write();
    }
  }

  // plot eff cuts
  for(int i_eff = 0; i_eff < 2; ++i_eff)
  {
    if(i_eff == eff_default) continue;
    for(int i_mode = 0; i_mode < 4; ++i_mode)
    {
      g_rhoSys_Eff[i_eff][i_mode]->Write();
    }
  }
  File_OutPut->Close();
}

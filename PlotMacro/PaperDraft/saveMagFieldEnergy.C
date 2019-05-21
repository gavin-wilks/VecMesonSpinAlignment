#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "../../Utility/draw.h"
#include "../../Utility/StSpinAlignmentCons.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TMath.h"

using namespace std;

float ErrorAdd(float x, float y)
{
  return sqrt(x*x+y*y);
}

float ErrTimes(float x, float y, float dx, float dy)
{
  return x*y*ErrorAdd(dx/x,dy/y);
}

float ErrDiv(float x, float y, float dx, float dy)
{
  return x/y*ErrorAdd(dx/x,dy/y);
}

void saveMagFieldEnergy()
{
  const double mass = 95.0; // MeV | mass of strange quark
  const double err_mass = 5.0; // MeV
  const double eVsquare = 0.0169; // 1eV^2 = 0.0169T => 16.9mT
  const double temp = 120.0; // MeV => 10^6 eV

  gStyle->SetOptDate(0);
  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/rho00_stat_sys.root");
  TGraphAsymmErrors *g_rho_1st_stat = (TGraphAsymmErrors*)File_Input->Get("rho00_1stEP_energy_stat");;
  TGraphAsymmErrors *g_rho_1st_sys  = (TGraphAsymmErrors*)File_Input->Get("rho00_1stEP_energy_sys");
  TGraphAsymmErrors *g_eMagOverT_1st = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_eMag_1st = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_eMagTesla_1st = new TGraphAsymmErrors();

  const int NumOfEnergies_1st = g_rho_1st_stat->GetN();
  for(int i_energy = 0; i_energy < NumOfEnergies_1st; ++i_energy)
  {
    double energy, rho;
    g_rho_1st_stat->GetPoint(i_energy,energy,rho);
    double err_stat  = g_rho_1st_stat->GetErrorYhigh(i_energy);
    double err_sys   = g_rho_1st_sys->GetErrorYhigh(i_energy);
    double error_rho = TMath::Sqrt(err_stat*err_stat+err_sys*err_sys);
    cout << "1st: i_energy = " << i_energy << ", energy = " << energy << ", rho = " << rho << " +/- " << err_stat << " +/- " << err_sys << ", error_rho = " << error_rho << endl;

    double eMagOverT      = 3.0*mass*TMath::Sqrt(9.0*(rho-1.0/3.0)); // eB/T = 3m_{s}#sqrt{9(#rho_{00} - 1/3)}
    double diff_eMagOverT = 0.5/TMath::Sqrt(9.0*(rho-1.0/3.0));
    double err_eMagOverT  = 3.0*mass*diff_eMagOverT*9.0*error_rho;
    g_eMagOverT_1st->SetPoint(i_energy,energy,eMagOverT);
    g_eMagOverT_1st->SetPointError(i_energy,0.0,0.0,err_eMagOverT,err_eMagOverT);
    cout << "1st: eMagOverT = " << eMagOverT << " +/- " << err_eMagOverT << endl;

    double eMag     = eMagOverT*temp; // eB = T*3m_{s}#sqrt{9(#rho_{00} - 1/3)}
    double err_eMag = err_eMagOverT*temp;
    g_eMag_1st->SetPoint(i_energy,energy,eMag);
    g_eMag_1st->SetPointError(i_energy,0.0,0.0,err_eMag,err_eMag);
    cout << "1st: eMag = " << eMag << " +/- " << err_eMag << endl;

    double eMagTesla     = eMag*1e6*1e6*eVsquare; // 1eV^2 = 0.0169T
    double err_eMagTesla = err_eMag*1e6*1e6*eVsquare;
    g_eMagTesla_1st->SetPoint(i_energy,energy,eMagTesla);
    g_eMagTesla_1st->SetPointError(i_energy,0.0,0.0,err_eMagTesla,err_eMagTesla);
    cout << "1st: eMag in Tesla = " << eMagTesla << " +/- " << err_eMagTesla << endl;
  }

  TGraphAsymmErrors *g_rho_2nd_stat = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_energy_stat");;
  TGraphAsymmErrors *g_rho_2nd_sys  = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_energy_sys");
  TGraphAsymmErrors *g_eMagOverT_2nd = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_eMag_2nd = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_eMagTesla_2nd = new TGraphAsymmErrors();

  const int NumOfEnergies_2nd = g_rho_2nd_stat->GetN();
  for(int i_energy = 0; i_energy < NumOfEnergies_2nd; ++i_energy)
  {
    double energy, rho;
    g_rho_2nd_stat->GetPoint(i_energy,energy,rho);
    double err_stat  = g_rho_2nd_stat->GetErrorYhigh(i_energy);
    double err_sys   = g_rho_2nd_sys->GetErrorYhigh(i_energy);
    double error_rho = TMath::Sqrt(err_stat*err_stat+err_sys*err_sys);
    cout << "2nd: i_energy = " << i_energy << ", energy = " << energy << ", rho = " << rho << " +/- " << err_stat << " +/- " << err_sys << ", error_rho = " << error_rho << endl;

    double eMagOverT      = 3.0*mass*TMath::Sqrt(9.0*(rho-1.0/3.0)); // eB/T = 3m_{s}#sqrt{9(#rho_{00} - 1/3)}
    double diff_eMagOverT = 0.5/TMath::Sqrt(9.0*(rho-1.0/3.0));
    double err_eMagOverT  = 3.0*mass*diff_eMagOverT*9.0*error_rho;
    g_eMagOverT_2nd->SetPoint(i_energy,energy,eMagOverT);
    g_eMagOverT_2nd->SetPointError(i_energy,0.0,0.0,err_eMagOverT,err_eMagOverT);
    cout << "2nd: eMagOverT = " << eMagOverT << " +/- " << err_eMagOverT << endl;

    double eMag     = eMagOverT*temp; // eB = T*3m_{s}#sqrt{9(#rho_{00} - 1/3)}
    double err_eMag = err_eMagOverT*temp;
    g_eMag_2nd->SetPoint(i_energy,energy,eMag);
    g_eMag_2nd->SetPointError(i_energy,0.0,0.0,err_eMag,err_eMag);
    cout << "2nd: eMag = " << eMag << " +/- " << err_eMag << endl;

    double eMagTesla     = eMag*1e6*1e6*eVsquare; // 1eV^2 = 0.0169T
    double err_eMagTesla = err_eMag*1e6*1e6*eVsquare;
    g_eMagTesla_2nd->SetPoint(i_energy,energy,eMagTesla);
    g_eMagTesla_2nd->SetPointError(i_energy,0.0,0.0,err_eMagTesla,err_eMagTesla);
    cout << "2nd: eMag in Tesla = " << eMagTesla << " +/- " << err_eMagTesla << endl;
  }

  // for mean magnetic field
  TFile *File_Input_mean = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/rho00_mean_BES.root");
  TGraphAsymmErrors *g_mean_rho_1st_stat = (TGraphAsymmErrors*)File_Input_mean->Get("rho00_1stEP_mean_stat");;
  TGraphAsymmErrors *g_mean_rho_1st_sys  = (TGraphAsymmErrors*)File_Input_mean->Get("rho00_1stEP_mean_sys");
  TGraphAsymmErrors *g_eMeanMagOverT_1st = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_eMeanMag_1st = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_eMeanMagTesla_1st = new TGraphAsymmErrors();

  const int NumOfMeanEnergies_1st = g_mean_rho_1st_stat->GetN();
  for(int i_energy = 0; i_energy < NumOfMeanEnergies_1st; ++i_energy)
  {
    double energy, rho;
    g_mean_rho_1st_stat->GetPoint(i_energy,energy,rho);
    double err_stat  = g_mean_rho_1st_stat->GetErrorYhigh(i_energy);
    double err_sys   = g_mean_rho_1st_sys->GetErrorYhigh(i_energy);
    double error_rho = TMath::Sqrt(err_stat*err_stat+err_sys*err_sys);
    cout << "1st: mean energy = " << energy << ", rho = " << rho << " +/- " << err_stat << " +/- " << err_sys << ", error_rho = " << error_rho << endl;

    double eMagOverT      = 3.0*mass*TMath::Sqrt(9.0*(rho-1.0/3.0)); // eB/T = 3m_{s}#sqrt{9(#rho_{00} - 1/3)}
    double diff_eMagOverT = 0.5/TMath::Sqrt(9.0*(rho-1.0/3.0));
    double err_eMagOverT  = 3.0*mass*diff_eMagOverT*9.0*error_rho;
    g_eMeanMagOverT_1st->SetPoint(i_energy,energy,eMagOverT);
    g_eMeanMagOverT_1st->SetPointError(i_energy,0.0,0.0,err_eMagOverT,err_eMagOverT);
    cout << "1st: mean eMagOverT = " << eMagOverT << " +/- " << err_eMagOverT << endl;

    double eMag     = eMagOverT*temp; // eB = T*3m_{s}#sqrt{9(#rho_{00} - 1/3)}
    double err_eMag = err_eMagOverT*temp;
    g_eMeanMag_1st->SetPoint(i_energy,energy,eMag);
    g_eMeanMag_1st->SetPointError(i_energy,0.0,0.0,err_eMag,err_eMag);
    cout << "1st: mean eMag = " << eMag << " +/- " << err_eMag << endl;

    double eMagTesla     = eMag*1e6*1e6*eVsquare; // 1eV^2 = 0.0169T
    double err_eMagTesla = err_eMag*1e6*1e6*eVsquare;
    g_eMeanMagTesla_1st->SetPoint(i_energy,energy,eMagTesla);
    g_eMeanMagTesla_1st->SetPointError(i_energy,0.0,0.0,err_eMagTesla,err_eMagTesla);
    cout << "1st: eMag in Tesla = " << eMagTesla << " +/- " << err_eMagTesla << endl;
  }

  TGraphAsymmErrors *g_mean_rho_2nd_stat = (TGraphAsymmErrors*)File_Input_mean->Get("rho00_2ndEP_mean_stat");;
  TGraphAsymmErrors *g_mean_rho_2nd_sys  = (TGraphAsymmErrors*)File_Input_mean->Get("rho00_2ndEP_mean_sys");
  TGraphAsymmErrors *g_eMeanMagOverT_2nd = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_eMeanMag_2nd = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_eMeanMagTesla_2nd = new TGraphAsymmErrors();

  const int NumOfMeanEnergies_2nd = g_mean_rho_2nd_stat->GetN();
  for(int i_energy = 0; i_energy < NumOfMeanEnergies_2nd; ++i_energy)
  {
    double energy, rho;
    g_mean_rho_2nd_stat->GetPoint(i_energy,energy,rho);
    double err_stat  = g_mean_rho_2nd_stat->GetErrorYhigh(i_energy);
    double err_sys   = g_mean_rho_2nd_sys->GetErrorYhigh(i_energy);
    double error_rho = TMath::Sqrt(err_stat*err_stat+err_sys*err_sys);
    cout << "2nd: mean energy = " << energy << ", rho = " << rho << " +/- " << err_stat << " +/- " << err_sys << ", error_rho = " << error_rho << endl;

    double eMagOverT      = 3.0*mass*TMath::Sqrt(9.0*(rho-1.0/3.0)); // eB/T = 3m_{s}#sqrt{9(#rho_{00} - 1/3)}
    double diff_eMagOverT = 0.5/TMath::Sqrt(9.0*(rho-1.0/3.0));
    double err_eMagOverT  = 3.0*mass*diff_eMagOverT*9.0*error_rho;
    g_eMeanMagOverT_2nd->SetPoint(i_energy,energy,eMagOverT);
    g_eMeanMagOverT_2nd->SetPointError(i_energy,0.0,0.0,err_eMagOverT,err_eMagOverT);
    cout << "2nd: mean eMagOverT = " << eMagOverT << " +/- " << err_eMagOverT << endl;

    double eMag     = eMagOverT*temp; // eB = T*3m_{s}#sqrt{9(#rho_{00} - 1/3)}
    double err_eMag = err_eMagOverT*temp;
    g_eMeanMag_2nd->SetPoint(i_energy,energy,eMag);
    g_eMeanMag_2nd->SetPointError(i_energy,0.0,0.0,err_eMag,err_eMag);
    cout << "2nd: mean eMag = " << eMag << " +/- " << err_eMag << endl;

    double eMagTesla     = eMag*1e6*1e6*eVsquare; // 1eV^2 = 0.0169T
    double err_eMagTesla = err_eMag*1e6*1e6*eVsquare;
    g_eMeanMagTesla_2nd->SetPoint(i_energy,energy,eMagTesla);
    g_eMeanMagTesla_2nd->SetPointError(i_energy,0.0,0.0,err_eMagTesla,err_eMagTesla);
    cout << "2nd: eMag in Tesla = " << eMagTesla << " +/- " << err_eMagTesla << endl;
  }

  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/magFiled_energy_SysErrors.root","RECREATE");
  File_OutPut->cd();
  g_eMagOverT_1st->SetName("g_eMagOverT_1st");
  g_eMagOverT_1st->Write();
  g_eMag_1st->SetName("g_eMag_1st");
  g_eMag_1st->Write();
  g_eMagTesla_1st->SetName("g_eMagTesla_1st");
  g_eMagTesla_1st->Write();

  g_eMagOverT_2nd->SetName("g_eMagOverT_2nd");
  g_eMagOverT_2nd->Write();
  g_eMag_2nd->SetName("g_eMag_2nd");
  g_eMag_2nd->Write();
  g_eMagTesla_2nd->SetName("g_eMagTesla_2nd");
  g_eMagTesla_2nd->Write();

  g_eMeanMagOverT_1st->SetName("g_eMeanMagOverT_1st");
  g_eMeanMagOverT_1st->Write();
  g_eMeanMag_1st->SetName("g_eMeanMag_1st");
  g_eMeanMag_1st->Write();
  g_eMeanMagTesla_1st->SetName("g_eMeanMagTesla_1st");
  g_eMeanMagTesla_1st->Write();

  g_eMeanMagOverT_2nd->SetName("g_eMeanMagOverT_2nd");
  g_eMeanMagOverT_2nd->Write();
  g_eMeanMag_2nd->SetName("g_eMeanMag_2nd");
  g_eMeanMag_2nd->Write();
  g_eMeanMagTesla_2nd->SetName("g_eMeanMagTesla_2nd");
  g_eMeanMagTesla_2nd->Write();
  File_OutPut->Close();
}

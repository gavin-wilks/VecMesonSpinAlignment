#include <TGraphAsymmErrors.h>
#include "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_sNN_May21_2021.h"

void makeData_rho00EnergyKstar()
{
  // K* STAR
  TGraphAsymmErrors *g_rhoKstar_stat = new TGraphAsymmErrors();
  g_rhoKstar_stat->SetName("g_rhoKstar_stat");
  TGraphAsymmErrors *g_rhoKstar_sys = new TGraphAsymmErrors();
  g_rhoKstar_sys->SetName("g_rhoKstar_sys");
  for(int i_point = 0; i_point < kstar_NsNN; ++i_point)
  {
    g_rhoKstar_stat->SetPoint(i_point,kstar_sNN[i_point],kstar_rho00_sNN[i_point]);
    g_rhoKstar_stat->SetPointError(i_point,0.0,0.0,kstar_rho00_sNN_stat[i_point],kstar_rho00_sNN_stat[i_point]);
    g_rhoKstar_sys->SetPoint(i_point,kstar_sNN[i_point],kstar_rho00_sNN[i_point]);
    g_rhoKstar_sys->SetPointError(i_point,0.0,0.0,kstar_rho00_sNN_syst[i_point],kstar_rho00_sNN_syst[i_point]);
  }

  // phi-meson ALICE
  TGraphAsymmErrors *g_rhoKstar_ALICE_stat = new TGraphAsymmErrors();
  g_rhoKstar_ALICE_stat->SetName("g_rhoKstar_ALICE_stat");
  g_rhoKstar_ALICE_stat->SetPoint(0,sNN_alice2760_kstar[0],Rho00_pT12_50_alice2760_kstar[0]);
  g_rhoKstar_ALICE_stat->SetPointError(0,0.0,0.0,err_Rho00_pT12_50_alice2760_kstar[0], err_Rho00_pT12_50_alice2760_kstar[0]);
  TGraphAsymmErrors *g_rhoKstar_ALICE_sys = new TGraphAsymmErrors();
  g_rhoKstar_ALICE_sys->SetName("g_rhoKstar_ALICE_sys");
  g_rhoKstar_ALICE_sys->SetPoint(0,sNN_alice2760_kstar[0],Rho00_pT12_50_alice2760_kstar[0]);
  g_rhoKstar_ALICE_sys->SetPointError(0,0.0,0.0,syst_Rho00_pT12_50_alice2760_kstar[0], syst_Rho00_pT12_50_alice2760_kstar[0]);

  // K* ALICE
  TGraphAsymmErrors *g_rhoPhi_ALICE_stat = new TGraphAsymmErrors();
  g_rhoPhi_ALICE_stat->SetName("g_rhoPhi_ALICE_stat");
  g_rhoPhi_ALICE_stat->SetPoint(0,sNN_alice2760_phim[0]-300.0,Rho00_pT12_50_alice2760_phim[0]);
  g_rhoPhi_ALICE_stat->SetPointError(0,0.0,0.0,err_Rho00_pT12_50_alice2760_phim[0], err_Rho00_pT12_50_alice2760_phim[0]);
  TGraphAsymmErrors *g_rhoPhi_ALICE_sys = new TGraphAsymmErrors();
  g_rhoPhi_ALICE_sys->SetName("g_rhoPhi_ALICE_sys");
  g_rhoPhi_ALICE_sys->SetPoint(0,sNN_alice2760_phim[0]-300.0,Rho00_pT12_50_alice2760_phim[0]);
  g_rhoPhi_ALICE_sys->SetPointError(0,0.0,0.0,syst_Rho00_pT12_50_alice2760_phim[0], syst_Rho00_pT12_50_alice2760_phim[0]);

  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_sNN_May21_2021.root","RECREATE");
  File_OutPut->cd();
  g_rhoKstar_stat->Write();
  g_rhoKstar_sys->Write();
  g_rhoKstar_ALICE_stat->Write();
  g_rhoKstar_ALICE_sys->Write();
  g_rhoPhi_ALICE_stat->Write();
  g_rhoPhi_ALICE_sys->Write();
  File_OutPut->Close();
}

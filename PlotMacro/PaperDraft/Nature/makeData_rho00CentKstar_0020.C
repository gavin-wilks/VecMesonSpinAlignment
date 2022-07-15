#include <TGraphAsymmErrors.h>
#include "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_Cent_March26_2022.h"

void makeData_rho00CentKstar_0020()
{
  TGraphAsymmErrors *g_rhoCentKstar_0020_stat = new TGraphAsymmErrors();
  g_rhoCentKstar_0020_stat->SetName("g_rhoCentKstar_0020_stat");
  TGraphAsymmErrors *g_rhoCentKstar_0020_sys = new TGraphAsymmErrors();
  g_rhoCentKstar_0020_sys->SetName("g_rhoCentKstar_0020_sys");

  // 39 GeV
  g_rhoCentKstar_0020_stat->SetPoint(0,39,kstar_rho00_Cent_39_central[0]);
  g_rhoCentKstar_0020_stat->SetPointError(0,0.0,0.0,kstar_rho00_Cent_39_central_stat[0],kstar_rho00_Cent_39_central_stat[0]);
  g_rhoCentKstar_0020_sys->SetPoint(0,39,kstar_rho00_Cent_39_central[0]);
  g_rhoCentKstar_0020_sys->SetPointError(0,0.0,0.0,kstar_rho00_Cent_39_central_syst[0],kstar_rho00_Cent_39_central_syst[0]);

  // 54 GeV
  g_rhoCentKstar_0020_stat->SetPoint(1,54.4,kstar_rho00_Cent_54_central[0]);
  g_rhoCentKstar_0020_stat->SetPointError(1,0.0,0.0,kstar_rho00_Cent_54_central_stat[0],kstar_rho00_Cent_54_central_stat[0]);
  g_rhoCentKstar_0020_sys->SetPoint(1,54.4,kstar_rho00_Cent_54_central[0]);
  g_rhoCentKstar_0020_sys->SetPointError(1,0.0,0.0,kstar_rho00_Cent_54_central_syst[0],kstar_rho00_Cent_54_central_syst[0]);

  // 200 GeV
  g_rhoCentKstar_0020_stat->SetPoint(2,200,kstar_rho00_Cent_200_central[0]);
  g_rhoCentKstar_0020_stat->SetPointError(2,0.0,0.0,kstar_rho00_Cent_200_central_stat[0],kstar_rho00_Cent_200_central_stat[0]);
  g_rhoCentKstar_0020_sys->SetPoint(2,200,kstar_rho00_Cent_200_central[0]);
  g_rhoCentKstar_0020_sys->SetPointError(2,0.0,0.0,kstar_rho00_Cent_200_central_syst[0],kstar_rho00_Cent_200_central_syst[0]);

  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_Cent_0020.root","RECREATE");
  File_OutPut->cd();
  g_rhoCentKstar_0020_stat->Write();
  g_rhoCentKstar_0020_sys->Write();
  File_OutPut->Close();
}

#include <TGraphAsymmErrors.h>

void makeData_rho00EpPhi_2060()
{
  TGraphAsymmErrors *g_rho2nd_stat = new TGraphAsymmErrors();
  g_rho2nd_stat->SetName("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rho2nd_sys = new TGraphAsymmErrors();
  g_rho2nd_sys->SetName("rho00_2ndEP_energy_sys");

  const double energy[6]  = {11.5, 19.6, 27.0, 39.0, 62.4, 200.0};
  const double rho[6]     = {0.307902,0.33548,0.340277,0.328164,0.320273,0.32944};
  const double errStat[6] = {0.0140411,0.00755227,0.00254671,0.00288413,0.00491966,0.00167392};
  const double errSys[6]  = {0.0191246,0.00535029,0.00328893,0.000990129,0.00309557,0.00188936};

  for(int iEnergy = 0; iEnergy < 6; ++iEnergy)
  {
    g_rho2nd_stat->SetPoint(iEnergy,energy[iEnergy],rho[iEnergy]);
    g_rho2nd_stat->SetPointError(iEnergy,0.0,0.0,errStat[iEnergy],errStat[iEnergy]);
    g_rho2nd_sys->SetPoint(iEnergy,energy[iEnergy],rho[iEnergy]);
    g_rho2nd_sys->SetPointError(iEnergy,0.0,0.0,errSys[iEnergy],errSys[iEnergy]);
  }


  TFile *File_OutPut = new TFile("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/NewF_JHChen/rho00_stat_sys_Xaxis.root","RECREATE");
  File_OutPut->cd();
  g_rho2nd_stat->Write();
  g_rho2nd_sys->Write();
  File_OutPut->Close();
}

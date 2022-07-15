#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig4()
{
  TFile *File_Input = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *g_rhoPhi_1st_stat_Laxis = (TGraphAsymmErrors*)File_Input->Get("rho00_1stEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_1st_sys_Laxis  = (TGraphAsymmErrors*)File_Input->Get("rho00_1stEP_energy_sys");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat_Laxis = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys_Laxis  = (TGraphAsymmErrors*)File_Input->Get("rho00_2ndEP_energy_sys");

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig4.txt";
  ofstream file_output;
  file_output.open(outputfile.c_str());

  file_output << "$\\sqrt{s_{NN}} (GeV)$" << endl;
  file_output << "2" << endl;
  file_output << "$\\rho_{00}$" << endl;
  file_output << "$20-60\%$ centrality" << endl;
  file_output << "$\\rho_{00}$" << endl;
  file_output << "$20-60\%$ centrality" << endl;
  file_output << "no" << endl;
  file_output << "none" << endl;
  file_output << "none" << endl;
  file_output << "2" << endl;
  file_output << "stat." << endl;
  file_output << "symmetric" << endl;
  file_output << "syst." << endl;
  file_output << "symmetric" << endl;

  const double beamEnergy[6] = {11.5,19.6,27.0,39.0,62.4,200.0};

  // 1st EP
  for(int iEnergy = 0; iEnergy < g_rhoPhi_1st_stat_Laxis->GetN(); ++iEnergy)
  {
    double energy, rho;
    g_rhoPhi_1st_stat_Laxis->GetPoint(iEnergy,energy,rho);
    double errStat = g_rhoPhi_1st_stat_Laxis->GetErrorYhigh(iEnergy);
    double errSyst = g_rhoPhi_1st_sys_Laxis->GetErrorYhigh(iEnergy);

    if(errStat > 0)
    {
      cout << beamEnergy[iEnergy] << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << beamEnergy[iEnergy] << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;

  // 2nd EP
  for(int iEnergy = 0; iEnergy < g_rhoPhi_2nd_stat_Laxis->GetN(); ++iEnergy)
  {
    double energy, rho;
    g_rhoPhi_2nd_stat_Laxis->GetPoint(iEnergy,energy,rho);
    double errStat = g_rhoPhi_2nd_stat_Laxis->GetErrorYhigh(iEnergy);
    double errSyst = g_rhoPhi_2nd_sys_Laxis->GetErrorYhigh(iEnergy);

    if(errStat > 0)
    {
      cout << energy << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << energy << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}

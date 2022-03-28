#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig9_phi()
{
  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rhoCent_0020_LXaxis.root");
  TGraphAsymmErrors *g_rhoPhi_1st_stat = (TGraphAsymmErrors*)File_InputPhi->Get("g_rho1st_0020_stat");
  TGraphAsymmErrors *g_rhoPhi_1st_sys = (TGraphAsymmErrors*)File_InputPhi->Get("g_rho1st_0020_sys");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat = (TGraphAsymmErrors*)File_InputPhi->Get("g_rho2nd_0020_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys = (TGraphAsymmErrors*)File_InputPhi->Get("g_rho2nd_0020_sys");

  const double beamEnergy[4] = {27.0,39.0,62.4,200.0};

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig9_phi.txt";
  ofstream file_output;
  file_output.open(outputfile.c_str());

  file_output << "$\\sqrt{s_{NN}} (GeV)$" << endl;
  file_output << "2" << endl;
  file_output << "$\\rho_{00}$" << endl;
  file_output << "1st EP" << endl;
  file_output << "$\\rho_{00}$" << endl;
  file_output << "2nd EP" << endl;
  file_output << "no" << endl;
  file_output << "none" << endl;
  file_output << "none" << endl;
  file_output << "2" << endl;
  file_output << "stat." << endl;
  file_output << "symmetric" << endl;
  file_output << "syst." << endl;
  file_output << "symmetric" << endl;

  // 1st EP
  for(int iEnergy = 0; iEnergy < g_rhoPhi_1st_stat->GetN(); ++iEnergy)
  {
    double energy, rho;
    g_rhoPhi_1st_stat->GetPoint(iEnergy,energy,rho);
    double errStat = g_rhoPhi_1st_stat->GetErrorYhigh(iEnergy);
    double errSyst = g_rhoPhi_1st_sys->GetErrorYhigh(iEnergy);

    if(errStat > 0)
    {
      cout << beamEnergy[iEnergy] << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << beamEnergy[iEnergy] << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;

  // 2nd EP
  for(int iEnergy = 0; iEnergy < g_rhoPhi_2nd_stat->GetN(); ++iEnergy)
  {
    double energy, rho;
    g_rhoPhi_2nd_stat->GetPoint(iEnergy,energy,rho);
    double errStat = g_rhoPhi_2nd_stat->GetErrorYhigh(iEnergy);
    double errSyst = g_rhoPhi_2nd_sys->GetErrorYhigh(iEnergy);

    if(errStat > 0)
    {
      cout << energy << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << energy << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}

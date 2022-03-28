#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig8d()
{
  string inputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/rhoCent_200GeV_Run14_LXaxis.root";
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TGraphAsymmErrors *g_rhoPhiCent_1st_stat = (TGraphAsymmErrors*)File_InPut->Get("g_rhoCent_200GeV_Run14_1st_stat");
  TGraphAsymmErrors *g_rhoPhiCent_1st_sys  = (TGraphAsymmErrors*)File_InPut->Get("g_rhoCent_200GeV_Run14_1st_sys");
  TGraphAsymmErrors *g_rhoPhiCent_2nd_stat = (TGraphAsymmErrors*)File_InPut->Get("g_rhoCent_200GeV_Run14_2nd_stat");
  TGraphAsymmErrors *g_rhoPhiCent_2nd_sys  = (TGraphAsymmErrors*)File_InPut->Get("g_rhoCent_200GeV_Run14_2nd_sys");

  const int centStart[9] = {0,5,10,20,30,40,50,60,70};
  const int centStop[9]  = {5,10,20,30,40,50,60,70,80};

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig8d.txt";
  ofstream file_output;
  file_output.open(outputfile.c_str());

  file_output << "Centrality (\%)" << endl;
  file_output << "2" << endl;
  file_output << "$\\rho_{00}$" << endl;
  file_output << "$1^{st}$-order EP" << endl;
  file_output << "$\\rho_{00}$" << endl;
  file_output << "$2^{nd}$-order EP" << endl;
  file_output << "yes" << endl;
  file_output << "none" << endl;
  file_output << "none" << endl;
  file_output << "2" << endl;
  file_output << "stat." << endl;
  file_output << "symmetric" << endl;
  file_output << "syst." << endl;
  file_output << "symmetric" << endl;

  // 1st
  for(int iCent = 0; iCent < g_rhoPhiCent_1st_stat->GetN(); ++iCent)
  {
    double centrality, rho;
    g_rhoPhiCent_1st_stat->GetPoint(iCent,centrality,rho);
    double errStat = g_rhoPhiCent_1st_stat->GetErrorYhigh(iCent);
    double errSyst = g_rhoPhiCent_1st_sys->GetErrorYhigh(iCent);

    if(errStat > 0)
    {
      cout << centStart[iCent] << '\t' << centStop[iCent] << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << centStart[iCent] << '\t' << centStop[iCent] << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;

  // 2nd
  for(int iCent = 0; iCent < g_rhoPhiCent_2nd_stat->GetN(); ++iCent)
  {
    double centrality, rho;
    g_rhoPhiCent_2nd_stat->GetPoint(iCent,centrality,rho);
    double errStat = g_rhoPhiCent_2nd_stat->GetErrorYhigh(iCent);
    double errSyst = g_rhoPhiCent_2nd_sys->GetErrorYhigh(iCent);

    if(errStat > 0)
    {
      cout << centStart[iCent] << '\t' << centStop[iCent] << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << centStart[iCent] << '\t' << centStop[iCent] << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}

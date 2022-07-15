#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig8g()
{
  string inputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/data_Kstar_rho00_Cent_Nov22_2021.root";
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TGraphAsymmErrors *g_rhoCent_2nd_stat = (TGraphAsymmErrors*)File_InPut->Get("g_rhoCent_200GeV_2nd_stat");
  TGraphAsymmErrors *g_rhoCent_2nd_sys  = (TGraphAsymmErrors*)File_InPut->Get("g_rhoCent_200GeV_2nd_sys");

  const int centStart[5] = {60,40,20,10,0};
  const int centStop[5]  = {80,60,40,20,10};

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig8g.txt";
  ofstream file_output;
  file_output.open(outputfile.c_str());

  file_output << "Centrality (\%)" << endl;
  file_output << "1" << endl;
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

  // 2nd
  for(int iCent = 0; iCent < g_rhoCent_2nd_stat->GetN(); ++iCent)
  {
    double centrality, rho;
    g_rhoCent_2nd_stat->GetPoint(iCent,centrality,rho);
    double errStat = g_rhoCent_2nd_stat->GetErrorYhigh(iCent);
    double errSyst = g_rhoCent_2nd_sys->GetErrorYhigh(iCent);

    if(errStat > 0)
    {
      cout << centStart[iCent] << '\t' << centStop[iCent] << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
      file_output << centStart[iCent] << '\t' << centStop[iCent] << '\t' << rho << '\t' << errStat << '\t' << errSyst << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}

#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig1d()
{
  TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/Kstar_Fig_signal_54_rap0p1_Oct2021.root");
  TH1F *h_mYieldsKstar = (TH1F*)File_InputKstar->Get("hist_Raw_Yield_theta_54GeV_Pt_2.00_2.50")->Clone("h_mYieldsKstar");
  h_mYieldsKstar->Draw("pE");

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig1d.txt";
  ofstream file_output;
  file_output.open(outputfile.c_str());

  file_output << "$\\left|cos(\\theta*)\\right|$" << endl;
  file_output << "1" << endl;
  file_output << "Raw Yields" << endl;
  file_output << "$20-60\%$ centrality" << endl;
  file_output << "yes" << endl;
  file_output << "none" << endl;
  file_output << "none" << endl;
  file_output << "1" << endl;
  file_output << "stat." << endl;
  file_output << "symmetric" << endl;

  for(int iBin = 0; iBin < h_mYieldsKstar->GetNbinsX(); ++iBin)
  {
    double binLowEdge  = h_mYieldsKstar->GetBinLowEdge(iBin+1);
    double binWitdh    = h_mYieldsKstar->GetBinWidth(iBin+1);
    double binHighEdge = binLowEdge+binWitdh;
    double binContent  = h_mYieldsKstar->GetBinContent(iBin+1);
    double binError    = h_mYieldsKstar->GetBinError(iBin+1);
    if(binError > 0)
    {
      cout << binLowEdge << '\t' << binHighEdge << '\t' << binContent << '\t' << binError << endl;
      file_output << binLowEdge << '\t' << binHighEdge << '\t' << binContent << '\t' << binError << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}

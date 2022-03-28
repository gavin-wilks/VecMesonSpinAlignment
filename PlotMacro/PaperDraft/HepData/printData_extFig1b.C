#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig1b()
{
  TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Kstar/Kstar_Fig_signal_54_rap0p1_Oct2021.root");
  TH1D *h_mMassKstar = (TH1D*)File_InputKstar->Get("ThetaIntegratedSignal_54GeV_PtBin2")->Clone("h_mMassKstar");
  h_mMassKstar->Draw("pE");

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig1b.txt";
  ofstream file_output;
  file_output.open(outputfile.c_str());

  file_output << "$M(\\pi^{\\pm},K^{\\mp}) (GeV/c^{2})$" << endl;
  file_output << "1" << endl;
  file_output << "Raw Yields" << endl;
  file_output << "$20-60\%$ centrality" << endl;
  file_output << "yes" << endl;
  file_output << "none" << endl;
  file_output << "none" << endl;
  file_output << "1" << endl;
  file_output << "stat." << endl;
  file_output << "symmetric" << endl;

  for(int iBin = 0; iBin < h_mMassKstar->GetNbinsX(); ++iBin)
  {
    double binLowEdge  = h_mMassKstar->GetBinLowEdge(iBin+1);
    double binWitdh    = h_mMassKstar->GetBinWidth(iBin+1);
    double binHighEdge = binLowEdge+binWitdh;
    double binContent  = h_mMassKstar->GetBinContent(iBin+1);
    double binError    = h_mMassKstar->GetBinError(iBin+1);
    if(binError > 0 && binLowEdge > 0.73 && binLowEdge < 1.08)
    {
      cout << binLowEdge << '\t' << binHighEdge << '\t' << binContent << '\t' << binError << endl;
      file_output << binLowEdge << '\t' << binHighEdge << '\t' << binContent << '\t' << binError << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}

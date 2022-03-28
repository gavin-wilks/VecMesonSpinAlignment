#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>

void printData_extFig1a()
{
  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/Phi/Second_ActualFit_27GeV_pT_2.root");
  TH1D *h_mMassPhi = (TH1D*)File_InputPhi->Get("imass_all_pt2")->Clone("h_mMassPhi");
  h_mMassPhi->Draw("pE");

  const long numOfEvent = 229186355;
  h_mMassPhi->Scale(1.0/numOfEvent);

  string outputfile = "/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/PaperDraft/Nature/HepData/dataExtFig1a.txt";
  ofstream file_output;
  file_output.open(outputfile.c_str());

  file_output << "$M(K^{+},K^{-}) (GeV/c^{2})$" << endl;
  file_output << "1" << endl;
  file_output << "Raw Yields" << endl;
  file_output << "$20-60\%$ centrality" << endl;
  file_output << "yes" << endl;
  file_output << "none" << endl;
  file_output << "none" << endl;
  file_output << "1" << endl;
  file_output << "stat." << endl;
  file_output << "symmetric" << endl;

  for(int iBin = 0; iBin < h_mMassPhi->GetNbinsX(); ++iBin)
  {
    double binLowEdge  = h_mMassPhi->GetBinLowEdge(iBin+1);
    double binWitdh    = h_mMassPhi->GetBinWidth(iBin+1);
    double binHighEdge = binLowEdge+binWitdh;
    double binContent  = h_mMassPhi->GetBinContent(iBin+1);
    double binError    = h_mMassPhi->GetBinError(iBin+1);
    if(binError > 0)
    {
      // cout << binLowEdge << '\t' << binHighEdge << '\t' << binContent << '\t' << binError << endl;
      file_output << binLowEdge << '\t' << binHighEdge << '\t' << binContent << '\t' << binError << endl;
    }
  }
  file_output << "***" << endl;

  file_output.close();
}
